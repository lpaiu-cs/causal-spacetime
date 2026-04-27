"""Protocol-dependent cross-slice transport and gauge utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.cross_slice import CrossSlicePredicateReport


@dataclass(frozen=True)
class SliceTransportRule:
    """Per-slice affine/reflection identification rule."""

    name: str
    description: str
    slice_labels: NDArray[np.int_]
    scale_by_slice: dict[int, float]
    shift_by_slice: dict[int, float]
    reflection_by_slice: dict[int, float]


def _as_1d(values: ArrayLike, name: str) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=float)
    if array.ndim == 2 and array.shape[1] == 1:
        array = array[:, 0]
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional or shape (n, 1)")
    return array


def _rule_covers_slice(rule: SliceTransportRule, slice_label: int) -> bool:
    label = int(slice_label)
    return (
        label in rule.scale_by_slice
        and label in rule.shift_by_slice
        and label in rule.reflection_by_slice
    )


def apply_slice_transport_rule(
    slice_local_coords: ArrayLike,
    slice_labels: ArrayLike,
    rule: SliceTransportRule,
) -> NDArray[np.float64]:
    """Apply a per-slice affine/reflection transport rule."""

    coords = _as_1d(slice_local_coords, "slice_local_coords")
    labels = np.asarray(slice_labels, dtype=int)
    if labels.ndim != 1 or labels.shape[0] != coords.shape[0]:
        raise ValueError("slice_labels must have same length as coordinates")
    transported = np.full(coords.shape, np.nan, dtype=float)
    for label in sorted(set(int(value) for value in labels if value >= 0)):
        if not _rule_covers_slice(rule, label):
            continue
        mask = labels == label
        transported[mask] = (
            rule.reflection_by_slice[label]
            * rule.scale_by_slice[label]
            * coords[mask]
            + rule.shift_by_slice[label]
        )
    return transported


def apply_random_slice_gauge_transforms(
    slice_local_coords: ArrayLike,
    slice_labels: ArrayLike,
    seed: int | None = None,
    allow_reflection: bool = True,
    scale_range: tuple[float, float] = (0.5, 2.0),
    shift_range: tuple[float, float] = (-2.0, 2.0),
) -> NDArray[np.float64]:
    """Apply independent random affine/reflection gauges per slice."""

    coords = _as_1d(slice_local_coords, "slice_local_coords")
    labels = np.asarray(slice_labels, dtype=int)
    if labels.ndim != 1 or labels.shape[0] != coords.shape[0]:
        raise ValueError("slice_labels must have same length as coordinates")
    scale_min, scale_max = float(scale_range[0]), float(scale_range[1])
    shift_min, shift_max = float(shift_range[0]), float(shift_range[1])
    if scale_min <= 0.0 or scale_max < scale_min:
        raise ValueError("scale_range must be positive and ordered")
    if shift_max < shift_min:
        raise ValueError("shift_range must be ordered")
    rng = np.random.default_rng(seed)
    transformed = np.full(coords.shape, np.nan, dtype=float)
    for label in sorted(set(int(value) for value in labels if value >= 0)):
        mask = labels == label
        scale = rng.uniform(scale_min, scale_max)
        shift = rng.uniform(shift_min, shift_max)
        reflection = rng.choice([-1.0, 1.0]) if allow_reflection else 1.0
        transformed[mask] = reflection * scale * coords[mask] + shift
    return transformed


def fit_1d_affine_from_anchors(
    local_values: ArrayLike,
    reference_values: ArrayLike,
) -> tuple[float, float, float]:
    """Fit ``reference ~= scale * local + offset`` from anchor values."""

    local = _as_1d(local_values, "local_values")
    reference = _as_1d(reference_values, "reference_values")
    if local.shape != reference.shape:
        raise ValueError("local_values and reference_values must have equal shape")
    finite = np.isfinite(local) & np.isfinite(reference)
    local = local[finite]
    reference = reference[finite]
    if local.size == 0:
        raise ValueError("at least one finite anchor is required")
    if local.size == 1:
        scale = 1.0
        offset = float(reference[0] - local[0])
    else:
        design = np.column_stack((local, np.ones(local.size)))
        scale, offset = np.linalg.lstsq(design, reference, rcond=None)[0]
        scale = float(scale)
        offset = float(offset)
    predicted = scale * local + offset
    rmse = float(np.sqrt(np.mean((predicted - reference) ** 2)))
    return scale, offset, rmse


def align_slices_by_anchor_points(
    slice_local_coords: ArrayLike,
    slice_labels: ArrayLike,
    anchor_indices_by_slice: dict[int, NDArray[np.int_]],
    anchor_reference_positions_by_slice: dict[int, NDArray[np.float64]],
) -> NDArray[np.float64]:
    """Align each slice with supplied anchor correspondences."""

    coords = _as_1d(slice_local_coords, "slice_local_coords")
    labels = np.asarray(slice_labels, dtype=int)
    if labels.ndim != 1 or labels.shape[0] != coords.shape[0]:
        raise ValueError("slice_labels must have same length as coordinates")
    aligned = np.full(coords.shape, np.nan, dtype=float)
    for label in sorted(set(int(value) for value in labels if value >= 0)):
        mask = labels == label
        anchor_indices = np.asarray(
            anchor_indices_by_slice.get(label, np.empty(0, dtype=int)),
            dtype=int,
        )
        references = np.asarray(
            anchor_reference_positions_by_slice.get(label, np.empty(0)),
            dtype=float,
        )
        if anchor_indices.size == 0:
            continue
        if anchor_indices.shape[0] != references.shape[0]:
            raise ValueError("anchor index/reference size mismatch")
        finite_anchors = np.isfinite(coords[anchor_indices]) & np.isfinite(references)
        if not np.any(finite_anchors):
            continue
        scale, offset, _ = fit_1d_affine_from_anchors(
            coords[anchor_indices][finite_anchors],
            references[finite_anchors],
        )
        aligned[mask] = scale * coords[mask] + offset
    return aligned


def same_position_under_transport(
    position_i: float,
    slice_i: int,
    position_j: float,
    slice_j: int,
    rule: SliceTransportRule,
    tolerance: float = 1e-6,
) -> CrossSlicePredicateReport:
    """Evaluate same-position as a transport-relative predicate."""

    if not _rule_covers_slice(rule, int(slice_i)) or not _rule_covers_slice(
        rule,
        int(slice_j),
    ):
        return CrossSlicePredicateReport(
            predicate_name="same_position",
            defined=False,
            reason="The supplied transport rule does not cover both slices.",
            value=None,
        )
    transformed_i = (
        rule.reflection_by_slice[int(slice_i)]
        * rule.scale_by_slice[int(slice_i)]
        * float(position_i)
        + rule.shift_by_slice[int(slice_i)]
    )
    transformed_j = (
        rule.reflection_by_slice[int(slice_j)]
        * rule.scale_by_slice[int(slice_j)]
        * float(position_j)
        + rule.shift_by_slice[int(slice_j)]
    )
    same = abs(transformed_i - transformed_j) <= float(tolerance)
    return CrossSlicePredicateReport(
        predicate_name="same_position",
        defined=True,
        reason=f"Evaluated under transport rule '{rule.name}'.",
        value=bool(same),
    )


def velocity_under_transport(
    positions: ArrayLike,
    slice_labels: ArrayLike,
    object_ids: ArrayLike,
    time_values_by_slice: dict[int, float],
    rule: SliceTransportRule,
) -> dict[int, dict[str, float]]:
    """Estimate object velocities after a transport rule has been chosen."""

    pos = _as_1d(positions, "positions")
    slices = np.asarray(slice_labels, dtype=int)
    objects = np.asarray(object_ids, dtype=int)
    if slices.shape != pos.shape or objects.shape != pos.shape:
        raise ValueError("positions, slice_labels, and object_ids must align")
    transported = apply_slice_transport_rule(pos, slices, rule)
    result: dict[int, dict[str, float]] = {}
    for object_id in sorted(set(int(value) for value in objects)):
        mask = objects == object_id
        valid = mask & np.isfinite(transported)
        valid &= np.asarray([int(s) in time_values_by_slice for s in slices])
        if np.count_nonzero(valid) < 2:
            continue
        times = np.asarray([time_values_by_slice[int(s)] for s in slices[valid]])
        values = transported[valid]
        order = np.argsort(times)
        times = times[order]
        values = values[order]
        dt = np.diff(times)
        dx = np.diff(values)
        finite = dt != 0.0
        if not np.any(finite):
            continue
        velocities = dx[finite] / dt[finite]
        result[object_id] = {
            "observation_count": float(values.size),
            "velocity_mean": float(np.mean(velocities)),
            "velocity_std": float(np.std(velocities)),
            "velocity_min": float(np.min(velocities)),
            "velocity_max": float(np.max(velocities)),
        }
    return result
