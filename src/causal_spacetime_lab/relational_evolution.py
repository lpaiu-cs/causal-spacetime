"""Transport-gauge relational spatial-evolution utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class PairDistanceOrderSignature:
    """Ordinal pair-distance signature for one slice."""

    object_pairs: NDArray[np.int_]
    pair_distance_values: NDArray[np.float64]
    pair_order_ranks: NDArray[np.float64]


def unordered_object_pairs(object_ids: ArrayLike) -> NDArray[np.int_]:
    """Return all unordered object-id pairs."""

    ids = np.asarray(object_ids, dtype=int)
    if ids.ndim != 1:
        raise ValueError("object_ids must be one-dimensional")
    unique = np.unique(ids)
    if unique.size < 2:
        return np.empty((0, 2), dtype=int)
    left, right = np.triu_indices(unique.size, k=1)
    return np.column_stack((unique[left], unique[right])).astype(int)


def _ranks_from_values(values: NDArray[np.float64], tolerance: float) -> NDArray:
    ranks = np.full(values.shape, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(values))
    if finite_indices.size == 0:
        return ranks
    ordered = finite_indices[np.argsort(values[finite_indices], kind="mergesort")]
    rank = 0
    ranks[ordered[0]] = float(rank)
    last_value = values[ordered[0]]
    for index in ordered[1:]:
        if abs(values[index] - last_value) > tolerance:
            rank += 1
            last_value = values[index]
        ranks[index] = float(rank)
    return ranks


def pair_distance_order_signature(
    positions_by_object: dict[int, float],
    object_pairs: NDArray[np.int_] | None = None,
    tolerance: float = 0.0,
) -> PairDistanceOrderSignature:
    """Compute a slice-local persistent-object pair-distance order signature."""

    if not positions_by_object:
        pairs = np.empty((0, 2), dtype=int)
        values = np.empty(0, dtype=float)
        return PairDistanceOrderSignature(pairs, values, values.copy())
    if object_pairs is None:
        pairs = unordered_object_pairs(np.asarray(sorted(positions_by_object)))
    else:
        pairs = np.asarray(object_pairs, dtype=int)
        if pairs.ndim != 2 or pairs.shape[1] != 2:
            raise ValueError("object_pairs must have shape (m, 2)")
    values = np.empty(pairs.shape[0], dtype=float)
    for index, (left, right) in enumerate(pairs):
        try:
            values[index] = abs(
                float(positions_by_object[int(left)])
                - float(positions_by_object[int(right)])
            )
        except KeyError:
            values[index] = np.nan
    ranks = _ranks_from_values(values, float(tolerance))
    return PairDistanceOrderSignature(
        object_pairs=pairs.astype(int, copy=False),
        pair_distance_values=values,
        pair_order_ranks=ranks,
    )


def _pair_key(pair: ArrayLike) -> tuple[int, int]:
    left, right = np.asarray(pair, dtype=int).tolist()
    return (int(min(left, right)), int(max(left, right)))


def _signature_by_pair(
    signature: PairDistanceOrderSignature,
) -> dict[tuple[int, int], tuple[float, float]]:
    return {
        _pair_key(pair): (float(value), float(rank))
        for pair, value, rank in zip(
            signature.object_pairs,
            signature.pair_distance_values,
            signature.pair_order_ranks,
            strict=True,
        )
    }


def pair_order_comparison_matrix(
    signature: PairDistanceOrderSignature,
) -> NDArray[np.int_]:
    """Return pair-pair order signs for a slice signature."""

    ranks = np.asarray(signature.pair_order_ranks, dtype=float)
    matrix = np.zeros((ranks.size, ranks.size), dtype=int)
    for i in range(ranks.size):
        for j in range(ranks.size):
            if not np.isfinite(ranks[i]) or not np.isfinite(ranks[j]):
                continue
            if ranks[i] < ranks[j]:
                matrix[i, j] = -1
            elif ranks[i] > ranks[j]:
                matrix[i, j] = 1
    return matrix


def _aligned_pair_ranks(
    signature_a: PairDistanceOrderSignature,
    signature_b: PairDistanceOrderSignature,
) -> tuple[list[tuple[int, int]], NDArray[np.float64], NDArray[np.float64]]:
    by_a = _signature_by_pair(signature_a)
    by_b = _signature_by_pair(signature_b)
    common = sorted(set(by_a) & set(by_b))
    ranks_a = np.asarray([by_a[pair][1] for pair in common], dtype=float)
    ranks_b = np.asarray([by_b[pair][1] for pair in common], dtype=float)
    return common, ranks_a, ranks_b


def signature_order_disagreement(
    signature_a: PairDistanceOrderSignature,
    signature_b: PairDistanceOrderSignature,
    ignore_ties: bool = True,
) -> float:
    """Compare pair-pair order relations between two slice signatures."""

    _, ranks_a, ranks_b = _aligned_pair_ranks(signature_a, signature_b)
    if ranks_a.size < 2:
        return 0.0
    total = 0
    changed = 0
    for i in range(ranks_a.size):
        for j in range(i + 1, ranks_a.size):
            if not (
                np.isfinite(ranks_a[i])
                and np.isfinite(ranks_a[j])
                and np.isfinite(ranks_b[i])
                and np.isfinite(ranks_b[j])
            ):
                continue
            sign_a = int(np.sign(ranks_a[i] - ranks_a[j]))
            sign_b = int(np.sign(ranks_b[i] - ranks_b[j]))
            if ignore_ties and (sign_a == 0 or sign_b == 0):
                continue
            total += 1
            changed += int(sign_a != sign_b)
    return float(changed / total) if total else 0.0


def signature_change_events(
    signature_a: PairDistanceOrderSignature,
    signature_b: PairDistanceOrderSignature,
) -> list[dict[str, object]]:
    """Return pair-pair comparisons whose order changes between signatures."""

    pairs, ranks_a, ranks_b = _aligned_pair_ranks(signature_a, signature_b)
    rows: list[dict[str, object]] = []
    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            if not (
                np.isfinite(ranks_a[i])
                and np.isfinite(ranks_a[j])
                and np.isfinite(ranks_b[i])
                and np.isfinite(ranks_b[j])
            ):
                continue
            sign_a = int(np.sign(ranks_a[i] - ranks_a[j]))
            sign_b = int(np.sign(ranks_b[i] - ranks_b[j]))
            if sign_a != sign_b:
                rows.append(
                    {
                        "pair_a": pairs[i],
                        "pair_b": pairs[j],
                        "order_a": sign_a,
                        "order_b": sign_b,
                    }
                )
    return rows


def relational_shape_history(
    positions_by_slice: dict[int, dict[int, float]],
    tolerance: float = 0.0,
) -> dict[int, PairDistanceOrderSignature]:
    """Build pair-distance order signatures for every slice."""

    return {
        int(label): pair_distance_order_signature(positions, tolerance=tolerance)
        for label, positions in sorted(positions_by_slice.items())
    }


def relational_change_rate_between_slices(
    history: dict[int, PairDistanceOrderSignature],
) -> list[dict[str, float]]:
    """Return order-disagreement rates between consecutive slice labels."""

    labels = sorted(history)
    rows: list[dict[str, float]] = []
    for left, right in zip(labels[:-1], labels[1:], strict=False):
        changes = signature_change_events(history[left], history[right])
        rows.append(
            {
                "slice_a": float(left),
                "slice_b": float(right),
                "order_disagreement": signature_order_disagreement(
                    history[left],
                    history[right],
                ),
                "change_count": float(len(changes)),
            }
        )
    return rows


def _signature_from_values(
    object_pairs: NDArray[np.int_],
    values: NDArray[np.float64],
    tolerance: float = 0.0,
) -> PairDistanceOrderSignature:
    return PairDistanceOrderSignature(
        object_pairs=np.asarray(object_pairs, dtype=int),
        pair_distance_values=np.asarray(values, dtype=float),
        pair_order_ranks=_ranks_from_values(np.asarray(values, dtype=float), tolerance),
    )


def apply_per_slice_monotone_distance_transform_to_history(
    history: dict[int, PairDistanceOrderSignature],
    transform_name: str = "square",
) -> dict[int, PairDistanceOrderSignature]:
    """Apply a positive monotone transform to each slice's pair distances."""

    transformed: dict[int, PairDistanceOrderSignature] = {}
    for label, signature in history.items():
        values = np.asarray(signature.pair_distance_values, dtype=float)
        if transform_name == "square":
            new_values = values**2
        elif transform_name == "sqrt":
            new_values = np.sqrt(np.maximum(values, 0.0))
        elif transform_name == "log1p":
            new_values = np.log1p(np.maximum(values, 0.0))
        else:
            raise ValueError(f"unknown transform_name: {transform_name}")
        transformed[int(label)] = _signature_from_values(
            signature.object_pairs,
            new_values,
        )
    return transformed


def apply_per_slice_affine_position_gauge(
    positions_by_slice: dict[int, dict[int, float]],
    seed: int | None = None,
    allow_reflection: bool = True,
    scale_range: tuple[float, float] = (0.5, 2.0),
    shift_range: tuple[float, float] = (-1.0, 1.0),
) -> dict[int, dict[int, float]]:
    """Apply independent per-slice affine/reflection gauges to positions."""

    scale_min, scale_max = float(scale_range[0]), float(scale_range[1])
    shift_min, shift_max = float(shift_range[0]), float(shift_range[1])
    if scale_min <= 0.0 or scale_max < scale_min:
        raise ValueError("scale_range must be positive and ordered")
    if shift_max < shift_min:
        raise ValueError("shift_range must be ordered")
    rng = np.random.default_rng(seed)
    output: dict[int, dict[int, float]] = {}
    for label, positions in sorted(positions_by_slice.items()):
        scale = rng.uniform(scale_min, scale_max)
        shift = rng.uniform(shift_min, shift_max)
        reflection = rng.choice([-1.0, 1.0]) if allow_reflection else 1.0
        output[int(label)] = {
            int(object_id): float(reflection * scale * value + shift)
            for object_id, value in positions.items()
        }
    return output


def compare_histories_order_disagreement(
    history_a: dict[int, PairDistanceOrderSignature],
    history_b: dict[int, PairDistanceOrderSignature],
) -> float:
    """Compare matching slice signatures in two histories."""

    common = sorted(set(history_a) & set(history_b))
    if not common:
        return 0.0
    values = [
        signature_order_disagreement(history_a[label], history_b[label])
        for label in common
    ]
    return float(np.mean(values)) if values else 0.0
