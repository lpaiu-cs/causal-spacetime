"""Distance-order diagnostics restricted to observer-selected slices."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.ordinal import (
    order_inversion_rate,
    pair_distance_values_1d,
)


def _as_pairs(pairs: ArrayLike) -> NDArray[np.int_]:
    array = np.asarray(pairs, dtype=int)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("pairs must have shape (m, 2)")
    return array


def pair_slice_labels(
    pairs: ArrayLike,
    slice_labels: ArrayLike,
) -> NDArray[np.int_]:
    """Return each pair's shared slice label, or ``-1`` if not same-slice."""

    pair_array = _as_pairs(pairs)
    labels = np.asarray(slice_labels, dtype=int)
    if labels.ndim != 1:
        raise ValueError("slice_labels must be one-dimensional")
    if np.any((pair_array < 0) | (pair_array >= labels.size)):
        raise IndexError("pairs contain an index out of range")
    left = labels[pair_array[:, 0]]
    right = labels[pair_array[:, 1]]
    result = np.full(pair_array.shape[0], -1, dtype=int)
    same = (left >= 0) & (left == right)
    result[same] = left[same]
    return result


def sliced_pair_distance_order_inversion_rate(
    true_positions: ArrayLike,
    estimated_positions: ArrayLike,
    pairs: ArrayLike,
    *,
    tolerance: float = 0.0,
) -> float:
    """Compare pair-distance order among supplied same-slice pairs."""

    pair_array = _as_pairs(pairs)
    if pair_array.shape[0] < 2:
        return float("nan")
    true_distances = pair_distance_values_1d(true_positions, pair_array)
    estimated_distances = pair_distance_values_1d(estimated_positions, pair_array)
    return order_inversion_rate(
        true_distances,
        estimated_distances,
        ignore_ties=True,
        tolerance=tolerance,
    )


def quadruplet_constraints_from_sliced_pair_distances(
    pair_indices: ArrayLike,
    pair_distance_values: ArrayLike,
    slice_labels_for_pairs: ArrayLike,
    num_constraints: int,
    seed: int | None = None,
    tolerance: float = 0.0,
    allow_cross_slice_comparisons: bool = False,
) -> NDArray[np.int_]:
    """Sample constraints from pair distances, optionally within each slice only.

    With ``allow_cross_slice_comparisons=False``, both compared pairs must come
    from the same slice label. Cross-slice comparisons require an additional
    calibration assumption and are disabled by default.
    """

    pairs = _as_pairs(pair_indices)
    values = np.asarray(pair_distance_values, dtype=float)
    labels = np.asarray(slice_labels_for_pairs, dtype=int)
    if values.ndim != 1 or values.shape[0] != pairs.shape[0]:
        raise ValueError("pair_distance_values must have shape (m,)")
    if labels.ndim != 1 or labels.shape[0] != pairs.shape[0]:
        raise ValueError("slice_labels_for_pairs must have shape (m,)")
    count = int(num_constraints)
    if count < 0:
        raise ValueError("num_constraints must be nonnegative")
    tol = float(tolerance)
    if tol < 0.0:
        raise ValueError("tolerance must be nonnegative")
    if count == 0:
        return np.empty((0, 4), dtype=int)

    valid_indices = np.flatnonzero(labels >= 0)
    if valid_indices.size < 2:
        return np.empty((0, 4), dtype=int)
    rng = np.random.default_rng(seed)
    rows = np.empty((count, 4), dtype=int)
    accepted = 0
    attempts = 0
    max_attempts = max(100_000, 200 * count)
    while accepted < count and attempts < max_attempts:
        batch = min(max(512, 4 * (count - accepted)), max_attempts - attempts)
        left = rng.choice(valid_indices, size=batch, replace=True)
        right = rng.choice(valid_indices, size=batch, replace=True)
        different = left != right
        if not allow_cross_slice_comparisons:
            different &= labels[left] == labels[right]
        left_smaller = values[left] + tol < values[right]
        right_smaller = values[right] + tol < values[left]
        left_rows = np.column_stack((pairs[left], pairs[right]))
        right_rows = np.column_stack((pairs[right], pairs[left]))
        selected = np.vstack(
            (
                left_rows[different & left_smaller],
                right_rows[different & right_smaller],
            )
        )
        take = min(selected.shape[0], count - accepted)
        if take:
            rows[accepted : accepted + take] = selected[:take]
            accepted += take
        attempts += batch
    return rows[:accepted].astype(int, copy=False)
