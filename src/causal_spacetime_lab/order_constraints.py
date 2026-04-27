"""Order-facing construction of finite ordinal constraints."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.ordinal_embedding import (
    sample_quadruplet_constraints_from_distance_values,
)


def quadruplet_constraints_from_order_values(
    values: ArrayLike,
    num_constraints: int,
    seed: int | None = None,
    tolerance: float = 0.0,
) -> NDArray[np.int_]:
    """Sample pair-distance constraints from scalar order values.

    The values are treated as a one-dimensional coordinate proxy. The returned
    rows ``(i, j, k, l)`` mean ``|v_i - v_j| < |v_k - v_l|``.
    """

    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError("values must have shape (n,)")
    distance_matrix = np.abs(array[:, None] - array[None, :])
    return sample_quadruplet_constraints_from_distance_values(
        distance_matrix,
        num_constraints,
        seed=seed,
        tolerance=tolerance,
    )


def quadruplet_constraints_from_pair_distance_order(
    pair_indices: ArrayLike,
    pair_order_values: ArrayLike,
    num_constraints: int,
    seed: int | None = None,
    tolerance: float = 0.0,
) -> NDArray[np.int_]:
    """Sample quadruplet constraints from explicit pair-distance order values."""

    pairs = np.asarray(pair_indices, dtype=int)
    values = np.asarray(pair_order_values, dtype=float)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError("pair_indices must have shape (m, 2)")
    if values.ndim != 1 or values.shape[0] != pairs.shape[0]:
        raise ValueError("pair_order_values must have shape (m,)")
    if np.any(pairs[:, 0] == pairs[:, 1]):
        raise ValueError("pair_indices must not contain self-pairs")
    count = int(num_constraints)
    if count < 0:
        raise ValueError("num_constraints must be nonnegative")
    tol = float(tolerance)
    if tol < 0.0:
        raise ValueError("tolerance must be nonnegative")
    if count == 0:
        return np.empty((0, 4), dtype=int)
    rng = np.random.default_rng(seed)
    rows = np.empty((count, 4), dtype=int)
    accepted = 0
    attempts = 0
    max_attempts = max(100_000, 100 * count)
    while accepted < count and attempts < max_attempts:
        batch = min(max(256, 4 * (count - accepted)), max_attempts - attempts)
        left = rng.integers(0, pairs.shape[0], size=batch)
        right = rng.integers(0, pairs.shape[0], size=batch)
        different = left != right
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
    if accepted < count:
        return rows[:accepted]
    return rows


def scalar_order_constraints_from_values(
    values: ArrayLike,
    num_constraints: int,
    seed: int | None = None,
    tolerance: float = 0.0,
) -> NDArray[np.int_]:
    """Sample scalar order constraints ``(i, j)`` meaning ``value_i < value_j``."""

    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError("values must have shape (n,)")
    count = int(num_constraints)
    if count < 0:
        raise ValueError("num_constraints must be nonnegative")
    tol = float(tolerance)
    if tol < 0.0:
        raise ValueError("tolerance must be nonnegative")
    rng = np.random.default_rng(seed)
    rows = np.empty((count, 2), dtype=int)
    accepted = 0
    attempts = 0
    max_attempts = max(100_000, 100 * count)
    while accepted < count and attempts < max_attempts:
        batch = min(max(256, 4 * (count - accepted)), max_attempts - attempts)
        left = rng.integers(0, array.shape[0], size=batch)
        right = rng.integers(0, array.shape[0], size=batch)
        different = left != right
        left_smaller = array[left] + tol < array[right]
        right_smaller = array[right] + tol < array[left]
        selected = np.vstack(
            (
                np.column_stack((left, right))[different & left_smaller],
                np.column_stack((right, left))[different & right_smaller],
            )
        )
        take = min(selected.shape[0], count - accepted)
        if take:
            rows[accepted : accepted + take] = selected[:take]
            accepted += take
        attempts += batch
    return rows[:accepted]


def quadruplet_constraints_from_successor_ticks(
    successor_tick_positions: ArrayLike,
    accessible_mask: ArrayLike,
    num_constraints: int,
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Return same-emission radar-return scalar order constraints.

    Earlier successor ticks encode a target-rank order relative to one observer
    protocol. This is not a pair-pair distance order among targets, so this
    helper returns two-column scalar constraints ``(i, j)`` meaning target ``i``
    is closer than target ``j`` under the common return-tick protocol.
    """

    ticks = np.asarray(successor_tick_positions, dtype=float)
    mask = np.asarray(accessible_mask, dtype=bool)
    if ticks.ndim != 1 or mask.ndim != 1 or ticks.shape[0] != mask.shape[0]:
        raise ValueError("successor ticks and accessible mask must be 1D arrays")
    values = np.full(ticks.shape, np.nan, dtype=float)
    values[mask] = ticks[mask]
    finite = np.isfinite(values)
    if not np.any(finite):
        return np.empty((0, 2), dtype=int)
    compact_indices = np.flatnonzero(finite)
    compact_constraints = scalar_order_constraints_from_values(
        values[finite],
        num_constraints,
        seed=seed,
    )
    if compact_constraints.size == 0:
        return compact_constraints
    return compact_indices[compact_constraints]
