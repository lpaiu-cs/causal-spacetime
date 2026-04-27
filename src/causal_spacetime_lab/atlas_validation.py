"""Validation diagnostics for overlapping reconstructed observer charts."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def sample_event_pairs(
    indices: ArrayLike,
    num_pairs: int,
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Sample ordered index pairs without self-pairs."""

    index_array = np.asarray(indices, dtype=int)
    if index_array.ndim != 1:
        raise ValueError("indices must be one-dimensional")
    if index_array.size < 2:
        return np.empty((0, 2), dtype=int)
    pair_count = int(num_pairs)
    if pair_count < 0:
        raise ValueError("num_pairs must be non-negative")
    if pair_count == 0:
        return np.empty((0, 2), dtype=int)

    rng = np.random.default_rng(seed)
    first = rng.choice(index_array, size=pair_count, replace=True)
    second = rng.choice(index_array, size=pair_count, replace=True)
    same = first == second
    while np.any(same):
        second[same] = rng.choice(
            index_array,
            size=np.count_nonzero(same),
            replace=True,
        )
        same = first == second
    return np.column_stack((first, second)).astype(int, copy=False)


def _as_coord_array(coords: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(coords, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("coords must have shape (n, 2)")
    return array


def minkowski_interval_squared_from_coords(
    coords: ArrayLike,
    pairs: ArrayLike,
) -> NDArray[np.float64]:
    """Return ``s^2 = dt^2 - dx^2`` for coordinate pairs."""

    coord_array = _as_coord_array(coords)
    pair_array = np.asarray(pairs, dtype=int)
    if pair_array.ndim != 2 or pair_array.shape[1] != 2:
        raise ValueError("pairs must have shape (n_pairs, 2)")
    if pair_array.size == 0:
        return np.empty(0, dtype=np.float64)
    if np.any((pair_array < 0) | (pair_array >= coord_array.shape[0])):
        raise IndexError("pairs contain coordinate indices out of range")

    delta = coord_array[pair_array[:, 1]] - coord_array[pair_array[:, 0]]
    return (delta[:, 0] ** 2 - delta[:, 1] ** 2).astype(np.float64, copy=False)


def chart_interval_disagreement(
    coords_a: ArrayLike,
    coords_b: ArrayLike,
    accessible_mask: ArrayLike,
    num_pairs: int,
    seed: int | None = None,
) -> dict[str, float]:
    """Compare invariant intervals assigned by two reconstructed charts."""

    chart_a = _as_coord_array(coords_a)
    chart_b = _as_coord_array(coords_b)
    if chart_a.shape != chart_b.shape:
        raise ValueError("coords_a and coords_b must have the same shape")
    accessible = np.asarray(accessible_mask, dtype=bool)
    if accessible.ndim != 1 or accessible.shape[0] != chart_a.shape[0]:
        raise ValueError("accessible_mask must have length matching coordinates")

    valid = accessible & np.all(np.isfinite(chart_a), axis=1)
    valid &= np.all(np.isfinite(chart_b), axis=1)
    indices = np.flatnonzero(valid)
    pairs = sample_event_pairs(indices, num_pairs, seed=seed)
    if pairs.shape[0] == 0:
        return {
            "pair_count": 0.0,
            "interval_rmse": float("nan"),
            "interval_mae": float("nan"),
            "interval_bias": float("nan"),
        }

    intervals_a = minkowski_interval_squared_from_coords(chart_a, pairs)
    intervals_b = minkowski_interval_squared_from_coords(chart_b, pairs)
    diff = intervals_b - intervals_a
    return {
        "pair_count": float(pairs.shape[0]),
        "interval_rmse": float(np.sqrt(np.mean(diff**2))),
        "interval_mae": float(np.mean(np.abs(diff))),
        "interval_bias": float(np.mean(diff)),
    }
