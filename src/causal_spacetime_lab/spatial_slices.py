"""Observer-relative simultaneity slice utilities from radar tick brackets."""

from __future__ import annotations

from collections import Counter

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _as_1d_int(values: ArrayLike, name: str) -> NDArray[np.int_]:
    array = np.asarray(values, dtype=int)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return array


def _as_1d_bool(values: ArrayLike, name: str) -> NDArray[np.bool_]:
    array = np.asarray(values, dtype=bool)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return array


def _validate_brackets(
    predecessor_tick_positions: ArrayLike,
    successor_tick_positions: ArrayLike,
    accessible_mask: ArrayLike,
) -> tuple[NDArray[np.int_], NDArray[np.int_], NDArray[np.bool_]]:
    predecessor = _as_1d_int(predecessor_tick_positions, "predecessor_tick_positions")
    successor = _as_1d_int(successor_tick_positions, "successor_tick_positions")
    accessible = _as_1d_bool(accessible_mask, "accessible_mask")
    if predecessor.shape != successor.shape or predecessor.shape != accessible.shape:
        raise ValueError("tick bracket arrays and accessible_mask must have same shape")
    return predecessor, successor, accessible


def radar_time_rank_from_tick_brackets(
    predecessor_tick_positions: ArrayLike,
    successor_tick_positions: ArrayLike,
    accessible_mask: ArrayLike,
) -> NDArray[np.int_]:
    """Return order-level radar-time ranks from tick bracket positions.

    The rank is ``predecessor + successor`` for accessible targets and ``-1``
    otherwise. It is an ordinal proxy for radar time and uses no clock labels.
    """

    predecessor, successor, accessible = _validate_brackets(
        predecessor_tick_positions,
        successor_tick_positions,
        accessible_mask,
    )
    ranks = np.full(predecessor.shape, -1, dtype=int)
    ranks[accessible] = predecessor[accessible] + successor[accessible]
    return ranks


def radar_distance_rank_from_tick_brackets(
    predecessor_tick_positions: ArrayLike,
    successor_tick_positions: ArrayLike,
    accessible_mask: ArrayLike,
) -> NDArray[np.int_]:
    """Return order-level radar-distance ranks from tick bracket positions."""

    predecessor, successor, accessible = _validate_brackets(
        predecessor_tick_positions,
        successor_tick_positions,
        accessible_mask,
    )
    ranks = np.full(predecessor.shape, -1, dtype=int)
    ranks[accessible] = successor[accessible] - predecessor[accessible]
    return ranks


def assign_slices_from_radar_time_rank(
    radar_time_ranks: ArrayLike,
    accessible_mask: ArrayLike,
    bin_width: int = 2,
) -> NDArray[np.int_]:
    """Assign observer-relative radar-time slice labels from ordinal ranks."""

    ranks = _as_1d_int(radar_time_ranks, "radar_time_ranks")
    accessible = _as_1d_bool(accessible_mask, "accessible_mask")
    if ranks.shape != accessible.shape:
        raise ValueError("radar_time_ranks and accessible_mask must have same shape")
    width = int(bin_width)
    if width < 1:
        raise ValueError("bin_width must be at least 1")
    labels = np.full(ranks.shape, -1, dtype=int)
    valid = accessible & (ranks >= 0)
    labels[valid] = ranks[valid] // width
    return labels


def slice_sizes(slice_labels: ArrayLike) -> dict[int, int]:
    """Return counts per nonnegative slice label."""

    labels = _as_1d_int(slice_labels, "slice_labels")
    counts = Counter(int(label) for label in labels if label >= 0)
    return dict(sorted(counts.items()))


def same_slice_unordered_pairs(
    slice_labels: ArrayLike,
    max_pairs_per_slice: int | None = None,
    seed: int | None = None,
) -> NDArray[np.int_]:
    """Return unordered event pairs whose endpoints share a nonnegative slice."""

    labels = _as_1d_int(slice_labels, "slice_labels")
    if max_pairs_per_slice is not None and int(max_pairs_per_slice) < 1:
        raise ValueError("max_pairs_per_slice must be positive when supplied")
    rng = np.random.default_rng(seed)
    rows: list[NDArray[np.int_]] = []
    for label in sorted(set(int(value) for value in labels if value >= 0)):
        indices = np.flatnonzero(labels == label)
        if indices.size < 2:
            continue
        left, right = np.triu_indices(indices.size, k=1)
        pairs = np.column_stack((indices[left], indices[right])).astype(int)
        if max_pairs_per_slice is not None and pairs.shape[0] > max_pairs_per_slice:
            selected = rng.choice(
                pairs.shape[0],
                size=int(max_pairs_per_slice),
                replace=False,
            )
            pairs = pairs[selected]
        rows.append(pairs)
    if not rows:
        return np.empty((0, 2), dtype=int)
    return np.vstack(rows).astype(int, copy=False)


def filter_pairs_by_same_slice(
    pairs: ArrayLike,
    slice_labels: ArrayLike,
) -> NDArray[np.int_]:
    """Return pairs whose endpoints share the same nonnegative slice label."""

    pair_array = np.asarray(pairs, dtype=int)
    labels = _as_1d_int(slice_labels, "slice_labels")
    if pair_array.ndim != 2 or pair_array.shape[1] != 2:
        raise ValueError("pairs must have shape (m, 2)")
    if np.any((pair_array < 0) | (pair_array >= labels.size)):
        raise IndexError("pairs contain an index out of range")
    left = labels[pair_array[:, 0]]
    right = labels[pair_array[:, 1]]
    keep = (left >= 0) & (left == right)
    return pair_array[keep].astype(int, copy=False)
