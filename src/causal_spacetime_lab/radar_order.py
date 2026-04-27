"""Radar-return order utilities using causal order and observer tick order."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _as_square_bool_matrix(matrix: ArrayLike) -> NDArray[np.bool_]:
    array = np.asarray(matrix, dtype=bool)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("causal_matrix must be a square boolean matrix")
    return array


def _validate_indices(indices: ArrayLike, n: int, name: str) -> NDArray[np.int_]:
    array = np.asarray(indices, dtype=int)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if np.any((array < 0) | (array >= n)):
        raise IndexError(f"{name} contains an index out of range")
    return array


def latest_predecessor_tick_index(
    causal_matrix: ArrayLike,
    observer_indices: ArrayLike,
    target_index: int,
) -> int | None:
    """Return the latest observer tick position preceding the target."""

    causal_order = _as_square_bool_matrix(causal_matrix)
    observers = _validate_indices(observer_indices, causal_order.shape[0], "observers")
    target = int(target_index)
    if target < 0 or target >= causal_order.shape[0]:
        raise IndexError("target_index out of range")
    positions = np.flatnonzero(causal_order[observers, target])
    if positions.size == 0:
        return None
    return int(positions[-1])


def earliest_successor_tick_index(
    causal_matrix: ArrayLike,
    observer_indices: ArrayLike,
    target_index: int,
) -> int | None:
    """Return the earliest observer tick position following the target."""

    causal_order = _as_square_bool_matrix(causal_matrix)
    observers = _validate_indices(observer_indices, causal_order.shape[0], "observers")
    target = int(target_index)
    if target < 0 or target >= causal_order.shape[0]:
        raise IndexError("target_index out of range")
    positions = np.flatnonzero(causal_order[target, observers])
    if positions.size == 0:
        return None
    return int(positions[0])


def radar_tick_brackets_from_order(
    causal_matrix: ArrayLike,
    observer_indices: ArrayLike,
    target_indices: ArrayLike,
) -> tuple[NDArray[np.int_], NDArray[np.int_], NDArray[np.bool_]]:
    """Return predecessor/successor observer tick positions for targets.

    This function uses causal order and the order of ``observer_indices`` only.
    It does not use numeric observer clock labels.
    """

    causal_order = _as_square_bool_matrix(causal_matrix)
    observers = _validate_indices(observer_indices, causal_order.shape[0], "observers")
    targets = _validate_indices(target_indices, causal_order.shape[0], "targets")
    predecessor = np.full(targets.size, -1, dtype=int)
    successor = np.full(targets.size, -1, dtype=int)

    for row, target in enumerate(targets):
        previous = latest_predecessor_tick_index(causal_order, observers, int(target))
        following = earliest_successor_tick_index(causal_order, observers, int(target))
        if previous is not None:
            predecessor[row] = previous
        if following is not None:
            successor[row] = following

    accessible = (predecessor >= 0) & (successor >= 0)
    return predecessor, successor, accessible


def radar_return_order_from_successor_ticks(
    successor_tick_positions: ArrayLike,
    accessible_mask: ArrayLike,
) -> NDArray[np.bool_]:
    """Return same-protocol radar-return order from successor tick positions.

    For accessible targets, ``M[i, j]`` is true when target ``i`` returns to the
    observer at an earlier tick than target ``j``. This is a same-emission or
    common-protocol radar-return order. It is not an absolute distance.
    """

    successors = np.asarray(successor_tick_positions, dtype=int)
    accessible = np.asarray(accessible_mask, dtype=bool)
    if successors.ndim != 1 or accessible.ndim != 1:
        raise ValueError("successor_tick_positions and accessible_mask must be 1D")
    if successors.shape != accessible.shape:
        raise ValueError("successor_tick_positions and accessible_mask shapes differ")
    order = successors[:, None] < successors[None, :]
    return order & accessible[:, None] & accessible[None, :]
