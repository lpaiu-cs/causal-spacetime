"""Longest-chain calculations on finite causal orders."""

from __future__ import annotations

from collections import deque

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _as_square_bool_matrix(matrix: ArrayLike) -> NDArray[np.bool_]:
    array = np.asarray(matrix, dtype=bool)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("C must be a square boolean matrix")
    return array


def _topological_order(C: NDArray[np.bool_]) -> NDArray[np.int_]:
    n = C.shape[0]
    indegree = C.sum(axis=0).astype(int)
    queue: deque[int] = deque(np.flatnonzero(indegree == 0).tolist())
    order: list[int] = []

    while queue:
        node = queue.popleft()
        order.append(node)
        for successor in np.flatnonzero(C[node]):
            indegree[successor] -= 1
            if indegree[successor] == 0:
                queue.append(int(successor))

    if len(order) != n:
        raise ValueError("C must represent an acyclic causal relation")
    return np.asarray(order, dtype=int)


def _order_from_times(
    event_times: ArrayLike | None,
    C: NDArray[np.bool_],
) -> NDArray[np.int_]:
    if event_times is None:
        return _topological_order(C)

    times = np.asarray(event_times, dtype=float)
    if times.shape != (C.shape[0],):
        raise ValueError("event_times must have shape (n,)")
    return np.argsort(times, kind="stable")


def _validate_index(index: int, n: int, name: str) -> int:
    integer_index = int(index)
    if integer_index < 0 or integer_index >= n:
        raise IndexError(f"{name} index out of range")
    return integer_index


def longest_chain_length(
    C: ArrayLike,
    start: int | None = None,
    end: int | None = None,
    event_times: ArrayLike | None = None,
) -> int:
    """Compute the longest chain length in a finite DAG.

    Chain length is the number of events in the chain. If ``start`` and ``end``
    are supplied, the returned value includes both endpoints when a chain exists.
    If no chain connects them, ``0`` is returned.
    """

    causal_matrix = _as_square_bool_matrix(C)
    n = causal_matrix.shape[0]
    if n == 0:
        return 0

    if (start is None) != (end is None):
        raise ValueError("start and end must either both be supplied or both omitted")

    order = _order_from_times(event_times, causal_matrix)

    if start is None and end is None:
        best_ending_at = np.ones(n, dtype=int)
        for node in order:
            next_length = best_ending_at[node] + 1
            for successor in np.flatnonzero(causal_matrix[node]):
                if next_length > best_ending_at[successor]:
                    best_ending_at[successor] = next_length
        return int(best_ending_at.max())

    start_index = _validate_index(start, n, "start")
    end_index = _validate_index(end, n, "end")
    if start_index == end_index:
        return 1

    best_ending_at = np.zeros(n, dtype=int)
    best_ending_at[start_index] = 1

    for node in order:
        if best_ending_at[node] == 0:
            continue
        next_length = best_ending_at[node] + 1
        for successor in np.flatnonzero(causal_matrix[node]):
            if next_length > best_ending_at[successor]:
                best_ending_at[successor] = next_length

    return int(best_ending_at[end_index])


def longest_chain_indices(
    C: ArrayLike,
    start: int,
    end: int,
    event_times: ArrayLike | None = None,
) -> NDArray[np.int_]:
    """Return one longest causal chain from ``start`` to ``end``.

    The returned chain includes both endpoints. If no chain exists, an empty
    array is returned.
    """

    causal_matrix = _as_square_bool_matrix(C)
    n = causal_matrix.shape[0]
    if n == 0:
        return np.empty(0, dtype=int)

    start_index = _validate_index(start, n, "start")
    end_index = _validate_index(end, n, "end")
    if start_index == end_index:
        return np.asarray([start_index], dtype=int)

    order = _order_from_times(event_times, causal_matrix)
    best_ending_at = np.zeros(n, dtype=int)
    predecessor = np.full(n, -1, dtype=int)
    best_ending_at[start_index] = 1

    for node in order:
        if best_ending_at[node] == 0:
            continue
        next_length = best_ending_at[node] + 1
        for successor in np.flatnonzero(causal_matrix[node]):
            if next_length > best_ending_at[successor]:
                best_ending_at[successor] = next_length
                predecessor[successor] = node

    if best_ending_at[end_index] == 0:
        return np.empty(0, dtype=int)

    chain: list[int] = []
    node = end_index
    while node != -1:
        chain.append(int(node))
        if node == start_index:
            break
        node = int(predecessor[node])

    if chain[-1] != start_index:
        return np.empty(0, dtype=int)
    return np.asarray(chain[::-1], dtype=int)
