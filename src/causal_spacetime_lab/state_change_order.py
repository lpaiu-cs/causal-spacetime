"""Order utilities for finite state-change trigger networks."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.state_change import StateChangeNetwork


def immediate_trigger_adjacency(network: StateChangeNetwork) -> NDArray[np.bool_]:
    """Return immediate trigger adjacency for event-to-event trigger edges."""

    n_events = len(network.events)
    adjacency = np.zeros((n_events, n_events), dtype=bool)
    for edge in network.trigger_edges:
        if edge.source_event_id < 0:
            continue
        if edge.source_event_id >= n_events or edge.target_event_id >= n_events:
            raise IndexError("trigger edge references an event outside the network")
        adjacency[edge.source_event_id, edge.target_event_id] = True
    return adjacency


def topological_order_from_adjacency(adjacency: NDArray[np.bool_]) -> NDArray[np.int_]:
    """Return one topological order, raising if a cycle is detected."""

    matrix = np.asarray(adjacency, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("adjacency must be a square matrix")
    n_items = matrix.shape[0]
    in_degree = np.sum(matrix, axis=0).astype(int)
    ready = [int(index) for index in np.flatnonzero(in_degree == 0)]
    order: list[int] = []
    cursor = 0
    while cursor < len(ready):
        node = ready[cursor]
        cursor += 1
        order.append(node)
        for successor in np.flatnonzero(matrix[node]):
            in_degree[int(successor)] -= 1
            if in_degree[int(successor)] == 0:
                ready.append(int(successor))
    if len(order) != n_items:
        raise ValueError("adjacency contains a directed cycle")
    return np.asarray(order, dtype=int)


def transitive_closure_dag(adjacency: NDArray[np.bool_]) -> NDArray[np.bool_]:
    """Return the transitive closure of a finite DAG adjacency matrix."""

    matrix = np.asarray(adjacency, dtype=bool)
    order = topological_order_from_adjacency(matrix)
    closure = matrix.copy()
    for node in reversed(order):
        for successor in np.flatnonzero(matrix[node]):
            closure[int(node)] |= closure[int(successor)]
    np.fill_diagonal(closure, False)
    return closure


def is_irreflexive(order_matrix: NDArray[np.bool_]) -> bool:
    """Return whether an order matrix has no self-relations."""

    matrix = np.asarray(order_matrix, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("order_matrix must be square")
    return bool(not np.any(np.diag(matrix)))


def is_transitive(order_matrix: NDArray[np.bool_]) -> bool:
    """Return whether all composed order relations are already present."""

    matrix = np.asarray(order_matrix, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("order_matrix must be square")
    composed = (matrix.astype(int) @ matrix.astype(int)) > 0
    return bool(np.all(~composed | matrix))


def causal_interval_indices(
    order_matrix: NDArray[np.bool_],
    source: int,
    target: int,
) -> NDArray[np.int_]:
    """Return indices ``k`` such that ``source ≺_T k ≺_T target``."""

    matrix = np.asarray(order_matrix, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("order_matrix must be square")
    n_items = matrix.shape[0]
    if source < 0 or target < 0 or source >= n_items or target >= n_items:
        raise IndexError("source and target must be valid event indices")
    return np.flatnonzero(matrix[int(source)] & matrix[:, int(target)]).astype(int)


def local_finiteness_report(order_matrix: NDArray[np.bool_]) -> dict[str, float]:
    """Report finite interval sizes for a finite order matrix."""

    matrix = np.asarray(order_matrix, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("order_matrix must be square")
    interval_sizes: list[int] = []
    for source, target in zip(*np.nonzero(matrix), strict=True):
        interval_sizes.append(
            int(causal_interval_indices(matrix, int(source), int(target)).size)
        )
    comparable_count = int(np.count_nonzero(matrix))
    return {
        "comparable_pair_count": float(comparable_count),
        "max_interval_size": float(max(interval_sizes, default=0)),
        "mean_interval_size": float(np.mean(interval_sizes)) if interval_sizes else 0.0,
        "finite_interval_fraction": 1.0,
    }
