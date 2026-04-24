"""Finite-speed regular lattice causal graph utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class LatticeCausalGraph:
    """A finite 1+1D lattice causal graph."""

    events: NDArray[np.int64]
    edges: NDArray[np.int64]


def _validate_t_steps(t_steps: int) -> int:
    steps = int(t_steps)
    if steps < 0:
        raise ValueError("t_steps must be non-negative")
    return steps


def regular_lattice_causal_graph_1p1(t_steps: int) -> LatticeCausalGraph:
    """Return the forward causal cone of a finite-speed 1+1D lattice graph.

    Nodes are reachable integer events ``(t, x)`` from ``(0, 0)`` through
    ``t_steps``. Each node at time ``t < t_steps`` has directed edges to
    ``(t + 1, x - 1)`` and ``(t + 1, x + 1)``.
    """

    steps = _validate_t_steps(t_steps)

    events: list[tuple[int, int]] = []
    index_by_event: dict[tuple[int, int], int] = {}
    for t in range(steps + 1):
        for x in range(-t, t + 1, 2):
            index_by_event[(t, x)] = len(events)
            events.append((t, x))

    edges: list[tuple[int, int]] = []
    for t, x in events:
        if t == steps:
            continue
        source = index_by_event[(t, x)]
        edges.append((source, index_by_event[(t + 1, x - 1)]))
        edges.append((source, index_by_event[(t + 1, x + 1)]))

    return LatticeCausalGraph(
        events=np.asarray(events, dtype=np.int64),
        edges=np.asarray(edges, dtype=np.int64).reshape((-1, 2)),
    )


def lattice_shell_counts(t_steps: int) -> NDArray[np.int64]:
    """Return reachable-node counts at each integer lattice time."""

    steps = _validate_t_steps(t_steps)
    return np.arange(1, steps + 2, dtype=np.int64)


def lattice_cumulative_counts(t_steps: int) -> NDArray[np.int64]:
    """Return cumulative reachable-node counts through each integer time."""

    return np.cumsum(lattice_shell_counts(t_steps))


def cumulative_counts_by_time(
    events: ArrayLike,
    times: ArrayLike,
) -> NDArray[np.int64]:
    """Count events with coordinate time less than or equal to each value."""

    event_array = np.asarray(events, dtype=float)
    if event_array.ndim != 2 or event_array.shape[1] != 2:
        raise ValueError("events must have shape (n, 2), with columns (t, x)")

    time_values = np.asarray(times, dtype=float)
    return np.asarray(
        [np.count_nonzero(event_array[:, 0] <= time) for time in time_values],
        dtype=np.int64,
    )


def edge_displacements(
    graph: LatticeCausalGraph,
) -> NDArray[np.int64]:
    """Return ``(dt, dx)`` for each directed lattice edge."""

    return graph.events[graph.edges[:, 1]] - graph.events[graph.edges[:, 0]]

