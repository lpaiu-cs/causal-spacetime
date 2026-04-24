from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.lattice import (
    cumulative_counts_by_time,
    edge_displacements,
    lattice_cumulative_counts,
    lattice_shell_counts,
    regular_lattice_causal_graph_1p1,
)


def test_regular_lattice_graph_has_expected_small_cone() -> None:
    graph = regular_lattice_causal_graph_1p1(2)

    assert graph.events.tolist() == [
        [0, 0],
        [1, -1],
        [1, 1],
        [2, -2],
        [2, 0],
        [2, 2],
    ]
    assert graph.edges.tolist() == [
        [0, 1],
        [0, 2],
        [1, 3],
        [1, 4],
        [2, 4],
        [2, 5],
    ]


def test_lattice_edge_displacements_are_finite_speed_diagonals() -> None:
    graph = regular_lattice_causal_graph_1p1(4)
    displacements = edge_displacements(graph)

    assert set(map(tuple, displacements.tolist())) == {(1, -1), (1, 1)}


def test_lattice_count_growth() -> None:
    assert lattice_shell_counts(3).tolist() == [1, 2, 3, 4]
    assert lattice_cumulative_counts(3).tolist() == [1, 3, 6, 10]


def test_cumulative_counts_by_time() -> None:
    events = np.array(
        [
            [0.1, 0.0],
            [1.2, 0.5],
            [1.9, -0.5],
            [3.0, 0.0],
        ]
    )
    counts = cumulative_counts_by_time(events, [0.0, 1.5, 3.0])

    assert counts.tolist() == [0, 2, 4]


def test_lattice_rejects_negative_t_steps() -> None:
    with pytest.raises(ValueError, match="t_steps must be non-negative"):
        regular_lattice_causal_graph_1p1(-1)

