from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_brackets import (
    assign_reference_rank_slices,
    bracket_width_rank_from_reference_brackets,
    earliest_successor_reference_position,
    latest_predecessor_reference_position,
    radar_time_rank_from_reference_brackets,
    reference_tick_brackets_from_order,
)
from causal_spacetime_lab.state_change_order import transitive_closure_dag


def _hand_closure() -> np.ndarray:
    adjacency = np.zeros((5, 5), dtype=bool)
    adjacency[0, 1] = True
    adjacency[1, 3] = True
    adjacency[0, 2] = True
    adjacency[2, 3] = True
    adjacency[3, 4] = True
    return transitive_closure_dag(adjacency)


def test_reference_position_helpers_deterministic_example() -> None:
    closure = _hand_closure()
    reference = np.asarray([0, 3, 4])

    assert latest_predecessor_reference_position(closure, reference, 2) == 0
    assert earliest_successor_reference_position(closure, reference, 2) == 1


def test_reference_tick_brackets_from_order_deterministic_example() -> None:
    closure = _hand_closure()
    reference = np.asarray([0, 3, 4])

    predecessors, successors, accessible = reference_tick_brackets_from_order(
        closure,
        reference,
    )

    assert np.array_equal(predecessors, np.asarray([-1, 0, 0, -1, -1]))
    assert np.array_equal(successors, np.asarray([-1, 1, 1, -1, -1]))
    assert np.array_equal(accessible, np.asarray([False, True, True, False, False]))


def test_reference_tick_brackets_handles_missing_sides() -> None:
    closure = _hand_closure()
    reference = np.asarray([3, 4])

    predecessors, successors, accessible = reference_tick_brackets_from_order(
        closure,
        reference,
        target_indices=np.asarray([0, 1, 2]),
    )

    assert np.array_equal(predecessors, np.asarray([-1, -1, -1]))
    assert np.array_equal(successors, np.asarray([0, 0, 0]))
    assert not np.any(accessible)


def test_rank_helpers_from_reference_brackets() -> None:
    predecessors = np.asarray([-1, 0, 0, -1])
    successors = np.asarray([-1, 2, 3, -1])
    accessible = np.asarray([False, True, True, False])

    time_ranks = radar_time_rank_from_reference_brackets(
        predecessors,
        successors,
        accessible,
    )
    width_ranks = bracket_width_rank_from_reference_brackets(
        predecessors,
        successors,
        accessible,
    )
    slices = assign_reference_rank_slices(time_ranks, accessible, bin_width=2)

    assert np.array_equal(time_ranks, np.asarray([-1, 2, 3, -1]))
    assert np.array_equal(width_ranks, np.asarray([-1, 2, 3, -1]))
    assert np.array_equal(slices, np.asarray([-1, 1, 1, -1]))
