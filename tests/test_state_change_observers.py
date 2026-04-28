from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change import generate_state_change_network
from causal_spacetime_lab.state_change_observers import (
    chain_bracketing_mask,
    chain_comparability_mask,
    earliest_successor_chain_position,
    greedy_chain_candidate_from_order,
    is_chain,
    latest_predecessor_chain_position,
    local_system_chain_candidates,
    longest_chain_candidate_from_order,
    random_chain_candidate_from_order,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)


def _hand_closure() -> np.ndarray:
    adjacency = np.zeros((5, 5), dtype=bool)
    adjacency[0, 1] = True
    adjacency[1, 3] = True
    adjacency[0, 2] = True
    adjacency[2, 3] = True
    adjacency[3, 4] = True
    return transitive_closure_dag(adjacency)


def test_is_chain_true_on_hand_coded_chain() -> None:
    closure = _hand_closure()

    assert is_chain(closure, np.asarray([0, 3, 4]))


def test_is_chain_false_on_unordered_indices() -> None:
    closure = _hand_closure()

    assert not is_chain(closure, np.asarray([1, 2]))


def test_chain_comparability_mask_deterministic_example() -> None:
    closure = _hand_closure()

    mask = chain_comparability_mask(closure, np.asarray([0, 3, 4]))

    assert np.array_equal(mask, np.asarray([True, True, True, True, True]))


def test_chain_bracketing_mask_deterministic_example() -> None:
    closure = _hand_closure()

    mask = chain_bracketing_mask(closure, np.asarray([0, 3, 4]))

    assert np.array_equal(mask, np.asarray([False, True, True, False, False]))


def test_chain_bracket_positions_deterministic_example() -> None:
    closure = _hand_closure()
    chain = np.asarray([0, 3, 4])

    assert latest_predecessor_chain_position(closure, chain, 2) == 0
    assert earliest_successor_chain_position(closure, chain, 2) == 1


def test_local_system_chain_candidates_returns_expected_candidates() -> None:
    network = generate_state_change_network(3, 20, seed=101)

    candidates = local_system_chain_candidates(network)

    assert len(candidates) == 3
    assert {candidate.source for candidate in candidates} == {"local_system"}
    assert all(candidate.chain_event_ids.size >= 2 for candidate in candidates)


def test_order_candidate_generators_return_valid_chains() -> None:
    network = generate_state_change_network(4, 30, seed=102)
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))

    candidates = [
        greedy_chain_candidate_from_order(closure),
        longest_chain_candidate_from_order(closure),
        random_chain_candidate_from_order(closure, seed=103),
    ]

    assert all(is_chain(closure, candidate.chain_event_ids) for candidate in candidates)
    assert {candidate.source for candidate in candidates} == {
        "greedy_order",
        "longest_chain",
        "random_order",
    }
