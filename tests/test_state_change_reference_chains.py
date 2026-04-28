from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change import generate_state_change_network
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_reference_chains import (
    chain_bracketing_mask,
    chain_comparability_mask,
    greedy_reference_chain_candidate_from_order,
    is_chain,
    local_system_reference_chain_candidates,
    longest_reference_chain_candidate_from_order,
    random_reference_chain_candidate_from_order,
)


def _hand_closure() -> np.ndarray:
    adjacency = np.zeros((5, 5), dtype=bool)
    adjacency[0, 1] = True
    adjacency[1, 3] = True
    adjacency[0, 2] = True
    adjacency[2, 3] = True
    adjacency[3, 4] = True
    return transitive_closure_dag(adjacency)


def test_is_chain_true_on_reference_chain() -> None:
    assert is_chain(_hand_closure(), np.asarray([0, 3, 4]))


def test_is_chain_false_on_unordered_reference_indices() -> None:
    assert not is_chain(_hand_closure(), np.asarray([1, 2]))


def test_chain_comparability_mask_deterministic_example() -> None:
    mask = chain_comparability_mask(_hand_closure(), np.asarray([0, 3, 4]))

    assert np.array_equal(mask, np.asarray([True, True, True, True, True]))


def test_chain_bracketing_mask_deterministic_example() -> None:
    mask = chain_bracketing_mask(_hand_closure(), np.asarray([0, 3, 4]))

    assert np.array_equal(mask, np.asarray([False, True, True, False, False]))


def test_local_system_reference_chain_candidates_returns_expected() -> None:
    network = generate_state_change_network(3, 20, seed=201)

    candidates = local_system_reference_chain_candidates(network)

    assert len(candidates) == 3
    assert {candidate.source for candidate in candidates} == {"local_system"}


def test_reference_chain_generators_return_valid_chains() -> None:
    network = generate_state_change_network(4, 30, seed=202)
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))

    candidates = [
        greedy_reference_chain_candidate_from_order(closure),
        longest_reference_chain_candidate_from_order(closure),
        random_reference_chain_candidate_from_order(closure, seed=203),
    ]

    assert all(is_chain(closure, candidate.chain_event_ids) for candidate in candidates)
