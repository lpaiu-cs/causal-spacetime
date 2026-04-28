from __future__ import annotations

import math

import numpy as np
import pytest

from causal_spacetime_lab.state_change import (
    StateChangeEvent,
    StateChangeNetwork,
    TriggerEdge,
)
from causal_spacetime_lab.state_change_observer_ambiguity import (
    candidate_overlap_fraction,
    chain_candidate_diversity,
    top_score_gap,
)
from causal_spacetime_lab.state_change_observer_quality import (
    evaluate_chain_candidate,
    rank_chain_candidates,
)
from causal_spacetime_lab.state_change_observers import ObserverChainCandidate
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)


def _network_and_order() -> tuple[StateChangeNetwork, np.ndarray]:
    network = StateChangeNetwork(
        events=[
            StateChangeEvent(0, 0, 0, 0, 1),
            StateChangeEvent(1, 1, 0, 0, 1),
            StateChangeEvent(2, 0, 1, 1, 2),
            StateChangeEvent(3, 1, 1, 1, 2),
            StateChangeEvent(4, 1, 2, 2, 3),
        ],
        trigger_edges=[
            TriggerEdge(-1, 0, "initial_seed"),
            TriggerEdge(-1, 1, "initial_seed"),
            TriggerEdge(0, 2, "local_successor"),
            TriggerEdge(1, 2, "external_trigger"),
            TriggerEdge(1, 3, "local_successor"),
            TriggerEdge(3, 4, "local_successor"),
            TriggerEdge(2, 4, "external_trigger"),
        ],
        system_event_ids={0: [0, 2], 1: [1, 3, 4]},
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    return network, closure


def test_evaluate_chain_candidate_returns_expected_fields() -> None:
    network, closure = _network_and_order()
    candidate = ObserverChainCandidate("manual", np.asarray([1, 3, 4]), "manual")

    report = evaluate_chain_candidate(network, closure, candidate)

    assert report.name == "manual"
    assert report.is_valid_chain
    assert report.chain_length == 3
    assert report.event_count == 5
    assert report.comparable_fraction == 1.0
    assert report.bracketed_fraction == 0.5
    assert report.local_system_purity == 1.0


def test_rank_chain_candidates_returns_sorted_scores() -> None:
    network, closure = _network_and_order()
    reports = [
        evaluate_chain_candidate(
            network,
            closure,
            ObserverChainCandidate("short", np.asarray([0, 2]), "manual"),
        ),
        evaluate_chain_candidate(
            network,
            closure,
            ObserverChainCandidate("long", np.asarray([1, 3, 4]), "manual"),
        ),
    ]

    ranked = rank_chain_candidates(reports)

    assert ranked[0]["score"] >= ranked[1]["score"]
    assert ranked[0]["rank"] == 1.0
    assert "regularity_score" in ranked[0]


def test_top_score_gap_handles_one_and_multiple_rows() -> None:
    assert math.isnan(top_score_gap([{"score": 0.5}]))
    assert top_score_gap([{"score": 0.7}, {"score": 0.4}]) == pytest.approx(0.3)


def test_candidate_overlap_fraction_deterministic_example() -> None:
    candidate_a = ObserverChainCandidate("a", np.asarray([0, 1, 2]), "manual")
    candidate_b = ObserverChainCandidate("b", np.asarray([1, 2, 3]), "manual")

    assert candidate_overlap_fraction(candidate_a, candidate_b) == 0.5


def test_chain_candidate_diversity_deterministic_example() -> None:
    candidates = [
        ObserverChainCandidate("a", np.asarray([0, 1]), "manual"),
        ObserverChainCandidate("b", np.asarray([1, 2]), "manual"),
        ObserverChainCandidate("c", np.asarray([3, 4]), "manual"),
    ]

    diversity = chain_candidate_diversity(candidates)

    assert diversity["candidate_count"] == 3.0
    assert diversity["min_pairwise_overlap"] == 0.0
    assert diversity["max_pairwise_overlap"] > 0.0
