from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change import (
    StateChangeEvent,
    StateChangeNetwork,
    TriggerEdge,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_reference_chains import (
    ReferenceChainCandidate,
)
from causal_spacetime_lab.state_change_reference_quality import (
    evaluate_reference_chain_candidate,
    rank_reference_chain_candidates,
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
    return network, transitive_closure_dag(immediate_trigger_adjacency(network))


def test_evaluate_reference_chain_candidate_returns_fields() -> None:
    network, closure = _network_and_order()
    candidate = ReferenceChainCandidate("manual", np.asarray([1, 3, 4]), "manual")

    report = evaluate_reference_chain_candidate(network, closure, candidate)

    assert report.is_valid_chain
    assert report.comparable_fraction == 1.0
    assert report.bracketed_fraction == 0.5
    assert report.local_system_purity == 1.0


def test_rank_reference_chain_candidates_returns_sorted_scores() -> None:
    network, closure = _network_and_order()
    reports = [
        evaluate_reference_chain_candidate(
            network,
            closure,
            ReferenceChainCandidate("short", np.asarray([0, 2]), "manual"),
        ),
        evaluate_reference_chain_candidate(
            network,
            closure,
            ReferenceChainCandidate("long", np.asarray([1, 3, 4]), "manual"),
        ),
    ]

    ranked = rank_reference_chain_candidates(reports)

    assert ranked[0]["score"] >= ranked[1]["score"]
    assert ranked[0]["rank"] == 1.0
