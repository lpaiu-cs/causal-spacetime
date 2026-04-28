"""Controlled state-change scenarios for echo-response motif diagnostics."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.state_change import (
    StateChangeEvent,
    StateChangeNetwork,
    TriggerEdge,
    generate_state_change_network,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_reference_chains import (
    ReferenceChainCandidate,
    local_system_reference_chain_candidates,
    longest_reference_chain_candidate_from_order,
)
from causal_spacetime_lab.state_change_reference_quality import (
    evaluate_reference_chain_candidate,
    rank_reference_chain_candidates,
)


def generate_reference_backbone_network(
    reference_length: int,
    *,
    system_id: int = 0,
) -> tuple[StateChangeNetwork, NDArray[np.int_]]:
    """Build a minimal finite network containing one reference backbone."""

    if reference_length < 1:
        raise ValueError("reference_length must be positive")
    events = [
        StateChangeEvent(
            event_id=index,
            system_id=int(system_id),
            local_index=index,
            previous_state=index,
            next_state=index + 1,
        )
        for index in range(reference_length)
    ]
    edges = [
        TriggerEdge(-1, 0, "initial_seed"),
    ]
    edges.extend(
        TriggerEdge(index - 1, index, "local_successor")
        for index in range(1, reference_length)
    )
    reference_chain = np.arange(reference_length, dtype=int)
    return (
        StateChangeNetwork(
            events=events,
            trigger_edges=edges,
            system_event_ids={int(system_id): reference_chain.tolist()},
        ),
        reference_chain,
    )


def _longest_local_reference(network: StateChangeNetwork) -> NDArray[np.int_]:
    candidates = local_system_reference_chain_candidates(network, min_length=1)
    if not candidates:
        return np.empty(0, dtype=int)
    candidates.sort(
        key=lambda candidate: (
            candidate.chain_event_ids.size,
            -candidate.chain_event_ids[0],
        ),
        reverse=True,
    )
    return np.asarray(candidates[0].chain_event_ids, dtype=int)


def _highest_utility_reference(
    network: StateChangeNetwork,
    closure: NDArray[np.bool_],
) -> NDArray[np.int_]:
    candidates: list[ReferenceChainCandidate] = local_system_reference_chain_candidates(
        network,
        min_length=1,
    )
    candidates.append(longest_reference_chain_candidate_from_order(closure))
    reports = [
        evaluate_reference_chain_candidate(network, closure, candidate)
        for candidate in candidates
    ]
    ranked = rank_reference_chain_candidates(reports)
    by_name = {candidate.name: candidate for candidate in candidates}
    return np.asarray(by_name[str(ranked[0]["name"])].chain_event_ids, dtype=int)


def generate_background_state_change_network_with_reference(
    num_systems: int,
    max_events: int,
    trigger_probability: float,
    *,
    reference_source: str = "longest_local",
    seed: int | None = None,
) -> tuple[StateChangeNetwork, NDArray[np.int_]]:
    """Generate a background network and select a reference protocol backbone.

    The selected reference chain is a protocol backbone for diagnostics, not a
    physical-observer claim.
    """

    network = generate_state_change_network(
        num_systems,
        max_events,
        trigger_probability=trigger_probability,
        seed=seed,
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    if reference_source == "longest_local":
        reference = _longest_local_reference(network)
    elif reference_source == "highest_utility":
        reference = _highest_utility_reference(network, closure)
    elif reference_source == "longest_order":
        reference = longest_reference_chain_candidate_from_order(
            closure,
        ).chain_event_ids
    else:
        raise ValueError(
            "reference_source must be one of: longest_local, highest_utility, "
            "longest_order"
        )
    return network, np.asarray(reference, dtype=int)
