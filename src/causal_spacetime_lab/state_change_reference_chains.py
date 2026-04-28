"""Reference-chain candidate utilities for state-change trigger networks."""

from __future__ import annotations

from causal_spacetime_lab.state_change_observers import (
    ReferenceChainCandidate,
    chain_bracketing_mask,
    chain_comparability_mask,
    earliest_successor_chain_position,
    greedy_reference_chain_candidate_from_order,
    is_chain,
    latest_predecessor_chain_position,
    local_system_reference_chain_candidates,
    longest_reference_chain_candidate_from_order,
    random_reference_chain_candidate_from_order,
)

__all__ = [
    "ReferenceChainCandidate",
    "chain_bracketing_mask",
    "chain_comparability_mask",
    "earliest_successor_chain_position",
    "greedy_reference_chain_candidate_from_order",
    "is_chain",
    "latest_predecessor_chain_position",
    "local_system_reference_chain_candidates",
    "longest_reference_chain_candidate_from_order",
    "random_reference_chain_candidate_from_order",
]
