"""Immediate-edge thinning utilities for state-change echo diagnostics."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from causal_spacetime_lab.state_change import StateChangeNetwork, TriggerEdge
from causal_spacetime_lab.state_change_echo import echo_delay_rank_for_emission
from causal_spacetime_lab.state_change_echo_motifs import EchoMotifRecord
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    topological_order_from_adjacency,
)


def edge_key(edge: TriggerEdge) -> tuple[int, int, str]:
    """Return a stable key for an immediate trigger edge."""

    return (
        int(edge.source_event_id),
        int(edge.target_event_id),
        str(edge.trigger_kind),
    )


def protected_source_target_pairs_for_reference_chain(
    reference_chain_event_ids: ArrayLike,
) -> set[tuple[int, int]]:
    """Return consecutive source-target pairs along a reference chain."""

    chain = np.asarray(reference_chain_event_ids, dtype=int)
    if chain.ndim != 1:
        raise ValueError("reference_chain_event_ids must be one-dimensional")
    return {
        (int(chain[index]), int(chain[index + 1]))
        for index in range(max(0, chain.size - 1))
    }


def protected_edge_keys_for_reference_chain(
    reference_chain_event_ids: ArrayLike,
) -> set[tuple[int, int, str]]:
    """Return local-successor edge keys along a reference chain."""

    return {
        (source, target, "local_successor")
        for source, target in protected_source_target_pairs_for_reference_chain(
            reference_chain_event_ids
        )
    }


def protected_edge_keys_for_motifs(
    motifs: list[EchoMotifRecord],
) -> set[tuple[int, int]]:
    """Return source-target pairs for planted motif paths."""

    pairs: set[tuple[int, int]] = set()
    for motif in motifs:
        path = [
            int(motif.reference_emission_event_id),
            *[int(value) for value in motif.outward_event_ids],
            int(motif.target_event_id),
            *[int(value) for value in motif.return_event_ids],
            int(motif.reference_return_event_id),
        ]
        pairs.update(
            (int(path[index]), int(path[index + 1]))
            for index in range(len(path) - 1)
        )
    return pairs


def thin_immediate_trigger_edges(
    network: StateChangeNetwork,
    removal_probability: float,
    *,
    protected_source_target_pairs: set[tuple[int, int]] | None = None,
    seed: int | None = None,
) -> tuple[StateChangeNetwork, int]:
    """Randomly remove immediate trigger edges while preserving acyclicity."""

    if not 0.0 <= removal_probability <= 1.0:
        raise ValueError("removal_probability must be in [0, 1]")
    rng = np.random.default_rng(seed)
    protected = protected_source_target_pairs or set()
    retained_edges: list[TriggerEdge] = []
    removed = 0
    for edge in network.trigger_edges:
        pair = (int(edge.source_event_id), int(edge.target_event_id))
        if edge.trigger_kind == "initial_seed" or pair in protected:
            retained_edges.append(edge)
            continue
        if rng.random() < removal_probability:
            removed += 1
        else:
            retained_edges.append(edge)
    thinned = StateChangeNetwork(
        list(network.events),
        retained_edges,
        {int(key): list(value) for key, value in network.system_event_ids.items()},
    )
    topological_order_from_adjacency(immediate_trigger_adjacency(thinned))
    return thinned, removed


def _classify_recovery_change(
    before_delay: int | None,
    after_delay: int | None,
    planted_delay: int,
) -> str:
    if before_delay is None and after_delay is None:
        return "unchanged"
    if before_delay is None and after_delay is not None:
        return "newly_reachable"
    if before_delay is not None and after_delay is None:
        return "became_missing"
    assert before_delay is not None and after_delay is not None
    before_shortcut = before_delay < planted_delay
    after_shortcut = after_delay < planted_delay
    before_exact = before_delay == planted_delay
    after_exact = after_delay == planted_delay
    if not before_shortcut and after_shortcut:
        return "shortcut_created"
    if before_shortcut and after_shortcut and after_delay < before_delay:
        return "shortcut_deepened"
    if before_shortcut and not after_shortcut:
        return "shortcut_removed"
    if after_delay > before_delay:
        return "delayed"
    if before_exact and not after_exact:
        return "exact_became_inexact"
    if not before_exact and after_exact:
        return "inexact_became_exact"
    return "unchanged"


def compare_recovery_before_after_edge_thinning(
    before_order_matrix: ArrayLike,
    after_order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    motifs: list[EchoMotifRecord],
) -> list[dict[str, float | str]]:
    """Compare motif echo recovery before and after edge thinning."""

    before = np.asarray(before_order_matrix, dtype=bool)
    after = np.asarray(after_order_matrix, dtype=bool)
    reference = np.asarray(reference_chain_event_ids, dtype=int)
    if before.shape != after.shape:
        raise ValueError("before and after order matrices must have equal shape")
    rows: list[dict[str, float | str]] = []
    for motif in motifs:
        before_delay = echo_delay_rank_for_emission(
            before,
            reference,
            motif.target_event_id,
            motif.emission_position,
        )
        after_delay = echo_delay_rank_for_emission(
            after,
            reference,
            motif.target_event_id,
            motif.emission_position,
        )
        effect = _classify_recovery_change(
            before_delay,
            after_delay,
            motif.planted_delay_rank,
        )
        rows.append(
            {
                "target_event_id": float(motif.target_event_id),
                "emission_position": float(motif.emission_position),
                "planted_delay_rank": float(motif.planted_delay_rank),
                "before_delay_rank": float(before_delay)
                if before_delay is not None
                else float("nan"),
                "after_delay_rank": float(after_delay)
                if after_delay is not None
                else float("nan"),
                "effect": effect,
            }
        )
    return rows
