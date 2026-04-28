"""Validation summaries for finite state-change trigger networks."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.state_change import StateChangeNetwork


def local_chain_lengths(network: StateChangeNetwork) -> dict[int, int]:
    """Return the number of state-changing events in each local system chain."""

    return {
        int(system_id): len(event_ids)
        for system_id, event_ids in sorted(network.system_event_ids.items())
    }


def state_change_network_summary(network: StateChangeNetwork) -> dict[str, float]:
    """Return basic finite state-change network counts."""

    chain_lengths = list(local_chain_lengths(network).values())
    local_edges = sum(
        edge.trigger_kind == "local_successor" for edge in network.trigger_edges
    )
    external_edges = sum(
        edge.trigger_kind == "external_trigger" for edge in network.trigger_edges
    )
    return {
        "event_count": float(len(network.events)),
        "system_count": float(len(network.system_event_ids)),
        "trigger_edge_count": float(len(network.trigger_edges)),
        "local_successor_edge_count": float(local_edges),
        "external_trigger_edge_count": float(external_edges),
        "mean_events_per_system": float(np.mean(chain_lengths))
        if chain_lengths
        else 0.0,
        "max_events_per_system": float(max(chain_lengths, default=0)),
    }


def trigger_graph_summary(
    adjacency: NDArray[np.bool_],
    closure: NDArray[np.bool_],
) -> dict[str, float]:
    """Return immediate and transitive trigger graph summaries."""

    immediate = np.asarray(adjacency, dtype=bool)
    order = np.asarray(closure, dtype=bool)
    if immediate.shape != order.shape or immediate.ndim != 2:
        raise ValueError("adjacency and closure must have the same matrix shape")
    n_events = immediate.shape[0]
    out_degree = np.sum(immediate, axis=1)
    in_degree = np.sum(immediate, axis=0)
    denominator = n_events * (n_events - 1)
    return {
        "immediate_edge_count": float(np.count_nonzero(immediate)),
        "causal_relation_count": float(np.count_nonzero(order)),
        "relation_density": float(np.count_nonzero(order) / denominator)
        if denominator
        else 0.0,
        "max_out_degree": float(np.max(out_degree)) if n_events else 0.0,
        "max_in_degree": float(np.max(in_degree)) if n_events else 0.0,
        "mean_out_degree": float(np.mean(out_degree)) if n_events else 0.0,
        "mean_in_degree": float(np.mean(in_degree)) if n_events else 0.0,
    }


def branching_statistics(network: StateChangeNetwork) -> dict[str, float]:
    """Report outgoing external trigger counts per event."""

    counts = np.zeros(len(network.events), dtype=int)
    for edge in network.trigger_edges:
        if edge.trigger_kind == "external_trigger" and edge.source_event_id >= 0:
            counts[edge.source_event_id] += 1
    return {
        "external_branching_event_count": float(np.count_nonzero(counts)),
        "max_external_triggers_from_event": float(np.max(counts))
        if counts.size
        else 0.0,
        "mean_external_triggers_from_event": float(np.mean(counts))
        if counts.size
        else 0.0,
    }
