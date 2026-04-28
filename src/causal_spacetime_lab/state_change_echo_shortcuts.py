"""Controlled shortcut injection and background edge perturbation utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change import (
    StateChangeEvent,
    StateChangeNetwork,
    TriggerEdge,
)
from causal_spacetime_lab.state_change_echo_interference import (
    earliest_echo_delay_from_spectrum,
    return_delay_spectrum,
)
from causal_spacetime_lab.state_change_echo_motifs import EchoMotifRecord
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    topological_order_from_adjacency,
)

SHORTCUT_MODES = {
    "target_to_early_reference",
    "return_path_to_early_reference",
    "decoy_path_to_early_reference",
}


@dataclass(frozen=True)
class ShortcutInjectionSpec:
    """Configuration for controlled shortcut-return stress tests."""

    probability: float
    min_depth: int = 1
    max_depth: int | None = None
    mode: str = "target_to_early_reference"

    def __post_init__(self) -> None:
        if not 0.0 <= self.probability <= 1.0:
            raise ValueError("probability must be in [0, 1]")
        if self.min_depth < 1:
            raise ValueError("min_depth must be at least 1")
        if self.max_depth is not None and self.max_depth < self.min_depth:
            raise ValueError("max_depth must be None or at least min_depth")
        if self.mode not in SHORTCUT_MODES:
            allowed = ", ".join(sorted(SHORTCUT_MODES))
            raise ValueError(f"mode must be one of: {allowed}")


@dataclass(frozen=True)
class ShortcutInjectionRecord:
    """Record of one controlled shortcut-return injection."""

    motif_target_event_id: int
    emission_position: int
    planted_delay_rank: int
    shortcut_delay_rank: int
    shortcut_reference_position: int
    shortcut_source_event_id: int
    shortcut_target_reference_event_id: int
    mode: str


def _validate_reference_chain(
    reference_chain_event_ids: ArrayLike,
    event_count: int,
) -> NDArray[np.int_]:
    chain = np.asarray(reference_chain_event_ids, dtype=int)
    if chain.ndim != 1:
        raise ValueError("reference_chain_event_ids must be one-dimensional")
    if chain.size == 0:
        raise ValueError("reference_chain_event_ids must be nonempty")
    if np.min(chain) < 0 or np.max(chain) >= event_count:
        raise IndexError("reference chain contains event ids outside the network")
    return chain


def _copy_network(
    network: StateChangeNetwork,
) -> tuple[list[StateChangeEvent], list[TriggerEdge], dict[int, list[int]]]:
    return (
        list(network.events),
        list(network.trigger_edges),
        {
            int(system_id): list(event_ids)
            for system_id, event_ids in network.system_event_ids.items()
        },
    )


def _next_system_id(system_event_ids: dict[int, list[int]]) -> int:
    return max(system_event_ids, default=-1) + 1


def _add_decoy_event(
    events: list[StateChangeEvent],
    system_event_ids: dict[int, list[int]],
) -> StateChangeEvent:
    system_id = _next_system_id(system_event_ids)
    event = StateChangeEvent(
        event_id=len(events),
        system_id=system_id,
        local_index=0,
        previous_state=0,
        next_state=1,
    )
    events.append(event)
    system_event_ids.setdefault(system_id, []).append(event.event_id)
    return event


def _is_acyclic(
    events: list[StateChangeEvent],
    edges: list[TriggerEdge],
    system_event_ids: dict[int, list[int]],
) -> bool:
    trial = StateChangeNetwork(events, edges, system_event_ids)
    try:
        topological_order_from_adjacency(immediate_trigger_adjacency(trial))
    except ValueError:
        return False
    return True


def _candidate_shortcut_delays(
    motif: EchoMotifRecord,
    spec: ShortcutInjectionSpec,
) -> NDArray[np.int_]:
    delays = np.arange(1, motif.planted_delay_rank, dtype=int)
    if delays.size == 0:
        return delays
    depths = motif.planted_delay_rank - delays
    mask = depths >= spec.min_depth
    if spec.max_depth is not None:
        mask &= depths <= spec.max_depth
    return delays[mask]


def _source_event_for_mode(
    motif: EchoMotifRecord,
    mode: str,
    rng: np.random.Generator,
) -> int:
    if mode == "return_path_to_early_reference" and motif.return_event_ids.size:
        return int(rng.choice(motif.return_event_ids))
    return int(motif.target_event_id)


def inject_shortcut_returns(
    network: StateChangeNetwork,
    reference_chain_event_ids: ArrayLike,
    motifs: list[EchoMotifRecord],
    spec: ShortcutInjectionSpec,
    seed: int | None = None,
) -> tuple[StateChangeNetwork, list[ShortcutInjectionRecord]]:
    """Inject controlled shortcut returns while preserving acyclicity.

    This is a controlled stress test for the earliest-return rule, distinct
    from generic background trigger perturbations.
    """

    chain = _validate_reference_chain(reference_chain_event_ids, len(network.events))
    rng = np.random.default_rng(seed)
    events, edges, system_event_ids = _copy_network(network)
    records: list[ShortcutInjectionRecord] = []
    for motif in motifs:
        if rng.random() > spec.probability:
            continue
        candidate_delays = _candidate_shortcut_delays(motif, spec)
        if candidate_delays.size == 0:
            continue
        shortcut_delay = int(rng.choice(candidate_delays))
        shortcut_reference_position = motif.emission_position + shortcut_delay
        if shortcut_reference_position >= chain.size:
            continue
        target_reference_event = int(chain[shortcut_reference_position])
        trial_events = list(events)
        trial_edges = list(edges)
        trial_system_ids = {
            int(system_id): list(event_ids)
            for system_id, event_ids in system_event_ids.items()
        }
        source_event_id = _source_event_for_mode(motif, spec.mode, rng)
        if spec.mode == "decoy_path_to_early_reference":
            decoy = _add_decoy_event(trial_events, trial_system_ids)
            trial_edges.append(
                TriggerEdge(motif.target_event_id, decoy.event_id, "external_trigger")
            )
            trial_edges.append(
                TriggerEdge(decoy.event_id, target_reference_event, "external_trigger")
            )
            source_event_id = decoy.event_id
        else:
            trial_edges.append(
                TriggerEdge(source_event_id, target_reference_event, "external_trigger")
            )
        if not _is_acyclic(trial_events, trial_edges, trial_system_ids):
            continue
        events = trial_events
        edges = trial_edges
        system_event_ids = trial_system_ids
        records.append(
            ShortcutInjectionRecord(
                motif_target_event_id=int(motif.target_event_id),
                emission_position=int(motif.emission_position),
                planted_delay_rank=int(motif.planted_delay_rank),
                shortcut_delay_rank=shortcut_delay,
                shortcut_reference_position=int(shortcut_reference_position),
                shortcut_source_event_id=int(source_event_id),
                shortcut_target_reference_event_id=target_reference_event,
                mode=spec.mode,
            )
        )
    return StateChangeNetwork(events, edges, system_event_ids), records


def add_random_acyclic_background_edges(
    network: StateChangeNetwork,
    edge_probability: float,
    seed: int | None = None,
    max_edges: int | None = None,
) -> tuple[StateChangeNetwork, int]:
    """Add generic acyclic background trigger edges."""

    if not 0.0 <= edge_probability <= 1.0:
        raise ValueError("edge_probability must be in [0, 1]")
    if max_edges is not None and max_edges < 0:
        raise ValueError("max_edges must be nonnegative when supplied")
    rng = np.random.default_rng(seed)
    events, edges, system_event_ids = _copy_network(network)
    adjacency = immediate_trigger_adjacency(network)
    topological_order_from_adjacency(adjacency)
    existing = {
        (int(edge.source_event_id), int(edge.target_event_id))
        for edge in edges
        if edge.source_event_id >= 0
    }
    added = 0
    n_events = len(events)
    candidates = [
        (source, target)
        for source in range(n_events)
        for target in range(n_events)
        if source != target and (source, target) not in existing
    ]
    rng.shuffle(candidates)
    for pair in candidates:
        if rng.random() > edge_probability:
            continue
        trial_edges = list(edges)
        trial_edges.append(TriggerEdge(pair[0], pair[1], "external_trigger"))
        if not _is_acyclic(events, trial_edges, system_event_ids):
            continue
        edges = trial_edges
        existing.add(pair)
        added += 1
        if max_edges is not None and added >= max_edges:
            return StateChangeNetwork(events, edges, system_event_ids), added
    return StateChangeNetwork(events, edges, system_event_ids), added


def _classify_delay_change(
    before_delay: int | None,
    after_delay: int | None,
    planted_delay: int,
) -> str:
    if before_delay is None and after_delay is not None:
        return "newly_reachable"
    if before_delay is not None and after_delay is None:
        return "became_missing"
    if before_delay is None and after_delay is None:
        return "unchanged"
    assert before_delay is not None and after_delay is not None
    before_shortcut = before_delay < planted_delay
    after_shortcut = after_delay < planted_delay
    if not before_shortcut and after_shortcut:
        return "new_shortcut"
    if before_shortcut and after_shortcut and after_delay < before_delay:
        return "shortcut_deepened"
    if before_shortcut and not after_shortcut:
        return "shortcut_removed"
    return "unchanged"


def classify_added_edge_effects_on_motifs(
    before_order_matrix: ArrayLike,
    after_order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    motifs: list[EchoMotifRecord],
) -> list[dict[str, float | str]]:
    """Classify how added edges change motif echo-delay recovery."""

    before = np.asarray(before_order_matrix, dtype=bool)
    after = np.asarray(after_order_matrix, dtype=bool)
    chain = _validate_reference_chain(reference_chain_event_ids, before.shape[0])
    if before.shape != after.shape:
        raise ValueError("before and after order matrices must have equal shape")
    rows: list[dict[str, float | str]] = []
    for motif in motifs:
        before_delay = earliest_echo_delay_from_spectrum(
            return_delay_spectrum(
                before,
                chain,
                motif.target_event_id,
                motif.emission_position,
            )
        )
        after_delay = earliest_echo_delay_from_spectrum(
            return_delay_spectrum(
                after,
                chain,
                motif.target_event_id,
                motif.emission_position,
            )
        )
        classification = _classify_delay_change(
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
                "effect": classification,
            }
        )
    return rows
