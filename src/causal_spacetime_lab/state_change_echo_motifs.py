"""Controlled echo-response motifs for state-change trigger networks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change import (
    StateChangeEvent,
    StateChangeNetwork,
    TriggerEdge,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    topological_order_from_adjacency,
)


@dataclass(frozen=True)
class EchoMotifSpec:
    """Specification for a controlled echo-response motif."""

    emission_position: int
    planted_delay_rank: int
    target_system_id: int | None = None
    outward_steps: int = 1
    return_steps: int = 1
    motif_kind: str = "single_target"

    def __post_init__(self) -> None:
        if self.emission_position < 0:
            raise ValueError("emission_position must be nonnegative")
        if self.planted_delay_rank < 1:
            raise ValueError("planted_delay_rank must be at least 1")
        if self.target_system_id is not None and self.target_system_id < 0:
            raise ValueError("target_system_id must be nonnegative when supplied")
        if self.outward_steps < 0:
            raise ValueError("outward_steps must be nonnegative")
        if self.return_steps < 0:
            raise ValueError("return_steps must be nonnegative")
        if not self.motif_kind:
            raise ValueError("motif_kind must be nonempty")


@dataclass(frozen=True)
class EchoMotifRecord:
    """Record of a planted causal trigger motif.

    A planted echo motif records controlled causal trigger structure. It is not
    a metric-distance object.
    """

    target_event_id: int
    emission_position: int
    planted_return_position: int
    planted_delay_rank: int
    outward_event_ids: NDArray[np.int_]
    return_event_ids: NDArray[np.int_]
    reference_emission_event_id: int
    reference_return_event_id: int
    motif_kind: str


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


def _add_storage_event(
    events: list[StateChangeEvent],
    system_event_ids: dict[int, list[int]],
    system_id: int,
) -> StateChangeEvent:
    local_index = len(system_event_ids.get(int(system_id), []))
    event = StateChangeEvent(
        event_id=len(events),
        system_id=int(system_id),
        local_index=int(local_index),
        previous_state=int(local_index),
        next_state=int(local_index + 1),
    )
    events.append(event)
    system_event_ids.setdefault(int(system_id), []).append(event.event_id)
    return event


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


def insert_echo_motif(
    network: StateChangeNetwork,
    reference_chain_event_ids: ArrayLike,
    spec: EchoMotifSpec,
    *,
    system_id_start: int | None = None,
) -> tuple[StateChangeNetwork, EchoMotifRecord]:
    """Insert one controlled echo-response motif into a finite trigger network.

    The planted delay rank is an order-level validation label. Newly inserted
    event ids are finite-storage identifiers, not physical time labels.
    """

    reference_chain = _validate_reference_chain(
        reference_chain_event_ids,
        len(network.events),
    )
    planted_return_position = spec.emission_position + spec.planted_delay_rank
    if planted_return_position >= reference_chain.size:
        raise ValueError("planted return position is outside the reference chain")

    events, trigger_edges, system_event_ids = _copy_network(network)
    system_cursor = (
        int(system_id_start)
        if system_id_start is not None
        else _next_system_id(system_event_ids)
    )

    def allocate_system_id(preferred: int | None = None) -> int:
        nonlocal system_cursor
        if preferred is not None:
            return int(preferred)
        system_id = system_cursor
        system_cursor += 1
        return system_id

    outward_event_ids: list[int] = []
    return_event_ids: list[int] = []
    previous_event_id = int(reference_chain[spec.emission_position])

    for _ in range(spec.outward_steps):
        event = _add_storage_event(
            events,
            system_event_ids,
            allocate_system_id(),
        )
        outward_event_ids.append(event.event_id)
        trigger_edges.append(
            TriggerEdge(previous_event_id, event.event_id, "external_trigger")
        )
        previous_event_id = event.event_id

    target_event = _add_storage_event(
        events,
        system_event_ids,
        allocate_system_id(spec.target_system_id),
    )
    trigger_edges.append(
        TriggerEdge(previous_event_id, target_event.event_id, "external_trigger")
    )
    previous_event_id = target_event.event_id

    for _ in range(spec.return_steps):
        event = _add_storage_event(
            events,
            system_event_ids,
            allocate_system_id(),
        )
        return_event_ids.append(event.event_id)
        trigger_edges.append(
            TriggerEdge(previous_event_id, event.event_id, "external_trigger")
        )
        previous_event_id = event.event_id

    reference_return_event_id = int(reference_chain[planted_return_position])
    trigger_edges.append(
        TriggerEdge(previous_event_id, reference_return_event_id, "external_trigger")
    )

    updated = StateChangeNetwork(events, trigger_edges, system_event_ids)
    try:
        topological_order_from_adjacency(immediate_trigger_adjacency(updated))
    except ValueError as exc:
        raise ValueError("inserted echo motif creates a cycle") from exc

    record = EchoMotifRecord(
        target_event_id=target_event.event_id,
        emission_position=int(spec.emission_position),
        planted_return_position=int(planted_return_position),
        planted_delay_rank=int(spec.planted_delay_rank),
        outward_event_ids=np.asarray(outward_event_ids, dtype=int),
        return_event_ids=np.asarray(return_event_ids, dtype=int),
        reference_emission_event_id=int(reference_chain[spec.emission_position]),
        reference_return_event_id=reference_return_event_id,
        motif_kind=spec.motif_kind,
    )
    return updated, record


def insert_multiple_echo_motifs(
    network: StateChangeNetwork,
    reference_chain_event_ids: ArrayLike,
    specs: list[EchoMotifSpec],
) -> tuple[StateChangeNetwork, list[EchoMotifRecord]]:
    """Insert multiple controlled echo-response motifs sequentially."""

    current = network
    records: list[EchoMotifRecord] = []
    for spec in specs:
        current, record = insert_echo_motif(
            current,
            reference_chain_event_ids,
            spec,
        )
        records.append(record)
    return current, records


def build_echo_motif_specs_for_reference_chain(
    reference_chain_event_ids: ArrayLike,
    emission_positions: ArrayLike,
    delay_ranks: ArrayLike,
    targets_per_emission: int = 1,
    outward_steps: int = 1,
    return_steps: int = 1,
    seed: int | None = None,
) -> list[EchoMotifSpec]:
    """Build valid motif specs for a reference chain and delay-rank set."""

    chain = np.asarray(reference_chain_event_ids, dtype=int)
    emissions = np.asarray(emission_positions, dtype=int)
    delays = np.asarray(delay_ranks, dtype=int)
    if chain.ndim != 1:
        raise ValueError("reference_chain_event_ids must be one-dimensional")
    if emissions.ndim != 1 or delays.ndim != 1:
        raise ValueError("emission_positions and delay_ranks must be vectors")
    if targets_per_emission < 1:
        raise ValueError("targets_per_emission must be positive")
    if outward_steps < 0 or return_steps < 0:
        raise ValueError("outward_steps and return_steps must be nonnegative")
    if delays.size == 0 or np.any(delays < 1):
        raise ValueError("delay_ranks must contain positive values")

    rng = np.random.default_rng(seed)
    specs: list[EchoMotifSpec] = []
    for emission_position in emissions:
        if emission_position < 0 or emission_position >= chain.size:
            continue
        valid_delays = delays[emission_position + delays < chain.size]
        if valid_delays.size == 0:
            continue
        chosen = rng.choice(
            valid_delays,
            size=targets_per_emission,
            replace=True,
        )
        for delay in np.atleast_1d(chosen):
            specs.append(
                EchoMotifSpec(
                    emission_position=int(emission_position),
                    planted_delay_rank=int(delay),
                    outward_steps=int(outward_steps),
                    return_steps=int(return_steps),
                )
            )
    return specs
