"""Finite state-change causal trigger networks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

TRIGGER_KINDS = {"local_successor", "external_trigger", "initial_seed"}


@dataclass(frozen=True)
class StateChangeEvent:
    """A transition of one local system from one state to another.

    A state-changing event is not a spacetime point with metric coordinates.
    ``event_id`` is implementation bookkeeping for finite storage and
    topological construction. It is not a global physical time.
    """

    event_id: int
    system_id: int
    local_index: int
    previous_state: int
    next_state: int


@dataclass(frozen=True)
class TriggerEdge:
    """Immediate causal trigger edge between state-changing events."""

    source_event_id: int
    target_event_id: int
    trigger_kind: str

    def __post_init__(self) -> None:
        if self.trigger_kind not in TRIGGER_KINDS:
            allowed = ", ".join(sorted(TRIGGER_KINDS))
            raise ValueError(f"trigger_kind must be one of: {allowed}")


@dataclass(frozen=True)
class StateChangeNetwork:
    """Finite state-change trigger graph grouped by local system chains."""

    events: list[StateChangeEvent]
    trigger_edges: list[TriggerEdge]
    system_event_ids: dict[int, list[int]]


def _add_event(
    events: list[StateChangeEvent],
    system_event_ids: dict[int, list[int]],
    system_id: int,
    local_index: int,
    previous_state: int,
    next_state: int,
) -> StateChangeEvent:
    event = StateChangeEvent(
        event_id=len(events),
        system_id=int(system_id),
        local_index=int(local_index),
        previous_state=int(previous_state),
        next_state=int(next_state),
    )
    events.append(event)
    system_event_ids.setdefault(int(system_id), []).append(event.event_id)
    return event


def generate_state_change_network(
    num_systems: int,
    max_events: int,
    *,
    initial_state: int = 0,
    trigger_probability: float = 0.25,
    max_triggers_per_event: int = 2,
    seed: int | None = None,
) -> StateChangeNetwork:
    """Generate a finite state-change causal trigger network.

    The construction loop is algorithmic bookkeeping for a finite DAG. It is
    not a global physical clock or synchronous update process.
    """

    if num_systems <= 0:
        raise ValueError("num_systems must be positive")
    if max_events < num_systems:
        raise ValueError("max_events must be at least num_systems")
    if not 0.0 <= trigger_probability <= 1.0:
        raise ValueError("trigger_probability must be in [0, 1]")
    if max_triggers_per_event < 0:
        raise ValueError("max_triggers_per_event must be nonnegative")

    rng = np.random.default_rng(seed)
    events: list[StateChangeEvent] = []
    trigger_edges: list[TriggerEdge] = []
    system_event_ids: dict[int, list[int]] = {
        system_id: [] for system_id in range(num_systems)
    }
    current_states = np.full(num_systems, int(initial_state), dtype=int)
    local_counts = np.zeros(num_systems, dtype=int)

    for system_id in range(num_systems):
        event = _add_event(
            events,
            system_event_ids,
            system_id,
            int(local_counts[system_id]),
            int(current_states[system_id]),
            int(current_states[system_id] + 1),
        )
        trigger_edges.append(
            TriggerEdge(
                source_event_id=-1,
                target_event_id=event.event_id,
                trigger_kind="initial_seed",
            )
        )
        current_states[system_id] += 1
        local_counts[system_id] += 1

    while len(events) < max_events:
        target_system = int(rng.integers(0, num_systems))
        previous_local_event_id = system_event_ids[target_system][-1]
        event = _add_event(
            events,
            system_event_ids,
            target_system,
            int(local_counts[target_system]),
            int(current_states[target_system]),
            int(current_states[target_system] + 1),
        )
        current_states[target_system] += 1
        local_counts[target_system] += 1
        trigger_edges.append(
            TriggerEdge(
                source_event_id=previous_local_event_id,
                target_event_id=event.event_id,
                trigger_kind="local_successor",
            )
        )

        if max_triggers_per_event == 0 or rng.random() > trigger_probability:
            continue
        candidate_sources = [
            older.event_id
            for older in events[: event.event_id]
            if older.system_id != target_system
        ]
        if not candidate_sources:
            continue
        trigger_count = min(
            int(rng.integers(1, max_triggers_per_event + 1)),
            len(candidate_sources),
        )
        sources = rng.choice(candidate_sources, size=trigger_count, replace=False)
        for source in np.sort(sources):
            trigger_edges.append(
                TriggerEdge(
                    source_event_id=int(source),
                    target_event_id=event.event_id,
                    trigger_kind="external_trigger",
                )
            )

    return StateChangeNetwork(
        events=events,
        trigger_edges=trigger_edges,
        system_event_ids=system_event_ids,
    )
