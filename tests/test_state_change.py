from __future__ import annotations

from causal_spacetime_lab.state_change import (
    StateChangeEvent,
    generate_state_change_network,
)
from causal_spacetime_lab.state_change_validation import (
    branching_statistics,
    state_change_network_summary,
)


def test_state_change_event_stores_transition_fields() -> None:
    event = StateChangeEvent(
        event_id=7,
        system_id=2,
        local_index=3,
        previous_state=4,
        next_state=5,
    )

    assert event.event_id == 7
    assert event.system_id == 2
    assert event.local_index == 3
    assert event.previous_state == 4
    assert event.next_state == 5


def test_generate_state_change_network_event_count_and_unique_ids() -> None:
    network = generate_state_change_network(4, 30, seed=11)

    ids = [event.event_id for event in network.events]
    assert len(network.events) <= 30
    assert len(ids) == len(set(ids))
    assert sorted(ids) == list(range(len(network.events)))


def test_generate_state_change_network_local_indices_increase() -> None:
    network = generate_state_change_network(3, 25, seed=12)

    for event_ids in network.system_event_ids.values():
        local_indices = [network.events[event_id].local_index for event_id in event_ids]
        assert local_indices == list(range(len(event_ids)))


def test_trigger_edges_point_from_older_ids_to_newer_ids() -> None:
    network = generate_state_change_network(4, 40, seed=13)

    for edge in network.trigger_edges:
        assert edge.source_event_id < edge.target_event_id


def test_state_change_network_summary_has_expected_keys() -> None:
    network = generate_state_change_network(3, 20, seed=14)

    summary = state_change_network_summary(network)
    branching = branching_statistics(network)

    assert summary["event_count"] == 20.0
    assert summary["system_count"] == 3.0
    assert "local_successor_edge_count" in summary
    assert "external_branching_event_count" in branching
