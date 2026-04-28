from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    build_echo_motif_specs_for_reference_chain,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_background_state_change_network_with_reference,
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    is_irreflexive,
    topological_order_from_adjacency,
    transitive_closure_dag,
)


def test_echo_motif_spec_validation_rejects_invalid_values() -> None:
    with pytest.raises(ValueError):
        EchoMotifSpec(emission_position=-1, planted_delay_rank=1)
    with pytest.raises(ValueError):
        EchoMotifSpec(emission_position=0, planted_delay_rank=0)
    with pytest.raises(ValueError):
        EchoMotifSpec(emission_position=0, planted_delay_rank=1, outward_steps=-1)
    with pytest.raises(ValueError):
        EchoMotifSpec(emission_position=0, planted_delay_rank=1, return_steps=-1)


def test_generate_reference_backbone_network_creates_valid_chain() -> None:
    network, reference = generate_reference_backbone_network(5)
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))

    assert np.array_equal(reference, np.arange(5))
    assert all(closure[reference[i], reference[i + 1]] for i in range(4))
    assert is_irreflexive(closure)


def test_insert_echo_motif_adds_target_and_trigger_edges() -> None:
    network, reference = generate_reference_backbone_network(8)
    updated, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 3, outward_steps=1, return_steps=1),
    )

    assert len(updated.events) == len(network.events) + 3
    assert motif.target_event_id >= len(network.events)
    assert motif.reference_emission_event_id == reference[2]
    assert motif.reference_return_event_id == reference[5]
    assert len(updated.trigger_edges) > len(network.trigger_edges)


def test_inserted_motif_is_acyclic() -> None:
    network, reference = generate_reference_backbone_network(8)
    updated, _ = insert_echo_motif(network, reference, EchoMotifSpec(2, 3))

    order = topological_order_from_adjacency(immediate_trigger_adjacency(updated))

    assert order.size == len(updated.events)


def test_build_echo_motif_specs_for_reference_chain_generates_valid_specs() -> None:
    reference = np.arange(8)
    specs = build_echo_motif_specs_for_reference_chain(
        reference,
        emission_positions=np.asarray([1, 2, 6]),
        delay_ranks=np.asarray([2, 3]),
        targets_per_emission=2,
        seed=0,
    )

    assert specs
    assert all(
        spec.emission_position + spec.planted_delay_rank < reference.size
        for spec in specs
    )


def test_background_generator_returns_network_and_reference_chain() -> None:
    network, reference = generate_background_state_change_network_with_reference(
        4,
        30,
        0.2,
        seed=0,
    )

    assert len(network.events) <= 30
    assert reference.ndim == 1
    assert reference.size >= 1
