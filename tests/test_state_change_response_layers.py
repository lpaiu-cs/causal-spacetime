from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_echo_motifs import insert_multiple_echo_motifs
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)
from causal_spacetime_lab.state_change_response_layers import (
    EchoResponseLayerSpec,
    build_layered_echo_motif_specs,
    planted_layer_labels_for_motifs,
    planted_response_signature_from_motifs,
)
from causal_spacetime_lab.state_change_response_signature import (
    echo_response_signature_from_motifs,
)


def test_build_layered_echo_motif_specs() -> None:
    _, reference = generate_reference_backbone_network(12)

    specs = build_layered_echo_motif_specs(
        reference,
        2,
        [EchoResponseLayerSpec(3, 2), EchoResponseLayerSpec(5, 1)],
        seed=0,
    )

    assert len(specs) == 3
    assert sorted(spec.planted_delay_rank for spec in specs) == [3, 3, 5]
    with pytest.raises(ValueError):
        build_layered_echo_motif_specs(
            reference,
            9,
            [EchoResponseLayerSpec(5, 1)],
        )


def test_planted_response_signature_from_motifs() -> None:
    network, reference = generate_reference_backbone_network(12)
    specs = build_layered_echo_motif_specs(
        reference,
        2,
        [EchoResponseLayerSpec(3, 2), EchoResponseLayerSpec(5, 1)],
        seed=0,
    )
    network, motifs = insert_multiple_echo_motifs(network, reference, specs)

    planted = planted_response_signature_from_motifs(motifs)
    labels = planted_layer_labels_for_motifs(motifs)

    assert planted.target_event_ids.size == 3
    assert sorted(planted.delay_ranks.tolist()) == [3, 3, 5]
    assert set(labels) == set(planted.target_event_ids.tolist())


def test_recovered_signature_matches_clean_layered_motifs() -> None:
    network, reference = generate_reference_backbone_network(12)
    specs = build_layered_echo_motif_specs(
        reference,
        2,
        [EchoResponseLayerSpec(3, 2), EchoResponseLayerSpec(5, 1)],
        seed=0,
    )
    network, motifs = insert_multiple_echo_motifs(network, reference, specs)
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))

    recovered = echo_response_signature_from_motifs(closure, reference, motifs)

    assert np.array_equal(np.sort(recovered.delay_ranks), np.asarray([3, 3, 5]))
    assert np.all(recovered.reachable_mask)

