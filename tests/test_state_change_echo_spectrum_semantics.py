from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_echo_spectrum_semantics import (
    compress_suffix_spectrum,
    full_transitive_return_spectrum,
    gated_echo_delay_from_spectrum,
    immediate_edge_return_spectrum,
    is_suffix_spectrum,
    retained_reference_return_spectrum,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)


def _case():
    network, reference = generate_reference_backbone_network(8)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 4, outward_steps=1, return_steps=0),
    )
    adjacency = immediate_trigger_adjacency(network)
    closure = transitive_closure_dag(adjacency)
    return reference, motif, adjacency, closure


def test_full_transitive_return_spectrum_returns_suffix() -> None:
    reference, motif, _, closure = _case()

    spectrum = full_transitive_return_spectrum(
        closure,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )

    assert np.array_equal(spectrum, np.asarray([4, 5]))
    assert is_suffix_spectrum(spectrum)


def test_retained_reference_return_spectrum_can_be_sparse() -> None:
    reference, motif, _, closure = _case()
    retained = reference[np.asarray([0, 2, 4, 7])]

    spectrum = retained_reference_return_spectrum(
        closure,
        retained,
        motif.target_event_id,
        emission_position=1,
    )

    assert np.array_equal(spectrum, np.asarray([2]))


def test_immediate_edge_return_spectrum_is_not_full_suffix() -> None:
    reference, motif, adjacency, _ = _case()

    spectrum = immediate_edge_return_spectrum(
        adjacency,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )

    assert np.array_equal(spectrum, np.asarray([4]))


def test_is_suffix_spectrum() -> None:
    assert is_suffix_spectrum([2, 3, 4])
    assert not is_suffix_spectrum([2, 4])
    assert is_suffix_spectrum([])


def test_compress_suffix_spectrum() -> None:
    report = compress_suffix_spectrum([2, 3, 4])

    assert report["is_empty"] == 0.0
    assert report["is_suffix"] == 1.0
    assert report["earliest_delay"] == 2.0
    assert report["latest_delay"] == 4.0


def test_gated_echo_delay_from_spectrum() -> None:
    assert gated_echo_delay_from_spectrum([2, 4, 5], 3) == 4
    assert gated_echo_delay_from_spectrum([2], 3) is None
