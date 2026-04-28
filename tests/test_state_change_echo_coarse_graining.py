from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_echo_coarse_graining import (
    expected_coarse_delay_rank_for_motif,
    protected_indices_for_reference_and_motifs,
    remap_echo_motif_record_for_event_thinning,
    remap_reference_chain,
    restrict_transitive_order_to_retained_events,
    return_spectrum_stability_report,
    sample_retained_indices,
    spectrum_jaccard,
    subsample_reference_chain_positions,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)


def _motif_case():
    network, reference = generate_reference_backbone_network(8)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 3, outward_steps=1, return_steps=1),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    return network, closure, reference, motif


def test_protected_indices_for_reference_and_motifs() -> None:
    _, _, reference, motif = _motif_case()

    protected = protected_indices_for_reference_and_motifs(reference, [motif])
    protected_with_paths = protected_indices_for_reference_and_motifs(
        reference,
        [motif],
        include_motif_paths=True,
    )

    assert set(reference).issubset(set(protected))
    assert motif.target_event_id in protected
    assert set(motif.outward_event_ids).issubset(set(protected_with_paths))
    assert set(motif.return_event_ids).issubset(set(protected_with_paths))


def test_sample_retained_indices_keeps_protected_events() -> None:
    protected = np.asarray([1, 3, 5])

    retained = sample_retained_indices(8, 0.0, protected, seed=0)

    assert np.array_equal(retained, protected)
    with pytest.raises(ValueError):
        sample_retained_indices(8, -0.1, protected)


def test_restrict_transitive_order_to_retained_events_preserves_closure() -> None:
    _, closure, reference, motif = _motif_case()
    retained = np.asarray([reference[2], motif.target_event_id, reference[5]])

    result = restrict_transitive_order_to_retained_events(closure, retained)

    assert result.restricted_order_matrix[
        result.old_to_new[int(reference[2])],
        result.old_to_new[int(motif.target_event_id)],
    ]
    assert result.restricted_order_matrix[
        result.old_to_new[int(motif.target_event_id)],
        result.old_to_new[int(reference[5])],
    ]


def test_remap_reference_chain() -> None:
    reference = np.asarray([10, 20, 30])
    mapping = {10: 0, 20: 1, 30: 2}

    remapped = remap_reference_chain(reference, mapping)

    assert np.array_equal(remapped, np.asarray([0, 1, 2]))
    with pytest.raises(KeyError):
        remap_reference_chain(reference, {10: 0})


def test_remap_echo_motif_record_for_event_thinning() -> None:
    _, _, reference, motif = _motif_case()
    protected = protected_indices_for_reference_and_motifs(reference, [motif])
    mapping = {int(old): index for index, old in enumerate(protected)}

    remapped = remap_echo_motif_record_for_event_thinning(motif, mapping)

    assert remapped is not None
    assert remapped.target_event_id == mapping[motif.target_event_id]
    assert remapped.outward_event_ids.size == 0


def test_spectrum_jaccard() -> None:
    assert spectrum_jaccard([1, 2, 3], [2, 3, 4]) == pytest.approx(0.5)
    assert spectrum_jaccard([], []) == 1.0


def test_return_spectrum_stability_report() -> None:
    report = return_spectrum_stability_report(np.asarray([2, 3]), np.asarray([3]))

    assert report["baseline_size"] == 2.0
    assert report["coarse_size"] == 1.0
    assert report["earliest_delay_shift"] == 1.0
    assert report["exact_earliest_preserved"] == 0.0


def test_subsample_reference_chain_positions() -> None:
    reference = np.arange(8)

    result = subsample_reference_chain_positions(
        reference,
        stride=3,
        protected_positions=np.asarray([2]),
    )

    assert 2 in result.retained_reference_positions
    assert 7 in result.retained_reference_positions
    assert np.array_equal(
        result.subsampled_reference_chain,
        reference[result.retained_reference_positions],
    )


def test_expected_coarse_delay_rank_for_motif() -> None:
    _, _, reference, motif = _motif_case()
    result = subsample_reference_chain_positions(
        reference,
        stride=2,
        protected_positions=np.asarray([motif.emission_position]),
    )

    assert expected_coarse_delay_rank_for_motif(motif, result) == 2
