from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_echo_interference import (
    earliest_echo_delay_from_spectrum,
    return_delay_spectrum,
    return_positions_for_target,
    return_spectrum_report_for_motif,
    shortcut_depth,
    shortcut_positions_before_planted,
    summarize_return_spectrum_reports,
)
from causal_spacetime_lab.state_change_echo_motifs import (
    EchoMotifSpec,
    insert_echo_motif,
)
from causal_spacetime_lab.state_change_echo_scenarios import (
    generate_reference_backbone_network,
)
from causal_spacetime_lab.state_change_echo_shortcuts import (
    ShortcutInjectionSpec,
    inject_shortcut_returns,
)
from causal_spacetime_lab.state_change_order import (
    immediate_trigger_adjacency,
    transitive_closure_dag,
)


def _clean_motif_case():
    network, reference = generate_reference_backbone_network(8)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 3, outward_steps=1, return_steps=1),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    return network, closure, reference, motif


def test_return_positions_for_target_deterministic_example() -> None:
    _, closure, reference, motif = _clean_motif_case()

    positions = return_positions_for_target(
        closure,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )

    assert np.array_equal(positions, np.asarray([5, 6, 7]))


def test_return_delay_spectrum_deterministic_example() -> None:
    _, closure, reference, motif = _clean_motif_case()

    spectrum = return_delay_spectrum(
        closure,
        reference,
        motif.target_event_id,
        motif.emission_position,
    )

    assert np.array_equal(spectrum, np.asarray([3, 4, 5]))


def test_earliest_echo_delay_from_spectrum() -> None:
    assert earliest_echo_delay_from_spectrum(np.asarray([5, 3, 4])) == 3
    assert earliest_echo_delay_from_spectrum(np.asarray([], dtype=int)) is None


def test_shortcut_positions_before_planted() -> None:
    shortcuts = shortcut_positions_before_planted(np.asarray([2, 3, 5]), 3)

    assert np.array_equal(shortcuts, np.asarray([2]))


def test_shortcut_depth_exact_early_late_missing() -> None:
    assert shortcut_depth(3, 3) == 0.0
    assert shortcut_depth(2, 3) == 1.0
    assert shortcut_depth(4, 3) == -1.0
    assert np.isnan(shortcut_depth(None, 3))


def test_return_spectrum_report_for_motif_exact_case() -> None:
    _, closure, reference, motif = _clean_motif_case()

    report = return_spectrum_report_for_motif(closure, reference, motif)

    assert report["planted_delay_rank"] == 3.0
    assert report["recovered_delay_rank"] == 3.0
    assert report["exact_recovery"] == 1.0
    assert report["early_shortcut"] == 0.0
    assert report["return_spectrum_string"] == "3;4;5"


def test_return_spectrum_report_for_motif_shortcut_case() -> None:
    network, _, reference, motif = _clean_motif_case()
    network, _ = inject_shortcut_returns(
        network,
        reference,
        [motif],
        ShortcutInjectionSpec(probability=1.0, min_depth=1, max_depth=1),
        seed=0,
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))

    report = return_spectrum_report_for_motif(closure, reference, motif)

    assert report["recovered_delay_rank"] == 2.0
    assert report["early_shortcut"] == 1.0
    assert report["shortcut_depth"] == 1.0


def test_summarize_return_spectrum_reports_fields() -> None:
    rows = [
        {
            "exact_recovery": 1.0,
            "early_shortcut": 0.0,
            "missing_return": 0.0,
            "late_recovery": 0.0,
            "shortcut_count": 0.0,
            "shortcut_depth": 0.0,
            "spectrum_size": 3.0,
            "earliest_delay": 3.0,
            "latest_delay": 5.0,
        },
        {
            "exact_recovery": 0.0,
            "early_shortcut": 1.0,
            "missing_return": 0.0,
            "late_recovery": 0.0,
            "shortcut_count": 1.0,
            "shortcut_depth": 2.0,
            "spectrum_size": 4.0,
            "earliest_delay": 2.0,
            "latest_delay": 5.0,
        },
    ]

    summary = summarize_return_spectrum_reports(rows)

    assert summary["motif_count"] == 2.0
    assert summary["exact_recovery_fraction"] == pytest.approx(0.5)
    assert summary["shortcut_fraction"] == pytest.approx(0.5)
    assert summary["mean_shortcut_depth"] == pytest.approx(2.0)
