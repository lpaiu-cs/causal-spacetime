from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_echo_motif_validation import (
    motif_order_recovery_rate,
    motif_recovery_report,
    recovered_delay_for_motif,
    summarize_motif_recovery,
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


def _clean_motif_case():
    network, reference = generate_reference_backbone_network(8)
    network, motif = insert_echo_motif(
        network,
        reference,
        EchoMotifSpec(2, 3, outward_steps=1, return_steps=1),
    )
    closure = transitive_closure_dag(immediate_trigger_adjacency(network))
    return closure, reference, motif


def test_recovered_delay_for_motif_equals_planted_delay_in_clean_case() -> None:
    closure, reference, motif = _clean_motif_case()

    assert recovered_delay_for_motif(closure, reference, motif) == 3


def test_motif_recovery_report_exact_case() -> None:
    closure, reference, motif = _clean_motif_case()

    report = motif_recovery_report(closure, reference, motif)

    assert report["planted_delay_rank"] == 3.0
    assert report["recovered_delay_rank"] == 3.0
    assert report["exact_recovery"] == 1.0
    assert report["early_shortcut"] == 0.0


def test_summarize_motif_recovery_fields() -> None:
    rows = [
        {
            "exact_recovery": 1.0,
            "early_shortcut": 0.0,
            "late_or_missing": 0.0,
            "recovery_error": 0.0,
            "recovered_delay_rank": 3.0,
            "planted_delay_rank": 3.0,
        },
        {
            "exact_recovery": 0.0,
            "early_shortcut": 1.0,
            "late_or_missing": 0.0,
            "recovery_error": -1.0,
            "recovered_delay_rank": 2.0,
            "planted_delay_rank": 3.0,
        },
    ]

    summary = summarize_motif_recovery(rows)

    assert summary["motif_count"] == 2.0
    assert summary["exact_recovery_fraction"] == pytest.approx(0.5)
    assert summary["early_shortcut_fraction"] == pytest.approx(0.5)
    assert summary["mean_absolute_recovery_error"] == pytest.approx(0.5)


def test_motif_order_recovery_rate_deterministic_rows() -> None:
    rows = [
        {"planted_delay_rank": 2.0, "recovered_delay_rank": 2.0},
        {"planted_delay_rank": 3.0, "recovered_delay_rank": 3.0},
        {"planted_delay_rank": 5.0, "recovered_delay_rank": 5.0},
    ]

    order = motif_order_recovery_rate(rows)

    assert order["comparable_motif_count"] == 3.0
    assert order["order_inversion_rate"] == 0.0
    assert order["order_agreement_rate"] == 1.0


def test_motif_order_recovery_rate_ignores_missing_recovery() -> None:
    rows = [
        {"planted_delay_rank": 2.0, "recovered_delay_rank": 2.0},
        {"planted_delay_rank": 3.0, "recovered_delay_rank": np.nan},
    ]

    order = motif_order_recovery_rate(rows)

    assert order["comparable_motif_count"] == 1.0
    assert np.isnan(order["order_inversion_rate"])
