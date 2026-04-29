from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_comparison import (
    compare_manifest_family_against_thresholds,
    family_report_card,
    summarize_nulls_by_taxonomy,
)
from causal_spacetime_lab.state_change_manifest_family_config import (
    default_family_comparison_config,
)
from causal_spacetime_lab.state_change_manifest_family_diagnostics import (
    ManifestFitDiagnosticRow,
)


def _row(heldout: float, gap: float) -> ManifestFitDiagnosticRow:
    return ManifestFitDiagnosticRow(
        manifest_id="m1",
        family_name="eligible_rank_gap",
        family_kind="structured",
        eligible=True,
        fitted=True,
        reason_not_fit="",
        embedding_dim=1,
        train_violation_rate=heldout - gap,
        heldout_violation_rate=heldout,
        generalization_gap=gap,
        train_hinge_loss=0.0,
        heldout_hinge_loss=0.0,
        target_count=4,
        train_constraint_count=10,
        heldout_constraint_count=5,
    )


def test_default_family_comparison_config_returns_fixed_thresholds() -> None:
    config = default_family_comparison_config()

    assert config.dims == [1, 2, 3]
    assert config.heldout_violation_threshold == 0.20


def test_compare_manifest_family_against_thresholds_pass_fail() -> None:
    config = default_family_comparison_config()

    passed = compare_manifest_family_against_thresholds([_row(0.1, 0.02)], config)
    failed = compare_manifest_family_against_thresholds([_row(0.5, 0.02)], config)

    assert passed[0]["threshold_pass"] == 1.0
    assert failed[0]["threshold_pass"] == 0.0


def test_summarize_nulls_by_taxonomy() -> None:
    rows = summarize_nulls_by_taxonomy(
        [
            {
                "null_type": "shuffled_sides",
                "mean_heldout_violation_rate": 0.5,
                "structured_minus_null_mean": -0.3,
            }
        ]
    )

    assert rows[0]["taxonomy_class"] == "destructive_null"


def test_family_report_card_deterministic_case() -> None:
    cards = family_report_card(
        [
            {
                "family_name": "eligible_rank_gap",
                "family_kind": "structured",
                "fitted_count": 1.0,
                "no_fit_count": 0.0,
                "mean_heldout_violation": 0.1,
            }
        ],
        [],
        [],
    )

    assert cards[0]["family_name"] == "eligible_rank_gap"
