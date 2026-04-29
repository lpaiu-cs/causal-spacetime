from __future__ import annotations

import pytest

from causal_spacetime_lab.state_change_manifest_family_robustness_metrics import (
    aggregate_family_robustness_metrics,
    destructive_null_gap_from_rows,
    symmetry_control_gap_from_rows,
)


def test_destructive_null_gap_from_rows_deterministic_case() -> None:
    rows = [
        {
            "null_type": "shuffled_sides",
            "mean_heldout_violation_rate": "0.5",
        }
    ]

    assert destructive_null_gap_from_rows("family", rows, 0.2) == pytest.approx(0.3)


def test_symmetry_control_gap_from_rows_deterministic_case() -> None:
    rows = [
        {
            "null_type": "permuted_targets",
            "mean_heldout_violation_rate": "0.22",
        }
    ]

    assert symmetry_control_gap_from_rows("family", rows, 0.2) == pytest.approx(0.02)


def test_aggregate_family_robustness_metrics_expected_keys() -> None:
    bundle = {
        "fit_summary": [
            {
                "family_name": "eligible_rank_gap",
                "family_kind": "structured",
                "manifest_count": "4",
                "fitted_count": "4",
                "no_fit_count": "0",
                "mean_heldout_violation": "0.1",
                "mean_generalization_gap": "0.02",
            }
        ],
        "null_taxonomy": [
            {
                "null_type": "shuffled_sides",
                "mean_heldout_violation_rate": "0.4",
            }
        ],
        "stricter_criteria": [
            {
                "family_name": "eligible_rank_gap",
                "threshold_pass": "1",
            }
        ],
        "failed_accounting": [{"row_type": "family_count"}],
        "no_retuning_audit": [{"passed": "1"}],
        "report_card": [],
    }

    rows = aggregate_family_robustness_metrics(bundle)

    assert rows[0]["family_name"] == "eligible_rank_gap"
    assert "destructive_null_gap" in rows[0]
