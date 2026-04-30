from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_blocking_analysis import (
    decompose_v2_blocking_by_family,
)
from causal_spacetime_lab.state_change_manifest_v2_counterfactuals import (
    structural_count_counterfactual_report,
    would_remain_blocked_ignoring_structural_count,
)


def _row() -> dict[str, float | str]:
    return {
        "family_name": "f",
        "family_kind": "structured",
        "manifest_count": 1.0,
        "fitted_fraction": 1.0,
        "no_fit_fraction": 0.0,
        "mean_heldout_violation": 0.5,
        "mean_generalization_gap": 0.0,
        "stricter_threshold_pass_fraction": 1.0,
        "destructive_null_gap": 0.2,
        "symmetry_control_gap": 0.01,
        "target_coverage_fraction": 1.0,
        "pair_node_coverage_fraction": 1.0,
        "restart_std": 0.01,
        "latent_order_disagreement": 0.01,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
    }


def test_structural_count_counterfactual_is_report_only() -> None:
    rows = structural_count_counterfactual_report(
        [_row()],
        default_cross_family_robustness_criteria(),
        hypothetical_manifest_count=3,
    )

    assert rows[0]["manifest_count_would_pass"] == 1.0
    assert rows[0]["analysis_label"] == "report_only_counterfactual_not_decision_change"


def test_remain_blocked_ignoring_structural_count_detects_measured() -> None:
    criteria = default_cross_family_robustness_criteria()
    records = decompose_v2_blocking_by_family([_row()], criteria)
    rows = would_remain_blocked_ignoring_structural_count(records)

    assert rows[0]["would_remain_blocked"] == 1.0
    assert "mean_heldout_violation" in str(rows[0]["measured_blockers"])

