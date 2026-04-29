from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    decide_family_robustness,
    default_cross_family_robustness_criteria,
)


def _strong_metrics() -> dict[str, float | str]:
    return {
        "family_name": "eligible_rank_gap",
        "family_kind": "structured",
        "manifest_count": 4.0,
        "fitted_fraction": 1.0,
        "no_fit_fraction": 0.0,
        "mean_heldout_violation": 0.08,
        "mean_generalization_gap": 0.02,
        "stricter_threshold_pass_fraction": 1.0,
        "destructive_null_gap": 0.2,
        "symmetry_control_gap": 0.01,
        "target_coverage_fraction": 0.9,
        "pair_node_coverage_fraction": 0.7,
        "restart_std": 0.02,
        "latent_order_disagreement": 0.1,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
    }


def test_default_cross_family_robustness_criteria_fixed_thresholds() -> None:
    criteria = default_cross_family_robustness_criteria()

    assert criteria.max_mean_heldout_violation == 0.20
    assert criteria.min_destructive_null_gap == 0.10


def test_decide_family_robustness_carry_forward_case() -> None:
    decision = decide_family_robustness(
        _strong_metrics(),
        default_cross_family_robustness_criteria(),
    )

    assert decision.decision == "carry_forward"
    assert decision.passed


def test_decide_family_robustness_blocked_case() -> None:
    metrics = _strong_metrics()
    metrics["mean_heldout_violation"] = 0.4

    decision = decide_family_robustness(
        metrics,
        default_cross_family_robustness_criteria(),
    )

    assert decision.decision == "blocked"
    assert "high_heldout_violation" in decision.failed_reasons


def test_decide_family_robustness_report_only_for_ineligible_control() -> None:
    metrics = _strong_metrics()
    metrics["family_kind"] = "ineligible_control"

    decision = decide_family_robustness(
        metrics,
        default_cross_family_robustness_criteria(),
    )

    assert decision.decision == "report_only"
