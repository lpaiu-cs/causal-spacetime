from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_carry_forward import (
    decide_v2_family_robustness,
    v2_decision_to_row,
    v2_diagnostic_completeness_by_family,
    v2_metrics_rows_from_bundle,
)


def _strong_row() -> dict[str, float | str]:
    return {
        "family_name": "strong",
        "family_kind": "structured",
        "manifest_count": 3.0,
        "fitted_fraction": 1.0,
        "no_fit_fraction": 0.0,
        "mean_heldout_violation": 0.05,
        "mean_generalization_gap": 0.01,
        "stricter_threshold_pass_fraction": 1.0,
        "destructive_null_gap": 0.2,
        "symmetry_control_gap": 0.01,
        "target_coverage_fraction": 1.0,
        "pair_node_coverage_fraction": 0.8,
        "restart_std": 0.01,
        "latent_order_disagreement": 0.1,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
    }


def test_v2_metrics_rows_from_bundle_reads_metric_rows() -> None:
    rows = v2_metrics_rows_from_bundle({"metrics": [{"family_name": "a"}]})

    assert rows == [{"family_name": "a"}]


def test_decide_v2_family_robustness_strong_synthetic_row() -> None:
    decisions = decide_v2_family_robustness(
        [_strong_row()],
        default_cross_family_robustness_criteria(),
    )

    assert decisions[0].decision == "carry_forward"


def test_v2_decision_to_row_includes_diagnostic_complete() -> None:
    decision = decide_v2_family_robustness(
        [_strong_row()],
        default_cross_family_robustness_criteria(),
    )[0]

    row = v2_decision_to_row(decision, diagnostic_complete=True)

    assert row["diagnostic_complete"] == 1.0


def test_v2_diagnostic_completeness_by_family() -> None:
    status = v2_diagnostic_completeness_by_family(
        {"completeness": [{"family_name": "a", "missing_metric_count": "0"}]}
    )

    assert status["a"]
