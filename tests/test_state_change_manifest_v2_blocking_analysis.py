from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_blocking_analysis import (
    classify_v2_blocking_criterion,
    decompose_v2_blocking_by_family,
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


def test_classify_manifest_count_structural_failure() -> None:
    record = classify_v2_blocking_criterion(
        _row(),
        "manifest_count",
        1.0,
        3.0,
        "min_required",
    )

    assert record.blocking_type == "structural_blocking"
    assert record.status == "fail"


def test_classify_heldout_measured_failure() -> None:
    row = _row()
    row["manifest_count"] = 5.0
    record = classify_v2_blocking_criterion(
        row,
        "mean_heldout_violation",
        0.5,
        0.2,
        "max_allowed",
    )

    assert record.blocking_type == "measured_blocking"
    assert record.status == "fail"


def test_decompose_returns_all_criteria() -> None:
    records = decompose_v2_blocking_by_family(
        [_row()],
        default_cross_family_robustness_criteria(),
    )

    assert len(records) == 14

