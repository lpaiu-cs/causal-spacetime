from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_failure_decomposition import (
    CriterionFailureRecord,
    decompose_family_metric_failures,
    summarize_failure_records,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
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
        "mean_generalization_gap": 0.03,
        "stricter_threshold_pass_fraction": 0.75,
        "destructive_null_gap": 0.2,
        "symmetry_control_gap": 0.02,
        "target_coverage_fraction": 0.9,
        "pair_node_coverage_fraction": 0.7,
        "restart_std": 0.03,
        "latent_order_disagreement": 0.1,
        "no_retuning_audit_pass": 1.0,
        "failed_accounting_present": 1.0,
    }


def test_decompose_family_metric_failures_measured_heldout_failure() -> None:
    metrics = _strong_metrics()
    metrics["mean_heldout_violation"] = 0.4

    records = decompose_family_metric_failures(
        metrics,
        default_cross_family_robustness_criteria(),
    )

    heldout = next(
        record
        for record in records
        if record.criterion_name == "mean_heldout_violation"
    )
    assert heldout.status == "measured_failure"
    assert heldout.root_cause_category == "heldout_failure"


def test_decompose_family_metric_failures_missing_metric() -> None:
    metrics = _strong_metrics()
    metrics.pop("target_coverage_fraction")

    records = decompose_family_metric_failures(
        metrics,
        default_cross_family_robustness_criteria(),
    )

    coverage = next(
        record
        for record in records
        if record.criterion_name == "target_coverage_fraction"
    )
    assert coverage.status == "missing_metric"
    assert coverage.root_cause_category == "missing_metric"


def test_summarize_failure_records_groups_rows() -> None:
    records = [
        CriterionFailureRecord(
            "family",
            "structured",
            "mean_heldout_violation",
            "hard",
            0.4,
            0.2,
            "measured_failure",
            "heldout_failure",
            "failed",
        ),
        CriterionFailureRecord(
            "family",
            "structured",
            "target_coverage_fraction",
            "warning",
            float("nan"),
            0.8,
            "missing_metric",
            "missing_metric",
            "missing",
        ),
    ]

    rows = summarize_failure_records(records)

    assert {row["status"] for row in rows} == {"measured_failure", "missing_metric"}
