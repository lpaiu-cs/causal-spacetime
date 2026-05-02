from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v4_blocking_analysis import (
    decompose_v4_blocking_by_family,
)


def _metric_row(**updates: float | str) -> dict[str, float | str]:
    row: dict[str, float | str] = {
        "family_name": "structured_family",
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
        "preconditions_passed": 1.0,
    }
    row.update(updates)
    return row


def _failed_roots(row: dict[str, float | str]) -> set[str]:
    records = decompose_v4_blocking_by_family(
        [row],
        [],
        [],
        default_cross_family_robustness_criteria(),
    )
    return {record.root_cause_category for record in records if record.status == "fail"}


def test_decompose_v4_blocking_by_family_detects_measured_failures() -> None:
    assert "heldout_failure" in _failed_roots(_metric_row(mean_heldout_violation=0.4))
    assert "stricter_pass_failure" in _failed_roots(
        _metric_row(stricter_threshold_pass_fraction=0.0)
    )
    assert "null_separation_failure" in _failed_roots(
        _metric_row(destructive_null_gap=0.0)
    )
    assert "latent_order_instability" in _failed_roots(
        _metric_row(latent_order_disagreement=0.6)
    )
    assert "coverage_failure" in _failed_roots(
        _metric_row(pair_node_coverage_fraction=0.0)
    )


def test_failed_control_becomes_control_family_blocking() -> None:
    records = decompose_v4_blocking_by_family(
        [_metric_row(family_name="failed", family_kind="failed_control")],
        [],
        [],
        default_cross_family_robustness_criteria(),
    )

    assert all(record.blocking_type == "control_family_blocking" for record in records)
