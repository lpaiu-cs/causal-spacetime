from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_blocking_analysis import (
    V3ProtocolBlockingRecord,
    decompose_v3_protocol_blocking_by_family,
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


def test_decompose_v3_protocol_blocking_detects_heldout_failure() -> None:
    records = decompose_v3_protocol_blocking_by_family(
        [_metric_row(mean_heldout_violation=0.4)],
        [],
        [],
        default_cross_family_robustness_criteria(),
    )

    assert any(
        record.root_cause_category == "heldout_failure"
        and record.status == "fail"
        for record in records
    )


def test_decompose_v3_protocol_blocking_detects_latent_order_failure() -> None:
    records = decompose_v3_protocol_blocking_by_family(
        [_metric_row(latent_order_disagreement=0.6)],
        [],
        [],
        default_cross_family_robustness_criteria(),
    )

    assert any(
        record.root_cause_category == "latent_order_instability"
        and record.status == "fail"
        for record in records
    )


def test_failed_control_becomes_control_family_blocking() -> None:
    records = decompose_v3_protocol_blocking_by_family(
        [_metric_row(family_name="failed", family_kind="failed_control")],
        [],
        [],
        default_cross_family_robustness_criteria(),
    )

    assert all(record.blocking_type == "control_family_blocking" for record in records)


def test_v3_protocol_blocking_record_validates_allowed_values() -> None:
    record = V3ProtocolBlockingRecord(
        family_name="family",
        family_kind="structured",
        criterion_name="mean_heldout_violation",
        observed_value=0.1,
        threshold_value=0.2,
        criterion_direction="max_allowed",
        blocking_type="not_blocking",
        status="pass",
        root_cause_category="heldout_failure",
        explanation="passes",
    )

    assert record.status == "pass"
