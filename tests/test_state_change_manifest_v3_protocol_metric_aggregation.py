from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v3_protocol_metric_aggregation import (
    aggregate_v3_protocol_required_metrics,
    v3_protocol_diagnostic_completeness_check,
)


def _complete_inputs() -> dict[str, list[dict[str, float | str]]]:
    fit = [
        {
            "family_name": "rank_gap_earliest_full_reference_v3",
            "family_kind": "structured",
            "manifest_count": 3.0,
            "fitted_count": 3.0,
            "no_fit_count": 0.0,
            "mean_heldout_violation": 0.1,
            "mean_generalization_gap": 0.02,
        }
    ]
    null = [
        {
            "family_name": "rank_gap_earliest_full_reference_v3",
            "family_kind": "structured",
            "taxonomy_class": "destructive_null",
            "mean_heldout_violation_rate": 0.3,
            "structured_heldout_violation_rate": 0.1,
        },
        {
            "family_name": "rank_gap_earliest_full_reference_v3",
            "family_kind": "structured",
            "taxonomy_class": "symmetry_control",
            "mean_heldout_violation_rate": 0.12,
            "structured_heldout_violation_rate": 0.1,
        },
    ]
    stricter = [
        {
            "family_name": "rank_gap_earliest_full_reference_v3",
            "family_kind": "structured",
            "threshold_pass": 1.0,
        }
    ]
    failed = [
        {
            "family_name": "rank_gap_earliest_full_reference_v3",
            "family_kind": "structured",
            "handoff_provenance_type": "hybrid_template_instantiated_from_profile",
        }
    ]
    coverage = [
        {
            "family_name": "rank_gap_earliest_full_reference_v3",
            "family_kind": "structured",
            "target_coverage_fraction": 1.0,
            "pair_node_coverage_fraction": 0.8,
        }
    ]
    restart = [
        {
            "family_name": "rank_gap_earliest_full_reference_v3",
            "family_kind": "structured",
            "restart_std": 0.01,
        }
    ]
    latent = [
        {
            "family_name": "rank_gap_earliest_full_reference_v3",
            "family_kind": "structured",
            "latent_order_disagreement": 0.1,
        }
    ]
    return {
        "fit_rows": fit,
        "null_rows": null,
        "stricter_rows": stricter,
        "failed_rows": failed,
        "coverage_rows": coverage,
        "restart_rows": restart,
        "latent_order_rows": latent,
    }


def test_aggregate_v3_protocol_required_metrics_returns_all_required_metrics() -> None:
    inputs = _complete_inputs()
    rows = aggregate_v3_protocol_required_metrics(
        **inputs,
        no_retuning_audit_pass=True,
    )
    row = rows[0]

    for key in (
        "manifest_count",
        "fitted_fraction",
        "no_fit_fraction",
        "mean_heldout_violation",
        "mean_generalization_gap",
        "stricter_threshold_pass_fraction",
        "destructive_null_gap",
        "symmetry_control_gap",
        "target_coverage_fraction",
        "pair_node_coverage_fraction",
        "restart_std",
        "latent_order_disagreement",
        "no_retuning_audit_pass",
        "failed_accounting_present",
    ):
        assert key in row
    assert row["hybrid_manifest_count"] == 1.0


def test_v3_protocol_diagnostic_completeness_check_detects_complete_case() -> None:
    inputs = _complete_inputs()
    rows = aggregate_v3_protocol_required_metrics(
        **inputs,
        no_retuning_audit_pass=True,
    )
    completeness = v3_protocol_diagnostic_completeness_check(rows)

    assert completeness[0]["diagnostic_complete"] == 1.0
    assert completeness[0]["missing_metrics"] == ""

