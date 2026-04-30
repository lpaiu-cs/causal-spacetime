from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v2_metric_aggregation import (
    aggregate_v2_required_metrics,
    v2_diagnostic_completeness_check,
)


def test_aggregate_v2_required_metrics_returns_all_14_metrics() -> None:
    fit_rows = [
        {
            "family_name": "rank_gap_more_protocol_columns_v2",
            "family_kind": "structured",
            "embedding_dim": "2",
            "manifest_count": "1",
            "fitted_count": "1",
            "no_fit_count": "0",
            "mean_heldout_violation": "0.1",
            "mean_generalization_gap": "0.02",
        }
    ]
    null_rows = [
        {
            "family_name": "rank_gap_more_protocol_columns_v2",
            "null_type": "shuffled_sides",
            "taxonomy_class": "destructive_null",
            "mean_heldout_violation_rate": "0.4",
        },
        {
            "family_name": "rank_gap_more_protocol_columns_v2",
            "null_type": "permuted_targets",
            "taxonomy_class": "symmetry_control",
            "mean_heldout_violation_rate": "0.11",
        },
    ]
    rows = aggregate_v2_required_metrics(
        fit_rows=fit_rows,
        null_rows=null_rows,
        stricter_rows=[
            {
                "family_name": "rank_gap_more_protocol_columns_v2",
                "threshold_pass": "1",
            }
        ],
        failed_rows=[{"family_name": "rank_gap_more_protocol_columns_v2"}],
        restart_rows=[
            {"family_name": "rank_gap_more_protocol_columns_v2", "restart_std": "0.01"}
        ],
        latent_order_rows=[
            {
                "family_name": "rank_gap_more_protocol_columns_v2",
                "latent_order_disagreement": "0.05",
            }
        ],
        coverage_rows=[
            {
                "family_name": "rank_gap_more_protocol_columns_v2",
                "target_coverage_fraction": "1.0",
                "pair_node_coverage_fraction": "0.8",
            }
        ],
        no_retuning_audit_pass=True,
    )

    row = rows[0]
    required = [
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
    ]
    assert all(metric in row for metric in required)


def test_v2_diagnostic_completeness_check_detects_no_missing_complete_case() -> None:
    metric_rows = [
        {
            "family_name": "family",
            "family_kind": "structured",
            "manifest_count": 1.0,
            "fitted_fraction": 1.0,
            "no_fit_fraction": 0.0,
            "mean_heldout_violation": 0.1,
            "mean_generalization_gap": 0.0,
            "stricter_threshold_pass_fraction": 1.0,
            "destructive_null_gap": 0.2,
            "symmetry_control_gap": 0.01,
            "target_coverage_fraction": 1.0,
            "pair_node_coverage_fraction": 0.8,
            "restart_std": 0.02,
            "latent_order_disagreement": 0.1,
            "no_retuning_audit_pass": 1.0,
            "failed_accounting_present": 1.0,
        }
    ]

    rows = v2_diagnostic_completeness_check(metric_rows)

    assert rows[0]["missing_metric_count"] == 0.0
