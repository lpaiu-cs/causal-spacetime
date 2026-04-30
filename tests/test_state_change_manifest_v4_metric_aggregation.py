from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    default_diagnostic_metric_requirements,
)
from causal_spacetime_lab.state_change_manifest_v4_metric_aggregation import (
    aggregate_v4_protocol_required_metrics,
    v4_protocol_diagnostic_completeness_check,
)


def test_aggregate_v4_protocol_required_metrics_returns_all_14_metrics() -> None:
    fit_rows = [
        {
            "family_name": "rank_gap_earliest_full_stability_v4",
            "family_kind": "structured",
            "manifest_id": "m1",
            "fit_status": "fit",
            "heldout_violation_rate": 0.1,
            "train_violation_rate": 0.05,
            "generalization_gap": 0.05,
        }
    ]
    common = {
        "family_name": "rank_gap_earliest_full_stability_v4",
        "family_kind": "structured",
    }
    rows = aggregate_v4_protocol_required_metrics(
        fit_rows=fit_rows,
        null_rows=[
            {**common, "null_type": "destructive_null", "mean_heldout": 0.4},
            {**common, "null_type": "symmetry_control", "mean_heldout": 0.15},
        ],
        stricter_rows=[{**common, "stricter_pass": 1.0}],
        failed_rows=[
            {
                **common,
                "manifest_count": 1.0,
                "failed_accounting_present": 1.0,
            }
        ],
        coverage_rows=[
            {
                **common,
                "target_coverage_fraction": 1.0,
                "pair_node_coverage_fraction": 1.0,
            }
        ],
        restart_rows=[{**common, "restart_std": 0.01}],
        latent_order_rows=[{**common, "latent_order_disagreement": 0.02}],
        no_retuning_audit_pass=True,
    )
    required = {item.metric_name for item in default_diagnostic_metric_requirements()}

    assert rows
    assert required.issubset(rows[0])
    assert "dominant_handoff_provenance_type" in rows[0]
    assert "all_structured_protocol_invariant" in rows[0]


def test_v4_protocol_diagnostic_completeness_check_detects_complete_case() -> None:
    required = {
        item.metric_name: 1.0 for item in default_diagnostic_metric_requirements()
    }
    rows = v4_protocol_diagnostic_completeness_check(
        [{"family_name": "family", **required}]
    )

    assert rows[0]["diagnostic_complete"] == 1.0
