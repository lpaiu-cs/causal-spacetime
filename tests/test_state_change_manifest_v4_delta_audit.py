from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v4_delta_audit import (
    planned_v3_to_v4_family_links,
    summarize_v3_to_v4_deltas,
    v3_to_v4_metric_delta_rows,
)


def test_v3_to_v4_delta_audit_uses_correct_improvement_direction() -> None:
    link = planned_v3_to_v4_family_links()[0]
    v3_family = link["v3_family_name"]
    v4_family = link["v4_family_name"]
    records = v3_to_v4_metric_delta_rows(
        [
            {
                "family_name": v3_family,
                "mean_heldout_violation": 0.4,
                "stricter_threshold_pass_fraction": 0.2,
                "destructive_null_gap": 0.1,
                "symmetry_control_gap": 0.2,
                "target_coverage_fraction": 0.7,
                "pair_node_coverage_fraction": 0.5,
                "restart_std": 0.2,
                "latent_order_disagreement": 0.5,
                "mean_generalization_gap": 0.2,
            }
        ],
        [
            {
                "family_name": v4_family,
                "mean_heldout_violation": 0.3,
                "stricter_threshold_pass_fraction": 0.3,
                "destructive_null_gap": 0.2,
                "symmetry_control_gap": 0.1,
                "target_coverage_fraction": 0.8,
                "pair_node_coverage_fraction": 0.6,
                "restart_std": 0.1,
                "latent_order_disagreement": 0.4,
                "mean_generalization_gap": 0.1,
            }
        ],
        [link],
    )
    summary = summarize_v3_to_v4_deltas(records)

    assert all(record.improved for record in records)
    assert summary[0]["worsened_metric_count"] == 0.0
