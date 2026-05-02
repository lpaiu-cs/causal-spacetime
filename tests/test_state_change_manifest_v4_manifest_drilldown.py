from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v4_manifest_drilldown import (
    load_v4_manifest_drilldown_rows,
    summarize_v4_manifest_drilldown_by_family,
)
from tests.v4_protocol_test_helpers import write_csv


def test_manifest_level_drilldown_groups_by_family(tmp_path: Path) -> None:
    data = tmp_path / "data"
    write_csv(
        data / "v4_protocol_manifest_family_fit_comparison.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "family",
                "family_kind": "structured",
                "train_violation_rate": "0.1",
                "heldout_violation_rate": "0.4",
            }
        ],
    )
    write_csv(
        data / "v4_protocol_manifest_family_coverage_metrics.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "family",
                "target_coverage_fraction": "1",
                "pair_node_coverage_fraction": "0.5",
            }
        ],
    )
    write_csv(
        data / "v4_protocol_manifest_family_restart_stability.csv",
        [{"manifest_id": "m1", "family_name": "family", "restart_std": "0.1"}],
    )
    write_csv(
        data / "v4_protocol_manifest_family_latent_order_stability.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "family",
                "latent_order_disagreement": "0.5",
            }
        ],
    )
    write_csv(
        data / "v4_protocol_manifest_metadata_audit.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "family",
                "measurement_protocol_id": "p",
                "profile_invariance_status": "protocol_invariant",
            }
        ],
    )
    write_csv(
        data / "v4_protocol_manifest_family_failed_accounting.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "family",
                "eligible": "1",
                "failed_reasons": "",
                "handoff_provenance_type": "hybrid_template_instantiated_from_profile",
            }
        ],
    )

    rows = load_v4_manifest_drilldown_rows(tmp_path)
    summary = summarize_v4_manifest_drilldown_by_family(rows)

    assert rows[0].family_name == "family"
    assert summary[0]["manifest_count"] == 1.0
