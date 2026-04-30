from __future__ import annotations

import csv
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v3_protocol_manifest_drilldown import (
    load_v3_protocol_manifest_drilldown_rows,
    summarize_manifest_drilldown_by_family,
)


def _write(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_manifest_level_drilldown_groups_by_family(tmp_path: Path) -> None:
    data = tmp_path / "data"
    _write(
        data / "v3_protocol_manifest_family_fit_comparison.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "family",
                "family_kind": "structured",
                "train_violation_rate": "0.1",
                "heldout_violation_rate": "0.3",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_family_coverage_metrics.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "family",
                "target_coverage_fraction": "1.0",
                "pair_node_coverage_fraction": "1.0",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_family_restart_stability.csv",
        [{"manifest_id": "m1", "family_name": "family", "restart_std": "0.1"}],
    )
    _write(
        data / "v3_protocol_manifest_family_latent_order_stability.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "family",
                "latent_order_disagreement": "0.4",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_metadata_audit.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "family",
                "measurement_protocol_id": "p",
                "profile_invariance_status": "protocol_invariant",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_family_failed_accounting.csv",
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

    rows = load_v3_protocol_manifest_drilldown_rows(tmp_path)
    summary = summarize_manifest_drilldown_by_family(rows)

    assert rows[0].family_name == "family"
    assert summary[0]["manifest_count"] == 1.0
