from __future__ import annotations

import csv
import json
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v4_outputs import (
    V4_PROTOCOL_OUTPUT_FILES,
)


def family_from_manifest_dir(manifest_dir: Path) -> str:
    manifest_path = sorted(manifest_dir.glob("*.json"))[0]
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    label = str(payload["profile_label"])
    return label.rsplit("_m", maxsplit=1)[0]


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def strong_metric_row(
    family_name: str = "strong_v4",
    family_kind: str = "structured",
) -> dict[str, str]:
    return {
        "family_name": family_name,
        "family_kind": family_kind,
        "manifest_count": "3",
        "fitted_fraction": "1",
        "no_fit_fraction": "0",
        "mean_heldout_violation": "0.05",
        "mean_generalization_gap": "0.01",
        "stricter_threshold_pass_fraction": "1",
        "destructive_null_gap": "0.2",
        "symmetry_control_gap": "0.01",
        "target_coverage_fraction": "1",
        "pair_node_coverage_fraction": "0.8",
        "restart_std": "0.01",
        "latent_order_disagreement": "0.1",
        "no_retuning_audit_pass": "1",
        "failed_accounting_present": "1",
        "dominant_handoff_provenance_type": "hybrid_template_instantiated_from_profile",
        "top_down_template_count": "0",
        "hybrid_manifest_count": "3",
        "bottom_up_manifest_count": "0",
        "report_only_manifest_count": "0",
        "all_manifests_have_provenance": "1",
    }


def write_v4_protocol_bundle(output_dir: Path, family_name: str) -> None:
    data_dir = output_dir / "data"
    metrics = [strong_metric_row(family_name)]
    write_csv(
        data_dir / "v4_protocol_cross_family_robustness_metrics.csv",
        metrics,
    )
    write_csv(
        data_dir / "v4_protocol_diagnostic_completeness_check.csv",
        [
            {
                "family_name": family_name,
                "family_kind": "structured",
                "required_metric_count": "14",
                "available_metric_count": "14",
                "missing_metric_count": "0",
                "completeness_fraction": "1",
                "diagnostic_complete": "1",
                "missing_metrics": "",
            }
        ],
    )
    write_csv(
        data_dir / "v4_protocol_manifest_generation.csv",
        [
            {
                "family_name": family_name,
                "family_kind": "structured",
                "manifest_count": "2",
            }
        ],
    )
    for _key, filename in V4_PROTOCOL_OUTPUT_FILES.items():
        path = data_dir / filename
        if path.exists():
            continue
        write_csv(
            path,
            [
                {
                    "family_name": family_name,
                    "family_kind": "structured",
                    "ok": "1",
                }
            ],
        )
