"""Final sanity checks for M44 v4 protocol manifest generation."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    default_diagnostic_metric_requirements,
)
from causal_spacetime_lab.state_change_manifest_v4_diagnostics import (
    manifest_v4_metadata_audit_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 protocol final sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    manifest_paths = sorted((config.output_dir / "manifests_v4").glob("*.json"))
    manifests = [
        json.loads(path.read_text(encoding="utf-8")) for path in manifest_paths
    ]
    audit_rows = manifest_v4_metadata_audit_rows(config.output_dir / "manifests_v4")
    structured = [row for row in audit_rows if row.get("family_kind") == "structured"]
    metrics = read_csv_rows(
        data_path(config.output_dir, "v4_protocol_cross_family_robustness_metrics.csv")
    )
    required = {item.metric_name for item in default_diagnostic_metric_requirements()}
    metrics_have_required = bool(metrics) and all(
        required.issubset(set(row)) for row in metrics
    )
    return [
        {"check": "v4_manifests_exist", "passed": float(bool(manifest_paths))},
        {
            "check": "v4_manifests_not_written_to_older_dirs",
            "passed": float(
                not list((config.output_dir / "manifests_v3").glob("m44_*.json"))
                and not list((config.output_dir / "manifests").glob("m44_*.json"))
            ),
        },
        {
            "check": "manifests_include_measurement_protocol_metadata",
            "passed": float(
                bool(manifests)
                and all(
                    isinstance(item.get("measurement_protocol"), dict)
                    for item in manifests
                )
            ),
        },
        {
            "check": "manifests_include_profile_metadata",
            "passed": float(
                bool(manifests)
                and all(
                    isinstance(item.get("profile_metadata"), dict)
                    for item in manifests
                )
            ),
        },
        {
            "check": "manifests_include_handoff_provenance_metadata",
            "passed": float(
                bool(manifests)
                and all(
                    isinstance(item.get("handoff_provenance"), dict)
                    for item in manifests
                )
            ),
        },
        {
            "check": "structured_manifests_protocol_invariant",
            "passed": float(
                bool(structured)
                and all(
                    row["profile_invariance_status"] == "protocol_invariant"
                    for row in structured
                )
            ),
        },
        {
            "check": "structured_manifests_parameter_complete",
            "passed": float(
                bool(structured)
                and all(float(row["parameter_complete"]) == 1.0 for row in structured)
            ),
        },
        {"check": "all_14_metrics_present", "passed": float(metrics_have_required)},
        {
            "check": "diagnostic_completeness_check_exists",
            "passed": float(
                data_path(
                    config.output_dir, "v4_protocol_diagnostic_completeness_check.csv"
                ).exists()
            ),
        },
        {
            "check": "no_v4_carry_forward_decision_file",
            "passed": float(
                not data_path(
                    config.output_dir,
                    "v4_protocol_cross_family_robustness_decisions.csv",
                ).exists()
            ),
        },
        {
            "check": "no_v4_stress_test_output",
            "passed": float(
                not any(
                    (config.output_dir / "data").glob("v4_*stress_test_result*.csv")
                )
            ),
        },
    ]


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(
            config.output_dir,
            "v4_protocol_manifest_generation_final_sanity.csv",
        ),
        ["check", "passed"],
    )
    print(f"Wrote v4 protocol manifest generation final sanity: {path}")


if __name__ == "__main__":
    main()
