"""Final sanity checks for M41 patched v3 manifest generation."""

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
from causal_spacetime_lab.state_change_manifest_v3_protocol_diagnostics import (
    manifest_metadata_audit_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol final sanity.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    manifest_paths = sorted((config.output_dir / "manifests_v3").glob("*.json"))
    manifests = [
        json.loads(path.read_text(encoding="utf-8")) for path in manifest_paths
    ]
    audit_rows = manifest_metadata_audit_rows(config.output_dir / "manifests_v3")
    structured = [row for row in audit_rows if row.get("family_kind") == "structured"]
    metrics = read_csv_rows(
        data_path(config.output_dir, "v3_protocol_cross_family_robustness_metrics.csv")
    )
    required = {item.metric_name for item in default_diagnostic_metric_requirements()}
    metrics_have_required = bool(metrics) and all(
        required.issubset(set(row)) for row in metrics
    )
    return [
        {"check": "v3_manifests_exist", "passed": float(bool(manifest_paths))},
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
            "check": "manifests_include_response_profile_metadata",
            "passed": float(
                bool(manifests)
                and all(
                    isinstance(item.get("profile_metadata"), dict) for item in manifests
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
        {
            "check": "structured_manifests_valid_provenance",
            "passed": float(
                bool(structured)
                and all(float(row["valid_provenance"]) == 1.0 for row in structured)
            ),
        },
        {"check": "all_14_metrics_present", "passed": float(metrics_have_required)},
        {
            "check": "diagnostic_completeness_check_exists",
            "passed": float(
                data_path(
                    config.output_dir, "v3_protocol_diagnostic_completeness_check.csv"
                ).exists()
            ),
        },
        {
            "check": "no_v3_carry_forward_decision_file",
            "passed": float(
                not data_path(
                    config.output_dir, "v3_cross_family_robustness_decisions.csv"
                ).exists()
            ),
        },
        {
            "check": "no_v3_stress_test_output",
            "passed": float(
                not bool(list((config.output_dir / "data").glob("v3_*stress*.csv")))
            ),
        },
    ]


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(
            config.output_dir, "v3_protocol_manifest_generation_final_sanity.csv"
        ),
        ["check", "passed"],
    )
    print(f"Wrote v3 protocol manifest generation final sanity: {path}")


if __name__ == "__main__":
    main()
