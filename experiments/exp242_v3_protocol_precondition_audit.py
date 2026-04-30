"""Audit v3 protocol metadata/provenance carry-forward preconditions."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v3_protocol_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v3_protocol_outputs import (
    load_v3_protocol_output_bundle,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preconditions import (
    evaluate_v3_protocol_preconditions,
    v3_protocol_precondition_report_to_row,
    v3_protocol_precondition_summary,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")
    manifest_dir: Path = Path("outputs/manifests_v3")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol precondition audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--manifest-dir", type=Path, default=Path("outputs/manifests_v3")
    )
    args = parser.parse_args()
    return ExperimentConfig(args.output_dir, args.manifest_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    bundle = load_v3_protocol_output_bundle(config.output_dir)
    reports = evaluate_v3_protocol_preconditions(config.manifest_dir, bundle)
    rows = [v3_protocol_precondition_report_to_row(report) for report in reports]
    return rows, v3_protocol_precondition_summary(rows)


def main() -> None:
    config = parse_args()
    rows, summary = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "v3_protocol_precondition_audit.csv"),
        ["family_name", "family_kind", "preconditions_passed"],
    )
    summary_path = write_csv(
        summary,
        data_path(config.output_dir, "v3_protocol_precondition_summary.csv"),
        ["family_kind", "status", "count"],
    )
    print(f"Wrote v3 protocol precondition audit: {path}, {summary_path}")


if __name__ == "__main__":
    main()
