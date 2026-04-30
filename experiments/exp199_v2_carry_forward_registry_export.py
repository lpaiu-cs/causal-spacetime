"""Export the v2 carry-forward registry."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, decisions_from_csv

from causal_spacetime_lab.state_change_manifest_v2_registry import (
    build_v2_carry_forward_registry,
    v2_registry_summary_rows,
    write_v2_carry_forward_registry,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 registry export."""

    output_dir: Path = DEFAULT_OUTPUT_DIR
    manifest_dir: Path = Path("outputs/manifests_v2")


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="Export v2 carry-forward registry.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v2"),
    )
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir, manifest_dir=args.manifest_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[Path | None, list[dict[str, float | str]]]:
    """Build and write the v2 carry-forward registry."""

    decisions = decisions_from_csv(
        data_path(config.output_dir, "v2_cross_family_robustness_decisions.csv")
    )
    if not decisions:
        return None, [
            {
                "registry_id": "",
                "criteria_name": "default_v2_fixed_m34",
                "family_name": "",
                "family_kind": "",
                "decision": "missing_decisions",
                "manifest_count": 0.0,
                "eligible_manifest_count": 0.0,
                "ineligible_manifest_count": 0.0,
                "failed_reasons": "missing_v2_decisions",
                "warning_reasons": "",
                "allowed_future_use": "none",
            }
        ]
    registry = build_v2_carry_forward_registry(
        decisions,
        config.manifest_dir,
        criteria_name="default_v2_fixed_m34",
    )
    path = write_v2_carry_forward_registry(registry, config.output_dir)
    return path, v2_registry_summary_rows(registry)


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write v2 registry export summary rows."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_carry_forward_registry_export.csv"),
        [
            "registry_id",
            "criteria_name",
            "family_name",
            "family_kind",
            "decision",
            "manifest_count",
            "eligible_manifest_count",
            "ineligible_manifest_count",
            "failed_reasons",
            "warning_reasons",
            "allowed_future_use",
        ],
    )


def main() -> None:
    config = parse_args()
    registry_path, rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    print(f"Wrote v2 registry summary: {output_path}")
    if registry_path is not None:
        print(f"Wrote v2 registry: {registry_path}")


if __name__ == "__main__":
    main()
