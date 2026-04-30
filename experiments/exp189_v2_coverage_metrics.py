"""Produce v2 target and pair-node coverage metrics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import DEFAULT_OUTPUT_DIR, data_path

from causal_spacetime_lab.state_change_manifest_v2_diagnostics import (
    compute_v2_coverage_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 coverage metrics."""

    manifest_dir: Path = Path("outputs/manifests_v2")
    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 coverage metrics.")
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v2"),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(manifest_dir=args.manifest_dir, output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Compute v2 coverage rows."""

    return compute_v2_coverage_rows(config.manifest_dir)


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write v2 coverage metrics."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_manifest_family_coverage_metrics.csv"),
        [
            "manifest_id",
            "family_name",
            "family_kind",
            "target_coverage_fraction",
            "pair_node_coverage_fraction",
            "constraint_count",
        ],
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    print(f"Wrote v2 coverage metrics: {output_path}")


if __name__ == "__main__":
    main()
