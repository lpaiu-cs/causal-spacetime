"""Audit required v2 bundle inputs before carry-forward evaluation."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_v2_outputs import (
    load_v2_output_bundle,
    v2_bundle_input_report,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 bundle input audit."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="Audit v2 bundle inputs.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Run v2 input audit."""

    return v2_bundle_input_report(load_v2_output_bundle(config.output_dir))


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write v2 input audit rows."""

    return write_csv(
        rows,
        output_dir / "data" / "v2_bundle_input_audit.csv",
        ["input_name", "filename", "row_count", "present"],
    )


def main() -> None:
    config = parse_args()
    path = write_outputs(run_experiment(config), config.output_dir)
    print(f"Wrote v2 bundle input audit: {path}")


if __name__ == "__main__":
    main()
