"""Audit required M44 v4 protocol bundle inputs."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v4_protocol_carry_forward_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v4_outputs import (
    load_v4_protocol_output_bundle,
    v4_protocol_bundle_input_report,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 protocol bundle input audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    bundle = load_v4_protocol_output_bundle(config.output_dir)
    return v4_protocol_bundle_input_report(bundle)


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v4_protocol_bundle_input_audit.csv"),
        ["input_name", "filename", "row_count", "present"],
    )
    print(f"Wrote v4 protocol bundle input audit: {path}")


if __name__ == "__main__":
    main()
