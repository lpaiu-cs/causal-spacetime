"""Coverage metrics for v4 protocol manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_v4_diagnostics import (
    compute_v4_protocol_coverage_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    manifest_dir: Path = Path("outputs/manifests_v4")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 coverage metrics.")
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v4"),
    )
    return ExperimentConfig(parser.parse_args().manifest_dir)


def main() -> None:
    config = parse_args()
    path = write_csv(
        compute_v4_protocol_coverage_rows(config.manifest_dir),
        Path("outputs/data/v4_protocol_manifest_family_coverage_metrics.csv"),
        ["manifest_id", "family_name", "family_kind"],
    )
    print(f"Wrote v4 protocol coverage metrics: {path}")


if __name__ == "__main__":
    main()
