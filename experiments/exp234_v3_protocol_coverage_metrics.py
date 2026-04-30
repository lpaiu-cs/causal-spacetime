"""Coverage metrics for protocol-invariant v3 manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_v3_protocol_diagnostics import (
    compute_v3_protocol_coverage_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    manifest_dir: Path = Path("outputs/manifests_v3")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol coverage metrics.")
    parser.add_argument(
        "--manifest-dir", type=Path, default=Path("outputs/manifests_v3")
    )
    args = parser.parse_args()
    return ExperimentConfig(args.manifest_dir)


def main() -> None:
    config = parse_args()
    rows = compute_v3_protocol_coverage_rows(config.manifest_dir)
    path = write_csv(
        rows,
        Path("outputs/data/v3_protocol_manifest_family_coverage_metrics.csv"),
        ["manifest_id", "family_name", "target_coverage_fraction"],
    )
    print(f"Wrote v3 protocol coverage metrics: {path}")


if __name__ == "__main__":
    main()
