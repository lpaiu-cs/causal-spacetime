"""Failed-control and report-only accounting for v4 manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import save_family_metric_figure, write_csv

from causal_spacetime_lab.state_change_manifest_v4_diagnostics import (
    compute_v4_protocol_failed_accounting_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    manifest_dir: Path = Path("outputs/manifests_v4")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 failed accounting.")
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path("outputs/manifests_v4"),
    )
    return ExperimentConfig(parser.parse_args().manifest_dir)


def main() -> None:
    config = parse_args()
    rows = compute_v4_protocol_failed_accounting_rows(config.manifest_dir)
    path = write_csv(
        rows,
        Path("outputs/data/v4_protocol_manifest_family_failed_accounting.csv"),
        ["manifest_id", "family_name", "family_kind"],
    )
    save_family_metric_figure(
        rows,
        metric="eligible",
        path=Path("outputs/figures/v4_protocol_manifest_family_failed_reason_counts.png"),
        ylabel="Eligible flag",
    )
    print(f"Wrote v4 protocol failed accounting: {path}")


if __name__ == "__main__":
    main()
