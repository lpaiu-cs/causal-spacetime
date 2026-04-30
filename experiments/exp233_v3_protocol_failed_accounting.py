"""Failed and report-only accounting for protocol-invariant v3 manifests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import save_family_metric_figure, write_csv

from causal_spacetime_lab.state_change_manifest_v3_protocol_diagnostics import (
    compute_v3_protocol_failed_accounting_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    manifest_dir: Path = Path("outputs/manifests_v3")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol failed accounting.")
    parser.add_argument(
        "--manifest-dir", type=Path, default=Path("outputs/manifests_v3")
    )
    args = parser.parse_args()
    return ExperimentConfig(args.manifest_dir)


def main() -> None:
    config = parse_args()
    rows = compute_v3_protocol_failed_accounting_rows(config.manifest_dir)
    path = write_csv(
        rows,
        Path("outputs/data/v3_protocol_manifest_family_failed_accounting.csv"),
        ["manifest_id", "family_name", "family_kind", "eligible"],
    )
    figure_rows = [
        {**row, "failed": 1.0 - float(row.get("eligible", 0.0))} for row in rows
    ]
    save_family_metric_figure(
        figure_rows,
        metric="failed",
        path=Path(
            "outputs/figures/v3_protocol_manifest_family_failed_reason_counts.png"
        ),
        ylabel="Failed or report-only rows",
    )
    print(f"Wrote v3 protocol failed accounting: {path}")


if __name__ == "__main__":
    main()
