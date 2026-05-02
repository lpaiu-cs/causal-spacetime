"""Drill down v4 blocked results to manifest-level variance."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar
from v4_blocked_v5_experiment_helpers import data_path

from causal_spacetime_lab.state_change_manifest_v4_manifest_drilldown import (
    load_v4_manifest_drilldown_rows,
    manifest_drilldown_row_to_csv,
    summarize_v4_manifest_drilldown_by_family,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 manifest-level drilldown.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    rows = load_v4_manifest_drilldown_rows(config.output_dir)
    return (
        [manifest_drilldown_row_to_csv(row) for row in rows],
        summarize_v4_manifest_drilldown_by_family(rows),
    )


def main() -> None:
    config = parse_args()
    rows, summary = run_experiment(config)
    paths = [
        write_csv(
            rows,
            data_path(config.output_dir, "v4_manifest_level_drilldown.csv"),
            ["manifest_id", "family_name", "heldout_violation_rate"],
        ),
        write_csv(
            summary,
            data_path(config.output_dir, "v4_manifest_level_summary.csv"),
            ["family_name", "mean_heldout", "std_heldout"],
        ),
    ]
    figure = save_metric_bar(
        summary,
        label_key="family_name",
        value_key="std_heldout",
        path=config.output_dir / "figures" / "v4_manifest_heldout_variance.png",
        ylabel="held-out std",
    )
    print("Wrote v4 manifest drilldown: " + ", ".join(map(str, paths)))
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
