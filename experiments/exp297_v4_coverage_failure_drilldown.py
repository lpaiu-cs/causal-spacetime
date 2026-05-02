"""Drill down v4 target and pair-node coverage failures."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import save_metric_bar
from v4_blocked_v5_experiment_helpers import (
    data_path,
    read_csv_rows,
    to_float_rows,
)

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v4_coverage_drilldown import (
    summarize_v4_coverage_failures,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 coverage drilldown.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    return summarize_v4_coverage_failures(
        to_float_rows(
            read_csv_rows(
                data_path(
                    config.output_dir,
                    "v4_protocol_manifest_family_coverage_metrics.csv",
                )
            )
        ),
        default_cross_family_robustness_criteria(),
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_csv(
        rows,
        data_path(config.output_dir, "v4_coverage_failure_drilldown.csv"),
        ["family_name", "target_coverage_fraction", "pair_node_coverage_fraction"],
    )
    figure = save_metric_bar(
        rows,
        label_key="family_name",
        value_key="pair_node_coverage_fraction",
        path=config.output_dir / "figures" / "v4_coverage_failures.png",
        ylabel="pair-node coverage",
    )
    print(f"Wrote v4 coverage drilldown: {path}")
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
