"""Drill down v4 destructive-null and symmetry-control weaknesses."""

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
from causal_spacetime_lab.state_change_manifest_v4_null_drilldown import (
    summarize_v4_null_taxonomy_failures,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V4 null drilldown.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    return summarize_v4_null_taxonomy_failures(
        to_float_rows(
            read_csv_rows(
                data_path(
                    config.output_dir,
                    "v4_protocol_manifest_family_null_taxonomy.csv",
                )
            )
        ),
        to_float_rows(
            read_csv_rows(
                data_path(
                    config.output_dir,
                    "v4_protocol_cross_family_robustness_metrics.csv",
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
        data_path(config.output_dir, "v4_null_failure_drilldown.csv"),
        ["family_name", "destructive_null_gap", "symmetry_control_gap"],
    )
    figure = save_metric_bar(
        rows,
        label_key="family_name",
        value_key="destructive_null_gap",
        path=config.output_dir / "figures" / "v4_null_failure_gaps.png",
        ylabel="destructive null gap",
    )
    print(f"Wrote v4 null drilldown: {path}")
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
