"""Compare v3-to-v4 remediation effects."""

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

from causal_spacetime_lab.state_change_manifest_v4_delta_audit import (
    planned_v3_to_v4_family_links,
    summarize_v3_to_v4_deltas,
    v3_to_v4_delta_record_to_row,
    v3_to_v4_metric_delta_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3-to-v4 metric delta audit.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    records = v3_to_v4_metric_delta_rows(
        to_float_rows(
            read_csv_rows(
                data_path(
                    config.output_dir,
                    "v3_protocol_cross_family_robustness_metrics.csv",
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
        planned_v3_to_v4_family_links(),
    )
    return (
        [v3_to_v4_delta_record_to_row(record) for record in records],
        summarize_v3_to_v4_deltas(records),
    )


def main() -> None:
    config = parse_args()
    rows, summary = run_experiment(config)
    paths = [
        write_csv(
            rows,
            data_path(config.output_dir, "v3_to_v4_metric_delta_audit.csv"),
            ["v3_family_name", "v4_family_name", "metric_name", "delta"],
        ),
        write_csv(
            summary,
            data_path(config.output_dir, "v3_to_v4_metric_delta_summary.csv"),
            ["v3_family_name", "v4_family_name", "improved_metric_count"],
        ),
    ]
    figure = save_metric_bar(
        summary,
        label_key="v4_family_name",
        value_key="improved_metric_count",
        path=config.output_dir / "figures" / "v3_to_v4_metric_delta_counts.png",
        ylabel="improved metric count",
    )
    print("Wrote v3-to-v4 delta audit: " + ", ".join(map(str, paths)))
    print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
