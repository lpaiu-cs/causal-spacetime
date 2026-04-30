"""Decompose v2 carry-forward failures under fixed criteria."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows
from v2_manifest_experiment_helpers import save_metric_bar

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v2_failure_decomposition import (
    decompose_v2_family_failures,
    summarize_v2_failure_records,
    v2_failure_record_rows,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 failure decomposition."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 failure decomposition.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _metric_rows(output_dir: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for row in read_csv_rows(
        data_path(output_dir, "v2_cross_family_robustness_metrics.csv")
    ):
        rows.append({key: value for key, value in row.items()})
    return rows


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Run criterion-level v2 failure decomposition."""

    records = decompose_v2_family_failures(
        _metric_rows(config.output_dir),
        default_cross_family_robustness_criteria(),
    )
    return v2_failure_record_rows(records), summarize_v2_failure_records(records)


def write_outputs(
    rows: list[dict[str, float | str]],
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write v2 failure decomposition outputs."""

    return (
        write_csv(
            rows,
            data_path(output_dir, "v2_carry_forward_failure_decomposition.csv"),
            [
                "family_name",
                "family_kind",
                "criterion_name",
                "criterion_type",
                "observed_value",
                "threshold_value",
                "status",
                "root_cause_category",
                "explanation",
            ],
        ),
        write_csv(
            summary_rows,
            data_path(output_dir, "v2_carry_forward_failure_summary.csv"),
            [
                "family_name",
                "family_kind",
                "root_cause_category",
                "status",
                "count",
            ],
        ),
    )


def save_figures(
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save v2 failure root cause figure."""

    root_rows: dict[str, float] = {}
    for row in summary_rows:
        if row["status"] == "measured_failure":
            root = str(row["root_cause_category"])
            root_rows[root] = root_rows.get(root, 0.0) + float(row["count"])
    rows = [
        {"root_cause_category": root, "count": count}
        for root, count in sorted(root_rows.items())
    ]
    return [
        save_metric_bar(
            rows,
            label_key="root_cause_category",
            value_key="count",
            path=output_dir
            / "figures"
            / "v2_carry_forward_failure_root_causes.png",
            ylabel="measured failure count",
        )
    ]


def main() -> None:
    config = parse_args()
    rows, summary = run_experiment(config)
    paths = write_outputs(rows, summary, config.output_dir)
    figures = save_figures(summary, config.output_dir)
    print("Wrote v2 failure decomposition: " + ", ".join(str(path) for path in paths))
    for figure in figures:
        print(f"Wrote figure: {figure}")


if __name__ == "__main__":
    main()
