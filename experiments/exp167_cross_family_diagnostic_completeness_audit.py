"""Audit diagnostic completeness for cross-family robustness metrics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from carry_forward_failure_experiment_helpers import (
    data_path,
    read_csv_rows,
)
from manifest_family_experiment_helpers import save_bar_figure, write_csv

from causal_spacetime_lab.state_change_manifest_diagnostic_completeness import (
    diagnostic_completeness_table,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for diagnostic completeness audit."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Diagnostic completeness audit.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Compute family-level diagnostic completeness."""

    metrics_rows = read_csv_rows(
        data_path(config.output_dir, "cross_family_robustness_metrics.csv")
    )
    if not metrics_rows:
        return [
            {
                "family_name": "__missing_inputs__",
                "required_metric_count": 14.0,
                "available_metric_count": 0.0,
                "missing_metric_count": 14.0,
                "completeness_fraction": 0.0,
                "missing_metrics": "cross_family_robustness_metrics",
            }
        ]
    return diagnostic_completeness_table(metrics_rows)


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write completeness audit CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "cross_family_diagnostic_completeness_audit.csv",
        [
            "family_name",
            "required_metric_count",
            "available_metric_count",
            "missing_metric_count",
            "completeness_fraction",
            "missing_metrics",
        ],
    )


def save_figures(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save completeness fraction figure."""

    labels = [str(row["family_name"]) for row in rows]
    values = [float(row["completeness_fraction"]) for row in rows]
    return [
        save_bar_figure(
            labels,
            values,
            output_dir / "figures" / "cross_family_diagnostic_completeness.png",
            ylabel="completeness fraction",
        )
    ]


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    figure_paths = save_figures(rows, config.output_dir)
    print(f"Wrote diagnostic completeness audit: {path}")
    for figure_path in figure_paths:
        print(f"Wrote figure: {figure_path}")


if __name__ == "__main__":
    main()
