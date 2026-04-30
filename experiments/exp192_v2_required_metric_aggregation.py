"""Aggregate all required diagnostic-complete v2 family metrics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import (
    DEFAULT_OUTPUT_DIR,
    REQUIRED_V2_METRICS,
    data_path,
    figure_path,
    read_csv_rows,
    save_metric_bar,
)

from causal_spacetime_lab.state_change_manifest_v2_metric_aggregation import (
    aggregate_v2_required_metrics,
    v2_diagnostic_completeness_check,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 metric aggregation."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 required metric aggregation.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _audit_pass(output_dir: Path) -> bool:
    rows = read_csv_rows(data_path(output_dir, "v2_no_retuning_audit.csv"))
    return bool(rows) and all(float(row.get("passed", 0.0)) == 1.0 for row in rows)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    """Aggregate v2 diagnostic-complete family metrics."""

    output_dir = config.output_dir
    metric_rows = aggregate_v2_required_metrics(
        fit_rows=read_csv_rows(
            data_path(output_dir, "v2_manifest_family_fit_summary.csv")
        ),
        null_rows=read_csv_rows(
            data_path(output_dir, "v2_manifest_family_null_taxonomy.csv")
        ),
        stricter_rows=read_csv_rows(
            data_path(output_dir, "v2_manifest_family_stricter_criteria.csv")
        ),
        failed_rows=read_csv_rows(
            data_path(output_dir, "v2_manifest_family_failed_accounting.csv")
        ),
        restart_rows=read_csv_rows(
            data_path(output_dir, "v2_manifest_family_restart_stability.csv")
        ),
        latent_order_rows=read_csv_rows(
            data_path(output_dir, "v2_manifest_family_latent_order_stability.csv")
        ),
        coverage_rows=read_csv_rows(
            data_path(output_dir, "v2_manifest_family_coverage_metrics.csv")
        ),
        no_retuning_audit_pass=_audit_pass(output_dir),
    )
    return metric_rows, v2_diagnostic_completeness_check(metric_rows)


def write_outputs(
    metric_rows: list[dict[str, float | str]],
    completeness_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write v2 metric aggregation outputs."""

    metric_fields = ["family_name", "family_kind", *REQUIRED_V2_METRICS]
    return (
        write_csv(
            metric_rows,
            data_path(output_dir, "v2_cross_family_robustness_metrics.csv"),
            metric_fields,
        ),
        write_csv(
            completeness_rows,
            data_path(output_dir, "v2_diagnostic_completeness_check.csv"),
            [
                "family_name",
                "family_kind",
                "required_metric_count",
                "available_metric_count",
                "missing_metric_count",
                "completeness_fraction",
                "missing_metrics",
                "not_carry_forward_evaluated",
            ],
        ),
    )


def save_figures(
    completeness_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save v2 diagnostic completeness figure."""

    return [
        save_metric_bar(
            completeness_rows,
            label_key="family_name",
            value_key="completeness_fraction",
            path=figure_path(output_dir, "v2_diagnostic_completeness.png"),
            ylabel="diagnostic completeness",
        )
    ]


def main() -> None:
    config = parse_args()
    metric_rows, completeness_rows = run_experiment(config)
    paths = write_outputs(metric_rows, completeness_rows, config.output_dir)
    figure_paths = save_figures(completeness_rows, config.output_dir)
    print("Wrote v2 metric aggregation: " + ", ".join(str(path) for path in paths))
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
