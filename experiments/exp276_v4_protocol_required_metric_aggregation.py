"""Aggregate diagnostic-complete metrics for preregistered v4 families."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import save_family_metric_figure, write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_v4_metric_aggregation import (
    aggregate_v4_protocol_required_metrics,
    v4_protocol_diagnostic_completeness_check,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="Aggregate v4 protocol metrics.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    data = config.output_dir / "data"
    audit_rows = read_csv_rows(data / "v4_protocol_no_retuning_audit.csv")
    audit_pass = bool(audit_rows) and all(
        float(row.get("passed", 0.0)) == 1.0 for row in audit_rows
    )
    metrics = aggregate_v4_protocol_required_metrics(
        fit_rows=read_csv_rows(data / "v4_protocol_manifest_family_fit_summary.csv"),
        null_rows=read_csv_rows(data / "v4_protocol_manifest_family_null_taxonomy.csv"),
        stricter_rows=read_csv_rows(
            data / "v4_protocol_manifest_family_stricter_criteria.csv"
        ),
        failed_rows=read_csv_rows(
            data / "v4_protocol_manifest_family_failed_accounting.csv"
        ),
        coverage_rows=read_csv_rows(
            data / "v4_protocol_manifest_family_coverage_metrics.csv"
        ),
        restart_rows=read_csv_rows(
            data / "v4_protocol_manifest_family_restart_stability.csv"
        ),
        latent_order_rows=read_csv_rows(
            data / "v4_protocol_manifest_family_latent_order_stability.csv"
        ),
        no_retuning_audit_pass=audit_pass,
    )
    return metrics, v4_protocol_diagnostic_completeness_check(metrics)


def main() -> None:
    config = parse_args()
    metrics, completeness = run_experiment(config)
    metrics_path = write_csv(
        metrics,
        data_path(config.output_dir, "v4_protocol_cross_family_robustness_metrics.csv"),
        ["family_name", "family_kind"],
    )
    completeness_path = write_csv(
        completeness,
        data_path(config.output_dir, "v4_protocol_diagnostic_completeness_check.csv"),
        ["family_name", "diagnostic_complete", "missing_metrics"],
    )
    save_family_metric_figure(
        completeness,
        metric="completeness_fraction",
        path=Path("outputs/figures/v4_protocol_diagnostic_completeness.png"),
        ylabel="Diagnostic completeness",
    )
    print(f"Wrote v4 protocol metric aggregation: {metrics_path}, {completeness_path}")


if __name__ == "__main__":
    main()
