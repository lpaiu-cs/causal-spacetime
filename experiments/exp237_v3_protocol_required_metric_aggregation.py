"""Aggregate required diagnostic-complete metrics for patched v3 families."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import save_family_metric_figure, write_csv
from v2_carry_forward_experiment_helpers import data_path, read_csv_rows

from causal_spacetime_lab.state_change_manifest_v3_protocol_metric_aggregation import (
    aggregate_v3_protocol_required_metrics,
    v3_protocol_diagnostic_completeness_check,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="Aggregate v3 protocol metrics.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    return ExperimentConfig(args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    data = config.output_dir / "data"
    fit_rows = read_csv_rows(data / "v3_protocol_manifest_family_fit_summary.csv")
    null_rows = read_csv_rows(data / "v3_protocol_manifest_family_null_taxonomy.csv")
    stricter_rows = read_csv_rows(
        data / "v3_protocol_manifest_family_stricter_criteria.csv"
    )
    failed_rows = read_csv_rows(
        data / "v3_protocol_manifest_family_failed_accounting.csv"
    )
    coverage_rows = read_csv_rows(
        data / "v3_protocol_manifest_family_coverage_metrics.csv"
    )
    restart_rows = read_csv_rows(
        data / "v3_protocol_manifest_family_restart_stability.csv"
    )
    latent_rows = read_csv_rows(
        data / "v3_protocol_manifest_family_latent_order_stability.csv"
    )
    audit_rows = read_csv_rows(data / "v3_protocol_no_retuning_audit.csv")
    audit_pass = bool(audit_rows) and all(
        float(row.get("passed", 0.0)) == 1.0 for row in audit_rows
    )
    metrics = aggregate_v3_protocol_required_metrics(
        fit_rows=fit_rows,
        null_rows=null_rows,
        stricter_rows=stricter_rows,
        failed_rows=failed_rows,
        coverage_rows=coverage_rows,
        restart_rows=restart_rows,
        latent_order_rows=latent_rows,
        no_retuning_audit_pass=audit_pass,
    )
    completeness = v3_protocol_diagnostic_completeness_check(metrics)
    return metrics, completeness


def main() -> None:
    config = parse_args()
    metrics, completeness = run_experiment(config)
    metrics_path = write_csv(
        metrics,
        data_path(config.output_dir, "v3_protocol_cross_family_robustness_metrics.csv"),
        ["family_name", "family_kind"],
    )
    completeness_path = write_csv(
        completeness,
        data_path(config.output_dir, "v3_protocol_diagnostic_completeness_check.csv"),
        ["family_name", "diagnostic_complete", "missing_metrics"],
    )
    save_family_metric_figure(
        completeness,
        metric="completeness_fraction",
        path=Path("outputs/figures/v3_protocol_diagnostic_completeness.png"),
        ylabel="Diagnostic completeness",
    )
    print(f"Wrote v3 protocol metric aggregation: {metrics_path}, {completeness_path}")


if __name__ == "__main__":
    main()
