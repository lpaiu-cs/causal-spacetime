"""Decompose carry-forward decisions into measured failures and missing metrics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from carry_forward_failure_experiment_helpers import (
    data_path,
    read_csv_rows,
    save_count_figure,
    write_failure_records,
)
from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_failure_decomposition import (
    CriterionFailureRecord,
    decompose_family_metric_failures,
    summarize_failure_records,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for failure decomposition."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Carry-forward failure decomposition.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[list[CriterionFailureRecord], list[dict[str, float | str]]]:
    """Load robustness metrics and decompose family-level criteria."""

    metrics_rows = read_csv_rows(
        data_path(config.output_dir, "cross_family_robustness_metrics.csv")
    )
    decisions_rows = read_csv_rows(
        data_path(config.output_dir, "cross_family_robustness_decisions.csv")
    )
    if not metrics_rows:
        missing = CriterionFailureRecord(
            family_name="__missing_inputs__",
            family_kind="report_only",
            criterion_name="cross_family_robustness_metrics",
            criterion_type="accounting",
            observed_value=float("nan"),
            threshold_value=1.0,
            status="missing_metric",
            root_cause_category="missing_output",
            explanation="Milestone 34 robustness metrics are unavailable.",
        )
        return [missing], summarize_failure_records([missing])
    criteria = default_cross_family_robustness_criteria()
    decision_by_family = {
        row.get("family_name", ""): row.get("decision", "")
        for row in decisions_rows
    }
    records: list[CriterionFailureRecord] = []
    for metrics in metrics_rows:
        enriched: dict[str, float | str] = dict(metrics)
        family_name = metrics.get("family_name", "")
        enriched["decision"] = decision_by_family.get(family_name, "")
        records.extend(decompose_family_metric_failures(enriched, criteria))
    return records, summarize_failure_records(records)


def write_outputs(
    records: list[CriterionFailureRecord],
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write detailed and summarized decomposition rows."""

    detail_path = write_failure_records(
        records,
        output_dir / "data" / "carry_forward_failure_decomposition.csv",
    )
    summary_path = write_csv(
        summary_rows,
        output_dir / "data" / "carry_forward_failure_summary.csv",
        ["family_name", "family_kind", "root_cause_category", "status", "count"],
    )
    return detail_path, summary_path


def save_figures(
    summary_rows: list[dict[str, float | str]],
    output_dir: Path,
) -> list[Path]:
    """Save failure-decomposition summary figures."""

    return [
        save_count_figure(
            summary_rows,
            "status",
            output_dir / "figures" / "carry_forward_failure_status_counts.png",
            ylabel="criterion group count",
        ),
        save_count_figure(
            summary_rows,
            "root_cause_category",
            output_dir / "figures" / "carry_forward_failure_root_causes.png",
            ylabel="criterion group count",
        ),
    ]


def main() -> None:
    config = parse_args()
    records, summary_rows = run_experiment(config)
    paths = write_outputs(records, summary_rows, config.output_dir)
    figure_paths = save_figures(summary_rows, config.output_dir)
    print(f"Wrote carry-forward failure decomposition: {', '.join(map(str, paths))}")
    for path in figure_paths:
        print(f"Wrote figure: {path}")


if __name__ == "__main__":
    main()
