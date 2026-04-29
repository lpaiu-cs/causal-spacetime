"""Report how missing metrics affect carry-forward decisions without imputation."""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from carry_forward_failure_experiment_helpers import (
    data_path,
    failure_record_from_row,
    read_csv_rows,
)
from manifest_family_experiment_helpers import write_csv

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for missing-metric impact report."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Missing-metric impact report.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Summarize measured failures and missing metrics by family."""

    failure_rows = read_csv_rows(
        data_path(config.output_dir, "carry_forward_failure_decomposition.csv")
    )
    completeness_rows = {
        row.get("family_name", ""): row
        for row in read_csv_rows(
            data_path(
                config.output_dir,
                "cross_family_diagnostic_completeness_audit.csv",
            )
        )
    }
    if not failure_rows:
        return [
            {
                "family_name": "__missing_inputs__",
                "hard_measured_failures": "",
                "hard_missing_metrics": "carry_forward_failure_decomposition",
                "warning_missing_metrics": "",
                "completeness_fraction": 0.0,
                "decision_changed": 0.0,
                "would_remain_blocked_ignoring_missing_metrics": 1.0,
            }
        ]

    grouped: dict[str, list] = defaultdict(list)
    for row in failure_rows:
        record = failure_record_from_row(row)
        grouped[record.family_name].append(record)

    rows: list[dict[str, float | str]] = []
    for family_name, records in sorted(grouped.items()):
        hard_measured = [
            record.criterion_name
            for record in records
            if record.criterion_type in {"hard", "accounting"}
            and record.status == "measured_failure"
        ]
        hard_missing = [
            record.criterion_name
            for record in records
            if record.criterion_type in {"hard", "accounting"}
            and record.status == "missing_metric"
        ]
        warning_missing = [
            record.criterion_name
            for record in records
            if record.criterion_type == "warning"
            and record.status == "missing_metric"
        ]
        completeness = completeness_rows.get(family_name, {})
        rows.append(
            {
                "family_name": family_name,
                "hard_measured_failures": ";".join(hard_measured),
                "hard_missing_metrics": ";".join(hard_missing),
                "warning_missing_metrics": ";".join(warning_missing),
                "completeness_fraction": float(
                    completeness.get("completeness_fraction", "nan")
                ),
                "decision_changed": 0.0,
                "would_remain_blocked_ignoring_missing_metrics": float(
                    bool(hard_measured)
                ),
            }
        )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write missing-metric impact report."""

    return write_csv(
        rows,
        output_dir / "data" / "missing_metric_impact_report.csv",
        [
            "family_name",
            "hard_measured_failures",
            "hard_missing_metrics",
            "warning_missing_metrics",
            "completeness_fraction",
            "decision_changed",
            "would_remain_blocked_ignoring_missing_metrics",
        ],
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    print(f"Wrote missing-metric impact report: {path}")


if __name__ == "__main__":
    main()
