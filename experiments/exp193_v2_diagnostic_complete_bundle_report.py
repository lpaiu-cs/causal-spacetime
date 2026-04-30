"""Report whether the v2 output bundle is diagnostic-complete."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_manifest_experiment_helpers import DEFAULT_OUTPUT_DIR, data_path, read_csv_rows


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 bundle report."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 diagnostic bundle report.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Build v2 diagnostic-complete bundle report rows."""

    metrics = read_csv_rows(
        data_path(config.output_dir, "v2_cross_family_robustness_metrics.csv")
    )
    completeness = {
        row["family_name"]: row
        for row in read_csv_rows(
            data_path(config.output_dir, "v2_diagnostic_completeness_check.csv")
        )
    }
    rows: list[dict[str, float | str]] = []
    for row in metrics:
        family = row["family_name"]
        complete = completeness.get(family, {})
        rows.append(
            {
                **row,
                "completeness_fraction": float(
                    complete.get("completeness_fraction", "nan")
                ),
                "missing_metrics": complete.get("missing_metrics", ""),
                "bundle_status": (
                    "diagnostic_complete"
                    if float(complete.get("missing_metric_count", "nan")) == 0.0
                    else "diagnostic_incomplete"
                ),
                "carry_forward_status": "not carry-forward evaluated",
                "interpretation_warning": (
                    "diagnostic-complete does not mean successful"
                ),
            }
        )
    if rows:
        complete_count = sum(
            1 for row in rows if row["bundle_status"] == "diagnostic_complete"
        )
        rows.append(
            {
                "family_name": "__bundle__",
                "family_kind": "bundle_summary",
                "completeness_fraction": float(complete_count / (len(rows))),
                "missing_metrics": "",
                "bundle_status": (
                    "diagnostic_complete"
                    if complete_count == len(metrics)
                    else "diagnostic_incomplete"
                ),
                "carry_forward_status": "not carry-forward evaluated",
                "interpretation_warning": (
                    "diagnostic completeness is not carry-forward success"
                ),
            }
        )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write v2 bundle report."""

    fallback = [
        "family_name",
        "family_kind",
        "completeness_fraction",
        "missing_metrics",
        "bundle_status",
        "carry_forward_status",
        "interpretation_warning",
    ]
    return write_csv(
        rows,
        data_path(output_dir, "v2_diagnostic_complete_bundle_report.csv"),
        fallback,
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    output_path = write_outputs(rows, config.output_dir)
    print(f"Wrote v2 diagnostic-complete bundle report: {output_path}")


if __name__ == "__main__":
    main()
