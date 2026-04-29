"""Make the carry-forward stress-test stop condition explicit."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from carry_forward_failure_experiment_helpers import data_path, read_csv_rows
from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_stop_condition import (
    evaluate_stress_test_stop_condition,
    stop_condition_to_row,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for stop-condition audit."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Stress-test stop-condition audit.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Compute stop-condition row from robustness decisions."""

    decision_rows = read_csv_rows(
        data_path(config.output_dir, "cross_family_robustness_decisions.csv")
    )
    report = evaluate_stress_test_stop_condition(decision_rows)
    row = stop_condition_to_row(report)
    if not report.stress_tests_allowed:
        row["stop_condition_statement"] = (
            "future stress tests are not allowed without carry_forward or "
            "provisional families"
        )
    else:
        row["stop_condition_statement"] = report.reason
    return [row]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write stop-condition audit CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "stress_test_stop_condition_audit.csv",
        [
            "carry_forward_count",
            "provisional_count",
            "blocked_count",
            "report_only_count",
            "failed_control_count",
            "stress_tests_allowed",
            "reason",
            "allowed_mode",
            "stop_condition_statement",
        ],
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    print(f"Wrote stress-test stop-condition audit: {path}")


if __name__ == "__main__":
    main()
