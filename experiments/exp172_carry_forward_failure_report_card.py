"""Aggregate Milestone 35 outputs into a family-level report card."""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from carry_forward_failure_experiment_helpers import data_path, read_csv_rows
from manifest_family_experiment_helpers import write_csv

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for the failure report card."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Carry-forward failure report card.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _split(value: str) -> list[str]:
    return [item for item in value.split(";") if item]


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Build one family-level failure report card row per family."""

    decisions = {
        row.get("family_name", ""): row
        for row in read_csv_rows(
            data_path(config.output_dir, "cross_family_robustness_decisions.csv")
        )
    }
    impact_rows = {
        row.get("family_name", ""): row
        for row in read_csv_rows(
            data_path(config.output_dir, "missing_metric_impact_report.csv")
        )
    }
    completeness_rows = {
        row.get("family_name", ""): row
        for row in read_csv_rows(
            data_path(
                config.output_dir,
                "cross_family_diagnostic_completeness_audit.csv",
            )
        )
    }
    stop_rows = read_csv_rows(
        data_path(config.output_dir, "stress_test_stop_condition_audit.csv")
    )
    stop_mode = stop_rows[0].get("allowed_mode", "unknown") if stop_rows else "unknown"
    remediation_by_cause: dict[str, list[str]] = defaultdict(list)
    for row in read_csv_rows(
        data_path(config.output_dir, "upstream_remediation_design_table.csv")
    ):
        remediation_by_cause[row.get("target_root_cause", "")].append(
            row.get("proposal_name", "")
        )

    families = sorted(set(decisions) | set(impact_rows) | set(completeness_rows))
    if not families:
        return [
            {
                "family_name": "__missing_inputs__",
                "decision": "blocked",
                "hard_measured_failures": "",
                "hard_missing_metrics": "milestone35_inputs",
                "warning_missing_metrics": "",
                "completeness_fraction": 0.0,
                "stop_condition_status": "stop",
                "recommended_remediation_categories": "missing_output",
                "interpretation_warning": (
                    "report-only continuation; no thresholds changed"
                ),
            }
        ]

    rows: list[dict[str, float | str]] = []
    for family in families:
        decision = decisions.get(family, {})
        impact = impact_rows.get(family, {})
        completeness = completeness_rows.get(family, {})
        measured = str(impact.get("hard_measured_failures", ""))
        hard_missing = str(impact.get("hard_missing_metrics", ""))
        warning_missing = str(impact.get("warning_missing_metrics", ""))
        categories = set()
        if measured:
            categories.add("heldout_failure")
        if hard_missing or warning_missing:
            categories.add("missing_metric")
        proposals = sorted(
            {
                proposal
                for category in categories
                for proposal in remediation_by_cause.get(category, [])
                if proposal
            }
        )
        rows.append(
            {
                "family_name": family,
                "decision": decision.get("decision", ""),
                "hard_measured_failures": measured,
                "hard_missing_metrics": hard_missing,
                "warning_missing_metrics": warning_missing,
                "completeness_fraction": float(
                    completeness.get("completeness_fraction", "nan")
                ),
                "stop_condition_status": stop_mode,
                "recommended_remediation_categories": ";".join(proposals),
                "interpretation_warning": (
                    "failure decomposition is report-only and does not change decisions"
                ),
            }
        )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write Milestone 35 report card."""

    return write_csv(
        rows,
        output_dir / "data" / "carry_forward_failure_report_card.csv",
        [
            "family_name",
            "decision",
            "hard_measured_failures",
            "hard_missing_metrics",
            "warning_missing_metrics",
            "completeness_fraction",
            "stop_condition_status",
            "recommended_remediation_categories",
            "interpretation_warning",
        ],
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    print(f"Wrote carry-forward failure report card: {path}")


if __name__ == "__main__":
    main()
