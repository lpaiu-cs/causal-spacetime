"""Aggregate Milestone 34 outputs into a robustness report card."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for aggregate robustness report."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Cross-family robustness report card.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as csv_file:
        return list(csv.DictReader(csv_file))


def _plan_by_family(output_dir: Path) -> dict[str, dict[str, str]]:
    rows = _read_csv(output_dir / "data" / "stress_test_handoff_plan.csv")
    return {row["family_name"]: row for row in rows}


def _registry_allowed_use(output_dir: Path) -> dict[str, str]:
    path = output_dir / "carry_forward" / "carry_forward_registry.json"
    if not path.exists():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {
        str(record["family_name"]): str(record.get("allowed_future_use", ""))
        for record in payload.get("records", [])
    }


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Build aggregate robustness report rows."""

    decisions = _read_csv(
        config.output_dir / "data" / "cross_family_robustness_decisions.csv"
    )
    plan = _plan_by_family(config.output_dir)
    allowed_use = _registry_allowed_use(config.output_dir)
    rows: list[dict[str, float | str]] = []
    for row in decisions:
        family_name = row.get("family_name", "")
        plan_row = plan.get(family_name, {})
        rows.append(
            {
                "family_name": family_name,
                "decision": row.get("decision", ""),
                "mean_heldout_violation": row.get("mean_heldout_violation", "nan"),
                "destructive_null_gap": row.get("destructive_null_gap", "nan"),
                "symmetry_control_gap": row.get("symmetry_control_gap", "nan"),
                "stricter_threshold_pass_fraction": row.get(
                    "stricter_threshold_pass_fraction",
                    "nan",
                ),
                "failed_reasons": row.get("failed_reasons", ""),
                "allowed_future_use": allowed_use.get(family_name, ""),
                "allowed_tests": plan_row.get("allowed_tests", ""),
                "interpretation_warning": (
                    "carry_forward_is_stress_test_eligibility_not_geometry"
                ),
            }
        )
    if not rows:
        rows.append(
            {
                "family_name": "__missing_inputs__",
                "decision": "blocked",
                "mean_heldout_violation": "nan",
                "destructive_null_gap": "nan",
                "symmetry_control_gap": "nan",
                "stricter_threshold_pass_fraction": "nan",
                "failed_reasons": "missing_robustness_decisions",
                "allowed_future_use": "report_only",
                "allowed_tests": "",
                "interpretation_warning": (
                    "carry_forward_is_stress_test_eligibility_not_geometry"
                ),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write aggregate robustness report card."""

    return write_csv(
        rows,
        output_dir / "data" / "cross_family_robustness_report_card.csv",
        [
            "family_name",
            "decision",
            "mean_heldout_violation",
            "destructive_null_gap",
            "symmetry_control_gap",
            "stricter_threshold_pass_fraction",
            "failed_reasons",
            "allowed_future_use",
            "allowed_tests",
            "interpretation_warning",
        ],
    )


def main() -> None:
    config = parse_args()
    path = write_outputs(run_experiment(config), config.output_dir)
    print(f"Wrote cross-family robustness report card: {path}")


if __name__ == "__main__":
    main()
