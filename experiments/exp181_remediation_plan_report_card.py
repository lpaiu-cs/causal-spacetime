"""Aggregate Milestone 36 remediation-plan outputs into a report card."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from carry_forward_failure_experiment_helpers import data_path, read_csv_rows
from remediation_plan_experiment_helpers import write_table

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for remediation report card."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Remediation report card.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Build action/family rows plus a global stop-condition row."""

    actions = read_csv_rows(
        data_path(config.output_dir, "failure_to_remediation_mapping.csv")
    )
    families = read_csv_rows(
        data_path(config.output_dir, "new_manifest_family_design_v2.csv")
    )
    spec_rows = read_csv_rows(
        data_path(config.output_dir, "future_manifest_run_spec.csv")
    )
    allowed_now = spec_rows[0].get("allowed_to_execute_now", "0") if spec_rows else "0"
    rows: list[dict[str, float | str]] = []
    family_names = [row.get("family_name", "") for row in families]
    for action in actions:
        rows.append(
            {
                "row_type": "action",
                "action_name": action.get("action_name", ""),
                "target_root_cause": action.get("target_root_cause", ""),
                "action_type": action.get("action_type", ""),
                "requires_new_preregistration": action.get(
                    "requires_new_preregistration",
                    "",
                ),
                "allowed_in_current_milestone": action.get(
                    "allowed_in_current_milestone",
                    "",
                ),
                "linked_v2_family": ";".join(family_names),
                "interpretation_warning": (
                    "report-only design; future execution requires preregistration"
                ),
            }
        )
    for family in families:
        rows.append(
            {
                "row_type": "family_design",
                "action_name": "",
                "target_root_cause": "",
                "action_type": "manifest_generation_design",
                "requires_new_preregistration": 1.0,
                "allowed_in_current_milestone": 0.0,
                "linked_v2_family": family.get("family_name", ""),
                "interpretation_warning": "planned v2 family is not a current result",
            }
        )
    rows.append(
        {
            "row_type": "global_stop_condition",
            "action_name": "future_execution_boundary",
            "target_root_cause": "stop_condition",
            "action_type": "no_retuning_guard",
            "requires_new_preregistration": 1.0,
            "allowed_in_current_milestone": allowed_now,
            "linked_v2_family": "",
            "interpretation_warning": (
                "no stress tests are allowed until carry_forward or provisional "
                "families exist"
            ),
        }
    )
    return rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write remediation report card."""

    return write_table(
        rows,
        output_dir,
        "remediation_plan_report_card.csv",
        [
            "row_type",
            "action_name",
            "target_root_cause",
            "action_type",
            "requires_new_preregistration",
            "allowed_in_current_milestone",
            "linked_v2_family",
            "interpretation_warning",
        ],
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    print(f"Wrote remediation plan report card: {path}")


if __name__ == "__main__":
    main()
