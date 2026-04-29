"""Export the Milestone 36 preregistered remediation plan."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from remediation_plan_experiment_helpers import build_plan_from_outputs, write_table

from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    remediation_plan_digest,
    remediation_plan_to_jsonable,
    write_remediation_plan,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for remediation plan export."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Preregistered remediation plan export."
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(
    config: ExperimentConfig,
) -> tuple[Path, list[dict[str, float | str]]]:
    """Build and export the remediation plan JSON."""

    plan = build_plan_from_outputs(config.output_dir)
    plan_path = config.output_dir / "remediation" / "remediation_plan_m36.json"
    write_remediation_plan(plan, plan_path)
    payload = remediation_plan_to_jsonable(plan)
    rows = [
        {
            "plan_id": plan.plan_id,
            "plan_digest": remediation_plan_digest(payload),
            "action_count": float(len(plan.actions)),
            "diagnostic_requirement_count": float(len(plan.diagnostic_requirements)),
            "new_manifest_family_count": float(len(plan.new_manifest_family_specs)),
            "execution_allowed_in_current_milestone": float(
                plan.execution_allowed_in_current_milestone
            ),
            "json_path": str(plan_path),
        }
    ]
    return plan_path, rows


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write remediation plan export summary."""

    return write_table(
        rows,
        output_dir,
        "preregistered_remediation_plan_export.csv",
        [
            "plan_id",
            "plan_digest",
            "action_count",
            "diagnostic_requirement_count",
            "new_manifest_family_count",
            "execution_allowed_in_current_milestone",
            "json_path",
        ],
    )


def main() -> None:
    config = parse_args()
    plan_path, rows = run_experiment(config)
    summary_path = write_outputs(rows, config.output_dir)
    print(f"Wrote remediation plan export: {summary_path} and {plan_path}")


if __name__ == "__main__":
    main()
