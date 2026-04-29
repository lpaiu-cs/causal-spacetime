"""Map Milestone 35 failure rows to planned remediation actions."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

from remediation_plan_experiment_helpers import (
    load_m35_failure_inputs,
    write_table,
)

from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    build_remediation_actions_from_failure_summary,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for failure-to-remediation mapping."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Failure-to-remediation mapping.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Build remediation action rows from Milestone 35 outputs."""

    failure_rows, completeness_rows = load_m35_failure_inputs(config.output_dir)
    actions = build_remediation_actions_from_failure_summary(
        failure_rows,
        completeness_rows,
    )
    return [asdict(action) for action in actions]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path,
) -> Path:
    """Write remediation mapping CSV."""

    return write_table(
        rows,
        output_dir,
        "failure_to_remediation_mapping.csv",
        [
            "action_name",
            "target_root_cause",
            "action_type",
            "description",
            "expected_effect",
            "requires_new_preregistration",
            "allowed_in_current_milestone",
        ],
    )


def main() -> None:
    config = parse_args()
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    print(f"Wrote failure-to-remediation mapping: {path}")


if __name__ == "__main__":
    main()
