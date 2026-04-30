"""Generate the v2 future stress-test handoff plan."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v2_carry_forward_experiment_helpers import data_path, read_v2_registry

from causal_spacetime_lab.state_change_manifest_stop_condition import (
    evaluate_stress_test_stop_condition,
)
from causal_spacetime_lab.state_change_manifest_v2_stress_handoff import (
    build_v2_stress_test_plan,
    v2_stress_test_plan_rows,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for v2 stress-test handoff planning."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse CLI arguments."""

    parser = argparse.ArgumentParser(description="V2 stress-test handoff plan.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Build v2 stress-test handoff plan rows."""

    registry = read_v2_registry(config.output_dir)
    if registry is None:
        return [
            {
                "row_type": "missing_registry",
                "family_name": "",
                "decision": "stop",
                "allowed": 0.0,
                "allowed_tests": "",
                "required_caveats": "missing_v2_registry",
                "stop_if_failed": 1.0,
            }
        ]
    entries = build_v2_stress_test_plan(registry)
    rows = v2_stress_test_plan_rows(entries)
    for row in rows:
        row["row_type"] = "family_plan"
    decision_rows = [
        {
            "family_name": record.family_name,
            "decision": record.decision,
        }
        for record in registry.records
    ]
    stop = evaluate_stress_test_stop_condition(decision_rows)
    rows.append(
        {
            "row_type": "stop_condition",
            "family_name": "__stop_condition__",
            "decision": stop.allowed_mode,
            "allowed": float(stop.stress_tests_allowed),
            "allowed_tests": "",
            "required_caveats": stop.reason,
            "stop_if_failed": float(not stop.stress_tests_allowed),
        }
    )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write v2 stress-test handoff plan."""

    return write_csv(
        rows,
        data_path(output_dir, "v2_stress_test_handoff_plan.csv"),
        [
            "row_type",
            "family_name",
            "decision",
            "allowed",
            "allowed_tests",
            "required_caveats",
            "stop_if_failed",
        ],
    )


def main() -> None:
    config = parse_args()
    path = write_outputs(run_experiment(config), config.output_dir)
    print(f"Wrote v2 stress-test handoff plan: {path}")


if __name__ == "__main__":
    main()
