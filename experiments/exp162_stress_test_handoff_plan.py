"""Generate future stress-test handoff plan from carry-forward registry."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from cross_family_robustness_experiment_helpers import (
    build_registry_from_output,
    read_registry_json,
)
from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_carry_forward import (
    registry_from_jsonable,
    write_carry_forward_registry,
)
from causal_spacetime_lab.state_change_manifest_stress_handoff import (
    build_stress_test_plan,
    stress_test_plan_table,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for future stress-test handoff plan."""

    output_dir: Path = DEFAULT_OUTPUT_DIR


def parse_args() -> ExperimentConfig:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Future stress-test handoff plan.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    return ExperimentConfig(output_dir=args.output_dir)


def _load_or_build_registry(output_dir: Path):
    path = output_dir / "carry_forward" / "carry_forward_registry.json"
    if path.exists():
        return registry_from_jsonable(read_registry_json(path))
    registry = build_registry_from_output(output_dir)
    write_carry_forward_registry(registry, path)
    return registry


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    """Build stress-test handoff plan rows."""

    registry = _load_or_build_registry(config.output_dir)
    entries = build_stress_test_plan(registry)
    rows = stress_test_plan_table(entries)
    if not any(row["decision"] == "carry_forward" for row in rows):
        rows.append(
            {
                "family_name": "__stop_rule__",
                "decision": "stop",
                "allowed": 0.0,
                "allowed_tests": "",
                "required_caveats": "no_carry_forward_family",
                "stop_if_failed": 1.0,
            }
        )
    return rows


def write_outputs(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write stress-test handoff plan CSV."""

    return write_csv(
        rows,
        output_dir / "data" / "stress_test_handoff_plan.csv",
        [
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
    rows = run_experiment(config)
    path = write_outputs(rows, config.output_dir)
    print(f"Wrote stress-test handoff plan: {path}")


if __name__ == "__main__":
    main()
