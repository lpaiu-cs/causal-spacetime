"""Generate the v3 protocol future stress-test handoff plan."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from manifest_family_experiment_helpers import write_csv
from v3_protocol_carry_forward_experiment_helpers import (
    data_path,
    read_v3_protocol_registry,
)

from causal_spacetime_lab.state_change_manifest_stop_condition import (
    evaluate_stress_test_stop_condition,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_stress_handoff import (
    build_v3_protocol_stress_test_plan,
    v3_protocol_stress_test_plan_rows,
)


@dataclass(frozen=True)
class ExperimentConfig:
    output_dir: Path = Path("outputs")


def parse_args() -> ExperimentConfig:
    parser = argparse.ArgumentParser(description="V3 protocol stress-test plan.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return ExperimentConfig(parser.parse_args().output_dir)


def run_experiment(config: ExperimentConfig) -> list[dict[str, float | str]]:
    registry = read_v3_protocol_registry(config.output_dir)
    if registry is None:
        return [
            {
                "row_type": "missing_registry",
                "family_name": "",
                "decision": "stop",
                "allowed": 0.0,
                "allowed_tests": "",
                "required_caveats": "missing_v3_protocol_registry",
                "stop_if_failed": 1.0,
            }
        ]
    entries = build_v3_protocol_stress_test_plan(registry)
    rows = v3_protocol_stress_test_plan_rows(entries)
    for row in rows:
        row["row_type"] = "family_plan"
    stop = evaluate_stress_test_stop_condition(
        [
            {"family_name": record.family_name, "decision": record.decision}
            for record in registry.records
        ]
    )
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


def main() -> None:
    config = parse_args()
    path = write_csv(
        run_experiment(config),
        data_path(config.output_dir, "v3_protocol_stress_test_handoff_plan.csv"),
        ["row_type", "family_name", "decision", "allowed"],
    )
    print(f"Wrote v3 protocol stress-test handoff plan: {path}")


if __name__ == "__main__":
    main()

