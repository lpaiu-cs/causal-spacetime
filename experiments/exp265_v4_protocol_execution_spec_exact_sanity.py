"""Exact sanity checks for v4 protocol execution spec loading."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_v4_execution_spec import (
    load_v4_protocol_execution_specs,
    mark_v4_protocol_execution_specs_executed,
)


def run_experiment() -> list[dict[str, float | str]]:
    path = Path("outputs/remediation/v4_protocol_preregistration_spec_m43.json")
    specs = load_v4_protocol_execution_specs(path)
    if not specs:
        return [
            {
                "check": "v4_preregistration_exists",
                "passed": 0.0,
                "detail": str(path),
            }
        ]
    executed = mark_v4_protocol_execution_specs_executed(specs)
    structured = [spec for spec in specs if spec.family_kind == "structured"]
    controls = [spec for spec in specs if spec.family_kind != "structured"]
    return [
        {"check": "structured_v4_specs_load", "passed": float(bool(structured))},
        {
            "check": "structured_planned_count_at_least_3",
            "passed": float(
                all(spec.planned_manifest_count >= 3 for spec in structured)
            ),
        },
        {"check": "failed_report_only_controls_load", "passed": float(bool(controls))},
        {
            "check": "planned_only_before_execution",
            "passed": float(
                all(spec.execution_status == "planned_only" for spec in specs)
            ),
        },
        {
            "check": "mark_executed_works",
            "passed": float(
                all(spec.execution_status == "executed" for spec in executed)
            ),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/v4_protocol_execution_spec_exact_sanity.csv"),
        ["check", "passed", "detail"],
    )
    print(f"Wrote v4 protocol execution spec sanity: {path}")


if __name__ == "__main__":
    main()
