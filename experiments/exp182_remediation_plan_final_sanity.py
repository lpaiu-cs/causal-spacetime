"""Final exact sanity checks for Milestone 36 remediation planning."""

from __future__ import annotations

import csv
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    default_diagnostic_metric_requirements,
    missing_metric_remediation_priority,
)
from causal_spacetime_lab.state_change_manifest_future_run_spec import (
    build_future_manifest_run_spec_from_plan,
)
from causal_spacetime_lab.state_change_manifest_remediation_audit import (
    future_run_spec_audit,
    remediation_plan_execution_audit,
)
from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    build_preregistered_remediation_plan,
    default_new_manifest_family_specs_v2,
    forbidden_remediation_interpretations,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment() -> list[dict[str, float | str]]:
    """Run final deterministic remediation-plan sanity checks."""

    requirements = default_diagnostic_metric_requirements()
    plan = build_preregistered_remediation_plan([], [])
    spec = build_future_manifest_run_spec_from_plan(plan)
    plan_audit = remediation_plan_execution_audit(plan)
    spec_audit = future_run_spec_audit(spec)
    checks = [
        ("diagnostic_requirements_count_14", len(requirements) == 14),
        (
            "hard_metric_priority_high",
            missing_metric_remediation_priority("mean_heldout_violation") == "high",
        ),
        (
            "v2_families_planned_only",
            all(
                spec_row["execution_status"] == "planned_only"
                for spec_row in default_new_manifest_family_specs_v2()
            ),
        ),
        ("report_only_plan_audit_passes", plan_audit["passed"] == 1.0),
        ("future_spec_disallows_execution", spec_audit["passed"] == 1.0),
        (
            "forbidden_interpretations_nonempty",
            bool(forbidden_remediation_interpretations()),
        ),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write final sanity CSV."""

    path = output_dir / "data" / "remediation_plan_final_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["check", "passed"])
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    path = write_outputs(run_experiment())
    print(f"Wrote remediation plan final sanity: {path}")


if __name__ == "__main__":
    main()
