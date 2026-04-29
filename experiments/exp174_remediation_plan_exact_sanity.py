"""Exact sanity checks for preregistered remediation-plan construction."""

from __future__ import annotations

import csv
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    build_preregistered_remediation_plan,
    remediation_plan_digest,
    remediation_plan_to_jsonable,
)

DEFAULT_OUTPUT_DIR = Path("outputs")


def run_experiment() -> list[dict[str, float | str]]:
    """Run deterministic remediation-plan checks."""

    failure_summary = [
        {
            "family_name": "family",
            "family_kind": "structured",
            "root_cause_category": "heldout_failure",
            "status": "measured_failure",
            "count": 1.0,
        },
        {
            "family_name": "family",
            "family_kind": "structured",
            "root_cause_category": "missing_metric",
            "status": "missing_metric",
            "count": 1.0,
        },
        {
            "family_name": "family",
            "family_kind": "structured",
            "root_cause_category": "coverage_failure",
            "status": "measured_failure",
            "count": 1.0,
        },
    ]
    completeness = [
        {
            "family_name": "family",
            "missing_metrics": "target_coverage_fraction;restart_std",
        }
    ]
    plan = build_preregistered_remediation_plan(failure_summary, completeness)
    payload = remediation_plan_to_jsonable(plan)
    digest_a = remediation_plan_digest(payload)
    digest_b = remediation_plan_digest(payload)
    action_names = {action.action_name for action in plan.actions}
    checks = [
        ("plan_execution_not_allowed", not plan.execution_allowed_in_current_milestone),
        ("forbidden_interpretations_nonempty", bool(plan.forbidden_interpretations)),
        (
            "missing_metric_action_exists",
            "add_missing_metric_collection" in action_names,
        ),
        (
            "new_manifest_families_planned_only",
            all(
                spec["execution_status"] == "planned_only"
                for spec in plan.new_manifest_family_specs
            ),
        ),
        ("digest_stable", digest_a == digest_b and bool(digest_a)),
    ]
    return [{"check": name, "passed": float(passed)} for name, passed in checks]


def write_outputs(
    rows: list[dict[str, float | str]],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> Path:
    """Write exact sanity CSV."""

    path = output_dir / "data" / "remediation_plan_exact_sanity.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["check", "passed"])
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    path = write_outputs(run_experiment())
    print(f"Wrote remediation plan exact sanity: {path}")


if __name__ == "__main__":
    main()
