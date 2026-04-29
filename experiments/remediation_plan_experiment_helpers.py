"""Shared helpers for Milestone 36 remediation-plan experiments."""

from __future__ import annotations

import json
from pathlib import Path

from carry_forward_failure_experiment_helpers import data_path, read_csv_rows
from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    DiagnosticMetricRequirement,
)
from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    RemediationAction,
    RemediationPlan,
    build_preregistered_remediation_plan,
)


def load_m35_failure_inputs(
    output_dir: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Load Milestone 35 failure and completeness rows."""

    return (
        read_csv_rows(data_path(output_dir, "carry_forward_failure_summary.csv")),
        read_csv_rows(
            data_path(output_dir, "cross_family_diagnostic_completeness_audit.csv")
        ),
    )


def build_plan_from_outputs(output_dir: Path) -> RemediationPlan:
    """Build a remediation plan from current Milestone 35 output rows."""

    failure_rows, completeness_rows = load_m35_failure_inputs(output_dir)
    return build_preregistered_remediation_plan(failure_rows, completeness_rows)


def remediation_plan_from_jsonable(payload: dict[str, object]) -> RemediationPlan:
    """Reconstruct a remediation plan from JSON-compatible data."""

    actions = [
        RemediationAction(**item)
        for item in payload.get("actions", [])  # type: ignore[arg-type]
    ]
    requirements = [
        DiagnosticMetricRequirement(**item)
        for item in payload.get("diagnostic_requirements", [])  # type: ignore[arg-type]
    ]
    return RemediationPlan(
        plan_id=str(payload.get("plan_id", "")),
        created_by_milestone=str(payload.get("created_by_milestone", "")),
        input_failure_report_files=[
            str(item)
            for item in payload.get("input_failure_report_files", [])  # type: ignore[arg-type]
        ],
        actions=actions,
        diagnostic_requirements=requirements,
        new_manifest_family_specs=[
            dict(item)
            for item in payload.get("new_manifest_family_specs", [])  # type: ignore[arg-type]
        ],
        forbidden_interpretations=[
            str(item)
            for item in payload.get("forbidden_interpretations", [])  # type: ignore[arg-type]
        ],
        execution_allowed_in_current_milestone=bool(
            payload.get("execution_allowed_in_current_milestone", False)
        ),
    )


def read_or_build_remediation_plan(output_dir: Path) -> RemediationPlan:
    """Read the exported remediation plan or build it from Milestone 35 rows."""

    plan_path = output_dir / "remediation" / "remediation_plan_m36.json"
    if plan_path.exists():
        payload = json.loads(plan_path.read_text(encoding="utf-8"))
        return remediation_plan_from_jsonable(payload)
    return build_plan_from_outputs(output_dir)


def write_table(
    rows: list[dict[str, float | str]],
    output_dir: Path,
    filename: str,
    fallback_fields: list[str],
) -> Path:
    """Write a data table under outputs/data."""

    return write_csv(rows, output_dir / "data" / filename, fallback_fields)
