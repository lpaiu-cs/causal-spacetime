"""Failure decomposition for v4 protocol carry-forward evaluation."""

from __future__ import annotations

from dataclasses import asdict, replace

from causal_spacetime_lab.state_change_manifest_failure_decomposition import (
    CriterionFailureRecord,
    decompose_family_metric_failures,
    summarize_failure_records,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    CrossFamilyRobustnessCriteria,
)
from causal_spacetime_lab.state_change_manifest_v4_preconditions import (
    V4ProtocolPreconditionReport,
)

_PRECONDITION_ROOTS = {
    "missing_measurement_protocol": "protocol_metadata_failure",
    "missing_profile_metadata": "profile_metadata_failure",
    "missing_handoff_provenance": "provenance_failure",
    "protocol_mixed": "protocol_invariance_failure",
    "underspecified_profile": "protocol_invariance_failure",
    "parameter_incomplete": "parameter_completeness_failure",
    "invalid_provenance": "provenance_failure",
    "structured_not_admissible": "profile_metadata_failure",
    "report_only_marked_eligible": "provenance_failure",
    "failed_control_marked_eligible": "provenance_failure",
    "diagnostic_incomplete": "diagnostic_completeness_failure",
}


def decompose_v4_protocol_family_failures(
    metric_rows: list[dict[str, float | str]],
    criteria: CrossFamilyRobustnessCriteria,
    preconditions: list[V4ProtocolPreconditionReport],
) -> list[CriterionFailureRecord]:
    """Decompose metric and precondition failures for v4 protocol families."""

    records: list[CriterionFailureRecord] = []
    for row in metric_rows:
        family_kind = str(row.get("family_kind", ""))
        for record in decompose_family_metric_failures(row, criteria):
            if record.criterion_name == "symmetry_control_gap":
                record = replace(record, root_cause_category="symmetry_control_failure")
            if family_kind in {"report_only", "failed_control"}:
                record = replace(record, root_cause_category="control_family_blocking")
            records.append(record)
    for report in preconditions:
        for failed in report.failed_preconditions:
            records.append(
                CriterionFailureRecord(
                    family_name=report.family_name,
                    family_kind=report.family_kind,
                    criterion_name=failed,
                    criterion_type="accounting",
                    observed_value=0.0,
                    threshold_value=1.0,
                    status="measured_failure",
                    root_cause_category=_PRECONDITION_ROOTS.get(failed, "unknown"),
                    explanation=(
                        f"{failed} blocks production carry-forward eligibility."
                    ),
                )
            )
    return records


def v4_protocol_failure_record_rows(
    records: list[CriterionFailureRecord],
) -> list[dict[str, float | str]]:
    """Convert v4 protocol failure records to CSV rows."""

    return [asdict(record) for record in records]


def summarize_v4_protocol_failure_records(
    records: list[CriterionFailureRecord],
) -> list[dict[str, float | str]]:
    """Summarize v4 protocol failure records."""

    return summarize_failure_records(records)


def v4_protocol_missing_or_failed_precondition_rows(
    preconditions: list[V4ProtocolPreconditionReport],
) -> list[dict[str, float | str]]:
    """Return explicit rows for v4 failed preconditions."""

    rows: list[dict[str, float | str]] = []
    for report in preconditions:
        if not report.failed_preconditions:
            rows.append(
                {
                    "family_name": report.family_name,
                    "family_kind": report.family_kind,
                    "precondition": "none",
                    "status": "pass",
                }
            )
        for failed in report.failed_preconditions:
            rows.append(
                {
                    "family_name": report.family_name,
                    "family_kind": report.family_kind,
                    "precondition": failed,
                    "status": "failed",
                }
            )
    return rows
