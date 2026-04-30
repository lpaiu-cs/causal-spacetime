"""V2 carry-forward failure decomposition utilities."""

from __future__ import annotations

from collections import Counter

from causal_spacetime_lab.state_change_manifest_failure_decomposition import (
    CriterionFailureRecord,
    decompose_family_metric_failures,
    failure_record_to_dict,
    summarize_failure_records,
)
from causal_spacetime_lab.state_change_manifest_family_robustness import (
    CrossFamilyRobustnessCriteria,
)


def decompose_v2_family_failures(
    metric_rows: list[dict[str, float | str]],
    criteria: CrossFamilyRobustnessCriteria,
) -> list[CriterionFailureRecord]:
    """Decompose v2 family metric failures under fixed criteria."""

    records: list[CriterionFailureRecord] = []
    for row in metric_rows:
        records.extend(decompose_family_metric_failures(row, criteria))
    return records


def summarize_v2_failure_records(
    records: list[CriterionFailureRecord],
) -> list[dict[str, float | str]]:
    """Summarize v2 failure records by family, root cause, and status."""

    return summarize_failure_records(records)


def v2_failure_record_rows(
    records: list[CriterionFailureRecord],
) -> list[dict[str, float | str]]:
    """Convert v2 failure records to rows."""

    return [failure_record_to_dict(record) for record in records]


def v2_missing_metric_impact_rows(
    records: list[CriterionFailureRecord],
) -> list[dict[str, float | str]]:
    """Separate measured hard failures from missing hard and warning metrics."""

    grouped: Counter[tuple[str, str, str]] = Counter()
    for record in records:
        category = ""
        if record.criterion_type in {"hard", "accounting"}:
            if record.status == "measured_failure":
                category = "hard_measured_failure"
            elif record.status == "missing_metric":
                category = "hard_missing_metric"
        elif record.criterion_type == "warning" and record.status == "missing_metric":
            category = "warning_missing_metric"
        if category:
            grouped[(record.family_name, record.family_kind, category)] += 1
    return [
        {
            "family_name": family_name,
            "family_kind": family_kind,
            "impact_category": category,
            "count": float(count),
        }
        for (family_name, family_kind, category), count in sorted(grouped.items())
    ]
