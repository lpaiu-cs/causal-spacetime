"""Report-only counterfactuals for v2 structural blocking."""

from __future__ import annotations

from collections import defaultdict

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    CrossFamilyRobustnessCriteria,
)
from causal_spacetime_lab.state_change_manifest_v2_blocking_analysis import (
    V2BlockingCriterionRecord,
    decompose_v2_blocking_by_family,
)


def _to_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def structural_count_counterfactual_report(
    metric_rows: list[dict[str, float | str]],
    criteria: CrossFamilyRobustnessCriteria,
    hypothetical_manifest_count: int = 3,
) -> list[dict[str, float | str]]:
    """Report manifest-count counterfactuals without changing M38 decisions.

    These counterfactuals do not change Milestone 38 decisions and do not
    justify threshold retuning.
    """

    rows: list[dict[str, float | str]] = []
    for row in metric_rows:
        family_kind = str(row.get("family_kind", ""))
        if family_kind != "structured":
            continue
        original_count = _to_float(row.get("manifest_count"))
        hypothetical_pass = hypothetical_manifest_count >= criteria.min_manifest_count
        recomputed = dict(row)
        recomputed["manifest_count"] = float(hypothetical_manifest_count)
        remaining_records = decompose_v2_blocking_by_family([recomputed], criteria)
        remaining_measured = [
            record.criterion_name
            for record in remaining_records
            if record.status == "fail" and record.blocking_type == "measured_blocking"
        ]
        rows.append(
            {
                "family_name": str(row.get("family_name", "")),
                "family_kind": family_kind,
                "original_manifest_count": original_count,
                "hypothetical_manifest_count": float(hypothetical_manifest_count),
                "manifest_count_would_pass": float(hypothetical_pass),
                "remaining_measured_blocker_count": float(len(remaining_measured)),
                "remaining_measured_blockers": ";".join(remaining_measured),
                "analysis_label": (
                    "report_only_counterfactual_not_decision_change"
                ),
            }
        )
    return rows


def would_remain_blocked_ignoring_structural_count(
    blocking_records: list[V2BlockingCriterionRecord],
) -> list[dict[str, float | str]]:
    """Report whether measured failures remain after ignoring count structure.

    These counterfactuals do not change Milestone 38 decisions and do not
    justify threshold retuning.
    """

    grouped: dict[str, list[V2BlockingCriterionRecord]] = defaultdict(list)
    family_kind_by_name: dict[str, str] = {}
    for record in blocking_records:
        grouped[record.family_name].append(record)
        family_kind_by_name[record.family_name] = record.family_kind

    rows: list[dict[str, float | str]] = []
    for family_name, records in sorted(grouped.items()):
        measured_blockers = [
            record.criterion_name
            for record in records
            if record.status == "fail" and record.blocking_type == "measured_blocking"
        ]
        structural_count_blockers = [
            record.criterion_name
            for record in records
            if record.status == "fail"
            and record.blocking_type == "structural_blocking"
            and record.criterion_name == "manifest_count"
        ]
        rows.append(
            {
                "family_name": family_name,
                "family_kind": family_kind_by_name.get(family_name, ""),
                "ignored_structural_count_failure": float(
                    bool(structural_count_blockers)
                ),
                "measured_blocker_count": float(len(measured_blockers)),
                "measured_blockers": ";".join(measured_blockers),
                "would_remain_blocked": float(bool(measured_blockers)),
                "analysis_label": (
                    "report_only_counterfactual_not_decision_change"
                ),
            }
        )
    return rows
