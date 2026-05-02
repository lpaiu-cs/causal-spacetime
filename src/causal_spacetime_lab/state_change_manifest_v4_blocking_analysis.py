"""Blocked-decision root-cause audit utilities for v4 protocol families."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass

import numpy as np

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    CrossFamilyRobustnessCriteria,
)

CRITERION_DIRECTIONS = {
    "min_required",
    "max_allowed",
    "boolean_required",
    "control_only",
}
BLOCKING_TYPES = {
    "measured_blocking",
    "precondition_blocking",
    "control_family_blocking",
    "not_blocking",
}
STATUSES = {"pass", "fail", "not_applicable", "missing", "warning"}
ROOT_CAUSES = {
    "heldout_failure",
    "generalization_gap_failure",
    "stricter_pass_failure",
    "null_separation_failure",
    "symmetry_control_failure",
    "restart_instability",
    "latent_order_instability",
    "coverage_failure",
    "precondition_failure",
    "control_family_blocking",
    "unknown",
}


@dataclass(frozen=True)
class V4ProtocolBlockingRecord:
    """One criterion-level v4 protocol blocking record."""

    family_name: str
    family_kind: str
    criterion_name: str
    observed_value: float
    threshold_value: float
    criterion_direction: str
    blocking_type: str
    status: str
    root_cause_category: str
    explanation: str

    def __post_init__(self) -> None:
        if self.criterion_direction not in CRITERION_DIRECTIONS:
            raise ValueError("criterion_direction is not allowed")
        if self.blocking_type not in BLOCKING_TYPES:
            raise ValueError("blocking_type is not allowed")
        if self.status not in STATUSES:
            raise ValueError("status is not allowed")
        if self.root_cause_category not in ROOT_CAUSES:
            raise ValueError("root_cause_category is not allowed")


@dataclass(frozen=True)
class _CriterionSpec:
    name: str
    metric_key: str
    threshold: float
    direction: str
    root_cause: str


def _to_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _is_control_family(family_kind: str) -> bool:
    return family_kind in {"failed_control", "ineligible_control", "report_only"}


def _passes(value: float, threshold: float, direction: str) -> bool:
    if direction == "min_required":
        return value >= threshold
    if direction == "max_allowed":
        return value <= threshold
    if direction == "boolean_required":
        return value >= 1.0
    return False


def _criterion_specs(
    criteria: CrossFamilyRobustnessCriteria,
) -> list[_CriterionSpec]:
    return [
        _CriterionSpec(
            "manifest_count",
            "manifest_count",
            float(criteria.min_manifest_count),
            "min_required",
            "unknown",
        ),
        _CriterionSpec(
            "fitted_fraction",
            "fitted_fraction",
            criteria.min_fitted_fraction,
            "min_required",
            "unknown",
        ),
        _CriterionSpec(
            "no_fit_fraction",
            "no_fit_fraction",
            criteria.max_no_fit_fraction,
            "max_allowed",
            "unknown",
        ),
        _CriterionSpec(
            "mean_heldout_violation",
            "mean_heldout_violation",
            criteria.max_mean_heldout_violation,
            "max_allowed",
            "heldout_failure",
        ),
        _CriterionSpec(
            "mean_generalization_gap",
            "mean_generalization_gap",
            criteria.max_mean_generalization_gap,
            "max_allowed",
            "generalization_gap_failure",
        ),
        _CriterionSpec(
            "stricter_threshold_pass_fraction",
            "stricter_threshold_pass_fraction",
            criteria.min_stricter_threshold_pass_fraction,
            "min_required",
            "stricter_pass_failure",
        ),
        _CriterionSpec(
            "destructive_null_gap",
            "destructive_null_gap",
            criteria.min_destructive_null_gap,
            "min_required",
            "null_separation_failure",
        ),
        _CriterionSpec(
            "symmetry_control_gap",
            "symmetry_control_gap",
            criteria.max_symmetry_control_gap,
            "max_allowed",
            "symmetry_control_failure",
        ),
        _CriterionSpec(
            "target_coverage_fraction",
            "target_coverage_fraction",
            criteria.min_target_coverage_fraction,
            "min_required",
            "coverage_failure",
        ),
        _CriterionSpec(
            "pair_node_coverage_fraction",
            "pair_node_coverage_fraction",
            criteria.min_pair_node_coverage_fraction,
            "min_required",
            "coverage_failure",
        ),
        _CriterionSpec(
            "restart_std",
            "restart_std",
            criteria.max_restart_std,
            "max_allowed",
            "restart_instability",
        ),
        _CriterionSpec(
            "latent_order_disagreement",
            "latent_order_disagreement",
            criteria.max_latent_order_disagreement,
            "max_allowed",
            "latent_order_instability",
        ),
        _CriterionSpec(
            "no_retuning_audit_pass",
            "no_retuning_audit_pass",
            1.0,
            "boolean_required",
            "unknown",
        ),
        _CriterionSpec(
            "failed_accounting_present",
            "failed_accounting_present",
            1.0,
            "boolean_required",
            "unknown",
        ),
        _CriterionSpec(
            "preconditions_passed",
            "preconditions_passed",
            1.0,
            "boolean_required",
            "precondition_failure",
        ),
    ]


def _record_for_spec(
    row: dict[str, float | str],
    spec: _CriterionSpec,
) -> V4ProtocolBlockingRecord:
    family_name = str(row.get("family_name", ""))
    family_kind = str(row.get("family_kind", ""))
    observed = _to_float(row.get(spec.metric_key))
    if spec.name == "symmetry_control_gap":
        observed = abs(observed)
    if _is_control_family(family_kind):
        return V4ProtocolBlockingRecord(
            family_name=family_name,
            family_kind=family_kind,
            criterion_name=spec.name,
            observed_value=observed,
            threshold_value=spec.threshold,
            criterion_direction="control_only",
            blocking_type="control_family_blocking",
            status="not_applicable",
            root_cause_category="control_family_blocking",
            explanation="Control or report-only family is not eligible.",
        )
    if not np.isfinite(observed):
        return V4ProtocolBlockingRecord(
            family_name=family_name,
            family_kind=family_kind,
            criterion_name=spec.name,
            observed_value=observed,
            threshold_value=spec.threshold,
            criterion_direction=spec.direction,
            blocking_type="precondition_blocking",
            status="missing",
            root_cause_category="precondition_failure",
            explanation=f"{spec.name} is missing or non-finite.",
        )
    passed = _passes(observed, spec.threshold, spec.direction)
    if passed:
        return V4ProtocolBlockingRecord(
            family_name=family_name,
            family_kind=family_kind,
            criterion_name=spec.name,
            observed_value=observed,
            threshold_value=spec.threshold,
            criterion_direction=spec.direction,
            blocking_type="not_blocking",
            status="pass",
            root_cause_category=spec.root_cause,
            explanation=f"{spec.name} satisfies the fixed M34 criterion.",
        )
    blocking_type = (
        "precondition_blocking"
        if spec.root_cause == "precondition_failure"
        else "measured_blocking"
    )
    return V4ProtocolBlockingRecord(
        family_name=family_name,
        family_kind=family_kind,
        criterion_name=spec.name,
        observed_value=observed,
        threshold_value=spec.threshold,
        criterion_direction=spec.direction,
        blocking_type=blocking_type,
        status="fail",
        root_cause_category=spec.root_cause,
        explanation=f"{spec.name} fails the fixed M34 criterion.",
    )


def _merge_rows_by_family(
    *row_lists: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    merged: dict[str, dict[str, float | str]] = {}
    for rows in row_lists:
        for row in rows:
            family_name = str(row.get("family_name", ""))
            if not family_name:
                continue
            merged.setdefault(family_name, {}).update(row)
    return list(merged.values())


def decompose_v4_blocking_by_family(
    metric_rows: list[dict[str, float | str]],
    decision_rows: list[dict[str, float | str]],
    precondition_rows: list[dict[str, float | str]],
    criteria: CrossFamilyRobustnessCriteria,
) -> list[V4ProtocolBlockingRecord]:
    """Decompose fixed-criteria v4 protocol blocking by family."""

    rows = _merge_rows_by_family(metric_rows, decision_rows, precondition_rows)
    records: list[V4ProtocolBlockingRecord] = []
    for row in rows:
        for spec in _criterion_specs(criteria):
            records.append(_record_for_spec(row, spec))
    return records


def blocking_record_to_row(
    record: V4ProtocolBlockingRecord,
) -> dict[str, float | str]:
    """Convert one blocking record to a CSV row."""

    return asdict(record)


def summarize_v4_blocking_records(
    records: list[V4ProtocolBlockingRecord],
) -> list[dict[str, float | str]]:
    """Summarize v4 protocol blocking records by root cause and status."""

    counts: Counter[tuple[str, str, str, str, str]] = Counter(
        (
            record.family_name,
            record.blocking_type,
            record.root_cause_category,
            record.status,
            record.family_kind,
        )
        for record in records
    )
    return [
        {
            "family_name": family_name,
            "family_kind": family_kind,
            "blocking_type": blocking_type,
            "root_cause_category": root_cause,
            "status": status,
            "count": float(count),
        }
        for (
            family_name,
            blocking_type,
            root_cause,
            status,
            family_kind,
        ), count in sorted(counts.items())
    ]
