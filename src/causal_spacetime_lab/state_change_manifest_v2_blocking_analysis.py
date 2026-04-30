"""Root-cause blocking analysis for v2 carry-forward decisions."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

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
    "structural_blocking",
    "measured_blocking",
    "diagnostic_blocking",
    "control_family_blocking",
    "not_blocking",
}
STATUSES = {"pass", "fail", "not_applicable", "missing"}


@dataclass(frozen=True)
class V2BlockingCriterionRecord:
    """One v2 criterion-level blocking classification."""

    family_name: str
    family_kind: str
    criterion_name: str
    observed_value: float
    threshold_value: float
    criterion_direction: str
    blocking_type: str
    status: str
    explanation: str

    def __post_init__(self) -> None:
        if self.criterion_direction not in CRITERION_DIRECTIONS:
            allowed = ", ".join(sorted(CRITERION_DIRECTIONS))
            raise ValueError(f"criterion_direction must be one of: {allowed}")
        if self.blocking_type not in BLOCKING_TYPES:
            allowed = ", ".join(sorted(BLOCKING_TYPES))
            raise ValueError(f"blocking_type must be one of: {allowed}")
        if self.status not in STATUSES:
            allowed = ", ".join(sorted(STATUSES))
            raise ValueError(f"status must be one of: {allowed}")


@dataclass(frozen=True)
class _CriterionSpec:
    name: str
    metric_key: str
    threshold: float
    direction: str


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
        ),
        _CriterionSpec(
            "fitted_fraction",
            "fitted_fraction",
            criteria.min_fitted_fraction,
            "min_required",
        ),
        _CriterionSpec(
            "no_fit_fraction",
            "no_fit_fraction",
            criteria.max_no_fit_fraction,
            "max_allowed",
        ),
        _CriterionSpec(
            "mean_heldout_violation",
            "mean_heldout_violation",
            criteria.max_mean_heldout_violation,
            "max_allowed",
        ),
        _CriterionSpec(
            "mean_generalization_gap",
            "mean_generalization_gap",
            criteria.max_mean_generalization_gap,
            "max_allowed",
        ),
        _CriterionSpec(
            "stricter_threshold_pass_fraction",
            "stricter_threshold_pass_fraction",
            criteria.min_stricter_threshold_pass_fraction,
            "min_required",
        ),
        _CriterionSpec(
            "destructive_null_gap",
            "destructive_null_gap",
            criteria.min_destructive_null_gap,
            "min_required",
        ),
        _CriterionSpec(
            "symmetry_control_gap",
            "symmetry_control_gap",
            criteria.max_symmetry_control_gap,
            "max_allowed",
        ),
        _CriterionSpec(
            "target_coverage_fraction",
            "target_coverage_fraction",
            criteria.min_target_coverage_fraction,
            "min_required",
        ),
        _CriterionSpec(
            "pair_node_coverage_fraction",
            "pair_node_coverage_fraction",
            criteria.min_pair_node_coverage_fraction,
            "min_required",
        ),
        _CriterionSpec(
            "restart_std",
            "restart_std",
            criteria.max_restart_std,
            "max_allowed",
        ),
        _CriterionSpec(
            "latent_order_disagreement",
            "latent_order_disagreement",
            criteria.max_latent_order_disagreement,
            "max_allowed",
        ),
        _CriterionSpec(
            "no_retuning_audit_pass",
            "no_retuning_audit_pass",
            1.0,
            "boolean_required",
        ),
        _CriterionSpec(
            "failed_accounting_present",
            "failed_accounting_present",
            1.0,
            "boolean_required",
        ),
    ]


def _failing_blocking_type(
    *,
    criterion_name: str,
    family_kind: str,
    family_metric_row: dict[str, float | str],
    direction: str,
) -> str:
    if _is_control_family(family_kind) or direction == "control_only":
        return "control_family_blocking"
    if direction == "boolean_required":
        return "diagnostic_blocking"
    if criterion_name == "manifest_count":
        return "structural_blocking"
    if criterion_name in {"fitted_fraction", "no_fit_fraction"}:
        manifest_count = _to_float(family_metric_row.get("manifest_count"))
        return (
            "structural_blocking"
            if np.isfinite(manifest_count) and manifest_count < 1.0
            else "measured_blocking"
        )
    return "measured_blocking"


def classify_v2_blocking_criterion(
    family_metric_row: dict[str, float | str],
    criterion_name: str,
    observed_value: float,
    threshold_value: float,
    criterion_direction: str,
) -> V2BlockingCriterionRecord:
    """Classify one v2 criterion as structural, measured, diagnostic, or control."""

    family_name = str(family_metric_row.get("family_name", ""))
    family_kind = str(family_metric_row.get("family_kind", ""))
    observed = _to_float(observed_value)
    threshold = _to_float(threshold_value)
    if _is_control_family(family_kind) or criterion_direction == "control_only":
        return V2BlockingCriterionRecord(
            family_name=family_name,
            family_kind=family_kind,
            criterion_name=criterion_name,
            observed_value=observed,
            threshold_value=threshold,
            criterion_direction=criterion_direction,
            blocking_type="control_family_blocking",
            status="not_applicable",
            explanation="Control family is intentionally not carry-forward eligible.",
        )
    if not np.isfinite(observed):
        return V2BlockingCriterionRecord(
            family_name=family_name,
            family_kind=family_kind,
            criterion_name=criterion_name,
            observed_value=observed,
            threshold_value=threshold,
            criterion_direction=criterion_direction,
            blocking_type="diagnostic_blocking",
            status="missing",
            explanation=f"{criterion_name} is missing or non-finite.",
        )
    passed = _passes(observed, threshold, criterion_direction)
    if passed:
        return V2BlockingCriterionRecord(
            family_name=family_name,
            family_kind=family_kind,
            criterion_name=criterion_name,
            observed_value=observed,
            threshold_value=threshold,
            criterion_direction=criterion_direction,
            blocking_type="not_blocking",
            status="pass",
            explanation=f"{criterion_name} satisfies the fixed criterion.",
        )
    blocking_type = _failing_blocking_type(
        criterion_name=criterion_name,
        family_kind=family_kind,
        family_metric_row=family_metric_row,
        direction=criterion_direction,
    )
    return V2BlockingCriterionRecord(
        family_name=family_name,
        family_kind=family_kind,
        criterion_name=criterion_name,
        observed_value=observed,
        threshold_value=threshold,
        criterion_direction=criterion_direction,
        blocking_type=blocking_type,
        status="fail",
        explanation=(
            f"{criterion_name} fails the fixed {criterion_direction} criterion."
        ),
    )


def decompose_v2_blocking_by_family(
    metric_rows: list[dict[str, float | str]],
    criteria: CrossFamilyRobustnessCriteria,
) -> list[V2BlockingCriterionRecord]:
    """Evaluate all v2 blocking criteria for every family metric row."""

    records: list[V2BlockingCriterionRecord] = []
    for row in metric_rows:
        for spec in _criterion_specs(criteria):
            observed = _to_float(row.get(spec.metric_key))
            records.append(
                classify_v2_blocking_criterion(
                    row,
                    spec.name,
                    observed,
                    spec.threshold,
                    spec.direction,
                )
            )
    return records


def summarize_v2_blocking_records(
    records: list[V2BlockingCriterionRecord],
) -> list[dict[str, float | str]]:
    """Summarize blocking records by family, blocking type, and status."""

    counts: Counter[tuple[str, str, str, str]] = Counter(
        (
            record.family_name,
            record.family_kind,
            record.blocking_type,
            record.status,
        )
        for record in records
    )
    return [
        {
            "family_name": family_name,
            "family_kind": family_kind,
            "blocking_type": blocking_type,
            "status": status,
            "count": float(count),
        }
        for (family_name, family_kind, blocking_type, status), count in sorted(
            counts.items()
        )
    ]


def blocking_type_priority(blocking_type: str) -> int:
    """Return priority order for blocking categories."""

    priorities = {
        "structural_blocking": 0,
        "measured_blocking": 1,
        "diagnostic_blocking": 2,
        "control_family_blocking": 3,
        "not_blocking": 4,
    }
    return priorities.get(blocking_type, 99)
