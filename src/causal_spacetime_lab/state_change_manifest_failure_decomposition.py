"""Failure-mode decomposition for carry-forward robustness decisions."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass

import numpy as np

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    CrossFamilyRobustnessCriteria,
)

CRITERION_TYPES = {"hard", "non_hard", "warning", "accounting"}
STATUSES = {
    "pass",
    "measured_failure",
    "missing_metric",
    "not_applicable",
    "warning_only",
}
ROOT_CAUSE_CATEGORIES = {
    "heldout_failure",
    "generalization_gap_failure",
    "null_separation_failure",
    "coverage_failure",
    "restart_instability",
    "latent_order_instability",
    "missing_output",
    "missing_metric",
    "accounting_control",
    "ineligible_control",
    "unknown",
}


@dataclass(frozen=True)
class CriterionFailureRecord:
    """One criterion-level carry-forward failure decomposition row."""

    family_name: str
    family_kind: str
    criterion_name: str
    criterion_type: str
    observed_value: float
    threshold_value: float
    status: str
    root_cause_category: str
    explanation: str

    def __post_init__(self) -> None:
        if self.criterion_type not in CRITERION_TYPES:
            allowed = ", ".join(sorted(CRITERION_TYPES))
            raise ValueError(f"criterion_type must be one of: {allowed}")
        if self.status not in STATUSES:
            allowed = ", ".join(sorted(STATUSES))
            raise ValueError(f"status must be one of: {allowed}")
        if self.root_cause_category not in ROOT_CAUSE_CATEGORIES:
            allowed = ", ".join(sorted(ROOT_CAUSE_CATEGORIES))
            raise ValueError(f"root_cause_category must be one of: {allowed}")


@dataclass(frozen=True)
class _CriterionSpec:
    name: str
    metric_key: str
    criterion_type: str
    relation: str
    threshold: float
    root_cause: str


def _to_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _criterion_specs(
    criteria: CrossFamilyRobustnessCriteria,
) -> list[_CriterionSpec]:
    return [
        _CriterionSpec(
            "manifest_count",
            "manifest_count",
            "hard",
            "min",
            float(criteria.min_manifest_count),
            "missing_output",
        ),
        _CriterionSpec(
            "fitted_fraction",
            "fitted_fraction",
            "non_hard",
            "min",
            criteria.min_fitted_fraction,
            "unknown",
        ),
        _CriterionSpec(
            "no_fit_fraction",
            "no_fit_fraction",
            "non_hard",
            "max",
            criteria.max_no_fit_fraction,
            "unknown",
        ),
        _CriterionSpec(
            "mean_heldout_violation",
            "mean_heldout_violation",
            "hard",
            "max",
            criteria.max_mean_heldout_violation,
            "heldout_failure",
        ),
        _CriterionSpec(
            "mean_generalization_gap",
            "mean_generalization_gap",
            "non_hard",
            "max",
            criteria.max_mean_generalization_gap,
            "generalization_gap_failure",
        ),
        _CriterionSpec(
            "stricter_threshold_pass_fraction",
            "stricter_threshold_pass_fraction",
            "non_hard",
            "min",
            criteria.min_stricter_threshold_pass_fraction,
            "unknown",
        ),
        _CriterionSpec(
            "destructive_null_gap",
            "destructive_null_gap",
            "non_hard",
            "min",
            criteria.min_destructive_null_gap,
            "null_separation_failure",
        ),
        _CriterionSpec(
            "symmetry_control_gap",
            "symmetry_control_gap",
            "non_hard",
            "abs_max",
            criteria.max_symmetry_control_gap,
            "null_separation_failure",
        ),
        _CriterionSpec(
            "target_coverage_fraction",
            "target_coverage_fraction",
            "warning",
            "min",
            criteria.min_target_coverage_fraction,
            "coverage_failure",
        ),
        _CriterionSpec(
            "pair_node_coverage_fraction",
            "pair_node_coverage_fraction",
            "warning",
            "min",
            criteria.min_pair_node_coverage_fraction,
            "coverage_failure",
        ),
        _CriterionSpec(
            "restart_std",
            "restart_std",
            "warning",
            "max",
            criteria.max_restart_std,
            "restart_instability",
        ),
        _CriterionSpec(
            "latent_order_disagreement",
            "latent_order_disagreement",
            "warning",
            "max",
            criteria.max_latent_order_disagreement,
            "latent_order_instability",
        ),
        _CriterionSpec(
            "no_retuning_audit_pass",
            "no_retuning_audit_pass",
            "accounting",
            "min",
            1.0,
            "missing_output",
        ),
        _CriterionSpec(
            "failed_accounting_present",
            "failed_accounting_present",
            "accounting",
            "min",
            1.0,
            "missing_output",
        ),
    ]


def _status_for_spec(value: float, spec: _CriterionSpec) -> str:
    if not np.isfinite(value):
        return "missing_metric"
    comparable = abs(value) if spec.relation == "abs_max" else value
    if spec.relation in {"max", "abs_max"}:
        return "measured_failure" if comparable > spec.threshold else "pass"
    if spec.relation == "min":
        return "measured_failure" if comparable < spec.threshold else "pass"
    return "not_applicable"


def _explanation(status: str, spec: _CriterionSpec) -> str:
    if status == "missing_metric":
        return f"{spec.metric_key} is unavailable and is not treated as success."
    if status == "measured_failure":
        return f"{spec.metric_key} failed the fixed {spec.relation} threshold."
    if status == "pass":
        return f"{spec.metric_key} satisfied the fixed criterion."
    return f"{spec.metric_key} is not applicable to this family."


def decompose_family_metric_failures(
    family_metrics: dict[str, float | str],
    criteria: CrossFamilyRobustnessCriteria,
) -> list[CriterionFailureRecord]:
    """Decompose one family's metrics into measured failures and missing metrics."""

    family_name = str(family_metrics.get("family_name", ""))
    family_kind = str(family_metrics.get("family_kind", ""))
    records: list[CriterionFailureRecord] = []
    if family_kind in {"failed_control", "ineligible_control"}:
        records.append(
            CriterionFailureRecord(
                family_name=family_name,
                family_kind=family_kind,
                criterion_name="family_control_status",
                criterion_type="accounting",
                observed_value=float("nan"),
                threshold_value=float("nan"),
                status="warning_only",
                root_cause_category=(
                    "accounting_control"
                    if family_kind == "failed_control"
                    else "ineligible_control"
                ),
                explanation="Control family remains visible in report-only accounting.",
            )
        )
    for spec in _criterion_specs(criteria):
        observed = _to_float(family_metrics.get(spec.metric_key, float("nan")))
        status = _status_for_spec(observed, spec)
        root = "missing_metric" if status == "missing_metric" else spec.root_cause
        records.append(
            CriterionFailureRecord(
                family_name=family_name,
                family_kind=family_kind,
                criterion_name=spec.name,
                criterion_type=spec.criterion_type,
                observed_value=observed,
                threshold_value=spec.threshold,
                status=status,
                root_cause_category=root,
                explanation=_explanation(status, spec),
            )
        )
    return records


def summarize_failure_records(
    records: list[CriterionFailureRecord],
) -> list[dict[str, float | str]]:
    """Group failure records by family, root cause, and status."""

    counts: Counter[tuple[str, str, str, str]] = Counter(
        (
            record.family_name,
            record.family_kind,
            record.root_cause_category,
            record.status,
        )
        for record in records
    )
    return [
        {
            "family_name": family_name,
            "family_kind": family_kind,
            "root_cause_category": root,
            "status": status,
            "count": float(count),
        }
        for (family_name, family_kind, root, status), count in sorted(counts.items())
    ]


def family_blocking_reasons(
    records: list[CriterionFailureRecord],
) -> list[str]:
    """Return hard/accounting measured failures and hard/accounting missing metrics."""

    blocking_types = {"hard", "accounting"}
    blocking_statuses = {"measured_failure", "missing_metric"}
    return [
        record.criterion_name
        for record in records
        if record.criterion_type in blocking_types
        and record.status in blocking_statuses
    ]


def failure_record_to_dict(
    record: CriterionFailureRecord,
) -> dict[str, float | str]:
    """Convert one failure record to a CSV row."""

    return asdict(record)
