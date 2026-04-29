"""Cross-family robustness criteria for frozen-manifest family diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

DECISION_VALUES = {
    "carry_forward",
    "provisional",
    "blocked",
    "report_only",
    "failed_control",
}


@dataclass(frozen=True)
class CrossFamilyRobustnessCriteria:
    """Fixed diagnostic thresholds for carry-forward eligibility."""

    name: str
    min_manifest_count: int = 3
    min_fitted_fraction: float = 0.5
    max_no_fit_fraction: float = 0.5
    max_mean_heldout_violation: float = 0.20
    max_mean_generalization_gap: float = 0.10
    min_stricter_threshold_pass_fraction: float = 0.5
    min_destructive_null_gap: float = 0.10
    max_symmetry_control_gap: float = 0.10
    min_target_coverage_fraction: float = 0.8
    min_pair_node_coverage_fraction: float = 0.5
    max_restart_std: float = 0.10
    max_latent_order_disagreement: float = 0.30
    require_no_retuning_audit: bool = True
    require_failed_accounting: bool = True

    def __post_init__(self) -> None:
        if self.min_manifest_count < 0:
            raise ValueError("min_manifest_count must be nonnegative")


@dataclass(frozen=True)
class FamilyRobustnessDecision:
    """Carry-forward decision for one manifest family."""

    family_name: str
    family_kind: str
    decision: str
    passed: bool
    failed_reasons: list[str]
    warning_reasons: list[str]
    key_metrics: dict[str, float]

    def __post_init__(self) -> None:
        if self.decision not in DECISION_VALUES:
            allowed = ", ".join(sorted(DECISION_VALUES))
            raise ValueError(f"decision must be one of: {allowed}")


def default_cross_family_robustness_criteria() -> CrossFamilyRobustnessCriteria:
    """Return the preregistered default carry-forward criteria."""

    return CrossFamilyRobustnessCriteria(name="default_cross_family_robustness")


def _float_metric(metrics: dict[str, float | str], key: str) -> float:
    try:
        return float(metrics.get(key, float("nan")))
    except (TypeError, ValueError):
        return float("nan")


def _check_max(
    value: float,
    threshold: float,
    reason: str,
    *,
    warnings: list[str],
    missing_reasons: list[str],
) -> str | None:
    if not np.isfinite(value):
        warnings.append(f"missing_{reason}")
        missing_reasons.append(f"missing_{reason}")
        return f"missing_{reason}"
    if value > threshold:
        return reason
    return None


def _check_min(
    value: float,
    threshold: float,
    reason: str,
    *,
    warnings: list[str],
    missing_reasons: list[str],
) -> str | None:
    if not np.isfinite(value):
        warnings.append(f"missing_{reason}")
        missing_reasons.append(f"missing_{reason}")
        return f"missing_{reason}"
    if value < threshold:
        return reason
    return None


def decide_family_robustness(
    family_metrics: dict[str, float | str],
    criteria: CrossFamilyRobustnessCriteria,
) -> FamilyRobustnessDecision:
    """Apply fixed robustness criteria to one family metric row."""

    family_name = str(family_metrics.get("family_name", ""))
    family_kind = str(family_metrics.get("family_kind", ""))
    key_metrics = {
        key: _float_metric(family_metrics, key)
        for key in [
            "manifest_count",
            "fitted_fraction",
            "no_fit_fraction",
            "mean_heldout_violation",
            "mean_generalization_gap",
            "stricter_threshold_pass_fraction",
            "destructive_null_gap",
            "symmetry_control_gap",
            "target_coverage_fraction",
            "pair_node_coverage_fraction",
            "restart_std",
            "latent_order_disagreement",
            "no_retuning_audit_pass",
            "failed_accounting_present",
        ]
    }
    if family_kind == "failed_control":
        return FamilyRobustnessDecision(
            family_name=family_name,
            family_kind=family_kind,
            decision="failed_control",
            passed=False,
            failed_reasons=["failed_control_family"],
            warning_reasons=[],
            key_metrics=key_metrics,
        )
    if family_kind == "ineligible_control":
        return FamilyRobustnessDecision(
            family_name=family_name,
            family_kind=family_kind,
            decision="report_only",
            passed=False,
            failed_reasons=["ineligible_control_family"],
            warning_reasons=[],
            key_metrics=key_metrics,
        )

    hard_failures: list[str] = []
    nonhard_failures: list[str] = []
    warning_reasons: list[str] = []
    missing_reasons: list[str] = []

    hard_checks = [
        _check_min(
            key_metrics["manifest_count"],
            float(criteria.min_manifest_count),
            "low_manifest_count",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
        _check_max(
            key_metrics["mean_heldout_violation"],
            criteria.max_mean_heldout_violation,
            "high_heldout_violation",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
    ]
    if criteria.require_no_retuning_audit:
        hard_checks.append(
            _check_min(
                key_metrics["no_retuning_audit_pass"],
                1.0,
                "missing_or_failed_no_retuning_audit",
                warnings=warning_reasons,
                missing_reasons=missing_reasons,
            )
        )
    if criteria.require_failed_accounting:
        hard_checks.append(
            _check_min(
                key_metrics["failed_accounting_present"],
                1.0,
                "missing_failed_family_accounting",
                warnings=warning_reasons,
                missing_reasons=missing_reasons,
            )
        )
    hard_failures.extend(reason for reason in hard_checks if reason is not None)

    optional_checks = [
        _check_min(
            key_metrics["fitted_fraction"],
            criteria.min_fitted_fraction,
            "low_fitted_fraction",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
        _check_max(
            key_metrics["no_fit_fraction"],
            criteria.max_no_fit_fraction,
            "high_no_fit_fraction",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
        _check_max(
            key_metrics["mean_generalization_gap"],
            criteria.max_mean_generalization_gap,
            "high_generalization_gap",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
        _check_min(
            key_metrics["stricter_threshold_pass_fraction"],
            criteria.min_stricter_threshold_pass_fraction,
            "low_stricter_threshold_pass_fraction",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
        _check_min(
            key_metrics["destructive_null_gap"],
            criteria.min_destructive_null_gap,
            "low_destructive_null_gap",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
        _check_max(
            abs(key_metrics["symmetry_control_gap"]),
            criteria.max_symmetry_control_gap,
            "high_symmetry_control_gap",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
        _check_min(
            key_metrics["target_coverage_fraction"],
            criteria.min_target_coverage_fraction,
            "low_target_coverage",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
        _check_min(
            key_metrics["pair_node_coverage_fraction"],
            criteria.min_pair_node_coverage_fraction,
            "low_pair_node_coverage",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
        _check_max(
            key_metrics["restart_std"],
            criteria.max_restart_std,
            "high_restart_std",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
        _check_max(
            key_metrics["latent_order_disagreement"],
            criteria.max_latent_order_disagreement,
            "high_latent_order_disagreement",
            warnings=warning_reasons,
            missing_reasons=missing_reasons,
        ),
    ]
    nonhard_failures.extend(reason for reason in optional_checks if reason is not None)

    if hard_failures:
        decision = "blocked"
        passed = False
    elif nonhard_failures:
        decision = "provisional" if len(nonhard_failures) <= 2 else "blocked"
        passed = False
    else:
        decision = "carry_forward"
        passed = True

    failed_reasons = hard_failures + nonhard_failures
    warnings = sorted(set(warning_reasons + missing_reasons))
    return FamilyRobustnessDecision(
        family_name=family_name,
        family_kind=family_kind,
        decision=decision,
        passed=passed,
        failed_reasons=failed_reasons,
        warning_reasons=warnings,
        key_metrics=key_metrics,
    )
