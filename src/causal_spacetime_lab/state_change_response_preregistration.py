"""Preregistration rules for pre-embedding response-comparison experiments."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PreregistrationRule:
    """A leakage-prevention rule for future representability experiments."""

    rule_name: str
    description: str
    violation_example: str


def default_preregistration_rules() -> list[PreregistrationRule]:
    """Return default pre-embedding preregistration rules."""

    return [
        PreregistrationRule(
            "protocol_choice_before_embedding",
            "Choose pairwise response-comparison protocol before any fit.",
            "Selecting the protocol after seeing held-out performance.",
        ),
        PreregistrationRule(
            "thresholds_fixed_before_embedding",
            "Freeze margin and validation thresholds before any fit.",
            "Changing the margin threshold after null-baseline failure.",
        ),
        PreregistrationRule(
            "heldout_not_used_for_fitting",
            "Keep held-out constraints unavailable to fitting procedures.",
            "Tuning a representation after seeing held-out errors.",
        ),
        PreregistrationRule(
            "nulls_declared_before_evaluation",
            "Declare null baselines before evaluating the candidate pool.",
            "Adding easier nulls after difficult nulls fail.",
        ),
        PreregistrationRule(
            "failed_pools_not_silently_removed",
            "Report failed handoff manifests as negative results.",
            "Dropping failed pools from a summary table.",
        ),
        PreregistrationRule(
            "no_metric_interpretation_without_calibration",
            "Do not interpret response-comparison constraints metrically.",
            "Calling an eligible pool a recovered spatial geometry.",
        ),
        PreregistrationRule(
            "no_target_exclusion_after_failure",
            "Freeze target inclusion before validation.",
            "Removing difficult targets after a failed gate.",
        ),
    ]


def preregistration_rule_table() -> list[dict[str, str]]:
    """Return preregistration rules as CSV-ready rows."""

    return [
        {
            "rule_name": rule.rule_name,
            "description": rule.description,
            "violation_example": rule.violation_example,
        }
        for rule in default_preregistration_rules()
    ]
