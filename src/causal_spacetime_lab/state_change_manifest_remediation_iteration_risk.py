"""Remediation-iteration risk audit utilities."""

from __future__ import annotations

from dataclasses import asdict, dataclass

RISK_CATEGORIES = {
    "justified_protocol_design",
    "diagnostic_completeness",
    "profile_resolution",
    "stability_design",
    "coverage_design",
    "overfitting_risk",
    "threshold_tuning_risk",
    "stress_test_bypass_risk",
}
RISK_LEVELS = {"low", "medium", "high", "blocking"}


@dataclass(frozen=True)
class RemediationIterationRiskRecord:
    """One risk audit row for a proposed remediation action."""

    proposed_action: str
    target_root_cause: str
    evidence_from_m45: str
    risk_category: str
    risk_level: str
    allowed_as_v5_design: bool
    reason: str

    def __post_init__(self) -> None:
        if self.risk_category not in RISK_CATEGORIES:
            raise ValueError("risk_category is not allowed")
        if self.risk_level not in RISK_LEVELS:
            raise ValueError("risk_level is not allowed")


def _root_causes_present(
    blocking_summary_rows: list[dict[str, float | str]],
) -> set[str]:
    roots: set[str] = set()
    for row in blocking_summary_rows:
        if str(row.get("status", "")) == "fail":
            roots.add(str(row.get("root_cause_category", "")))
    return roots


def _delta_evidence(
    delta_summary_rows: list[dict[str, float | str]],
) -> str:
    improved = 0.0
    worsened = 0.0
    for row in delta_summary_rows:
        try:
            improved += float(row.get("improved_metric_count", 0.0))
            worsened += float(row.get("worsened_metric_count", 0.0))
        except (TypeError, ValueError):
            continue
    return f"v3_to_v4_improved={improved};v3_to_v4_worsened={worsened}"


def _classify_design(
    row: dict[str, float | str],
    roots: set[str],
    delta_evidence: str,
) -> RemediationIterationRiskRecord:
    action = str(row.get("family_name", row.get("proposed_action", "")))
    linked = str(row.get("linked_v4_root_causes", ""))
    rationale = str(row.get("remediation_rationale", ""))
    policy_blob = ";".join(
        str(row.get(key, ""))
        for key in [
            "profile_resolution_policy",
            "target_inclusion_policy",
            "null_taxonomy_policy",
            "stability_policy",
            "iteration_risk_category",
            "comparison_method",
            "margin_policy",
            "margin_value",
        ]
    ).lower()
    target = linked or "unspecified"
    evidence = f"{delta_evidence};linked_root_causes={linked};rationale={rationale}"
    if "threshold" in action.lower() or "threshold" in policy_blob:
        return RemediationIterationRiskRecord(
            proposed_action=action,
            target_root_cause=target,
            evidence_from_m45=evidence,
            risk_category="threshold_tuning_risk",
            risk_level="blocking",
            allowed_as_v5_design=False,
            reason="Changing thresholds is not an allowed remediation.",
        )
    if "stress" in action.lower() or "stress" in policy_blob:
        return RemediationIterationRiskRecord(
            proposed_action=action,
            target_root_cause=target,
            evidence_from_m45=evidence,
            risk_category="stress_test_bypass_risk",
            risk_level="blocking",
            allowed_as_v5_design=False,
            reason="Stress-test bypass is blocked while no family is eligible.",
        )
    if "try more families until one passes" in action.lower() or (
        "try more families until one passes" in rationale.lower()
    ):
        return RemediationIterationRiskRecord(
            proposed_action=action,
            target_root_cause=target,
            evidence_from_m45=evidence,
            risk_category="overfitting_risk",
            risk_level="blocking",
            allowed_as_v5_design=False,
            reason="Open-ended family search would be de facto tuning.",
        )
    removes_hard_cases = (
        "remove_hard_cases" in policy_blob
        or "remove hard cases" in policy_blob
        or (
            "post_hoc_hard_case_removal" in policy_blob
            and "no_post_hoc_hard_case_removal" not in policy_blob
        )
    )
    if removes_hard_cases:
        return RemediationIterationRiskRecord(
            proposed_action=action,
            target_root_cause=target,
            evidence_from_m45=evidence,
            risk_category="overfitting_risk",
            risk_level="high",
            allowed_as_v5_design=False,
            reason="Removing hard cases without protocol rationale is high risk.",
        )
    declared_category = str(row.get("iteration_risk_category", ""))
    if declared_category in RISK_CATEGORIES - {
        "overfitting_risk",
        "threshold_tuning_risk",
        "stress_test_bypass_risk",
    }:
        declared_levels = {
            "diagnostic_completeness": "low",
            "profile_resolution": "low",
            "coverage_design": "low" if "coverage_failure" in roots else "medium",
            "stability_design": "medium",
            "justified_protocol_design": "medium",
        }
        return RemediationIterationRiskRecord(
            proposed_action=action,
            target_root_cause=target,
            evidence_from_m45=evidence,
            risk_category=declared_category,
            risk_level=declared_levels.get(declared_category, "medium"),
            allowed_as_v5_design=True,
            reason="Allowed only as planned upstream remediation under fixed criteria.",
        )
    if "coverage" in policy_blob or "coverage_failure" in linked:
        category = "coverage_design"
        level = "low" if "coverage_failure" in roots else "medium"
    elif "resolution" in policy_blob or "tie" in policy_blob:
        category = "profile_resolution"
        level = "low"
    elif "stability" in policy_blob or "replicated" in action:
        category = "stability_design"
        level = "medium"
    elif "protocol" in policy_blob:
        category = "justified_protocol_design"
        level = "medium"
    else:
        category = "overfitting_risk"
        level = "medium"
    return RemediationIterationRiskRecord(
        proposed_action=action,
        target_root_cause=target,
        evidence_from_m45=evidence,
        risk_category=category,
        risk_level=level,
        allowed_as_v5_design=True,
        reason="Allowed only as planned upstream remediation under fixed criteria.",
    )


def audit_v5_remediation_iteration_risk(
    blocking_summary_rows: list[dict[str, float | str]],
    delta_summary_rows: list[dict[str, float | str]],
    proposed_v5_design_rows: list[dict[str, float | str]],
) -> list[RemediationIterationRiskRecord]:
    """Audit whether proposed v5 changes are justified or tuning risk."""

    roots = _root_causes_present(blocking_summary_rows)
    delta_evidence = _delta_evidence(delta_summary_rows)
    return [
        _classify_design(row, roots, delta_evidence)
        for row in proposed_v5_design_rows
    ]


def remediation_iteration_risk_table(
    records: list[RemediationIterationRiskRecord],
) -> list[dict[str, float | str]]:
    """Convert remediation risk records to CSV-safe rows."""

    rows: list[dict[str, float | str]] = []
    for record in records:
        row = asdict(record)
        row["allowed_as_v5_design"] = float(record.allowed_as_v5_design)
        rows.append(row)
    return rows
