"""Apply fixed M34 criteria to the diagnostic-complete v4 protocol bundle."""

from __future__ import annotations

from dataclasses import replace

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    CrossFamilyRobustnessCriteria,
    FamilyRobustnessDecision,
    decide_family_robustness,
)
from causal_spacetime_lab.state_change_manifest_v4_preconditions import (
    V4ProtocolPreconditionReport,
)


def _to_float_or_string(value: str) -> float | str:
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def v4_protocol_metrics_rows_from_bundle(
    bundle: dict[str, list[dict[str, str]]],
) -> list[dict[str, float | str]]:
    """Return v4 protocol metric rows from the M44 bundle."""

    return [
        {key: _to_float_or_string(value) for key, value in row.items()}
        for row in bundle.get("metrics", [])
    ]


def v4_protocol_diagnostic_completeness_by_family(
    bundle: dict[str, list[dict[str, str]]],
) -> dict[str, bool]:
    """Return diagnostic-complete status by v4 protocol family."""

    status: dict[str, bool] = {}
    for row in bundle.get("completeness", []):
        family = row.get("family_name", "")
        if family and not family.startswith("__"):
            try:
                missing = float(row.get("missing_metric_count", "nan"))
            except ValueError:
                missing = float("nan")
            status[family] = missing == 0.0
    return status


def decide_v4_protocol_family_robustness(
    metric_rows: list[dict[str, float | str]],
    criteria: CrossFamilyRobustnessCriteria,
    preconditions: list[V4ProtocolPreconditionReport],
) -> list[FamilyRobustnessDecision]:
    """Apply fixed M34 criteria with v4 provenance-aware precondition gating."""

    precondition_by_family = {item.family_name: item for item in preconditions}
    decisions: list[FamilyRobustnessDecision] = []
    for row in metric_rows:
        family_kind = str(row.get("family_kind", ""))
        if family_kind == "report_only":
            base = decide_family_robustness(
                {**row, "family_kind": "ineligible_control"}, criteria
            )
            decisions.append(
                replace(
                    base,
                    family_kind=family_kind,
                    decision="report_only",
                    failed_reasons=["report_only_family"],
                )
            )
            continue
        base = decide_family_robustness(row, criteria)
        precondition = precondition_by_family.get(base.family_name)
        if (
            base.family_kind == "structured"
            and precondition is not None
            and not precondition.preconditions_passed
        ):
            failed = sorted(
                set(base.failed_reasons + precondition.failed_preconditions)
            )
            warnings = sorted(
                set(base.warning_reasons + precondition.warning_preconditions)
            )
            decisions.append(
                replace(
                    base,
                    decision="blocked",
                    passed=False,
                    failed_reasons=failed,
                    warning_reasons=warnings,
                )
            )
        else:
            decisions.append(base)
    return decisions


def v4_protocol_decision_to_row(
    decision: FamilyRobustnessDecision,
    *,
    diagnostic_complete: bool,
    preconditions_passed: bool,
    failed_preconditions: list[str],
    missing_inputs: list[str] | None = None,
) -> dict[str, float | str]:
    """Convert one v4 protocol decision to a CSV row."""

    row: dict[str, float | str] = {
        "family_name": decision.family_name,
        "family_kind": decision.family_kind,
        "decision": decision.decision,
        "passed": float(decision.passed),
        "failed_reasons": ";".join(decision.failed_reasons),
        "warning_reasons": ";".join(decision.warning_reasons),
        "missing_inputs": ";".join(missing_inputs or []),
        "diagnostic_complete": float(diagnostic_complete),
        "preconditions_passed": float(preconditions_passed),
        "failed_preconditions": ";".join(failed_preconditions),
    }
    for key, value in decision.key_metrics.items():
        row[key] = float(value)
    return row
