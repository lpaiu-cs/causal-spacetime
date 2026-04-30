"""Apply fixed carry-forward criteria to v2 family metrics."""

from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    CrossFamilyRobustnessCriteria,
    FamilyRobustnessDecision,
    decide_family_robustness,
)


def _to_float_or_string(value: str) -> float | str:
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def v2_metrics_rows_from_bundle(
    bundle: dict[str, list[dict[str, str]]],
) -> list[dict[str, float | str]]:
    """Return v2 family metric rows from the diagnostic-complete bundle."""

    rows: list[dict[str, float | str]] = []
    for row in bundle.get("metrics", []):
        rows.append({key: _to_float_or_string(value) for key, value in row.items()})
    return rows


def decide_v2_family_robustness(
    metric_rows: list[dict[str, float | str]],
    criteria: CrossFamilyRobustnessCriteria,
) -> list[FamilyRobustnessDecision]:
    """Apply fixed M34 cross-family criteria to v2 metric rows."""

    return [decide_family_robustness(row, criteria) for row in metric_rows]


def v2_diagnostic_completeness_by_family(
    bundle: dict[str, list[dict[str, str]]],
) -> dict[str, bool]:
    """Return diagnostic-complete status by v2 family."""

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


def v2_decision_to_row(
    decision: FamilyRobustnessDecision,
    *,
    missing_inputs: list[str] | None = None,
    diagnostic_complete: bool | None = None,
) -> dict[str, float | str]:
    """Convert a v2 family decision to a CSV row."""

    row: dict[str, float | str] = {
        "family_name": decision.family_name,
        "family_kind": decision.family_kind,
        "decision": decision.decision,
        "passed": float(decision.passed),
        "failed_reasons": ";".join(decision.failed_reasons),
        "warning_reasons": ";".join(decision.warning_reasons),
        "missing_inputs": ";".join(missing_inputs or []),
        "diagnostic_complete": (
            float(diagnostic_complete) if diagnostic_complete is not None else 0.0
        ),
    }
    for key, value in decision.key_metrics.items():
        row[key] = float(value)
    return row
