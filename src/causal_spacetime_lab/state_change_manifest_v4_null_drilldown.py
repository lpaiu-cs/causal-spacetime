"""Null taxonomy drilldown utilities for v4 blocked families."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    CrossFamilyRobustnessCriteria,
)


def _to_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _mean(values: list[float]) -> float:
    finite = [value for value in values if np.isfinite(value)]
    return float(np.mean(finite)) if finite else float("nan")


def summarize_v4_null_taxonomy_failures(
    null_rows: list[dict[str, float | str]],
    metric_rows: list[dict[str, float | str]],
    criteria: CrossFamilyRobustnessCriteria,
) -> list[dict[str, float | str]]:
    """Summarize destructive-null and symmetry-control failures by family."""

    metric_by_family = {str(row.get("family_name", "")): row for row in metric_rows}
    grouped: dict[str, list[dict[str, float | str]]] = defaultdict(list)
    for row in null_rows:
        grouped[str(row.get("family_name", ""))].append(row)
    rows: list[dict[str, float | str]] = []
    for family_name, family_rows in sorted(grouped.items()):
        by_class: dict[str, list[float]] = defaultdict(list)
        structured_values: list[float] = []
        for row in family_rows:
            taxonomy = str(row.get("taxonomy_class", ""))
            by_class[taxonomy].append(_to_float(row.get("mean_heldout_violation_rate")))
            structured_values.append(
                _to_float(row.get("structured_heldout_violation_rate"))
            )
        metric_row = metric_by_family.get(family_name, {})
        structured = _to_float(metric_row.get("mean_heldout_violation"))
        if not np.isfinite(structured):
            structured = _mean(structured_values)
        destructive_mean = _mean(by_class["destructive_null"])
        symmetry_mean = _mean(by_class["symmetry_control"])
        marginal_mean = _mean(by_class["marginal_preserving_null"])
        destructive_gap = destructive_mean - structured
        symmetry_gap = symmetry_mean - structured
        reasons: list[str] = []
        if not np.isfinite(destructive_gap) or (
            destructive_gap < criteria.min_destructive_null_gap
        ):
            reasons.append("destructive_null_gap_failure")
        if np.isfinite(symmetry_gap) and (
            abs(symmetry_gap) > criteria.max_symmetry_control_gap
        ):
            reasons.append("symmetry_control_gap_failure")
        rows.append(
            {
                "family_name": family_name,
                "structured_heldout": structured,
                "destructive_null_mean_heldout": destructive_mean,
                "destructive_null_gap": destructive_gap,
                "destructive_null_gap_pass": float(
                    np.isfinite(destructive_gap)
                    and destructive_gap >= criteria.min_destructive_null_gap
                ),
                "symmetry_control_mean_heldout": symmetry_mean,
                "symmetry_control_gap": symmetry_gap,
                "symmetry_control_gap_pass": float(
                    np.isfinite(symmetry_gap)
                    and abs(symmetry_gap) <= criteria.max_symmetry_control_gap
                ),
                "marginal_preserving_null_mean_heldout": marginal_mean,
                "null_failure_reason": ";".join(reasons) or "none",
            }
        )
    return rows
