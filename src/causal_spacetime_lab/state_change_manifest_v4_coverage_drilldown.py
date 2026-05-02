"""Coverage failure drilldown for v4 protocol families."""

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


def summarize_v4_coverage_failures(
    coverage_rows: list[dict[str, float | str]],
    criteria: CrossFamilyRobustnessCriteria,
) -> list[dict[str, float | str]]:
    """Summarize target and pair-node coverage failures by family."""

    grouped: dict[str, list[dict[str, float | str]]] = defaultdict(list)
    for row in coverage_rows:
        grouped[str(row.get("family_name", ""))].append(row)
    rows: list[dict[str, float | str]] = []
    for family_name, family_rows in sorted(grouped.items()):
        target = _mean(
            [_to_float(row.get("target_coverage_fraction")) for row in family_rows]
        )
        pair_node = _mean(
            [_to_float(row.get("pair_node_coverage_fraction")) for row in family_rows]
        )
        target_pass = (
            np.isfinite(target)
            and target >= criteria.min_target_coverage_fraction
        )
        pair_pass = (
            np.isfinite(pair_node)
            and pair_node >= criteria.min_pair_node_coverage_fraction
        )
        failures = []
        if not target_pass:
            failures.append("target_coverage_failure")
        if not pair_pass:
            failures.append("pair_node_coverage_failure")
        rows.append(
            {
                "family_name": family_name,
                "target_coverage_fraction": target,
                "target_coverage_threshold": criteria.min_target_coverage_fraction,
                "target_coverage_pass": float(target_pass),
                "pair_node_coverage_fraction": pair_node,
                "pair_node_coverage_threshold": (
                    criteria.min_pair_node_coverage_fraction
                ),
                "pair_node_coverage_pass": float(pair_pass),
                "dominant_coverage_failure": ";".join(failures) or "none",
            }
        )
    return rows
