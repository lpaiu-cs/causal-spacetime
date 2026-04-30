"""Criterion-margin utilities for v2 blocked-decision audits."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np

from causal_spacetime_lab.state_change_manifest_v2_blocking_analysis import (
    V2BlockingCriterionRecord,
)


@dataclass(frozen=True)
class CriterionMarginRecord:
    """Margin between observed v2 criterion value and fixed threshold."""

    family_name: str
    criterion_name: str
    observed_value: float
    threshold_value: float
    margin: float
    normalized_margin: float
    passed: bool
    blocking_type: str


def _margin(record: V2BlockingCriterionRecord) -> float:
    if record.status == "not_applicable":
        return float("nan")
    if record.criterion_direction == "min_required":
        return record.observed_value - record.threshold_value
    if record.criterion_direction == "max_allowed":
        return record.threshold_value - record.observed_value
    if record.criterion_direction == "boolean_required":
        return 1.0 if record.status == "pass" else -1.0
    return float("nan")


def _normalized(margin: float, threshold: float) -> float:
    if not np.isfinite(margin):
        return float("nan")
    scale = max(abs(threshold), 1.0) if np.isfinite(threshold) else 1.0
    return float(margin / scale)


def criterion_margins_from_blocking_records(
    records: list[V2BlockingCriterionRecord],
) -> list[CriterionMarginRecord]:
    """Convert blocking records into signed criterion margins."""

    margins: list[CriterionMarginRecord] = []
    for record in records:
        margin = _margin(record)
        margins.append(
            CriterionMarginRecord(
                family_name=record.family_name,
                criterion_name=record.criterion_name,
                observed_value=record.observed_value,
                threshold_value=record.threshold_value,
                margin=margin,
                normalized_margin=_normalized(margin, record.threshold_value),
                passed=record.status in {"pass", "not_applicable"},
                blocking_type=record.blocking_type,
            )
        )
    return margins


def summarize_criterion_margins(
    margins: list[CriterionMarginRecord],
) -> list[dict[str, float | str]]:
    """Summarize worst criterion margins for each family."""

    grouped: dict[str, list[CriterionMarginRecord]] = defaultdict(list)
    for margin in margins:
        grouped[margin.family_name].append(margin)

    rows: list[dict[str, float | str]] = []
    for family_name, family_margins in sorted(grouped.items()):
        finite = [
            margin for margin in family_margins if np.isfinite(margin.normalized_margin)
        ]
        worst = min(finite, key=lambda item: item.normalized_margin) if finite else None
        rows.append(
            {
                "family_name": family_name,
                "worst_margin": float(worst.margin) if worst else float("nan"),
                "worst_normalized_margin": (
                    float(worst.normalized_margin) if worst else float("nan")
                ),
                "worst_criterion": worst.criterion_name if worst else "",
                "failed_criterion_count": float(
                    sum(not margin.passed for margin in family_margins)
                ),
                "structural_failure_count": float(
                    sum(
                        margin.blocking_type == "structural_blocking"
                        and not margin.passed
                        for margin in family_margins
                    )
                ),
                "measured_failure_count": float(
                    sum(
                        margin.blocking_type == "measured_blocking"
                        and not margin.passed
                        for margin in family_margins
                    )
                ),
                "diagnostic_failure_count": float(
                    sum(
                        margin.blocking_type == "diagnostic_blocking"
                        and not margin.passed
                        for margin in family_margins
                    )
                ),
            }
        )
    return rows
