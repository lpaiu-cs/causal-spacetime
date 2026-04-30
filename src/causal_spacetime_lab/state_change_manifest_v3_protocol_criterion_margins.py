"""Criterion-margin utilities for v3 protocol blocked-decision audits."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass

import numpy as np

from causal_spacetime_lab.state_change_manifest_v3_protocol_blocking_analysis import (
    V3ProtocolBlockingRecord,
)


@dataclass(frozen=True)
class V3ProtocolCriterionMargin:
    """Signed margin between observed criterion value and fixed threshold."""

    family_name: str
    criterion_name: str
    observed_value: float
    threshold_value: float
    margin: float
    normalized_margin: float
    passed: bool
    root_cause_category: str
    blocking_type: str


def _margin(record: V3ProtocolBlockingRecord) -> float:
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


def criterion_margins_from_v3_protocol_blocking_records(
    records: list[V3ProtocolBlockingRecord],
) -> list[V3ProtocolCriterionMargin]:
    """Convert blocking records to signed criterion margins."""

    margins: list[V3ProtocolCriterionMargin] = []
    for record in records:
        margin = _margin(record)
        margins.append(
            V3ProtocolCriterionMargin(
                family_name=record.family_name,
                criterion_name=record.criterion_name,
                observed_value=record.observed_value,
                threshold_value=record.threshold_value,
                margin=margin,
                normalized_margin=_normalized(margin, record.threshold_value),
                passed=record.status in {"pass", "not_applicable"},
                root_cause_category=record.root_cause_category,
                blocking_type=record.blocking_type,
            )
        )
    return margins


def criterion_margin_to_row(
    margin: V3ProtocolCriterionMargin,
) -> dict[str, float | str]:
    """Convert one margin to a CSV row."""

    row = asdict(margin)
    row["passed"] = float(margin.passed)
    return row


def summarize_v3_protocol_criterion_margins(
    margins: list[V3ProtocolCriterionMargin],
) -> list[dict[str, float | str]]:
    """Summarize worst criterion margin and failure counts by family."""

    grouped: dict[str, list[V3ProtocolCriterionMargin]] = defaultdict(list)
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
                "worst_criterion": worst.criterion_name if worst else "",
                "worst_margin": float(worst.margin) if worst else float("nan"),
                "failed_criterion_count": float(
                    sum(not margin.passed for margin in family_margins)
                ),
                "measured_failure_count": float(
                    sum(
                        not margin.passed
                        and margin.blocking_type == "measured_blocking"
                        for margin in family_margins
                    )
                ),
                "precondition_failure_count": float(
                    sum(
                        not margin.passed
                        and margin.blocking_type == "precondition_blocking"
                        for margin in family_margins
                    )
                ),
                "control_failure_count": float(
                    sum(
                        margin.blocking_type == "control_family_blocking"
                        for margin in family_margins
                    )
                ),
            }
        )
    return rows

