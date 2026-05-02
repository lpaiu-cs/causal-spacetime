"""V3-to-v4 remediation delta audit utilities."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class V3ToV4DeltaRecord:
    """One metric delta between a linked v3 and v4 family."""

    metric_name: str
    v3_family_name: str
    v4_family_name: str
    v3_value: float
    v4_value: float
    delta: float
    improved: bool
    worsened: bool
    interpretation: str


LOWER_IS_BETTER = {
    "mean_heldout_violation",
    "mean_generalization_gap",
    "symmetry_control_gap",
    "restart_std",
    "latent_order_disagreement",
}
HIGHER_IS_BETTER = {
    "stricter_threshold_pass_fraction",
    "destructive_null_gap",
    "target_coverage_fraction",
    "pair_node_coverage_fraction",
}
DELTA_METRICS = tuple(sorted(LOWER_IS_BETTER | HIGHER_IS_BETTER))


def planned_v3_to_v4_family_links() -> list[dict[str, str]]:
    """Return preregistered v3-to-v4 family comparison links."""

    return [
        {
            "v3_family_name": "rank_gap_earliest_full_reference_v3",
            "v4_family_name": "rank_gap_earliest_full_stability_v4",
        },
        {
            "v3_family_name": "rank_gap_earliest_retained_reference_v3",
            "v4_family_name": "rank_gap_earliest_retained_resolution_v4",
        },
        {
            "v3_family_name": "rank_gap_gated_full_reference_v3",
            "v4_family_name": "rank_gap_gated_full_stability_v4",
        },
        {
            "v3_family_name": "combined_earliest_full_reference_v3",
            "v4_family_name": "combined_earliest_full_stability_v4",
        },
        {
            "v3_family_name": "rank_gap_earliest_full_reference_v3",
            "v4_family_name": "rank_gap_high_resolution_reference_set_v4",
        },
        {
            "v3_family_name": "failed_controls_v3",
            "v4_family_name": "failed_controls_v4",
        },
        {
            "v3_family_name": "rank_gap_immediate_edge_report_only_v3",
            "v4_family_name": "report_only_immediate_edge_v4",
        },
    ]


def _to_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _row_by_family(
    rows: list[dict[str, float | str]],
) -> dict[str, dict[str, float | str]]:
    return {str(row.get("family_name", "")): row for row in rows}


def _improvement(metric: str, v3_value: float, v4_value: float) -> tuple[bool, bool]:
    if not np.isfinite(v3_value) or not np.isfinite(v4_value):
        return False, False
    if metric in LOWER_IS_BETTER:
        return v4_value < v3_value, v4_value > v3_value
    if metric in HIGHER_IS_BETTER:
        return v4_value > v3_value, v4_value < v3_value
    return False, False


def v3_to_v4_metric_delta_rows(
    v3_metric_rows: list[dict[str, float | str]],
    v4_metric_rows: list[dict[str, float | str]],
    links: list[dict[str, str]],
) -> list[V3ToV4DeltaRecord]:
    """Compare linked v3 and v4 family metrics with fixed directions."""

    v3_by_family = _row_by_family(v3_metric_rows)
    v4_by_family = _row_by_family(v4_metric_rows)
    records: list[V3ToV4DeltaRecord] = []
    for link in links:
        v3_family = link["v3_family_name"]
        v4_family = link["v4_family_name"]
        v3_row = v3_by_family.get(v3_family, {})
        v4_row = v4_by_family.get(v4_family, {})
        for metric in DELTA_METRICS:
            v3_value = _to_float(v3_row.get(metric))
            v4_value = _to_float(v4_row.get(metric))
            delta = v4_value - v3_value
            improved, worsened = _improvement(metric, v3_value, v4_value)
            if improved:
                interpretation = "improved"
            elif worsened:
                interpretation = "worsened"
            elif np.isfinite(delta):
                interpretation = "unchanged"
            else:
                interpretation = "missing_input"
            records.append(
                V3ToV4DeltaRecord(
                    metric_name=metric,
                    v3_family_name=v3_family,
                    v4_family_name=v4_family,
                    v3_value=v3_value,
                    v4_value=v4_value,
                    delta=delta,
                    improved=improved,
                    worsened=worsened,
                    interpretation=interpretation,
                )
            )
    return records


def v3_to_v4_delta_record_to_row(
    record: V3ToV4DeltaRecord,
) -> dict[str, float | str]:
    """Convert one delta record to a CSV-safe row."""

    row = asdict(record)
    row["improved"] = float(record.improved)
    row["worsened"] = float(record.worsened)
    return row


def summarize_v3_to_v4_deltas(
    records: list[V3ToV4DeltaRecord],
) -> list[dict[str, float | str]]:
    """Summarize v3-to-v4 metric changes by linked family pair."""

    grouped: dict[tuple[str, str], list[V3ToV4DeltaRecord]] = defaultdict(list)
    for record in records:
        grouped[(record.v3_family_name, record.v4_family_name)].append(record)
    rows: list[dict[str, float | str]] = []
    for (v3_family, v4_family), family_records in sorted(grouped.items()):
        improved = [record for record in family_records if record.improved]
        worsened = [record for record in family_records if record.worsened]
        unchanged = [
            record for record in family_records
            if record.interpretation == "unchanged"
        ]
        dominant_improvement = (
            max(improved, key=lambda item: abs(item.delta)).metric_name
            if improved
            else ""
        )
        dominant_regression = (
            max(worsened, key=lambda item: abs(item.delta)).metric_name
            if worsened
            else ""
        )
        rows.append(
            {
                "v3_family_name": v3_family,
                "v4_family_name": v4_family,
                "improved_metric_count": float(len(improved)),
                "worsened_metric_count": float(len(worsened)),
                "unchanged_metric_count": float(len(unchanged)),
                "dominant_improvement": dominant_improvement,
                "dominant_regression": dominant_regression,
            }
        )
    return rows
