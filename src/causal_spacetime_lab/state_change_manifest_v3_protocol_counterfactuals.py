"""Report-only counterfactual diagnostics for v3 protocol blocked decisions."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from causal_spacetime_lab.state_change_manifest_v3_protocol_blocking_analysis import (
    V3ProtocolBlockingRecord,
)


def _to_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _structured_rows(
    metric_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    return [row for row in metric_rows if row.get("family_kind") == "structured"]


def heldout_threshold_counterfactual_report(
    metric_rows: list[dict[str, float | str]],
    *,
    hypothetical_thresholds: list[float],
) -> list[dict[str, float | str]]:
    """These counterfactuals are report-only.

    They do not change M42 decisions or justify threshold retuning.
    """

    rows: list[dict[str, float | str]] = []
    for metric_row in _structured_rows(metric_rows):
        heldout = _to_float(metric_row.get("mean_heldout_violation"))
        for threshold in hypothetical_thresholds:
            rows.append(
                {
                    "family_name": str(metric_row.get("family_name", "")),
                    "mean_heldout_violation": heldout,
                    "hypothetical_threshold": float(threshold),
                    "would_pass_heldout_only": float(
                        np.isfinite(heldout) and heldout <= threshold
                    ),
                    "output_note": "report_only_not_decision_changing",
                }
            )
    return rows


def latent_order_threshold_counterfactual_report(
    metric_rows: list[dict[str, float | str]],
    *,
    hypothetical_thresholds: list[float],
) -> list[dict[str, float | str]]:
    """These counterfactuals are report-only.

    They do not change M42 decisions or justify threshold retuning.
    """

    rows: list[dict[str, float | str]] = []
    for metric_row in _structured_rows(metric_rows):
        latent = _to_float(metric_row.get("latent_order_disagreement"))
        for threshold in hypothetical_thresholds:
            rows.append(
                {
                    "family_name": str(metric_row.get("family_name", "")),
                    "latent_order_disagreement": latent,
                    "hypothetical_threshold": float(threshold),
                    "would_pass_latent_order_only": float(
                        np.isfinite(latent) and latent <= threshold
                    ),
                    "output_note": "report_only_not_decision_changing",
                }
            )
    return rows


def family_would_remain_blocked_after_single_fix(
    blocking_records: list[V3ProtocolBlockingRecord],
    *,
    ignored_root_cause: str,
) -> list[dict[str, float | str]]:
    """These counterfactuals are report-only.

    They do not change M42 decisions or justify threshold retuning.
    """

    grouped: dict[str, list[V3ProtocolBlockingRecord]] = defaultdict(list)
    for record in blocking_records:
        if record.family_kind == "structured":
            grouped[record.family_name].append(record)
    rows: list[dict[str, float | str]] = []
    for family_name, records in sorted(grouped.items()):
        remaining = [
            record
            for record in records
            if record.status == "fail"
            and record.blocking_type != "not_blocking"
            and record.root_cause_category != ignored_root_cause
        ]
        ignored = [
            record
            for record in records
            if record.status == "fail"
            and record.root_cause_category == ignored_root_cause
        ]
        rows.append(
            {
                "family_name": family_name,
                "ignored_root_cause": ignored_root_cause,
                "ignored_failure_count": float(len(ignored)),
                "remaining_blocker_count": float(len(remaining)),
                "would_remain_blocked": float(bool(remaining)),
                "remaining_root_causes": ";".join(
                    sorted({record.root_cause_category for record in remaining})
                ),
                "output_note": "report_only_not_decision_changing",
            }
        )
    return rows
