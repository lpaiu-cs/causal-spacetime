"""Aggregate required metrics for protocol-invariant v3 output bundles."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    default_diagnostic_metric_requirements,
)
from causal_spacetime_lab.state_change_manifest_v2_metric_aggregation import (
    aggregate_v2_required_metrics,
)


def aggregate_v3_protocol_required_metrics(
    *,
    fit_rows: list[dict[str, float | str]],
    null_rows: list[dict[str, float | str]],
    stricter_rows: list[dict[str, float | str]],
    failed_rows: list[dict[str, float | str]],
    coverage_rows: list[dict[str, float | str]],
    restart_rows: list[dict[str, float | str]],
    latent_order_rows: list[dict[str, float | str]],
    no_retuning_audit_pass: bool,
) -> list[dict[str, float | str]]:
    """Aggregate all 14 required metrics plus provenance summaries."""

    rows = aggregate_v2_required_metrics(
        fit_rows=fit_rows,
        null_rows=null_rows,
        stricter_rows=stricter_rows,
        failed_rows=failed_rows,
        restart_rows=restart_rows,
        latent_order_rows=latent_order_rows,
        coverage_rows=coverage_rows,
        no_retuning_audit_pass=no_retuning_audit_pass,
    )
    provenance_by_family = _provenance_counts(failed_rows)
    for row in rows:
        if str(row.get("family_kind", "")) in {"report_only", "failed_control"}:
            _fill_report_only_metric_sentinels(row)
        family = str(row.get("family_name", ""))
        counts = provenance_by_family.get(family, {})
        row["dominant_handoff_provenance_type"] = _dominant(counts)
        row["top_down_template_count"] = float(
            counts.get("top_down_preregistered_template", 0)
        )
        row["hybrid_manifest_count"] = float(
            counts.get("hybrid_template_instantiated_from_profile", 0)
        )
        row["bottom_up_manifest_count"] = float(
            counts.get("bottom_up_profile_derived", 0)
        )
        row["report_only_manifest_count"] = float(
            counts.get("report_only_control", 0)
        )
        row["all_manifests_have_provenance"] = float(
            bool(counts) and counts.get("missing", 0) == 0
        )
    return rows


def _fill_report_only_metric_sentinels(row: dict[str, float | str]) -> None:
    """Use finite accounting sentinels for report-only diagnostic rows."""

    sentinels = {
        "mean_heldout_violation": 1.0,
        "mean_generalization_gap": 0.0,
        "destructive_null_gap": -1.0,
        "symmetry_control_gap": 1.0,
        "restart_std": 0.0,
        "latent_order_disagreement": 0.0,
    }
    for key, value in sentinels.items():
        if not np.isfinite(_to_float(row.get(key, float("nan")))):
            row[key] = value


def v3_protocol_diagnostic_completeness_check(
    metric_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Check diagnostic completeness for protocol-invariant v3 metrics."""

    required = [item.metric_name for item in default_diagnostic_metric_requirements()]
    rows: list[dict[str, float | str]] = []
    for metric_row in metric_rows:
        missing: list[str] = []
        for metric in required:
            value = _to_float(metric_row.get(metric, float("nan")))
            if not np.isfinite(value):
                missing.append(metric)
        available = len(required) - len(missing)
        rows.append(
            {
                "family_name": str(metric_row.get("family_name", "")),
                "family_kind": str(metric_row.get("family_kind", "")),
                "required_metric_count": float(len(required)),
                "available_metric_count": float(available),
                "missing_metric_count": float(len(missing)),
                "completeness_fraction": float(available / len(required)),
                "diagnostic_complete": float(not missing),
                "missing_metrics": ";".join(missing),
                "not_carry_forward_evaluated": 1.0,
                "profile_invariance_status": "protocol_invariant",
                "admissible_for_pairwise_dissimilarity": float(
                    str(metric_row.get("family_kind", "")) == "structured"
                ),
                "dominant_handoff_provenance_type": str(
                    metric_row.get("dominant_handoff_provenance_type", "")
                ),
                "all_manifests_have_provenance": float(
                    metric_row.get("all_manifests_have_provenance", 0.0)
                ),
            }
        )
    return rows


def _provenance_counts(
    rows: list[dict[str, float | str]],
) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {}
    for row in rows:
        family = str(row.get("family_name", ""))
        provenance = str(row.get("handoff_provenance_type", "missing"))
        if not provenance:
            provenance = "missing"
        counts.setdefault(family, {})
        counts[family][provenance] = counts[family].get(provenance, 0) + 1
    return counts


def _dominant(counts: dict[str, int]) -> str:
    if not counts:
        return ""
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def _to_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")
