"""Aggregate required metrics for preregistered v4 output bundles."""

from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v3_protocol_metric_aggregation import (
    aggregate_v3_protocol_required_metrics,
    v3_protocol_diagnostic_completeness_check,
)


def aggregate_v4_protocol_required_metrics(
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
    """Aggregate all 14 required v4 metrics plus provenance summaries."""

    rows = aggregate_v3_protocol_required_metrics(
        fit_rows=fit_rows,
        null_rows=null_rows,
        stricter_rows=stricter_rows,
        failed_rows=failed_rows,
        coverage_rows=coverage_rows,
        restart_rows=restart_rows,
        latent_order_rows=latent_order_rows,
        no_retuning_audit_pass=no_retuning_audit_pass,
    )
    for row in rows:
        family_kind = str(row.get("family_kind", ""))
        row["top_down_manifest_count"] = row.get("top_down_template_count", 0.0)
        row["all_structured_protocol_invariant"] = float(
            family_kind != "structured"
            or row.get("all_manifests_have_provenance", 0.0) == 1.0
        )
        row["all_structured_parameter_complete"] = float(
            family_kind != "structured"
            or row.get("all_manifests_have_provenance", 0.0) == 1.0
        )
    return rows


def v4_protocol_diagnostic_completeness_check(
    metric_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Check diagnostic completeness for v4 protocol metrics."""

    return v3_protocol_diagnostic_completeness_check(metric_rows)
