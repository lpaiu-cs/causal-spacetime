"""Diagnostic completeness audit for cross-family robustness metrics."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class DiagnosticCompletenessReport:
    """Completeness report for one family metric row."""

    family_name: str
    required_metric_count: int
    available_metric_count: int
    missing_metric_count: int
    completeness_fraction: float
    missing_metrics: list[str]


def required_cross_family_metrics() -> list[str]:
    """Return metrics required by carry-forward criteria."""

    return [
        "manifest_count",
        "fitted_fraction",
        "no_fit_fraction",
        "mean_heldout_violation",
        "mean_generalization_gap",
        "stricter_threshold_pass_fraction",
        "destructive_null_gap",
        "symmetry_control_gap",
        "target_coverage_fraction",
        "pair_node_coverage_fraction",
        "restart_std",
        "latent_order_disagreement",
        "no_retuning_audit_pass",
        "failed_accounting_present",
    ]


def _available(value: object) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def diagnostic_completeness_for_family(
    family_metrics: dict[str, float | str],
) -> DiagnosticCompletenessReport:
    """Compute required-metric completeness for one family."""

    required = required_cross_family_metrics()
    missing = [
        metric
        for metric in required
        if metric not in family_metrics or not _available(family_metrics[metric])
    ]
    available_count = len(required) - len(missing)
    return DiagnosticCompletenessReport(
        family_name=str(family_metrics.get("family_name", "")),
        required_metric_count=len(required),
        available_metric_count=available_count,
        missing_metric_count=len(missing),
        completeness_fraction=available_count / len(required) if required else 1.0,
        missing_metrics=missing,
    )


def diagnostic_completeness_table(
    metrics_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Return family-level diagnostic completeness rows."""

    rows: list[dict[str, float | str]] = []
    for metrics in metrics_rows:
        report = diagnostic_completeness_for_family(metrics)
        row = asdict(report)
        row["missing_metrics"] = ";".join(report.missing_metrics)
        rows.append(row)
    return rows
