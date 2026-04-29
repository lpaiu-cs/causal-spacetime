"""Family-level comparison utilities for frozen-manifest diagnostics."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from causal_spacetime_lab.state_change_manifest_family_config import (
    FamilyComparisonConfig,
)
from causal_spacetime_lab.state_change_manifest_family_diagnostics import (
    ManifestFitDiagnosticRow,
    summarize_family_fit_diagnostics,
)
from causal_spacetime_lab.state_change_manifest_null_taxonomy import (
    classify_null_type,
)


def compare_manifest_family_against_thresholds(
    fit_rows: list[ManifestFitDiagnosticRow],
    config: FamilyComparisonConfig,
) -> list[dict[str, float | str]]:
    """Apply fixed family-level diagnostic thresholds."""

    summaries = summarize_family_fit_diagnostics(fit_rows)
    rows: list[dict[str, float | str]] = []
    for summary in summaries:
        heldout = float(summary["mean_heldout_violation"])
        gap = float(summary["mean_generalization_gap"])
        heldout_pass = bool(
            np.isfinite(heldout)
            and heldout <= config.heldout_violation_threshold
        )
        gap_pass = bool(np.isfinite(gap) and gap <= config.generalization_gap_threshold)
        rows.append(
            {
                **summary,
                "heldout_threshold": float(config.heldout_violation_threshold),
                "generalization_gap_threshold": float(
                    config.generalization_gap_threshold
                ),
                "heldout_pass": float(heldout_pass),
                "generalization_gap_pass": float(gap_pass),
                "threshold_pass": float(heldout_pass and gap_pass),
            }
        )
    return rows


def summarize_nulls_by_taxonomy(
    null_summary_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Summarize representation null results by taxonomy class."""

    grouped: dict[str, list[dict[str, float | str]]] = defaultdict(list)
    for row in null_summary_rows:
        taxonomy = classify_null_type(str(row.get("null_type", "")))
        grouped[taxonomy].append(row)
    rows: list[dict[str, float | str]] = []
    for taxonomy_class, group in sorted(grouped.items()):
        values = np.asarray(
            [float(row["mean_heldout_violation_rate"]) for row in group],
            dtype=float,
        )
        deltas = np.asarray(
            [float(row["structured_minus_null_mean"]) for row in group],
            dtype=float,
        )
        rows.append(
            {
                "taxonomy_class": taxonomy_class,
                "null_count": float(len(group)),
                "mean_null_heldout_violation": float(np.nanmean(values)),
                "mean_structured_minus_null": float(np.nanmean(deltas)),
            }
        )
    return rows


def family_report_card(
    family_summary_rows: list[dict[str, float | str]],
    null_taxonomy_rows: list[dict[str, float | str]],
    failure_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Build a compact family-level report card."""

    by_family: dict[str, list[dict[str, float | str]]] = defaultdict(list)
    for row in family_summary_rows:
        by_family[str(row["family_name"])].append(row)

    failure_summary = ";".join(
        f"{row.get('reason', '')}:{int(float(row['count']))}"
        for row in failure_rows
        if row.get("row_type") == "failure_reason"
        if float(row.get("count", 0.0)) > 0
    )
    null_summary = ";".join(
        f"{row['taxonomy_class']}={float(row['mean_null_heldout_violation']):.3f}"
        for row in null_taxonomy_rows
    )

    cards: list[dict[str, float | str]] = []
    for family_name, rows in sorted(by_family.items()):
        fitted = sum(float(row["fitted_count"]) for row in rows)
        no_fit = sum(float(row["no_fit_count"]) for row in rows)
        best_values = np.asarray(
            [
                float(row.get("best_heldout_violation", row["mean_heldout_violation"]))
                for row in rows
            ],
            dtype=float,
        )
        mean_values = np.asarray(
            [float(row["mean_heldout_violation"]) for row in rows],
            dtype=float,
        )
        finite_best = best_values[np.isfinite(best_values)]
        finite_mean = mean_values[np.isfinite(mean_values)]
        best_heldout = (
            float(np.min(finite_best)) if finite_best.size else float("nan")
        )
        mean_heldout = (
            float(np.mean(finite_mean)) if finite_mean.size else float("nan")
        )
        cards.append(
            {
                "family_name": family_name,
                "family_kind": str(rows[0]["family_kind"]),
                "fitted_count": float(fitted),
                "no_fit_count": float(no_fit),
                "best_heldout_violation": best_heldout,
                "mean_heldout_violation": mean_heldout,
                "null_comparison_summary": null_summary,
                "failure_reason_summary": failure_summary,
                "interpretation_warning": (
                    "family_level_latent_diagnostic_not_physical_geometry"
                ),
            }
        )
    return cards
