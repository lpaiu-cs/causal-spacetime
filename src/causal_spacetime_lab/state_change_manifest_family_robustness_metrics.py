"""Aggregate robustness metrics from manifest-family output bundles."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_manifest_family_outputs import (
    family_names_from_bundle,
)
from causal_spacetime_lab.state_change_manifest_null_taxonomy import (
    classify_null_type,
)


def _to_float(value: object, default: float = float("nan")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _finite_mean(values: list[float]) -> float:
    finite = [value for value in values if np.isfinite(value)]
    return float(np.mean(finite)) if finite else float("nan")


def _finite_max(values: list[float]) -> float:
    finite = [value for value in values if np.isfinite(value)]
    return float(np.max(finite)) if finite else float("nan")


def _family_rows(
    rows: list[dict[str, str]],
    family_name: str,
) -> list[dict[str, str]]:
    return [row for row in rows if row.get("family_name") == family_name]


def destructive_null_gap_from_rows(
    family_name: str,
    null_rows: list[dict[str, str]],
    structured_heldout: float,
) -> float:
    """Return destructive-null heldout minus structured heldout."""

    candidates = [
        row
        for row in null_rows
        if classify_null_type(row.get("null_type", "")) == "destructive_null"
        or row.get("taxonomy_class") == "destructive_null"
    ]
    family_specific = [
        row for row in candidates if row.get("family_name") == family_name
    ]
    selected = family_specific or candidates
    null_mean = _finite_mean(
        [_to_float(row.get("mean_heldout_violation_rate")) for row in selected]
    )
    if not np.isfinite(null_mean) or not np.isfinite(structured_heldout):
        return float("nan")
    return float(null_mean - structured_heldout)


def symmetry_control_gap_from_rows(
    family_name: str,
    null_rows: list[dict[str, str]],
    structured_heldout: float,
) -> float:
    """Return absolute symmetry-control heldout difference."""

    candidates = [
        row
        for row in null_rows
        if classify_null_type(row.get("null_type", "")) == "symmetry_control"
        or row.get("taxonomy_class") == "symmetry_control"
    ]
    family_specific = [
        row for row in candidates if row.get("family_name") == family_name
    ]
    selected = family_specific or candidates
    null_mean = _finite_mean(
        [_to_float(row.get("mean_heldout_violation_rate")) for row in selected]
    )
    if not np.isfinite(null_mean) or not np.isfinite(structured_heldout):
        return float("nan")
    return float(abs(null_mean - structured_heldout))


def stricter_pass_fraction_for_family(
    family_name: str,
    stricter_rows: list[dict[str, str]],
) -> float:
    """Return fraction of stricter-threshold rows passing for a family."""

    rows = _family_rows(stricter_rows, family_name)
    values = [_to_float(row.get("threshold_pass")) for row in rows]
    return _finite_mean(values)


def _audit_pass(audit_rows: list[dict[str, str]]) -> float:
    if not audit_rows:
        return 0.0
    values = [_to_float(row.get("passed")) for row in audit_rows]
    return float(all(value == 1.0 for value in values))


def _failed_accounting_present(rows: list[dict[str, str]]) -> float:
    return float(bool(rows))


def _family_kind(
    family_name: str,
    fit_rows: list[dict[str, str]],
    report_rows: list[dict[str, str]],
) -> str:
    for row in fit_rows + report_rows:
        if row.get("family_name") == family_name:
            return row.get("family_kind", "")
    return ""


def aggregate_family_robustness_metrics(
    bundle: dict[str, list[dict[str, str]]],
) -> list[dict[str, float | str]]:
    """Aggregate family-level metrics for robustness decisions."""

    fit_rows = bundle.get("fit_summary", [])
    null_rows = bundle.get("null_taxonomy", [])
    stricter_rows = bundle.get("stricter_criteria", [])
    audit_rows = bundle.get("no_retuning_audit", [])
    failed_rows = bundle.get("failed_accounting", [])
    report_rows = bundle.get("report_card", [])
    family_names = family_names_from_bundle(bundle)

    rows: list[dict[str, float | str]] = []
    for family_name in family_names:
        family_fit = _family_rows(fit_rows, family_name)
        manifest_count = _finite_max(
            [_to_float(row.get("manifest_count")) for row in family_fit]
        )
        fitted_count = float(
            sum(_to_float(row.get("fitted_count"), 0.0) for row in family_fit)
        )
        no_fit_count = float(
            sum(_to_float(row.get("no_fit_count"), 0.0) for row in family_fit)
        )
        denominator = fitted_count + no_fit_count
        fitted_fraction = fitted_count / denominator if denominator else float("nan")
        no_fit_fraction = no_fit_count / denominator if denominator else float("nan")
        structured_heldout = _finite_mean(
            [_to_float(row.get("mean_heldout_violation")) for row in family_fit]
        )
        row = {
            "family_name": family_name,
            "family_kind": _family_kind(family_name, fit_rows, report_rows),
            "manifest_count": manifest_count,
            "fitted_count": fitted_count,
            "no_fit_count": no_fit_count,
            "fitted_fraction": float(fitted_fraction),
            "no_fit_fraction": float(no_fit_fraction),
            "mean_heldout_violation": structured_heldout,
            "mean_generalization_gap": _finite_mean(
                [_to_float(item.get("mean_generalization_gap")) for item in family_fit]
            ),
            "stricter_threshold_pass_fraction": stricter_pass_fraction_for_family(
                family_name,
                stricter_rows,
            ),
            "destructive_null_gap": destructive_null_gap_from_rows(
                family_name,
                null_rows,
                structured_heldout,
            ),
            "symmetry_control_gap": symmetry_control_gap_from_rows(
                family_name,
                null_rows,
                structured_heldout,
            ),
            "target_coverage_fraction": float("nan"),
            "pair_node_coverage_fraction": float("nan"),
            "restart_std": float("nan"),
            "latent_order_disagreement": float("nan"),
            "no_retuning_audit_pass": _audit_pass(audit_rows),
            "failed_accounting_present": _failed_accounting_present(failed_rows),
        }
        rows.append(row)
    return rows
