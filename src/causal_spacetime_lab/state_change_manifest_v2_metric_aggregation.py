"""Aggregate diagnostic-complete v2 family metrics."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_manifest_diagnostic_schema import (
    default_diagnostic_metric_requirements,
)
from causal_spacetime_lab.state_change_manifest_family_robustness_metrics import (
    destructive_null_gap_from_rows,
    stricter_pass_fraction_for_family,
    symmetry_control_gap_from_rows,
)


def _to_float(value: object, default: float = float("nan")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _finite(values: list[float]) -> list[float]:
    return [value for value in values if np.isfinite(value)]


def _mean(values: list[float]) -> float:
    finite = _finite(values)
    return float(np.mean(finite)) if finite else float("nan")


def _family_rows(
    rows: list[dict[str, float | str]],
    family_name: str,
) -> list[dict[str, float | str]]:
    return [row for row in rows if row.get("family_name") == family_name]


def _family_names(*row_groups: list[dict[str, float | str]]) -> list[str]:
    names: set[str] = set()
    for rows in row_groups:
        names.update(str(row.get("family_name", "")) for row in rows)
    return sorted(name for name in names if name)


def _family_kind(
    family_name: str,
    *row_groups: list[dict[str, float | str]],
) -> str:
    for rows in row_groups:
        for row in rows:
            if row.get("family_name") == family_name:
                return str(row.get("family_kind", ""))
    return ""


def _control_placeholder(value: float, family_kind: str, placeholder: float) -> float:
    if np.isfinite(value):
        return value
    if family_kind in {"failed_control", "ineligible_control"}:
        return float(placeholder)
    return value


def _string_rows(
    rows: list[dict[str, float | str]],
) -> list[dict[str, str]]:
    return [{key: str(value) for key, value in row.items()} for row in rows]


def aggregate_v2_required_metrics(
    *,
    fit_rows: list[dict[str, float | str]],
    null_rows: list[dict[str, float | str]],
    stricter_rows: list[dict[str, float | str]],
    failed_rows: list[dict[str, float | str]],
    restart_rows: list[dict[str, float | str]],
    latent_order_rows: list[dict[str, float | str]],
    coverage_rows: list[dict[str, float | str]],
    no_retuning_audit_pass: bool,
) -> list[dict[str, float | str]]:
    """Aggregate all 14 required v2 family-level metrics."""

    family_names = _family_names(
        fit_rows,
        null_rows,
        stricter_rows,
        failed_rows,
        restart_rows,
        latent_order_rows,
        coverage_rows,
    )
    rows: list[dict[str, float | str]] = []
    null_string_rows = _string_rows(null_rows)
    stricter_string_rows = _string_rows(stricter_rows)
    for family_name in family_names:
        family_fit = _family_rows(fit_rows, family_name)
        family_kind = _family_kind(
            family_name,
            fit_rows,
            failed_rows,
            coverage_rows,
        )
        fitted_count = float(
            sum(_to_float(row.get("fitted_count"), 0.0) for row in family_fit)
        )
        no_fit_count = float(
            sum(_to_float(row.get("no_fit_count"), 0.0) for row in family_fit)
        )
        manifest_count = max(
            _finite([_to_float(row.get("manifest_count")) for row in family_fit])
            or [float(len(_family_rows(failed_rows, family_name)))]
        )
        denominator = fitted_count + no_fit_count
        structured_heldout = _mean(
            [_to_float(row.get("mean_heldout_violation")) for row in family_fit]
        )
        structured_heldout = _control_placeholder(
            structured_heldout,
            family_kind,
            1.0,
        )
        generalization_gap = _control_placeholder(
            _mean(
                [
                    _to_float(row.get("mean_generalization_gap"))
                    for row in family_fit
                ]
            ),
            family_kind,
            0.0,
        )
        coverage = _family_rows(coverage_rows, family_name)
        restart = _family_rows(restart_rows, family_name)
        latent = _family_rows(latent_order_rows, family_name)
        destructive_null_gap = _control_placeholder(
            destructive_null_gap_from_rows(
                family_name,
                null_string_rows,
                structured_heldout,
            ),
            family_kind,
            0.0,
        )
        symmetry_control_gap = _control_placeholder(
            symmetry_control_gap_from_rows(
                family_name,
                null_string_rows,
                structured_heldout,
            ),
            family_kind,
            0.0,
        )
        restart_std = _control_placeholder(
            _mean([_to_float(row.get("restart_std")) for row in restart]),
            family_kind,
            0.0,
        )
        latent_order_disagreement = _control_placeholder(
            _mean(
                [
                    _to_float(row.get("latent_order_disagreement"))
                    for row in latent
                ]
            ),
            family_kind,
            0.0,
        )
        rows.append(
            {
                "family_name": family_name,
                "family_kind": family_kind,
                "manifest_count": float(manifest_count),
                "fitted_fraction": float(fitted_count / denominator)
                if denominator
                else 0.0,
                "no_fit_fraction": float(no_fit_count / denominator)
                if denominator
                else 1.0,
                "mean_heldout_violation": structured_heldout,
                "mean_generalization_gap": generalization_gap,
                "stricter_threshold_pass_fraction": stricter_pass_fraction_for_family(
                    family_name,
                    stricter_string_rows,
                ),
                "destructive_null_gap": destructive_null_gap,
                "symmetry_control_gap": symmetry_control_gap,
                "target_coverage_fraction": _mean(
                    [_to_float(row.get("target_coverage_fraction")) for row in coverage]
                ),
                "pair_node_coverage_fraction": _mean(
                    [
                        _to_float(row.get("pair_node_coverage_fraction"))
                        for row in coverage
                    ]
                ),
                "restart_std": restart_std,
                "latent_order_disagreement": latent_order_disagreement,
                "no_retuning_audit_pass": float(no_retuning_audit_pass),
                "failed_accounting_present": float(bool(failed_rows)),
            }
        )
    return rows


def v2_diagnostic_completeness_check(
    metric_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Check diagnostic completeness for v2 family metric rows."""

    required = [item.metric_name for item in default_diagnostic_metric_requirements()]
    rows: list[dict[str, float | str]] = []
    for metric_row in metric_rows:
        missing: list[str] = []
        for metric in required:
            value = metric_row.get(metric, float("nan"))
            if not np.isfinite(_to_float(value)):
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
                "missing_metrics": ";".join(missing),
                "not_carry_forward_evaluated": 1.0,
            }
        )
    return rows
