"""Restart and latent-order instability drilldown for v4 protocol families."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    CrossFamilyRobustnessCriteria,
)


def _to_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _mean(values: list[float]) -> float:
    finite = [value for value in values if np.isfinite(value)]
    return float(np.mean(finite)) if finite else float("nan")


def summarize_v4_stability_failures(
    restart_rows: list[dict[str, float | str]],
    latent_order_rows: list[dict[str, float | str]],
    criteria: CrossFamilyRobustnessCriteria,
) -> list[dict[str, float | str]]:
    """Summarize restart and latent-order stability failures."""

    families = sorted(
        {str(row.get("family_name", "")) for row in restart_rows}
        | {str(row.get("family_name", "")) for row in latent_order_rows}
    )
    restart_by_family: dict[str, list[float]] = defaultdict(list)
    latent_by_family: dict[str, list[float]] = defaultdict(list)
    for row in restart_rows:
        restart_by_family[str(row.get("family_name", ""))].append(
            _to_float(row.get("restart_std"))
        )
    for row in latent_order_rows:
        latent_by_family[str(row.get("family_name", ""))].append(
            _to_float(row.get("latent_order_disagreement"))
        )
    rows: list[dict[str, float | str]] = []
    for family in families:
        restart_std = _mean(restart_by_family[family])
        latent = _mean(latent_by_family[family])
        restart_pass = (
            np.isfinite(restart_std)
            and restart_std <= criteria.max_restart_std
        )
        latent_pass = (
            np.isfinite(latent) and latent <= criteria.max_latent_order_disagreement
        )
        failures = []
        if not restart_pass:
            failures.append("restart_instability")
        if not latent_pass:
            failures.append("latent_order_instability")
        rows.append(
            {
                "family_name": family,
                "restart_std": restart_std,
                "restart_std_threshold": criteria.max_restart_std,
                "restart_std_pass": float(restart_pass),
                "latent_order_disagreement": latent,
                "latent_order_threshold": criteria.max_latent_order_disagreement,
                "latent_order_pass": float(latent_pass),
                "dominant_stability_failure": ";".join(failures) or "none",
            }
        )
    return rows
