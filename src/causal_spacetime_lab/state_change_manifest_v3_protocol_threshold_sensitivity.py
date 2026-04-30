"""Threshold sensitivity for v3 protocol carry-forward evaluation."""

from __future__ import annotations

from dataclasses import replace

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_carry_forward import (
    decide_v3_protocol_family_robustness,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preconditions import (
    V3ProtocolPreconditionReport,
)


def v3_protocol_threshold_sensitivity_rows(
    metric_rows: list[dict[str, float | str]],
    preconditions: list[V3ProtocolPreconditionReport],
    *,
    heldout_thresholds: list[float],
    null_gap_thresholds: list[float],
    stricter_pass_thresholds: list[float],
) -> list[dict[str, float | str]]:
    """Run sensitivity sweeps without mutating fixed default criteria."""

    rows: list[dict[str, float | str]] = []
    default = default_cross_family_robustness_criteria()
    for heldout in heldout_thresholds:
        for null_gap in null_gap_thresholds:
            for stricter in stricter_pass_thresholds:
                criteria = replace(
                    default,
                    max_mean_heldout_violation=float(heldout),
                    min_destructive_null_gap=float(null_gap),
                    min_stricter_threshold_pass_fraction=float(stricter),
                )
                decisions = decide_v3_protocol_family_robustness(
                    metric_rows,
                    criteria,
                    preconditions,
                )
                counts = {
                    "carry_forward": 0,
                    "provisional": 0,
                    "blocked": 0,
                    "report_only": 0,
                    "failed_control": 0,
                }
                for decision in decisions:
                    counts[decision.decision] += 1
                rows.append(
                    {
                        "heldout_threshold": float(heldout),
                        "null_gap_threshold": float(null_gap),
                        "stricter_pass_threshold": float(stricter),
                        **{
                            f"{key}_count": float(value)
                            for key, value in counts.items()
                        },
                        "analysis_label": "sensitivity_only_not_retuning",
                    }
                )
    return rows

