"""Validation helpers for controlled echo-response motifs."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from causal_spacetime_lab.ordinal import order_agreement_rate, order_inversion_rate
from causal_spacetime_lab.state_change_echo import echo_delay_rank_for_emission
from causal_spacetime_lab.state_change_echo_motifs import EchoMotifRecord


def recovered_delay_for_motif(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    motif: EchoMotifRecord,
) -> int | None:
    """Recover the echo-delay rank for one planted motif target."""

    return echo_delay_rank_for_emission(
        order_matrix,
        reference_chain_event_ids,
        motif.target_event_id,
        motif.emission_position,
    )


def motif_recovery_report(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    motif: EchoMotifRecord,
) -> dict[str, float | str]:
    """Return motif-level planted versus recovered delay diagnostics."""

    recovered = recovered_delay_for_motif(
        order_matrix,
        reference_chain_event_ids,
        motif,
    )
    planted = int(motif.planted_delay_rank)
    exact = recovered == planted
    early = recovered is not None and recovered < planted
    late_or_missing = recovered is None or recovered > planted
    return {
        "target_event_id": float(motif.target_event_id),
        "emission_position": float(motif.emission_position),
        "planted_return_position": float(motif.planted_return_position),
        "planted_delay_rank": float(planted),
        "recovered_delay_rank": float(recovered)
        if recovered is not None
        else float("nan"),
        "exact_recovery": float(exact),
        "early_shortcut": float(early),
        "late_or_missing": float(late_or_missing),
        "recovery_error": float(recovered - planted)
        if recovered is not None
        else float("nan"),
        "motif_kind": motif.motif_kind,
    }


def _finite_values(rows: list[dict[str, float | str]], key: str) -> np.ndarray:
    values: list[float] = []
    for row in rows:
        value = float(row[key])
        if np.isfinite(value):
            values.append(value)
    return np.asarray(values, dtype=float)


def summarize_motif_recovery(
    rows: list[dict[str, float | str]],
) -> dict[str, float]:
    """Summarize motif recovery fractions and delay errors."""

    motif_count = len(rows)
    if motif_count == 0:
        return {
            "motif_count": 0.0,
            "exact_recovery_fraction": float("nan"),
            "early_shortcut_fraction": float("nan"),
            "late_or_missing_fraction": float("nan"),
            "mean_absolute_recovery_error": float("nan"),
            "mean_recovered_delay_rank": float("nan"),
            "mean_planted_delay_rank": float("nan"),
        }
    errors = _finite_values(rows, "recovery_error")
    recovered = _finite_values(rows, "recovered_delay_rank")
    planted = _finite_values(rows, "planted_delay_rank")
    return {
        "motif_count": float(motif_count),
        "exact_recovery_fraction": float(
            np.mean([float(row["exact_recovery"]) for row in rows])
        ),
        "early_shortcut_fraction": float(
            np.mean([float(row["early_shortcut"]) for row in rows])
        ),
        "late_or_missing_fraction": float(
            np.mean([float(row["late_or_missing"]) for row in rows])
        ),
        "mean_absolute_recovery_error": float(np.mean(np.abs(errors)))
        if errors.size
        else float("nan"),
        "mean_recovered_delay_rank": float(np.mean(recovered))
        if recovered.size
        else float("nan"),
        "mean_planted_delay_rank": float(np.mean(planted))
        if planted.size
        else float("nan"),
    }


def motif_order_recovery_rate(
    motif_rows: list[dict[str, float | str]],
) -> dict[str, float]:
    """Compare planted and recovered delay-rank order across motifs."""

    planted: list[float] = []
    recovered: list[float] = []
    for row in motif_rows:
        recovered_value = float(row["recovered_delay_rank"])
        if not np.isfinite(recovered_value):
            continue
        planted.append(float(row["planted_delay_rank"]))
        recovered.append(recovered_value)
    if len(planted) < 2:
        inversion = float("nan")
        agreement = float("nan")
    else:
        inversion = order_inversion_rate(
            np.asarray(planted, dtype=float),
            np.asarray(recovered, dtype=float),
            ignore_ties=True,
        )
        agreement = order_agreement_rate(
            np.asarray(planted, dtype=float),
            np.asarray(recovered, dtype=float),
            ignore_ties=True,
        )
    return {
        "comparable_motif_count": float(len(planted)),
        "order_inversion_rate": inversion,
        "order_agreement_rate": agreement,
    }
