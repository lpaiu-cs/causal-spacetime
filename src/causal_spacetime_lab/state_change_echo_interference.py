"""Return-spectrum and shortcut-interference diagnostics for echo motifs."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change_echo import echo_delay_rank_for_emission
from causal_spacetime_lab.state_change_echo_motifs import EchoMotifRecord


def _validate_order_matrix(order_matrix: ArrayLike) -> NDArray[np.bool_]:
    matrix = np.asarray(order_matrix, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("order_matrix must be square")
    return matrix


def _validate_indices(indices: ArrayLike, n_events: int) -> NDArray[np.int_]:
    values = np.asarray(indices, dtype=int)
    if values.ndim != 1:
        raise ValueError("indices must be one-dimensional")
    if values.size and (np.min(values) < 0 or np.max(values) >= n_events):
        raise IndexError("indices are outside the order matrix")
    return values


def _validate_emission_position(chain: NDArray[np.int_], emission_position: int) -> int:
    position = int(emission_position)
    if position < 0 or position >= chain.size:
        raise IndexError("emission_position is outside the reference chain")
    return position


def return_positions_for_target(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    target_index: int,
    emission_position: int,
) -> NDArray[np.int_]:
    """Return all later reference positions reachable from a target.

    This only reports possible returns after the selected emission position. It
    does not require the target to be reachable after the emission.
    """

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(reference_chain_event_ids, matrix.shape[0])
    target = int(_validate_indices(np.asarray([target_index]), matrix.shape[0])[0])
    emission = _validate_emission_position(chain, emission_position)
    positions = np.arange(emission + 1, chain.size, dtype=int)
    if positions.size == 0:
        return np.empty(0, dtype=int)
    return positions[matrix[target, chain[positions]]]


def return_delay_spectrum(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    target_index: int,
    emission_position: int,
) -> NDArray[np.int_]:
    """Return sorted return-delay ranks for one target and emission."""

    positions = return_positions_for_target(
        order_matrix,
        reference_chain_event_ids,
        target_index,
        emission_position,
    )
    return np.sort(positions - int(emission_position)).astype(int)


def is_target_reachable_after_emission(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    target_index: int,
    emission_position: int,
) -> bool:
    """Return whether the selected reference emission precedes the target."""

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(reference_chain_event_ids, matrix.shape[0])
    target = int(_validate_indices(np.asarray([target_index]), matrix.shape[0])[0])
    emission = _validate_emission_position(chain, emission_position)
    return bool(matrix[int(chain[emission]), target])


def earliest_echo_delay_from_spectrum(delay_spectrum: ArrayLike) -> int | None:
    """Return the earliest return delay, or None for an empty spectrum."""

    spectrum = np.asarray(delay_spectrum, dtype=int)
    if spectrum.ndim != 1:
        raise ValueError("delay_spectrum must be one-dimensional")
    return int(np.min(spectrum)) if spectrum.size else None


def shortcut_positions_before_planted(
    delay_spectrum: ArrayLike,
    planted_delay_rank: int,
) -> NDArray[np.int_]:
    """Return all return delays earlier than the planted delay rank."""

    spectrum = np.asarray(delay_spectrum, dtype=int)
    if spectrum.ndim != 1:
        raise ValueError("delay_spectrum must be one-dimensional")
    if planted_delay_rank < 1:
        raise ValueError("planted_delay_rank must be positive")
    return spectrum[spectrum < int(planted_delay_rank)]


def shortcut_depth(
    recovered_delay_rank: int | None,
    planted_delay_rank: int,
) -> float:
    """Return planted minus recovered delay rank.

    Exact recovery gives 0. Early shortcut returns give positive values. Late
    recoveries give negative values. Missing recoveries return NaN.
    """

    if planted_delay_rank < 1:
        raise ValueError("planted_delay_rank must be positive")
    if recovered_delay_rank is None:
        return float("nan")
    return float(int(planted_delay_rank) - int(recovered_delay_rank))


def return_spectrum_report_for_motif(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    motif: EchoMotifRecord,
) -> dict[str, float | str]:
    """Return spectrum-level shortcut classification for one planted motif."""

    reachable = is_target_reachable_after_emission(
        order_matrix,
        reference_chain_event_ids,
        motif.target_event_id,
        motif.emission_position,
    )
    spectrum = return_delay_spectrum(
        order_matrix,
        reference_chain_event_ids,
        motif.target_event_id,
        motif.emission_position,
    )
    recovered = (
        echo_delay_rank_for_emission(
            order_matrix,
            reference_chain_event_ids,
            motif.target_event_id,
            motif.emission_position,
        )
        if reachable
        else None
    )
    planted = int(motif.planted_delay_rank)
    shortcuts = shortcut_positions_before_planted(spectrum, planted)
    exact = recovered == planted
    early = recovered is not None and recovered < planted
    missing = recovered is None
    late = recovered is not None and recovered > planted
    depth = shortcut_depth(recovered, planted)
    return {
        "target_event_id": float(motif.target_event_id),
        "emission_position": float(motif.emission_position),
        "planted_delay_rank": float(planted),
        "recovered_delay_rank": float(recovered)
        if recovered is not None
        else float("nan"),
        "target_reachable_after_emission": float(reachable),
        "spectrum_size": float(spectrum.size),
        "earliest_delay": float(np.min(spectrum)) if spectrum.size else float("nan"),
        "latest_delay": float(np.max(spectrum)) if spectrum.size else float("nan"),
        "shortcut_count": float(shortcuts.size),
        "shortcut_depth": depth,
        "exact_recovery": float(exact),
        "early_shortcut": float(early),
        "missing_return": float(missing),
        "late_recovery": float(late),
        "return_spectrum_string": ";".join(str(int(value)) for value in spectrum),
    }


def summarize_return_spectrum_reports(
    rows: list[dict[str, float | str]],
) -> dict[str, float]:
    """Summarize return-spectrum shortcut classifications."""

    if not rows:
        return {
            "motif_count": 0.0,
            "exact_recovery_fraction": float("nan"),
            "shortcut_fraction": float("nan"),
            "missing_fraction": float("nan"),
            "late_fraction": float("nan"),
            "mean_shortcut_count": float("nan"),
            "mean_shortcut_depth": float("nan"),
            "mean_spectrum_size": float("nan"),
            "mean_earliest_delay": float("nan"),
            "mean_latest_delay": float("nan"),
        }

    def mean_finite(key: str) -> float:
        values = np.asarray([float(row[key]) for row in rows], dtype=float)
        values = values[np.isfinite(values)]
        return float(np.mean(values)) if values.size else float("nan")

    positive_depths = np.asarray(
        [
            float(row["shortcut_depth"])
            for row in rows
            if np.isfinite(float(row["shortcut_depth"]))
            and float(row["shortcut_depth"]) > 0.0
        ],
        dtype=float,
    )
    return {
        "motif_count": float(len(rows)),
        "exact_recovery_fraction": mean_finite("exact_recovery"),
        "shortcut_fraction": mean_finite("early_shortcut"),
        "missing_fraction": mean_finite("missing_return"),
        "late_fraction": mean_finite("late_recovery"),
        "mean_shortcut_count": mean_finite("shortcut_count"),
        "mean_shortcut_depth": float(np.mean(positive_depths))
        if positive_depths.size
        else 0.0,
        "mean_spectrum_size": mean_finite("spectrum_size"),
        "mean_earliest_delay": mean_finite("earliest_delay"),
        "mean_latest_delay": mean_finite("latest_delay"),
    }
