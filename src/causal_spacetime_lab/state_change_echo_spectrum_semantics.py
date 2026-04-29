"""Return-spectrum semantics for state-change echo protocols."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change_echo_interference import (
    earliest_echo_delay_from_spectrum,
    return_delay_spectrum,
)


def _validate_square_matrix(matrix_like: ArrayLike, name: str) -> NDArray[np.bool_]:
    matrix = np.asarray(matrix_like, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{name} must be square")
    return matrix


def _validate_reference(
    reference_chain_event_ids: ArrayLike,
    n_events: int,
) -> NDArray[np.int_]:
    reference = np.asarray(reference_chain_event_ids, dtype=int)
    if reference.ndim != 1:
        raise ValueError("reference_chain_event_ids must be one-dimensional")
    if reference.size and (np.min(reference) < 0 or np.max(reference) >= n_events):
        raise IndexError("reference chain contains event ids outside the matrix")
    return reference


def _validate_target(target_index: int, n_events: int) -> int:
    target = int(target_index)
    if target < 0 or target >= n_events:
        raise IndexError("target_index is outside the matrix")
    return target


def _validate_emission(reference: NDArray[np.int_], emission_position: int) -> int:
    emission = int(emission_position)
    if emission < 0 or emission >= reference.size:
        raise IndexError("emission_position is outside the reference chain")
    return emission


def full_transitive_return_spectrum(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    target_index: int,
    emission_position: int,
) -> NDArray[np.int_]:
    """Return S_full for a target under a full transitive order.

    In a full transitive closure with a full reference chain, this spectrum is
    generally a suffix starting from the earliest return.
    """

    return return_delay_spectrum(
        order_matrix,
        reference_chain_event_ids,
        target_index,
        emission_position,
    )


def retained_reference_return_spectrum(
    order_matrix: ArrayLike,
    retained_reference_chain_event_ids: ArrayLike,
    target_index: int,
    emission_position: int,
) -> NDArray[np.int_]:
    """Return S_retained against retained or subsampled reference ticks.

    This reports only returns visible in the supplied retained reference chain,
    so it can be sparse in the retained reference-rank coordinates.
    """

    return return_delay_spectrum(
        order_matrix,
        retained_reference_chain_event_ids,
        target_index,
        emission_position,
    )


def immediate_edge_return_spectrum(
    immediate_adjacency: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    target_index: int,
    emission_position: int,
) -> NDArray[np.int_]:
    """Return S_immediate using only immediate target-to-reference edges."""

    adjacency = _validate_square_matrix(immediate_adjacency, "immediate_adjacency")
    reference = _validate_reference(reference_chain_event_ids, adjacency.shape[0])
    target = _validate_target(target_index, adjacency.shape[0])
    emission = _validate_emission(reference, emission_position)
    later_positions = np.arange(emission + 1, reference.size, dtype=int)
    if later_positions.size == 0:
        return np.empty(0, dtype=int)
    returning = later_positions[adjacency[target, reference[later_positions]]]
    return np.sort(returning - emission).astype(int)


def is_suffix_spectrum(
    spectrum: ArrayLike,
    *,
    earliest: int | None = None,
    latest: int | None = None,
) -> bool:
    """Return whether a finite spectrum is consecutive from earliest to latest."""

    values = np.asarray(spectrum, dtype=int)
    if values.ndim != 1:
        raise ValueError("spectrum must be one-dimensional")
    unique = np.unique(values)
    if unique.size == 0:
        return True
    start = int(np.min(unique)) if earliest is None else int(earliest)
    end = int(np.max(unique)) if latest is None else int(latest)
    if end < start:
        return unique.size == 0
    return np.array_equal(unique, np.arange(start, end + 1, dtype=int))


def compress_suffix_spectrum(spectrum: ArrayLike) -> dict[str, float]:
    """Compress a return spectrum into suffix summary fields."""

    values = np.asarray(spectrum, dtype=int)
    if values.ndim != 1:
        raise ValueError("spectrum must be one-dimensional")
    earliest = earliest_echo_delay_from_spectrum(values)
    latest = int(np.max(values)) if values.size else None
    return {
        "is_empty": float(values.size == 0),
        "is_suffix": float(is_suffix_spectrum(values)),
        "earliest_delay": float(earliest) if earliest is not None else float("nan"),
        "latest_delay": float(latest) if latest is not None else float("nan"),
        "spectrum_size": float(values.size),
    }


def gated_echo_delay_from_spectrum(
    spectrum: ArrayLike,
    gate_delay_rank: int,
) -> int | None:
    """Return earliest delay satisfying a predeclared gate.

    Gated echo is a separate protocol selected before evaluation. It is not a
    post-hoc correction to the earliest-return rule.
    """

    if gate_delay_rank < 1:
        raise ValueError("gate_delay_rank must be at least 1")
    values = np.asarray(spectrum, dtype=int)
    if values.ndim != 1:
        raise ValueError("spectrum must be one-dimensional")
    retained = np.sort(values[values >= int(gate_delay_rank)])
    return int(retained[0]) if retained.size else None

