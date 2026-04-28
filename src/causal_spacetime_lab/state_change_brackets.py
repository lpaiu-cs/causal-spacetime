"""Order-level brackets from selected reference chains."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


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


def latest_predecessor_reference_position(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    target_index: int,
) -> int | None:
    """Return latest reference-chain position preceding a target."""

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(reference_chain_event_ids, matrix.shape[0])
    _validate_indices(np.asarray([target_index]), matrix.shape[0])
    positions = np.flatnonzero(matrix[chain, int(target_index)])
    return int(positions[-1]) if positions.size else None


def earliest_successor_reference_position(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    target_index: int,
) -> int | None:
    """Return earliest reference-chain position succeeding a target."""

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(reference_chain_event_ids, matrix.shape[0])
    _validate_indices(np.asarray([target_index]), matrix.shape[0])
    positions = np.flatnonzero(matrix[int(target_index), chain])
    return int(positions[0]) if positions.size else None


def reference_tick_brackets_from_order(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    target_indices: ArrayLike | None = None,
    include_reference_events: bool = False,
) -> tuple[NDArray[np.int_], NDArray[np.int_], NDArray[np.bool_]]:
    """Return order-level predecessor/successor brackets for target events.

    By default, events on the reference chain are marked inaccessible because
    they are reference ticks rather than bracketed targets.
    """

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(reference_chain_event_ids, matrix.shape[0])
    targets = (
        np.arange(matrix.shape[0], dtype=int)
        if target_indices is None
        else _validate_indices(target_indices, matrix.shape[0])
    )
    chain_set = set(chain.tolist())
    predecessor_positions = np.full(targets.size, -1, dtype=int)
    successor_positions = np.full(targets.size, -1, dtype=int)
    accessible_mask = np.zeros(targets.size, dtype=bool)
    for row_index, target in enumerate(targets):
        if not include_reference_events and int(target) in chain_set:
            continue
        predecessor = latest_predecessor_reference_position(matrix, chain, int(target))
        successor = earliest_successor_reference_position(matrix, chain, int(target))
        if predecessor is not None:
            predecessor_positions[row_index] = predecessor
        if successor is not None:
            successor_positions[row_index] = successor
        accessible_mask[row_index] = predecessor is not None and successor is not None
    return predecessor_positions, successor_positions, accessible_mask


def radar_time_rank_from_reference_brackets(
    predecessor_positions: ArrayLike,
    successor_positions: ArrayLike,
    accessible_mask: ArrayLike,
) -> NDArray[np.int_]:
    """Return order-level radar-time ranks from reference brackets."""

    predecessors = np.asarray(predecessor_positions, dtype=int)
    successors = np.asarray(successor_positions, dtype=int)
    accessible = np.asarray(accessible_mask, dtype=bool)
    if predecessors.shape != successors.shape or predecessors.shape != accessible.shape:
        raise ValueError("predecessor, successor, and mask arrays must match")
    ranks = np.full(predecessors.shape, -1, dtype=int)
    ranks[accessible] = predecessors[accessible] + successors[accessible]
    return ranks


def bracket_width_rank_from_reference_brackets(
    predecessor_positions: ArrayLike,
    successor_positions: ArrayLike,
    accessible_mask: ArrayLike,
) -> NDArray[np.int_]:
    """Return order-level bracket-width rank, not metric radar distance."""

    predecessors = np.asarray(predecessor_positions, dtype=int)
    successors = np.asarray(successor_positions, dtype=int)
    accessible = np.asarray(accessible_mask, dtype=bool)
    if predecessors.shape != successors.shape or predecessors.shape != accessible.shape:
        raise ValueError("predecessor, successor, and mask arrays must match")
    ranks = np.full(predecessors.shape, -1, dtype=int)
    ranks[accessible] = successors[accessible] - predecessors[accessible]
    return ranks


def assign_reference_rank_slices(
    radar_time_ranks: ArrayLike,
    accessible_mask: ArrayLike,
    bin_width: int = 2,
) -> NDArray[np.int_]:
    """Assign rank-slice labels from radar-time ranks."""

    if bin_width < 1:
        raise ValueError("bin_width must be at least 1")
    ranks = np.asarray(radar_time_ranks, dtype=int)
    accessible = np.asarray(accessible_mask, dtype=bool)
    if ranks.shape != accessible.shape:
        raise ValueError("radar_time_ranks and accessible_mask must match")
    labels = np.full(ranks.shape, -1, dtype=int)
    labels[accessible] = ranks[accessible] // int(bin_width)
    return labels
