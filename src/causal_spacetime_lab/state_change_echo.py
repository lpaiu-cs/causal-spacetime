"""Same-emission echo-order utilities for state-change trigger networks."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.ordinal import order_agreement_rate, order_inversion_rate


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


def echo_return_position_for_emission(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    target_index: int,
    emission_position: int,
) -> int | None:
    """Return the first reference-chain return position for one emission.

    The return position is defined only when the selected reference event
    precedes the target and the target precedes a later reference event.
    """

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(reference_chain_event_ids, matrix.shape[0])
    target = int(_validate_indices(np.asarray([target_index]), matrix.shape[0])[0])
    emission = _validate_emission_position(chain, emission_position)
    emission_event = int(chain[emission])
    if not matrix[emission_event, target]:
        return None
    later_positions = np.arange(emission + 1, chain.size, dtype=int)
    if later_positions.size == 0:
        return None
    returning = later_positions[matrix[target, chain[later_positions]]]
    return int(returning[0]) if returning.size else None


def echo_delay_rank_for_emission(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    target_index: int,
    emission_position: int,
) -> int | None:
    """Return order-level echo-delay rank, not a metric distance."""

    return_position = echo_return_position_for_emission(
        order_matrix,
        reference_chain_event_ids,
        target_index,
        emission_position,
    )
    if return_position is None:
        return None
    return int(return_position - int(emission_position))


def echo_reachability_from_emission(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    emission_position: int,
    target_indices: ArrayLike | None = None,
    include_reference_events: bool = False,
) -> tuple[NDArray[np.int_], NDArray[np.int_], NDArray[np.bool_]]:
    """Return echo-return positions and delay ranks for target events.

    By default, events on the reference chain are marked unreachable because
    they are part of the chosen reference backbone rather than echo targets.
    """

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(reference_chain_event_ids, matrix.shape[0])
    emission = _validate_emission_position(chain, emission_position)
    targets = (
        np.arange(matrix.shape[0], dtype=int)
        if target_indices is None
        else _validate_indices(target_indices, matrix.shape[0])
    )
    reference_set = set(chain.tolist())
    return_positions = np.full(targets.size, -1, dtype=int)
    delay_ranks = np.full(targets.size, -1, dtype=int)
    reachable_mask = np.zeros(targets.size, dtype=bool)
    for row_index, target in enumerate(targets):
        if not include_reference_events and int(target) in reference_set:
            continue
        return_position = echo_return_position_for_emission(
            matrix,
            chain,
            int(target),
            emission,
        )
        if return_position is None:
            continue
        return_positions[row_index] = return_position
        delay_ranks[row_index] = return_position - emission
        reachable_mask[row_index] = True
    return return_positions, delay_ranks, reachable_mask


def echo_order_matrix_from_delay_ranks(
    delay_ranks: ArrayLike,
    reachable_mask: ArrayLike,
    strict: bool = True,
) -> NDArray[np.bool_]:
    """Return the pairwise echo-delay order induced by reachable targets."""

    delays = np.asarray(delay_ranks, dtype=int)
    reachable = np.asarray(reachable_mask, dtype=bool)
    if delays.ndim != 1 or reachable.ndim != 1 or delays.shape != reachable.shape:
        raise ValueError("delay_ranks and reachable_mask must be matching vectors")
    if strict:
        matrix = delays[:, None] < delays[None, :]
    else:
        matrix = delays[:, None] <= delays[None, :]
    return matrix & reachable[:, None] & reachable[None, :]


def compare_echo_orders(
    delay_ranks_a: ArrayLike,
    delay_ranks_b: ArrayLike,
    reachable_a: ArrayLike,
    reachable_b: ArrayLike,
) -> dict[str, float]:
    """Compare echo-delay orderings on targets reachable in both protocols."""

    delays_a = np.asarray(delay_ranks_a, dtype=float)
    delays_b = np.asarray(delay_ranks_b, dtype=float)
    mask_a = np.asarray(reachable_a, dtype=bool)
    mask_b = np.asarray(reachable_b, dtype=bool)
    if (
        delays_a.shape != delays_b.shape
        or delays_a.shape != mask_a.shape
        or delays_a.shape != mask_b.shape
    ):
        raise ValueError("delay and reachability arrays must match")
    common = mask_a & mask_b
    common_count = int(np.count_nonzero(common))
    if common_count < 2:
        inversion = float("nan")
        agreement = float("nan")
    else:
        inversion = order_inversion_rate(
            delays_a[common],
            delays_b[common],
            ignore_ties=True,
        )
        agreement = order_agreement_rate(
            delays_a[common],
            delays_b[common],
            ignore_ties=True,
        )
    return {
        "common_reachable_count": float(common_count),
        "common_reachable_fraction": float(common_count / delays_a.size)
        if delays_a.size
        else 0.0,
        "order_inversion_rate": inversion,
        "order_agreement_rate": agreement,
    }
