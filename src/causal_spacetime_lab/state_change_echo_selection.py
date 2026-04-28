"""Selection helpers for same-emission echo-order diagnostics."""

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


def _validate_emission_position(chain: NDArray[np.int_], emission_position: int) -> int:
    position = int(emission_position)
    if position < 0 or position >= chain.size:
        raise IndexError("emission_position is outside the reference chain")
    return position


def emission_positions_for_reference_chain(
    reference_chain_event_ids: ArrayLike,
    strategy: str = "interior_quantiles",
    count: int = 3,
) -> NDArray[np.int_]:
    """Select reference-chain positions for fixed-emission diagnostics.

    The first and last reference-chain positions are omitted. Chains with fewer
    than three events therefore have no interior emission position.
    """

    chain = np.asarray(reference_chain_event_ids, dtype=int)
    if chain.ndim != 1:
        raise ValueError("reference_chain_event_ids must be one-dimensional")
    if count < 0:
        raise ValueError("count must be nonnegative")
    interior = np.arange(1, max(chain.size - 1, 1), dtype=int)
    if chain.size < 3 or interior.size == 0 or count == 0:
        return np.empty(0, dtype=int)
    if strategy == "all_interior":
        return interior
    if strategy == "interior_quantiles":
        selected_count = min(int(count), interior.size)
        positions = np.rint(
            np.linspace(0, interior.size - 1, selected_count)
        ).astype(int)
        return np.unique(interior[positions])
    if strategy == "early_middle_late":
        selected_count = min(3, interior.size)
        positions = np.rint(
            np.linspace(0, interior.size - 1, selected_count)
        ).astype(int)
        return np.unique(interior[positions])
    raise ValueError(
        "strategy must be one of: all_interior, interior_quantiles, "
        "early_middle_late"
    )


def select_targets_after_emission(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    emission_position: int,
    target_indices: ArrayLike | None = None,
) -> NDArray[np.int_]:
    """Return targets causally after the selected reference emission."""

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(reference_chain_event_ids, matrix.shape[0])
    emission = _validate_emission_position(chain, emission_position)
    targets = (
        np.arange(matrix.shape[0], dtype=int)
        if target_indices is None
        else _validate_indices(target_indices, matrix.shape[0])
    )
    emission_event = int(chain[emission])
    return targets[matrix[emission_event, targets]]


def select_echo_reachable_targets(
    order_matrix: ArrayLike,
    reference_chain_event_ids: ArrayLike,
    emission_position: int,
    target_indices: ArrayLike | None = None,
) -> NDArray[np.int_]:
    """Return targets after an emission with a later reference return."""

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(reference_chain_event_ids, matrix.shape[0])
    emission = _validate_emission_position(chain, emission_position)
    targets = (
        np.arange(matrix.shape[0], dtype=int)
        if target_indices is None
        else _validate_indices(target_indices, matrix.shape[0])
    )
    after = select_targets_after_emission(matrix, chain, emission, targets)
    if after.size == 0 or emission + 1 >= chain.size:
        return np.empty(0, dtype=int)
    later_chain = chain[emission + 1 :]
    has_return = np.any(matrix[after[:, None], later_chain[None, :]], axis=1)
    return after[has_return]
