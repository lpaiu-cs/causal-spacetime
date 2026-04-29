"""Scalar ordinal representability diagnostics for response-order signs."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change_order import topological_order_from_adjacency
from causal_spacetime_lab.state_change_response_signature_comparison import (
    response_order_cycle_count,
)


def _validate_sign_matrix(order_sign_matrix: ArrayLike) -> NDArray[np.int_]:
    signs = np.asarray(order_sign_matrix, dtype=int)
    if signs.ndim != 2 or signs.shape[0] != signs.shape[1]:
        raise ValueError("order_sign_matrix must be square")
    return signs


def response_order_directed_edges(order_sign_matrix: ArrayLike) -> NDArray[np.int_]:
    """Return directed edges for nonzero response-order constraints.

    Sign -1 means target i has smaller response rank than target j, so the
    directed edge is i -> j. Sign 1 contributes the reverse edge.
    """

    signs = _validate_sign_matrix(order_sign_matrix)
    edges: set[tuple[int, int]] = set()
    n_targets = signs.shape[0]
    for i in range(n_targets):
        for j in range(i + 1, n_targets):
            if signs[i, j] < 0:
                edges.add((i, j))
            elif signs[i, j] > 0:
                edges.add((j, i))
    return np.asarray(sorted(edges), dtype=int).reshape((-1, 2))


def _adjacency_from_signs(order_sign_matrix: ArrayLike) -> NDArray[np.bool_]:
    signs = _validate_sign_matrix(order_sign_matrix)
    adjacency = np.zeros(signs.shape, dtype=bool)
    for source, target in response_order_directed_edges(signs):
        adjacency[int(source), int(target)] = True
    return adjacency


def has_response_order_cycle(order_sign_matrix: ArrayLike) -> bool:
    """Return whether the nonzero response-order relation has any cycle."""

    adjacency = _adjacency_from_signs(order_sign_matrix)
    try:
        topological_order_from_adjacency(adjacency)
    except ValueError:
        return True
    return False


def response_order_topological_ranks(order_sign_matrix: ArrayLike) -> NDArray[np.int_]:
    """Return topological scalar ranks when response order is acyclic.

    Tied or unresolved pairs impose no scalar-rank constraint. This is scalar
    ordinal representability, not metric representability.
    """

    adjacency = _adjacency_from_signs(order_sign_matrix)
    topo = topological_order_from_adjacency(adjacency)
    ranks = np.zeros(adjacency.shape[0], dtype=int)
    for node in topo:
        predecessors = np.flatnonzero(adjacency[:, int(node)])
        if predecessors.size:
            ranks[int(node)] = int(np.max(ranks[predecessors]) + 1)
    return ranks


def scalar_rank_representation_error(
    order_sign_matrix: ArrayLike,
    ranks: ArrayLike,
) -> dict[str, float]:
    """Report violations of response-order signs by scalar ranks."""

    signs = _validate_sign_matrix(order_sign_matrix)
    scalar = np.asarray(ranks, dtype=float)
    if scalar.ndim != 1 or scalar.size != signs.shape[0]:
        raise ValueError("ranks must be a vector matching the sign matrix")
    comparable = 0
    violations = 0
    for i in range(signs.shape[0] - 1):
        for j in range(i + 1, signs.shape[0]):
            sign = int(signs[i, j])
            if sign == 0:
                continue
            comparable += 1
            if sign < 0 and not scalar[i] < scalar[j]:
                violations += 1
            if sign > 0 and not scalar[i] > scalar[j]:
                violations += 1
    return {
        "comparable_pair_count": float(comparable),
        "violation_count": float(violations),
        "violation_fraction": float(violations / comparable) if comparable else 0.0,
    }


def scalar_representability_report(order_sign_matrix: ArrayLike) -> dict[str, float]:
    """Report scalar ordinal representability preconditions."""

    signs = _validate_sign_matrix(order_sign_matrix)
    target_count = signs.shape[0]
    total_pairs = target_count * (target_count - 1) // 2
    nonzero = int(
        sum(
            signs[i, j] != 0
            for i in range(max(0, target_count - 1))
            for j in range(i + 1, target_count)
        )
    )
    has_cycle = has_response_order_cycle(signs)
    if has_cycle:
        scalar_representable = 0.0
        rank_span = float("nan")
    else:
        ranks = response_order_topological_ranks(signs)
        scalar_representable = 1.0
        rank_span = float(np.max(ranks) - np.min(ranks)) if ranks.size else 0.0
    return {
        "target_count": float(target_count),
        "nonzero_pair_count": float(nonzero),
        "tie_or_unresolved_pair_count": float(total_pairs - nonzero),
        "has_cycle": float(has_cycle),
        "directed_3cycle_count": float(response_order_cycle_count(signs)),
        "scalar_representable": scalar_representable,
        "rank_span": rank_span,
    }

