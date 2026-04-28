"""Reference-chain candidates in finite state-change trigger networks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.state_change import StateChangeNetwork
from causal_spacetime_lab.state_change_order import topological_order_from_adjacency

REFERENCE_CANDIDATE_SOURCES = {
    "local_system",
    "greedy_order",
    "longest_chain",
    "manual",
    "random_order",
}

OBSERVER_CANDIDATE_SOURCES = REFERENCE_CANDIDATE_SOURCES


@dataclass(frozen=True)
class ReferenceChainCandidate:
    """A finite reference-chain candidate.

    A reference chain is a candidate reference backbone for order-level
    protocol diagnostics. It is not a calibrated clock, metric observer, or
    uniquely selected physical observer.
    """

    name: str
    chain_event_ids: NDArray[np.int_]
    source: str

    def __post_init__(self) -> None:
        if self.source not in REFERENCE_CANDIDATE_SOURCES:
            allowed = ", ".join(sorted(REFERENCE_CANDIDATE_SOURCES))
            raise ValueError(f"source must be one of: {allowed}")
        object.__setattr__(
            self,
            "chain_event_ids",
            np.asarray(self.chain_event_ids, dtype=int),
        )


ObserverChainCandidate = ReferenceChainCandidate


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


def is_chain(order_matrix: ArrayLike, chain_event_ids: ArrayLike) -> bool:
    """Return whether every earlier chain element precedes every later one."""

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(chain_event_ids, matrix.shape[0])
    for i, left in enumerate(chain):
        for right in chain[i + 1 :]:
            if not matrix[int(left), int(right)]:
                return False
    return True


def chain_comparability_mask(
    order_matrix: ArrayLike,
    chain_event_ids: ArrayLike,
    target_indices: ArrayLike | None = None,
) -> NDArray[np.bool_]:
    """Return whether each target is comparable to the chain or on it."""

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(chain_event_ids, matrix.shape[0])
    targets = (
        np.arange(matrix.shape[0], dtype=int)
        if target_indices is None
        else _validate_indices(target_indices, matrix.shape[0])
    )
    chain_set = set(chain.tolist())
    mask = np.zeros(targets.shape[0], dtype=bool)
    for index, target in enumerate(targets):
        if int(target) in chain_set:
            mask[index] = True
        elif chain.size:
            mask[index] = bool(
                np.any(matrix[chain, int(target)])
                or np.any(matrix[int(target), chain])
            )
    return mask


def chain_bracketing_mask(
    order_matrix: ArrayLike,
    chain_event_ids: ArrayLike,
    target_indices: ArrayLike | None = None,
) -> NDArray[np.bool_]:
    """Return whether each non-chain target has chain events before and after."""

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(chain_event_ids, matrix.shape[0])
    targets = (
        np.arange(matrix.shape[0], dtype=int)
        if target_indices is None
        else _validate_indices(target_indices, matrix.shape[0])
    )
    chain_set = set(chain.tolist())
    mask = np.zeros(targets.shape[0], dtype=bool)
    for index, target in enumerate(targets):
        if int(target) in chain_set or not chain.size:
            continue
        has_predecessor = bool(np.any(matrix[chain, int(target)]))
        has_successor = bool(np.any(matrix[int(target), chain]))
        mask[index] = has_predecessor and has_successor
    return mask


def latest_predecessor_chain_position(
    order_matrix: ArrayLike,
    chain_event_ids: ArrayLike,
    target_index: int,
) -> int | None:
    """Return latest reference-chain position preceding a target, if any."""

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(chain_event_ids, matrix.shape[0])
    _validate_indices(np.asarray([target_index]), matrix.shape[0])
    positions = np.flatnonzero(matrix[chain, int(target_index)])
    return int(positions[-1]) if positions.size else None


def earliest_successor_chain_position(
    order_matrix: ArrayLike,
    chain_event_ids: ArrayLike,
    target_index: int,
) -> int | None:
    """Return earliest reference-chain position succeeding a target, if any."""

    matrix = _validate_order_matrix(order_matrix)
    chain = _validate_indices(chain_event_ids, matrix.shape[0])
    _validate_indices(np.asarray([target_index]), matrix.shape[0])
    positions = np.flatnonzero(matrix[int(target_index), chain])
    return int(positions[0]) if positions.size else None


def local_system_chain_candidates(
    network: StateChangeNetwork,
    min_length: int = 2,
) -> list[ObserverChainCandidate]:
    """Return local-system chains as reference-chain candidates."""

    if min_length < 1:
        raise ValueError("min_length must be positive")
    candidates: list[ReferenceChainCandidate] = []
    for system_id, event_ids in sorted(network.system_event_ids.items()):
        if len(event_ids) >= min_length:
            candidates.append(
                ReferenceChainCandidate(
                    name=f"local_system_{system_id}",
                    chain_event_ids=np.asarray(event_ids, dtype=int),
                    source="local_system",
                )
            )
    return candidates


def greedy_chain_candidate_from_order(
    order_matrix: ArrayLike,
    *,
    start_event_id: int | None = None,
    name: str = "greedy_order_chain",
) -> ObserverChainCandidate:
    """Build a heuristic reference chain using high-reach successors."""

    matrix = _validate_order_matrix(order_matrix)
    n_events = matrix.shape[0]
    if n_events == 0:
        chain = np.empty(0, dtype=int)
    else:
        future_reach = np.sum(matrix, axis=1)
        past_reach = np.sum(matrix, axis=0)
        if start_event_id is None:
            min_past = np.min(past_reach)
            starts = np.flatnonzero(past_reach == min_past)
            current = int(starts[np.argmax(future_reach[starts])])
        else:
            _validate_indices(np.asarray([start_event_id]), n_events)
            current = int(start_event_id)
        chain_values = [current]
        while True:
            successors = np.flatnonzero(matrix[current])
            if not successors.size:
                break
            scores = future_reach[successors]
            current = int(successors[np.argmax(scores)])
            chain_values.append(current)
        chain = np.asarray(chain_values, dtype=int)
    return ReferenceChainCandidate(
        name=name,
        chain_event_ids=chain,
        source="greedy_order",
    )


def longest_chain_candidate_from_order(
    order_matrix: ArrayLike,
    name: str = "longest_order_chain",
) -> ObserverChainCandidate:
    """Return a longest reference chain using the transitive order DAG."""

    matrix = _validate_order_matrix(order_matrix)
    n_events = matrix.shape[0]
    if n_events == 0:
        chain = np.empty(0, dtype=int)
    else:
        topo = topological_order_from_adjacency(matrix)
        lengths = np.ones(n_events, dtype=int)
        next_node = np.full(n_events, -1, dtype=int)
        for node in reversed(topo):
            successors = np.flatnonzero(matrix[int(node)])
            if successors.size:
                best = int(successors[np.argmax(lengths[successors])])
                lengths[int(node)] = 1 + lengths[best]
                next_node[int(node)] = best
        current = int(np.argmax(lengths))
        values = [current]
        while next_node[current] >= 0:
            current = int(next_node[current])
            values.append(current)
        chain = np.asarray(values, dtype=int)
    return ReferenceChainCandidate(
        name=name,
        chain_event_ids=chain,
        source="longest_chain",
    )


def random_chain_candidate_from_order(
    order_matrix: ArrayLike,
    seed: int | None = None,
    name: str = "random_order_chain",
) -> ObserverChainCandidate:
    """Return a random valid reference chain by walking through successors."""

    matrix = _validate_order_matrix(order_matrix)
    n_events = matrix.shape[0]
    rng = np.random.default_rng(seed)
    if n_events == 0:
        chain = np.empty(0, dtype=int)
    else:
        current = int(rng.integers(0, n_events))
        values = [current]
        while True:
            successors = np.flatnonzero(matrix[current])
            if not successors.size:
                break
            current = int(rng.choice(successors))
            values.append(current)
        chain = np.asarray(values, dtype=int)
    return ReferenceChainCandidate(
        name=name,
        chain_event_ids=chain,
        source="random_order",
    )


local_system_reference_chain_candidates = local_system_chain_candidates
greedy_reference_chain_candidate_from_order = greedy_chain_candidate_from_order
longest_reference_chain_candidate_from_order = longest_chain_candidate_from_order
random_reference_chain_candidate_from_order = random_chain_candidate_from_order
