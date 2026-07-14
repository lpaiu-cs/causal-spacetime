"""Optional Numba replay of the validated 2D-order Metropolis trajectory.

This module is an analysis accelerator, not a distinct sampler. Proposal
indices and uniforms are generated in Python with the same conditional RNG
calls as :func:`two_orders.mcmc_2d_order_fast`. The compiled kernel applies
the same permutation swaps and action updates, while iterating only over the
causal rectangles affected by each swap. Tests require sampled permutations
to match the validated NumPy implementation exactly.

Numba is intentionally optional and is imported only when this module is
used by long-running replay analyses.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from .action import smearing_f2
from .two_orders import _IncrementalState, chain_observables


@njit
def _fenwick_add(tree: np.ndarray, position: int) -> None:
    position += 1
    while position < tree.size:
        tree[position] += 1
        position += position & -position


@njit
def _fenwick_query(tree: np.ndarray, position: int) -> int:
    position += 1
    total = 0
    while position > 0:
        total += tree[position]
        position -= position & -position
    return total


@njit
def _recompute_endpoint(
    permutation: np.ndarray,
    causal: np.ndarray,
    two_paths: np.ndarray,
    endpoint: int,
    tree: np.ndarray,
) -> None:
    """Recompute one row and column of C@C by permutation range counts."""
    n = permutation.size
    for index in range(n):
        two_paths[endpoint, index] = 0
        two_paths[index, endpoint] = 0
        tree[index + 1] = 0
    for later in range(endpoint + 1, n):
        if causal[endpoint, later]:
            two_paths[endpoint, later] = _fenwick_query(
                tree, permutation[later] - 1
            ) - _fenwick_query(tree, permutation[endpoint])
        _fenwick_add(tree, permutation[later])
    for index in range(1, n + 1):
        tree[index] = 0
    for earlier in range(endpoint - 1, -1, -1):
        if causal[earlier, endpoint]:
            two_paths[earlier, endpoint] = _fenwick_query(
                tree, permutation[endpoint] - 1
            ) - _fenwick_query(tree, permutation[earlier])
        _fenwick_add(tree, permutation[earlier])


@njit
def _total_weight(causal: np.ndarray, two_paths: np.ndarray, lut: np.ndarray) -> float:
    total = 0.0
    for first in range(causal.shape[0]):
        for second in range(causal.shape[1]):
            if causal[first, second]:
                total += lut[two_paths[first, second]]
    return total


@njit
def _endpoint_weight(
    causal: np.ndarray,
    two_paths: np.ndarray,
    lut: np.ndarray,
    first_endpoint: int,
    second_endpoint: int,
) -> float:
    """Weight of the union of rows/columns touching either endpoint."""
    total = 0.0
    n = causal.shape[0]
    for later in range(n):
        if causal[first_endpoint, later]:
            total += lut[two_paths[first_endpoint, later]]
        if causal[second_endpoint, later]:
            total += lut[two_paths[second_endpoint, later]]
    for earlier in range(n):
        if earlier != first_endpoint and earlier != second_endpoint:
            if causal[earlier, first_endpoint]:
                total += lut[two_paths[earlier, first_endpoint]]
            if causal[earlier, second_endpoint]:
                total += lut[two_paths[earlier, second_endpoint]]
    return total


@njit
def _update_middle_rectangle(
    causal: np.ndarray,
    two_paths: np.ndarray,
    lut: np.ndarray,
    predecessors: np.ndarray,
    successors: np.ndarray,
    delta: int,
    first_endpoint: int,
    second_endpoint: int,
    total: float,
) -> float:
    n = causal.shape[0]
    for earlier in range(n):
        if not predecessors[earlier]:
            continue
        for later in range(n):
            if not successors[later]:
                continue
            old_count = two_paths[earlier, later]
            two_paths[earlier, later] += delta
            if (
                earlier != first_endpoint
                and earlier != second_endpoint
                and later != first_endpoint
                and later != second_endpoint
                and causal[earlier, later]
            ):
                total += lut[two_paths[earlier, later]] - lut[old_count]
    return total


@njit
def _swap_state(
    permutation: np.ndarray,
    causal: np.ndarray,
    two_paths: np.ndarray,
    lut: np.ndarray,
    first_endpoint: int,
    second_endpoint: int,
    old_lines: np.ndarray,
    total: float,
    tree: np.ndarray,
) -> float:
    n = permutation.size
    for index in range(n):
        old_lines[0, index] = causal[index, first_endpoint]
        old_lines[1, index] = causal[first_endpoint, index]
        old_lines[2, index] = causal[index, second_endpoint]
        old_lines[3, index] = causal[second_endpoint, index]

    total -= _endpoint_weight(causal, two_paths, lut, first_endpoint, second_endpoint)
    total = _update_middle_rectangle(
        causal,
        two_paths,
        lut,
        old_lines[0],
        old_lines[1],
        -1,
        first_endpoint,
        second_endpoint,
        total,
    )
    total = _update_middle_rectangle(
        causal,
        two_paths,
        lut,
        old_lines[2],
        old_lines[3],
        -1,
        first_endpoint,
        second_endpoint,
        total,
    )

    temporary = permutation[first_endpoint]
    permutation[first_endpoint] = permutation[second_endpoint]
    permutation[second_endpoint] = temporary
    for later in range(n):
        causal[first_endpoint, later] = (
            later > first_endpoint and permutation[first_endpoint] < permutation[later]
        )
        causal[second_endpoint, later] = (
            later > second_endpoint
            and permutation[second_endpoint] < permutation[later]
        )
    for earlier in range(n):
        causal[earlier, first_endpoint] = (
            earlier < first_endpoint
            and permutation[earlier] < permutation[first_endpoint]
        )
        causal[earlier, second_endpoint] = (
            earlier < second_endpoint
            and permutation[earlier] < permutation[second_endpoint]
        )

    total = _update_middle_rectangle(
        causal,
        two_paths,
        lut,
        causal[:, first_endpoint],
        causal[first_endpoint, :],
        1,
        first_endpoint,
        second_endpoint,
        total,
    )
    total = _update_middle_rectangle(
        causal,
        two_paths,
        lut,
        causal[:, second_endpoint],
        causal[second_endpoint, :],
        1,
        first_endpoint,
        second_endpoint,
        total,
    )
    _recompute_endpoint(permutation, causal, two_paths, first_endpoint, tree)
    _recompute_endpoint(permutation, causal, two_paths, second_endpoint, tree)
    total += _endpoint_weight(causal, two_paths, lut, first_endpoint, second_endpoint)
    return total


@njit
def _replay_kernel(
    permutation: np.ndarray,
    causal: np.ndarray,
    two_paths: np.ndarray,
    lut: np.ndarray,
    first_indices: np.ndarray,
    second_indices: np.ndarray,
    uniforms: np.ndarray,
    beta: float,
    eps: float,
    burn: int,
    sample_every: int,
    resync_every: int,
) -> tuple[np.ndarray, int]:
    n = permutation.size
    total = _total_weight(causal, two_paths, lut)
    action = 2.0 * eps * n - 4.0 * eps * eps * total
    accepted = 0
    old_lines = np.empty((4, n), dtype=np.bool_)
    tree = np.zeros(n + 1, dtype=np.int32)
    sample_count = 0
    for step in range(first_indices.size):
        if (
            first_indices[step] != second_indices[step]
            and step >= burn
            and step % sample_every == 0
        ):
            sample_count += 1
    permutations = np.empty((sample_count, n), dtype=np.int64)

    sample_index = 0
    for step in range(first_indices.size):
        first = first_indices[step]
        second = second_indices[step]
        if first == second:
            continue
        total = _swap_state(
            permutation,
            causal,
            two_paths,
            lut,
            first,
            second,
            old_lines,
            total,
            tree,
        )
        proposed = 2.0 * eps * n - 4.0 * eps * eps * total
        if np.log(uniforms[step]) < -beta * (proposed - action):
            action = proposed
            accepted += 1
        else:
            total = _swap_state(
                permutation,
                causal,
                two_paths,
                lut,
                first,
                second,
                old_lines,
                total,
                tree,
            )
            action = 2.0 * eps * n - 4.0 * eps * eps * total
        if step % resync_every == resync_every - 1:
            total = _total_weight(causal, two_paths, lut)
            action = 2.0 * eps * n - 4.0 * eps * eps * total
        if step >= burn and step % sample_every == 0:
            permutations[sample_index] = permutation
            sample_index += 1
    return permutations, accepted


def proposal_stream(
    n_elements: int, steps: int, seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate exactly the conditional RNG stream used by the NumPy sampler."""
    rng = np.random.default_rng(seed)
    first = np.empty(steps, dtype=np.int32)
    second = np.empty(steps, dtype=np.int32)
    uniforms = np.ones(steps, dtype=np.float64)
    for step in range(steps):
        i, j = rng.integers(0, n_elements, 2)
        first[step] = i
        second[step] = j
        if i != j:
            uniforms[step] = rng.uniform()
    return first, second, uniforms


def mcmc_2d_order_replay_accelerated(
    pi0: np.ndarray,
    beta: float,
    eps: float,
    steps: int,
    seed: int,
    sample_every: int = 1000,
    burn_frac: float = 0.5,
    collect_perms: bool = False,
    resync_every: int = 1000,
) -> tuple[list[dict[str, float]], float, list[np.ndarray]]:
    """Replay ``mcmc_2d_order_fast`` with an identical proposal trajectory."""
    if sample_every < 1 or resync_every < 1:
        raise ValueError("sampling and resync intervals must be positive")
    state = _IncrementalState(pi0, eps)
    first, second, uniforms = proposal_stream(state.pi.size, steps, seed)
    permutations, accepted = _replay_kernel(
        state.pi.copy(),
        state.C.copy(),
        state.M.copy(),
        smearing_f2(np.arange(state.pi.size + 1), eps),
        first,
        second,
        uniforms,
        beta,
        eps,
        int(steps * burn_frac),
        sample_every,
        resync_every,
    )
    samples: list[dict[str, float]] = []
    retained: list[np.ndarray] = []
    for permutation in permutations:
        exact = _IncrementalState(permutation, eps)
        row = chain_observables(exact.C)
        row["S"] = exact.action()
        samples.append(row)
        if collect_perms:
            retained.append(permutation.copy())
    return samples, accepted / steps, retained
