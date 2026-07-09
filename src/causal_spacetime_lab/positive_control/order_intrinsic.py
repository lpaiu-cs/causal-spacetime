"""Order-intrinsic reference chains and profiles (no coordinates).

For a causal order with no embedding (a dynamics-generated causal set, P3),
reference chains cannot be supplied geometrically; they are selected from the
order itself as long disjoint chains, and their tick labels are positions along
the chain. Targets are elements two-sided bracketed by every selected chain.
The bracket-width profiles are then measured exactly as in PC-V1 and consumed
by the frozen dissimilarity/fit/gate pipeline.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.discrete_radar import find_radar_ticks_from_order
from causal_spacetime_lab.positive_control.echo_profiles import EchoProfileMatrix


def global_longest_chain(
    causal: NDArray[np.bool_],
    topo: NDArray[np.int_],
    allowed: NDArray[np.bool_],
) -> NDArray[np.int_]:
    """Longest chain among ``allowed`` nodes via topological longest-path DP."""

    n = causal.shape[0]
    depth = np.where(allowed, 1, 0).astype(int)
    pred = np.full(n, -1, dtype=int)
    for u in topo:
        if not allowed[u]:
            continue
        succ = np.flatnonzero(causal[u] & allowed)
        improve = succ[depth[u] + 1 > depth[succ]]
        depth[improve] = depth[u] + 1
        pred[improve] = u
    if depth.max() < 1:
        return np.array([], dtype=int)
    end = int(np.argmax(depth))
    chain: list[int] = []
    while end != -1:
        chain.append(end)
        end = int(pred[end])
    return np.array(chain[::-1], dtype=int)


def select_disjoint_chains(
    causal: NDArray[np.bool_],
    times: NDArray[np.float64],
    chain_count: int,
    min_length: int,
) -> list[NDArray[np.int_]]:
    """Greedily extract ``chain_count`` disjoint long chains from the order."""

    topo = np.argsort(times)
    allowed = np.ones(causal.shape[0], dtype=bool)
    chains: list[NDArray[np.int_]] = []
    for _ in range(chain_count):
        chain = global_longest_chain(causal, topo, allowed)
        if chain.size < min_length:
            break
        chains.append(chain)
        allowed[chain] = False
    return chains


def select_bracketed_targets(
    causal: NDArray[np.bool_],
    chains: list[NDArray[np.int_]],
    max_targets: int,
    seed: int,
) -> NDArray[np.int_]:
    """Return elements (outside the chains) two-sided bracketed by every chain."""

    used = np.concatenate(chains) if chains else np.array([], dtype=int)
    candidates = np.setdiff1d(np.arange(causal.shape[0]), used)
    eligible = [
        int(t)
        for t in candidates
        if all(
            find_radar_ticks_from_order(
                causal, chain, int(t), np.arange(chain.size, dtype=float)
            )
            is not None
            for chain in chains
        )
    ]
    eligible = np.array(eligible, dtype=int)
    if eligible.size > max_targets:
        rng = np.random.default_rng(seed)
        eligible = np.sort(rng.choice(eligible, max_targets, replace=False))
    return eligible


def measure_order_intrinsic_profiles(
    causal: NDArray[np.bool_],
    chains: list[NDArray[np.int_]],
    targets: NDArray[np.int_],
) -> EchoProfileMatrix:
    """Bracket-width profiles using per-chain position labels as clocks."""

    n, k = targets.size, len(chains)
    delays = np.full((n, k), np.nan, dtype=np.float64)
    reachable = np.zeros((n, k), dtype=bool)
    for row, target in enumerate(targets):
        for col, chain in enumerate(chains):
            ticks = find_radar_ticks_from_order(
                causal, chain, int(target), np.arange(chain.size, dtype=float)
            )
            if ticks is not None:
                delays[row, col] = ticks[1] - ticks[0]
                reachable[row, col] = True
    return EchoProfileMatrix(delays, reachable, targets.copy())
