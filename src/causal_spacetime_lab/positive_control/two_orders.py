"""2D orders (restricted causal-set ensemble) and Metropolis MCMC over them.

A 2D order on N labelled elements is the intersection of two total orders.
Canonical representation: a permutation ``pi``; element i sits at lightcone
coordinates (i, pi[i]) and ``i < j`` iff ``i < j`` and ``pi[i] < pi[j]``.
The uniform measure over permutations is exactly the causal structure of a
Poisson sprinkling (conditioned on N) into a flat 2D causal diamond in
lightcone coordinates -- so the beta = 0 ensemble IS sprinkled-2D-like.
The restriction does NOT exclude crystalline orders geometrically (the
complete bipartite order is the permutation (k..1, 2k..k+1)); it suppresses
them entropically only, so a continuum/crystal competition in beta remains.

MCMC: Metropolis over permutations with random-transposition proposals
(symmetric), stationary for ``exp(-beta * S_eps)`` with the smeared BD
action of :mod:`.action`. Validated in tests against exact enumeration at
small N (Gibbs distribution reproduced, all states visited).
"""

from __future__ import annotations

from math import gamma

import numpy as np

from .action import abundances, smeared_action_2d, smearing_f2


def perm_to_causal_matrix(pi: np.ndarray) -> np.ndarray:
    """Strict causal matrix of the 2D order represented by permutation pi."""
    idx = np.arange(pi.size)
    return (idx[:, None] < idx[None, :]) & (pi[:, None] < pi[None, :])


def bipartite_perm(n_elements: int) -> np.ndarray:
    """Permutation realizing the complete bipartite (2-layer crystal) order."""
    k = n_elements // 2
    return np.array(
        list(range(k - 1, -1, -1)) + list(range(n_elements - 1, k - 1, -1)),
        dtype=int,
    )


def balanced_layered_perm(
    n_elements: int,
    layer_count: int,
    seed: int,
    min_layer_size: int = 6,
) -> np.ndarray:
    """Return a randomized complete layered 2D order.

    Positions and values are partitioned into matching contiguous blocks.
    Reversing values inside each block makes elements within a layer
    incomparable, while every element in an earlier layer precedes every
    element in a later layer. Random layer sizes prevent the construction
    from depending on one specially balanced partition.
    """
    if layer_count < 2:
        raise ValueError("layer_count must be at least 2")
    if min_layer_size < 1:
        raise ValueError("min_layer_size must be positive")
    required = layer_count * min_layer_size
    if required > n_elements:
        raise ValueError("minimum layer sizes exceed n_elements")

    rng = np.random.default_rng(seed)
    sizes = np.full(layer_count, min_layer_size, dtype=int)
    sizes += rng.multinomial(
        n_elements - required, np.full(layer_count, 1 / layer_count)
    )
    parts: list[np.ndarray] = []
    start = 0
    for size in sizes:
        stop = start + int(size)
        parts.append(np.arange(start, stop, dtype=int)[::-1])
        start = stop
    return np.concatenate(parts)


def windowed_transpositions(
    permutation: np.ndarray,
    moves: int,
    window: int,
    seed: int,
) -> np.ndarray:
    """Apply random transpositions whose positions are at most ``window`` apart."""
    pi = np.asarray(permutation, dtype=int)
    if moves < 0:
        raise ValueError("moves must be non-negative")
    if window < 1:
        raise ValueError("window must be positive")
    if pi.ndim != 1 or not np.array_equal(np.sort(pi), np.arange(pi.size)):
        raise ValueError("permutation must contain each integer in [0, N) exactly once")

    result = pi.copy()
    if result.size < 2:
        return result
    rng = np.random.default_rng(seed)
    for _ in range(moves):
        first = int(rng.integers(0, result.size))
        low = max(0, first - window)
        high = min(result.size, first + window + 1)
        second = int(rng.integers(low, high))
        if first != second:
            result[first], result[second] = result[second], result[first]
    return result


def order_height(causal_matrix: np.ndarray) -> int:
    """Length (number of elements) of the longest chain.

    Works for any index labelling: ancestor counts strictly increase along
    relations, so sorting by them yields a valid topological order.
    """
    n = causal_matrix.shape[0]
    heights = np.ones(n, dtype=int)
    for j in np.argsort(causal_matrix.sum(axis=0)):
        below = np.where(causal_matrix[:, j])[0]
        if below.size:
            heights[j] = 1 + heights[below].max()
    return int(heights.max())


def myrheim_meyer_dimension(causal_matrix: np.ndarray) -> float:
    """Myrheim-Meyer dimension from the ordering fraction.

    Inverts f(d) = 1.5 Gamma(d/2+1) Gamma(d+1) / Gamma(3d/2+1)
    (f(1) = 1, f(2) = 0.5). Validated on sprinkled 2D/3D diamonds.
    NOTE (E1/E2 lesson): crystalline bipartite orders also have f ~ 0.5,
    i.e. MM dimension ~ 2 -- never use this alone as a manifoldlikeness
    judge; combine with abundances (n1, n2) and height.
    """
    n = causal_matrix.shape[0]
    frac = causal_matrix.sum() / (n * (n - 1) / 2)
    if frac <= 1e-9 or frac >= 1 - 1e-9:
        return float("nan")

    def g(d: float) -> float:
        return 1.5 * gamma(d / 2 + 1) * gamma(d + 1) / gamma(3 * d / 2 + 1) - frac

    lo, hi = 0.3, 12.0
    if g(lo) * g(hi) > 0:
        return float("nan")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if g(lo) * g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def chain_observables(causal_matrix: np.ndarray) -> dict[str, float]:
    relations, (n0, n1, n2) = abundances(causal_matrix, kmax=2)
    return {
        "R": float(relations),
        "n0": float(n0),
        "n1": float(n1),
        "n2": float(n2),
        "mm_dim": myrheim_meyer_dimension(causal_matrix),
        "height": float(order_height(causal_matrix)),
    }


class _IncrementalState:
    """C, M = C@C and T = sum f2 over related pairs, updated in O(N^2)/move.

    A value swap at positions (a, b) only changes relations with endpoint a
    or b. Middle-element contributions of a and b to M are rank-1 outer
    products; rows/columns a and b of M are recomputed exactly afterwards.
    """

    def __init__(self, pi: np.ndarray, eps: float):
        self.eps = eps
        self.pi = np.array(pi, dtype=int).copy()
        n = self.pi.size
        self.lut = smearing_f2(np.arange(n + 1), eps)
        self.C = perm_to_causal_matrix(self.pi)
        self.Cf = self.C.astype(np.float32)
        self.M = (self.Cf @ self.Cf).astype(np.int32)
        self._recompute_T()

    def _recompute_T(self) -> None:
        self.T = float(self.lut[self.M[self.C]].sum())

    def action(self) -> float:
        return 2.0 * self.eps * self.pi.size - 4.0 * self.eps * self.eps * self.T

    def swap(self, a: int, b: int) -> None:
        pi, C, M = self.pi, self.C, self.M
        idx = np.arange(pi.size)
        # remove old middle-element contributions of a and b
        M -= np.outer(C[:, a], C[a, :]).astype(np.int32)
        M -= np.outer(C[:, b], C[b, :]).astype(np.int32)
        pi[a], pi[b] = pi[b], pi[a]
        # rewrite the four changed lines of C
        C[a, :] = (idx > a) & (pi[a] < pi)
        C[:, a] = (idx < a) & (pi < pi[a])
        C[b, :] = (idx > b) & (pi[b] < pi)
        C[:, b] = (idx < b) & (pi < pi[b])
        self.Cf = C.astype(np.float32)
        # add new middle-element contributions
        M += np.outer(C[:, a], C[a, :]).astype(np.int32)
        M += np.outer(C[:, b], C[b, :]).astype(np.int32)
        # endpoint rows/columns of a and b: recompute exactly
        for r in (a, b):
            M[r, :] = (self.Cf[r, :] @ self.Cf).astype(np.int32)
            M[:, r] = (self.Cf @ self.Cf[:, r]).astype(np.int32)
        self._recompute_T()


def mcmc_2d_order_fast(
    pi0: np.ndarray,
    beta: float,
    eps: float,
    steps: int,
    seed: int,
    sample_every: int = 1000,
    burn_frac: float = 0.5,
    collect_perms: bool = False,
) -> tuple[list[dict[str, float]], float, list[np.ndarray]]:
    """Same chain as :func:`mcmc_2d_order` (identical RNG stream and
    trajectory) with O(N^2)/move incremental action updates, usable at
    N of several hundred. Rejected proposals are undone by re-swapping.

    Returns (samples, acceptance, perms) where perms holds sampled
    permutations if ``collect_perms`` (for downstream discriminator runs).
    """
    rng = np.random.default_rng(seed)
    state = _IncrementalState(pi0, eps)
    action = state.action()
    burn = int(steps * burn_frac)
    samples: list[dict[str, float]] = []
    perms: list[np.ndarray] = []
    accepted = 0
    for t in range(steps):
        i, j = rng.integers(0, state.pi.size, 2)
        if i == j:
            continue
        state.swap(i, j)
        proposed = state.action()
        if np.log(rng.uniform()) < -beta * (proposed - action):
            action = proposed
            accepted += 1
        else:
            state.swap(i, j)
        if t >= burn and t % sample_every == 0:
            obs = chain_observables(state.C)
            obs["S"] = action
            samples.append(obs)
            if collect_perms:
                perms.append(state.pi.copy())
    return samples, accepted / steps, perms


def mcmc_2d_order(
    pi0: np.ndarray,
    beta: float,
    eps: float,
    steps: int,
    seed: int,
    sample_every: int = 1000,
    burn_frac: float = 0.5,
) -> tuple[list[dict[str, float]], float]:
    """Metropolis chain over 2D orders; returns (samples, acceptance rate).

    Each sample dict holds S plus :func:`chain_observables` fields, taken
    every ``sample_every`` steps after a ``burn_frac`` burn-in.
    """
    rng = np.random.default_rng(seed)
    pi = np.array(pi0, dtype=int).copy()
    action = smeared_action_2d(perm_to_causal_matrix(pi), eps)
    burn = int(steps * burn_frac)
    samples: list[dict[str, float]] = []
    accepted = 0
    for t in range(steps):
        i, j = rng.integers(0, pi.size, 2)
        if i == j:
            continue
        pi[i], pi[j] = pi[j], pi[i]
        proposed = smeared_action_2d(perm_to_causal_matrix(pi), eps)
        if np.log(rng.uniform()) < -beta * (proposed - action):
            action = proposed
            accepted += 1
        else:
            pi[i], pi[j] = pi[j], pi[i]
        if t >= burn and t % sample_every == 0:
            obs = chain_observables(perm_to_causal_matrix(pi))
            obs["S"] = action
            samples.append(obs)
    return samples, accepted / steps
