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

from .action import abundances, smeared_action_2d


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
