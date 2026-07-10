"""Tests for the 2D-orders ensemble and its Metropolis MCMC."""

from __future__ import annotations

from itertools import permutations

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1, causal_matrix_minkowski
from causal_spacetime_lab.positive_control.action import smeared_action_2d
from causal_spacetime_lab.positive_control.two_orders import (
    bipartite_perm,
    chain_observables,
    mcmc_2d_order,
    myrheim_meyer_dimension,
    order_height,
    perm_to_causal_matrix,
)
from causal_spacetime_lab.sprinkling import sprinkle_minkowski_causal_diamond


def test_perm_special_cases():
    n = 12
    chain = perm_to_causal_matrix(np.arange(n))
    assert np.array_equal(chain, np.triu(np.ones((n, n), dtype=bool), 1))
    antichain = perm_to_causal_matrix(np.arange(n)[::-1])
    assert not antichain.any()
    bip = perm_to_causal_matrix(bipartite_perm(n))
    k = n // 2
    expected = np.zeros((n, n), dtype=bool)
    expected[:k, k:] = True
    assert np.array_equal(bip, expected)


def test_order_height():
    n = 12
    assert order_height(perm_to_causal_matrix(np.arange(n))) == n
    assert order_height(perm_to_causal_matrix(np.arange(n)[::-1])) == 1
    assert order_height(perm_to_causal_matrix(bipartite_perm(n))) == 2


def _sprinkled(n: int, dim: int, seed: int):
    events = sprinkle_minkowski_causal_diamond(n, spacetime_dim=dim, T=2.0, seed=seed)
    fn = causal_matrix_1p1 if dim == 2 else causal_matrix_minkowski
    C = np.array(fn(events), dtype=bool)
    np.fill_diagonal(C, False)
    return C


def test_myrheim_meyer_dimension_on_sprinklings():
    d2 = np.mean([myrheim_meyer_dimension(_sprinkled(200, 2, s)) for s in range(4)])
    assert abs(d2 - 2.0) < 0.25
    d3 = np.mean([myrheim_meyer_dimension(_sprinkled(200, 3, s)) for s in range(4)])
    assert abs(d3 - 3.0) < 0.35


def test_uniform_permutations_match_sprinkled_diamond():
    # Uniform 2D orders == sprinkled 2D diamond causal structure
    # (lightcone coordinates), so abundances must agree statistically.
    n = 60
    rng = np.random.default_rng(1)
    perm_n0 = np.mean(
        [chain_observables(perm_to_causal_matrix(rng.permutation(n)))["n0"] for _ in range(40)]
    )
    spr_n0 = np.mean([chain_observables(_sprinkled(n, 2, s))["n0"] for s in range(40)])
    assert abs(perm_n0 - spr_n0) / spr_n0 < 0.15


def test_mcmc_exact_gibbs_small_n():
    # Gold standard: N=5 (120 permutations), empirical distribution must
    # reproduce the exact Gibbs measure and visit every state.
    n, beta, eps = 5, 0.5, 0.3
    perms = [np.array(p) for p in permutations(range(n))]
    actions = np.array(
        [smeared_action_2d(perm_to_causal_matrix(p), eps) for p in perms]
    )
    weights = np.exp(-beta * (actions - actions.min()))
    target = weights / weights.sum()
    index = {p.tobytes(): k for k, p in enumerate(perms)}

    rng = np.random.default_rng(2)
    pi = perms[0].copy()
    action = smeared_action_2d(perm_to_causal_matrix(pi), eps)
    steps = 150_000
    counts = np.zeros(len(perms))
    for t in range(steps):
        i, j = rng.integers(0, n, 2)
        if i == j:
            continue
        pi[i], pi[j] = pi[j], pi[i]
        proposed = smeared_action_2d(perm_to_causal_matrix(pi), eps)
        if np.log(rng.uniform()) < -beta * (proposed - action):
            action = proposed
        else:
            pi[i], pi[j] = pi[j], pi[i]
        if t >= steps // 10:
            counts[index[pi.tobytes()]] += 1
    empirical = counts / counts.sum()
    assert int((counts > 0).sum()) == len(perms)
    assert 0.5 * np.abs(target - empirical).sum() < 0.06


def test_mcmc_beta0_reproduces_uniform_ensemble():
    n = 30
    rng = np.random.default_rng(3)
    samples, acceptance = mcmc_2d_order(
        rng.permutation(n), beta=0.0, eps=0.12, steps=30_000, seed=4,
        sample_every=300,
    )
    assert acceptance > 0.95
    mcmc_n0 = np.mean([s["n0"] for s in samples])
    direct_n0 = np.mean(
        [chain_observables(perm_to_causal_matrix(rng.permutation(n)))["n0"] for _ in range(150)]
    )
    assert abs(mcmc_n0 - direct_n0) / direct_n0 < 0.12
