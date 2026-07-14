"""Tests for the 2D-orders ensemble and its Metropolis MCMC."""

from __future__ import annotations

from itertools import permutations

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1, causal_matrix_minkowski
from causal_spacetime_lab.positive_control.action import smeared_action_2d
from causal_spacetime_lab.positive_control.two_orders import (
    _IncrementalState,
    balanced_layered_perm,
    bipartite_perm,
    chain_observables,
    mcmc_2d_order,
    mcmc_2d_order_fast,
    myrheim_meyer_dimension,
    order_height,
    perm_to_causal_matrix,
    windowed_transpositions,
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


def test_balanced_layered_permutation_is_seeded_and_has_requested_height():
    first = balanced_layered_perm(120, layer_count=12, seed=3, min_layer_size=6)
    repeat = balanced_layered_perm(120, layer_count=12, seed=3, min_layer_size=6)
    other = balanced_layered_perm(120, layer_count=12, seed=4, min_layer_size=6)
    assert np.array_equal(np.sort(first), np.arange(120))
    assert np.array_equal(first, repeat)
    assert not np.array_equal(first, other)
    assert order_height(perm_to_causal_matrix(first)) == 12


def test_windowed_transpositions_are_seeded_and_preserve_permutation():
    original = np.arange(80)
    first = windowed_transpositions(original, moves=30, window=5, seed=7)
    repeat = windowed_transpositions(original, moves=30, window=5, seed=7)
    assert np.array_equal(np.sort(first), original)
    assert np.array_equal(first, repeat)
    assert not np.array_equal(first, original)

    one_move = windowed_transpositions(original, moves=1, window=5, seed=8)
    changed = np.flatnonzero(one_move != original)
    assert changed.size in (0, 2)
    if changed.size:
        assert int(np.ptp(changed)) <= 5


def test_order_height_index_order_invariant():
    # Height must not depend on element labelling (sprinkled matrices are
    # not time-sorted). Relabel a chain adversarially and check.
    n = 12
    rng = np.random.default_rng(5)
    C = perm_to_causal_matrix(np.arange(n))
    relabel = rng.permutation(n)
    shuffled = C[np.ix_(relabel, relabel)]
    assert order_height(shuffled) == n
    for seed in range(3):
        C = _sprinkled(80, 2, seed)
        relabel = rng.permutation(80)
        assert order_height(C[np.ix_(relabel, relabel)]) == order_height(C)


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


def test_incremental_state_matches_full_recompute():
    rng = np.random.default_rng(6)
    state = _IncrementalState(rng.permutation(60), eps=0.12)
    for _ in range(300):
        i, j = rng.integers(0, 60, 2)
        if i == j:
            continue
        state.swap(i, j)
    Cf = state.C.astype(np.float32)
    assert np.array_equal(state.M, (Cf @ Cf).astype(np.int32))
    assert abs(state.T - float(state.lut[state.M[state.C]].sum())) < 1e-6


def test_fast_sampler_identical_trajectory():
    # Same seed => same RNG stream => the fast sampler must reproduce the
    # reference sampler's trajectory exactly.
    rng = np.random.default_rng(7)
    pi0 = rng.permutation(30)
    slow, acc_slow = mcmc_2d_order(
        pi0, beta=0.8, eps=0.2, steps=8_000, seed=8, sample_every=200
    )
    fast, acc_fast, _ = mcmc_2d_order_fast(
        pi0, beta=0.8, eps=0.2, steps=8_000, seed=8, sample_every=200
    )
    assert acc_slow == acc_fast
    assert len(slow) == len(fast)
    for s, f in zip(slow, fast):
        assert abs(s["S"] - f["S"]) < 1e-6
        assert s["n0"] == f["n0"] and s["height"] == f["height"]


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
