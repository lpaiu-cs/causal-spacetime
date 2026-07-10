"""Tests for the verified 2D Benincasa-Dowker action (raw and smeared)."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.positive_control.action import (
    bd_action_2d,
    smeared_action_2d,
)
from causal_spacetime_lab.positive_control.rewire import transitive_closure
from causal_spacetime_lab.sprinkling import sprinkle_minkowski_causal_diamond


def _random_closed_poset(n: int, density: float, rng: np.random.Generator):
    upper = np.triu(rng.uniform(size=(n, n)) < density, 1)
    closed = transitive_closure(upper)
    np.fill_diagonal(closed, False)
    return closed


def test_smeared_eps1_reduces_to_raw_exactly():
    rng = np.random.default_rng(0)
    for _ in range(20):
        C = _random_closed_poset(10, rng.uniform(0.05, 0.5), rng)
        assert smeared_action_2d(C, 1.0) == bd_action_2d(C)


def test_bipartite_closed_form():
    # Complete bipartite N/2 x N/2: every relation is a link, so
    # S_eps = 2 eps N - 4 eps^2 (N^2/4) = eps N (2 - eps N) exactly.
    for n in [10, 18, 100]:
        k = n // 2
        C = np.zeros((n, n), dtype=bool)
        C[:k, k : 2 * k] = True
        for eps in [1.0, 0.5, 0.12, 0.04]:
            expected = eps * n * (2.0 - eps * n)
            assert abs(smeared_action_2d(C, eps) - expected) < 1e-9


def test_chain_and_antichain_raw_action():
    n = 40
    chain = np.triu(np.ones((n, n), dtype=bool), 1)
    antichain = np.zeros((n, n), dtype=bool)
    # chain: n0 = n-1, n1 = n-2, n2 = n-3 -> S = 2n; antichain: S = 2n.
    assert bd_action_2d(chain) == 2 * n
    assert bd_action_2d(antichain) == 2 * n


def _sprinkled(n: int, seed: int):
    events = sprinkle_minkowski_causal_diamond(n, spacetime_dim=2, T=2.0, seed=seed)
    C = np.array(causal_matrix_1p1(events), dtype=bool)
    np.fill_diagonal(C, False)
    return C


def test_sprinkled_flat_smeared_action_near_zero_and_damped():
    # Flat-space value is 0; smearing damps per-realization fluctuations.
    raw = [smeared_action_2d(_sprinkled(300, s), 1.0) for s in range(6)]
    smeared = [smeared_action_2d(_sprinkled(300, s), 0.08) for s in range(6)]
    assert np.std(smeared) < 0.2 * np.std(raw)
    assert all(abs(v) < 8.0 for v in smeared)
