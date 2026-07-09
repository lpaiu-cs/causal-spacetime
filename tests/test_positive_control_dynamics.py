"""Tests for the P3 geometry-free dynamics and order-intrinsic selection."""

from __future__ import annotations

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.positive_control.dynamics import (
    relation_density,
    transitive_percolation,
)
from causal_spacetime_lab.positive_control.order_intrinsic import (
    global_longest_chain,
    measure_order_intrinsic_profiles,
    select_bracketed_targets,
    select_disjoint_chains,
)
from causal_spacetime_lab.positive_control.rewire import transitive_closure
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond


def test_transitive_percolation_is_transitive_acyclic_and_deterministic() -> None:
    causal, idx = transitive_percolation(200, 0.02, seed=1)
    assert np.array_equal(causal, transitive_closure(causal))  # transitive
    assert not np.any(np.diag(causal))  # irreflexive
    assert np.array_equal(idx, np.arange(200))
    again, _ = transitive_percolation(200, 0.02, seed=1)
    assert np.array_equal(causal, again)  # deterministic
    # growth order respects the relation (time-respecting): no i->j with i>j
    assert not np.any(np.tril(causal, k=-1))


def test_density_increases_with_p() -> None:
    idx = np.arange(400, dtype=float)
    lo, _ = transitive_percolation(400, 0.004, seed=0)
    hi, _ = transitive_percolation(400, 0.02, seed=0)
    assert relation_density(lo, idx) < relation_density(hi, idx)


def test_global_longest_chain_on_a_total_order() -> None:
    # a total order 0<1<...<n-1: the longest chain is everything
    n = 30
    causal = transitive_closure(np.triu(np.ones((n, n), bool), k=1))
    chain = global_longest_chain(causal, np.arange(n),
                                 np.ones(n, dtype=bool))
    assert chain.size == n
    assert list(chain) == list(range(n))


def test_disjoint_chains_are_disjoint_and_ordered() -> None:
    ev = sprinkle_1p1_causal_diamond(800, T=2.0, seed=2)
    causal = causal_matrix_1p1(ev)
    chains = select_disjoint_chains(causal, ev[:, 0], chain_count=4, min_length=10)
    assert len(chains) == 4
    seen = np.concatenate(chains)
    assert seen.size == len(set(seen.tolist()))  # disjoint
    for chain in chains:  # each is a genuine causal chain
        for a in range(chain.size - 1):
            assert causal[chain[a], chain[a + 1]]


def test_order_intrinsic_profiles_bracketed_and_measurable() -> None:
    ev = sprinkle_1p1_causal_diamond(1000, T=2.0, seed=3)
    causal = causal_matrix_1p1(ev)
    chains = select_disjoint_chains(causal, ev[:, 0], 5, 15)
    targets = select_bracketed_targets(causal, chains, max_targets=30, seed=3)
    assert targets.size >= 10
    profiles = measure_order_intrinsic_profiles(causal, chains, targets)
    assert profiles.delay_ranks.shape == (targets.size, len(chains))
    # every selected target is two-sided bracketed by every chain
    assert bool(np.all(profiles.reachable))
