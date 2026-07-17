"""Regression tests for the chain-selector robustness probe.

The probe's docs conclusion (docs/selector_robustness_probe.md) rests on
deterministic cells that must not silently drift: the anchored selector
produces valid disjoint chains, the uniform positive control passes under
both the frozen greedy selector and the non-greedy anchored one, the
frozen P6a layered hard negative still reaches the numerical gates and
blocks there under both, and the crystal never yields chains. Seeds are
fixed; every assertion is a deterministic replay of probe cells.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from selector_robustness_probe import (  # noqa: E402
    FROZEN_VARIANT,
    VARIANTS,
    _majority_verdict,
    evaluate,
    order_inputs,
    select_anchored_chains,
)

from causal_spacetime_lab.positive_control.two_orders import (  # noqa: E402
    balanced_layered_perm,
    bipartite_perm,
    windowed_transpositions,
)

GREEDY_FROZEN = VARIANTS[0]
ANCHORED = VARIANTS[5]


def test_variant_table_is_the_documented_grid():
    assert GREEDY_FROZEN[0] == FROZEN_VARIANT
    assert [v[1] for v in VARIANTS].count("anchored") == 3
    assert {(v[2], v[3]) for v in VARIANTS} == {
        (6, 25), (4, 25), (8, 25), (6, 20), (6, 30),
    }


def test_anchored_selector_returns_valid_disjoint_chains():
    causal, times, _ = order_inputs(np.random.default_rng(200).permutation(600))
    chains = select_anchored_chains(causal, times, 6, 25, seed=0)
    again = select_anchored_chains(causal, times, 6, 25, seed=0)
    assert len(chains) == 6
    assert all(np.array_equal(a, b) for a, b in zip(chains, again, strict=True))
    flat = np.concatenate(chains)
    assert len(set(flat.tolist())) == flat.size
    for chain in chains:
        assert chain.size >= 25
        assert all(
            causal[a, b] for a, b in zip(chain[:-1], chain[1:], strict=True)
        )


def test_uniform_positive_control_passes_under_both_selectors():
    causal, times, coords = order_inputs(
        np.random.default_rng(200).permutation(600)
    )
    for variant in (GREEDY_FROZEN, ANCHORED):
        result = evaluate(causal, times, coords, variant, 200)
        assert result["status"] == "ok", (variant[0], result)
        assert result["gate_pass"], (variant[0], result)


def test_layered_hard_negative_blocks_at_the_gates_under_both_selectors():
    """The frozen P6a (k = 25, moves = 100) cell, seed 100: chain-rich
    (reaches the numerical gates -- no structural block) yet blocked
    there, whichever selector found its chains."""

    base = balanced_layered_perm(600, layer_count=25, seed=100, min_layer_size=6)
    pi = windowed_transpositions(base, moves=100, window=60, seed=10_100)
    causal, times, coords = order_inputs(pi)
    for variant in (GREEDY_FROZEN, ANCHORED):
        result = evaluate(causal, times, coords, variant, 100)
        assert result["status"] == "ok", (variant[0], result)
        assert not result["gate_pass"], (variant[0], result)


def test_crystal_blocks_structurally_under_both_selectors():
    causal, times, coords = order_inputs(bipartite_perm(600))
    for variant in (GREEDY_FROZEN, ANCHORED):
        result = evaluate(causal, times, coords, variant, 0)
        assert result["status"].startswith("structural_block"), result


def test_majority_verdict_criterion():
    assert _majority_verdict({"evaluated": 0, "gate_passes": 0}) is None
    assert _majority_verdict({"evaluated": 4, "gate_passes": 2}) == "pass"
    assert _majority_verdict({"evaluated": 4, "gate_passes": 1}) == "block"
