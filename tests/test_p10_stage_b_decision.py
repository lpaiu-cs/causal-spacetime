"""Regressions for P10 Stage B's frozen decision block and its gates.

Stage B is closed (the B0 yardstick gate failed), so no real B1 data
exist and none may be generated. Precisely because of that, the frozen
8.3 decision logic can only be verified synthetically — the P9 pattern.
A review found the first aggregator claimed to compute these verdicts
and emitted a stub instead; these tests are what prevents that from
recurring silently.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p10_stage_b import (  # noqa: E402
    B0_SAMPLES,
    B1_CHAIN_SEEDS,
    LADDER,
    SAMPLES_PER_CHAIN,
    STARTS,
    TOST_MARGIN,
    _require_b0_gate,
    b0_gate_verdict,
    chain_is_complete,
    evaluate_frozen_hypotheses,
    require_all_chains,
    require_frozen_chain,
    scaled_constants,
)


def _cloud(rng, center, n=48, spread=0.01):
    return list(center + spread * rng.standard_normal(n))


def test_tracking_and_deepening_supports_the_conjunction():
    rng = np.random.default_rng(1)
    e = {600: _cloud(rng, 0.16), 900: _cloud(rng, 0.14), 1200: _cloud(rng, 0.12)}
    s = {600: _cloud(rng, 0.16), 900: _cloud(rng, 0.14), 1200: _cloud(rng, 0.12)}
    out = evaluate_frozen_hypotheses(e, s)
    assert out["h_track_supported"]
    assert out["h_deepen_supported"]
    assert out["consistent_with_continuum_limit"]


def test_a_real_offset_at_one_rung_fails_h_track():
    """An E-S offset outside the margin at a single rung is enough."""

    rng = np.random.default_rng(2)
    e = {600: _cloud(rng, 0.16), 900: _cloud(rng, 0.14),
         1200: _cloud(rng, 0.12 + TOST_MARGIN + 0.03)}
    s = {600: _cloud(rng, 0.16), 900: _cloud(rng, 0.14), 1200: _cloud(rng, 0.12)}
    out = evaluate_frozen_hypotheses(e, s)
    assert not out["per_rung"]["1200"]["tost_pass"]
    assert not out["h_track_supported"]
    assert not out["consistent_with_continuum_limit"]


def test_wide_uncertainty_fails_tost_even_with_zero_offset():
    """TOST is an equivalence test: a CI wider than the margin fails
    even when the point estimate is exactly zero. This is the distinction
    the P10-A record was corrected to respect, now enforced in code."""

    rng = np.random.default_rng(3)
    e = {600: _cloud(rng, 0.16, n=6, spread=0.08)}
    s = {600: _cloud(rng, 0.16, n=6, spread=0.08)}
    out = evaluate_frozen_hypotheses(e, s)
    assert not out["per_rung"]["600"]["tost_pass"]


def test_flat_e_curve_fails_h_deepen():
    rng = np.random.default_rng(4)
    e = {600: _cloud(rng, 0.15), 900: _cloud(rng, 0.15), 1200: _cloud(rng, 0.15)}
    s = {600: _cloud(rng, 0.15), 900: _cloud(rng, 0.15), 1200: _cloud(rng, 0.15)}
    out = evaluate_frozen_hypotheses(e, s)
    assert out["h_track_supported"]          # tracking alone can hold
    assert not out["h_deepen_supported"]     # but nothing deepened
    assert not out["consistent_with_continuum_limit"]


def test_a_rung_with_no_surviving_chain_blocks_h_track():
    """A missing rung is not skipped -- H-TRACK requires every rung."""

    rng = np.random.default_rng(5)
    e = {600: _cloud(rng, 0.15), 900: [], 1200: _cloud(rng, 0.13)}
    s = {600: _cloud(rng, 0.15), 900: _cloud(rng, 0.14), 1200: _cloud(rng, 0.13)}
    out = evaluate_frozen_hypotheses(e, s)
    assert not out["per_rung"]["900"]["evaluable"]
    assert not out["h_track_supported"]


def test_the_b0_gate_locks_b1_under_the_committed_record():
    """The frozen B0 record says the yardstick did not fall, so the
    guard must refuse -- this refusal is the preregistration operating."""

    with pytest.raises(SystemExit, match="gate FAILED"):
        _require_b0_gate()


def test_a_shortened_chain_is_incomplete():
    """The frozen sampler can silently drop a scheduled retention point
    (an i == j draw lands on it, ~1/N per point), and the sampler may
    not be touched. A 47-row chain must therefore fail completeness --
    its m = 48 diagnostic basis and its nominal step labels are gone."""

    full = [{"chain_complete": "True"} for _ in range(SAMPLES_PER_CHAIN)]
    assert chain_is_complete(full)
    assert not chain_is_complete(full[:-1])                # short
    flagged = full[:-1] + [{"chain_complete": "False"}]
    assert not chain_is_complete(flagged)                  # runner flagged it


def test_random_start_shards_alone_cannot_reach_the_hypotheses():
    """Aggregating three random-start chains would auto-satisfy the melt
    criterion and evaluate the frozen hypotheses from half the design.
    Presence of all six (rung, start) chains is required first; whether
    a present chain survives is then the screen's job."""

    all_six = {(n, s) for n in LADDER for s in STARTS}
    require_all_chains(all_six)                       # complete: no exit
    with pytest.raises(SystemExit, match="missing"):
        require_all_chains({(n, "random") for n in LADDER})
    with pytest.raises(SystemExit, match="1200"):
        require_all_chains(all_six - {(1200, "bipartite")})


def test_b1_seeds_do_not_collide_with_p6b_reference_orders():
    """The first table used 1000-1050; P6-B's uniform reference orders
    are exactly default_rng(1000..1019).permutation(600), and every B1
    random start is default_rng(seed + 1).permutation(N) -- so the
    (600, random) chain would have literally reused a P6-B reference
    order. Pin the whole derived-stream envelope clear of that range."""

    p6b_reference = set(range(1000, 1020))
    for (n, start), seed in B1_CHAIN_SEEDS.items():
        derived = {seed, seed + 1} | {seed + 7 * k for k in range(48)}
        assert not (derived & p6b_reference), (n, start)


def test_stale_or_wrong_seed_chains_are_rejected():
    """Grouping by (n, start) alone would let a correctly named CSV
    produced with a different seed replace the preregistered chain, or
    a seventh chain slip in. Both must exit."""

    good_seed = float(B1_CHAIN_SEEDS[(600, "random")])
    good = [{"chain_seed": good_seed} for _ in range(3)]
    require_frozen_chain(600, "random", good)          # no exit
    with pytest.raises(SystemExit, match="not the preregistered chain"):
        require_frozen_chain(600, "random",
                             [{"chain_seed": 1000.0} for _ in range(3)])
    with pytest.raises(SystemExit, match="not one of the six"):
        require_frozen_chain(600, "uniform", good)


def test_the_b0_gate_requires_full_rungs_not_just_a_falling_ci():
    """A run whose instrument structurally blocked most samples could
    otherwise authorize B1 from whatever it happened to measure -- one
    surviving sample per endpoint yields a degenerate 'passing' CI.
    The gate is completeness AND fall, never fall alone."""

    full = {600: B0_SAMPLES, 900: B0_SAMPLES, 1200: B0_SAMPLES}
    falling = (-0.03, (-0.05, -0.01))
    assert b0_gate_verdict(full, *falling)
    assert not b0_gate_verdict(full, -0.03, (-0.05, +0.01))   # CI touches 0
    short = dict(full, **{1200: B0_SAMPLES - 1})
    assert not b0_gate_verdict(short, *falling)               # one block
    assert not b0_gate_verdict({600: 1, 900: 1, 1200: 1}, *falling)


def test_the_anchor_operating_point_is_the_frozen_instrument():
    assert scaled_constants(600) == {
        "chain_count": 6, "min_chain_len": 25, "max_targets": 44,
        "min_targets": 20, "train_c": 3000, "heldout_c": 800,
    }
