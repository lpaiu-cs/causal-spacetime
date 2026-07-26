"""Regressions for the Stage B' scale-fixed observable.

The one that carries everything is the anchor: at ``delta = 0`` the
restricted scorer must reproduce the frozen scorer BIT FOR BIT on the
same seed and pair stream — that is what makes it a restriction of the
frozen observable rather than a second definition (the lesson this
programme has now paid for twice).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p10_stage_bprime import (  # noqa: E402
    B0PRIME_SEED_BASE,
    DELTA_MARGIN,
    FINE_BIN_EDGES,
    HIGH_EDGE,
    INSTRUMENT_DERIVED_OFFSETS,
    SEEDFIX_SCORE_SEED_BASE,
    SPACED_SCORE_OFFSET,
    SPACED_STRIDE,
    _seed_collision_map,
    _spaced_seed,
    _spaced_window_map,
    bprime_gate_verdict,
    margin_restricted_order_error,
)

from causal_spacetime_lab.ordinal_embedding import (  # noqa: E402
    embedding_distance_order_error,
)


def _scene(seed=7, n=40):
    rng = np.random.default_rng(seed)
    true_x = np.sort(rng.uniform(-1, 1, size=n))
    fitted = true_x + 0.05 * rng.standard_normal(n)
    return fitted, true_x


def test_delta_zero_reproduces_the_frozen_scorer_exactly():
    """Same seed, same count, delta = 0: the pair streams coincide and
    the two scorers must agree to the last bit."""

    fitted, true_x = _scene()
    frozen = embedding_distance_order_error(
        fitted.reshape(-1, 1), true_x.reshape(-1, 1),
        num_pair_comparisons=8000, seed=123,
    )
    restricted, n_eligible = margin_restricted_order_error(
        fitted, true_x, delta=0.0, num_pair_comparisons=8000, seed=123,
    )
    assert restricted == frozen
    assert n_eligible > 7000          # only the i==j mask rejects


def test_raising_the_margin_makes_the_question_easier():
    """Restricting to wider-margin comparisons must not increase the
    error: near-ties are where a noisy fit flips signs."""

    fitted, true_x = _scene(seed=11)
    errors = []
    for delta in (0.0, 0.05, DELTA_MARGIN, 0.4):
        err, n_el = margin_restricted_order_error(
            fitted, true_x, delta, num_pair_comparisons=32_000, seed=5,
        )
        errors.append(err)
        assert n_el > 500
    assert errors == sorted(errors, reverse=True)


def test_a_perfect_fit_scores_zero_at_any_margin():
    _, true_x = _scene(seed=13)
    for delta in (0.0, DELTA_MARGIN):
        err, _ = margin_restricted_order_error(
            true_x.copy(), true_x, delta,
            num_pair_comparisons=8000, seed=3,
        )
        assert err == 0.0


def test_eligible_fraction_matches_the_continuum_geometry():
    """Uniform points in the order's own coordinates: the a-priori
    threshold was the continuum 25th percentile, so about 75% of
    comparisons must survive it. Loose bounds; this pins the
    normalization convention (x in units of the region), where an error
    of a factor N would send the fraction to ~0 or ~1."""

    rng = np.random.default_rng(17)
    true_x = rng.uniform(0, 1, size=300) - rng.uniform(0, 1, size=300)
    _, n_eligible = margin_restricted_order_error(
        true_x + 0.01 * rng.standard_normal(300), true_x, DELTA_MARGIN,
        num_pair_comparisons=32_000, seed=19,
    )
    fraction = n_eligible / 32_000
    assert 0.60 < fraction < 0.85


def test_the_bprime_gate_requires_completeness_floors_and_fall():
    full = {600: 20, 900: 20, 1200: 20}
    falling = (-0.03, (-0.05, -0.01))
    assert bprime_gate_verdict(full, True, *falling)
    assert not bprime_gate_verdict(full, False, *falling)      # floor unmet
    assert not bprime_gate_verdict(full, True, -0.03, (-0.05, 0.01))
    assert not bprime_gate_verdict({**full, 1200: 19}, True, *falling)


def test_bprime_seeds_are_fresh():
    """Clear of every range used anywhere in the programme, including
    P6-B's 1000-1019 and the B1 envelope 30000-30379.

    This test checked derived seeds against EARLIER ranges only, and
    that was its hole: the base+k+9 scorer seeds also had to be checked
    against the B'0 chain seeds themselves, and were not — the round-six
    review finding. The two tests below close it."""

    used = (set(range(0, 10)) | set(range(100, 120)) | set(range(400, 420))
            | set(range(500, 520)) | set(range(820, 960))
            | set(range(1000, 1020)) | set(range(1100, 1160))
            | set(range(30000, 30380)))
    for base in B0PRIME_SEED_BASE.values():
        derived = {base + k for k in range(20)} | {base + k + 9 for k in range(20)}
        assert not (derived & used), base


def test_the_old_scorer_collision_is_exactly_the_diagnosed_one():
    """Pin the round-six diagnosis: 51 of 60 derived scorer seeds
    coincide with some sample's chain seed (values 40009-40059). If a
    seed table change ever alters this map, the prereg record's
    correction note goes stale and this fails."""

    diag = _seed_collision_map()
    assert diag["chain_seed_span"] == [40000, 40059]
    assert diag["old_scorer_seed_span"] == [40009, 40068]
    assert diag["n_old_scorer_seeds_colliding_with_chains"] == 51
    assert diag["n_scorer_seeds_total"] == 60


def test_the_seedfix_namespace_is_disjoint_from_everything():
    """The corrected scorer namespace must be fresh with respect to the
    chain seeds, the old scorer seeds, AND every documented earlier
    range — the full freshness rule, this time including ourselves."""

    used = (set(range(0, 10)) | set(range(100, 120)) | set(range(400, 420))
            | set(range(500, 520)) | set(range(820, 960))
            | set(range(1000, 1020)) | set(range(1100, 1160))
            | set(range(30000, 30380)))
    chain = {b + k for b in B0PRIME_SEED_BASE.values() for k in range(20)}
    old_scorer = {b + k + 9 for b in B0PRIME_SEED_BASE.values()
                  for k in range(20)}
    # the frozen instrument itself derives seed+9 (truth scorer) and
    # seed+100 (fit learning) from every chain seed -- the clean
    # namespace must clear those derived windows too
    instrument_derived = {b + k + 100 for b in B0PRIME_SEED_BASE.values()
                          for k in range(20)}
    clean = {b + k for b in SEEDFIX_SCORE_SEED_BASE.values()
             for k in range(20)}
    assert len(clean) == 60
    assert not (clean & (chain | old_scorer | instrument_derived | used))
    assert _seed_collision_map()["clean_disjoint_from_all"]


def test_the_high_edge_is_the_actual_bin_edge_not_a_rounded_literal():
    """A 0.478 literal once labelled a pool that begins at 0.4784; the
    threshold must be the bin edge itself."""

    assert HIGH_EDGE == FINE_BIN_EDGES[11]
    assert abs(HIGH_EDGE - 0.4784) < 1e-12


def test_the_spaced_windows_are_private_and_fresh():
    """Round eight: every stream the instrument derives (offsets 0, 3,
    5, 9, 61, 100) and the spaced scorer (+150) must live inside its
    own row's stride-200 window; the 60 windows must be pairwise
    disjoint and clear of every namespace ever used, including the
    B'0 spans and the seedfix namespace."""

    diag = _spaced_window_map()
    assert diag["pairwise_disjoint"]
    assert diag["offsets_inside_stride"]
    assert max(INSTRUMENT_DERIVED_OFFSETS) < SPACED_SCORE_OFFSET
    assert SPACED_SCORE_OFFSET < SPACED_STRIDE
    assert SPACED_SCORE_OFFSET not in INSTRUMENT_DERIVED_OFFSETS

    used = (set(range(0, 10)) | set(range(100, 120)) | set(range(400, 420))
            | set(range(500, 520)) | set(range(820, 960))
            | set(range(1000, 1020)) | set(range(1100, 1160))
            | set(range(9001, 9060)) | set(range(30000, 30380))
            | set(range(40000, 40169)) | set(range(41000, 41060)))
    spaced_all = {
        _spaced_seed(n, k) + off
        for n in B0PRIME_SEED_BASE
        for k in range(20)
        for off in (*INSTRUMENT_DERIVED_OFFSETS, SPACED_SCORE_OFFSET)
    }
    assert not (spaced_all & used)
