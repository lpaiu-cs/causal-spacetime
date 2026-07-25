"""Regressions for A1: the finite-resolution cost of the rigidity results.

Descriptive measurements, not preregistered and not gating anything. They
are pinned here because Section 8.7 now travels with operating conditions,
and a silent drift in those conditions would change what the manuscript is
entitled to say.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "theory"
sys.path.insert(0, str(EXPERIMENT_DIR))

from t1_g4b_unlabeled_2plus1d import (  # noqa: E402
    RANK_TOLERANCE,
    dissimilarity,
    flatten,
    jacobian,
    profile_to_dissimilarity,
    radial_profile,
    scene,
)
from t1_g4c_conditioning import (  # noqa: E402
    INSTRUMENT_TICKS,
    MIN_HEADROOM,
    amplification,
    check_error_model_sensitivity,
    check_frozen_instrument,
    check_minimum_observer_count_is_not_the_operating_point,
    diameter,
    instrument_tick_ladder,
    margin,
)

# ------------------------------------------------------------------
# the observable has one definition
# ------------------------------------------------------------------

def test_the_split_observable_reproduces_the_whole_one():
    """`profile_to_dissimilarity(radial_profile(...))` must BE
    `dissimilarity(...)`. The split exists so a perturbed-profile study
    can enter halfway through; if the halves ever stop composing, this
    module is measuring a different instrument than G4b proved about."""

    for d, R, n in ((2, 8, 12), (3, 5, 9), (4, 6, 8)):
        X, P = scene(n, R, seed=11 * n + R, d=d)
        theta = flatten(X, P)
        np.testing.assert_allclose(
            profile_to_dissimilarity(radial_profile(theta, n, R, d), R),
            dissimilarity(theta, n, R, d),
            rtol=0, atol=0,
        )


# ------------------------------------------------------------------
# the margin, and the trap next to it
# ------------------------------------------------------------------

def test_margin_reports_nullity_so_a_collapsed_rank_cannot_hide():
    """`sigma_min` alone is "smallest above the cutoff", not "smallest
    nonzero". Every margin therefore carries its nullity and its
    headroom, so a configuration that has quietly lost a dimension of
    rank cannot be read as a well-conditioned one."""

    X, P = scene(48, 8, seed=6100 + 8, d=2)
    reading = margin(flatten(X, P), 48, 8, 2)
    assert reading["nullity"] == reading["gauge"] == 3
    assert reading["rigid"]
    assert reading["margin_is_readable"]
    assert reading["headroom_over_rank_tolerance"] > MIN_HEADROOM


def test_the_rank_cutoff_trap_is_real_and_would_flatter_the_worse_scene():
    """At d = 2, R = 32 the margin sits about four times the cutoff.
    Nudging the cutoff one notch reclassifies a genuine flex direction as
    null and the reported margin jumps by more than an order of
    magnitude -- upward. This pins the trap so the guard is not removed
    as redundant."""

    X, P = scene(48, 32, seed=6100 + 32, d=2)
    spectrum = np.linalg.svd(jacobian(flatten(X, P), 48, 32, 2),
                             compute_uv=False)
    at_shipped = spectrum[spectrum > spectrum[0] * RANK_TOLERANCE].min()
    at_looser = spectrum[spectrum > spectrum[0] * (5 * RANK_TOLERANCE)].min()
    assert at_shipped / spectrum[0] < 10 * RANK_TOLERANCE   # close to the edge
    assert at_looser / at_shipped > 10                      # and it jumps up


# ------------------------------------------------------------------
# the headline
# ------------------------------------------------------------------

def test_the_frozen_instrument_survives_its_own_resolution():
    """At the configuration the pipeline actually uses, a delta/2 readout
    bound costs about two of itself in position, and the response is
    linear in delta -- so the 2+1D metric result is usable at the
    instrument's resolution rather than only in the exact model."""

    outcome = check_frozen_instrument()
    assert outcome["passed"]
    assert outcome["margin"]["rigid"]
    assert 1.5 < outcome["amp_at_the_instruments_own_resolution"] < 3.0
    assert outcome["amp_sd_there"] > 0.0          # dispersion is reported
    assert outcome["position_error_over_L_there"] < 0.05
    assert outcome["relative_spread_of_amp_over_the_fine_half"] < 0.1


def test_the_headline_is_taken_from_the_instruments_own_tick_count():
    """Looked up by tick count, never by position in the ladder. The
    manuscript quotes this number."""

    outcome = check_frozen_instrument()
    assert outcome["instrument_ticks"] == INSTRUMENT_TICKS == 96
    row = next(r for r in outcome["rows"] if r["ticks"] == INSTRUMENT_TICKS)
    assert row["amp"] == outcome["amp_at_the_instruments_own_resolution"]
    assert row["position_error_over_L"] == \
        outcome["position_error_over_L_there"]


def test_the_tick_ladder_agrees_with_g4a_rather_than_restating_it():
    """G4a's tracked table is the authority on the instrument's delta at
    a given tick count. The two extra rungs are labelled as extensions."""

    ladder = instrument_tick_ladder()
    assert ladder["shared_rungs_agree_with_g4a"]
    assert not ladder["mismatched_ticks"]
    sources = {r["ticks"]: r["source"] for r in ladder["rungs"]}
    assert sources[96] == "g4a table"
    assert sources[1536] == "extension of the same rule"


def test_position_error_falls_linearly_with_tick_count():
    """Halving delta halves the error. No floor, no threshold."""

    rows = {r["ticks"]: r for r in check_frozen_instrument()["rows"]}
    for coarse, fine in ((96, 192), (192, 384), (384, 768)):
        ratio = (rows[coarse]["position_error_over_L"]
                 / rows[fine]["position_error_over_L"])
        assert 1.8 < ratio < 2.2, (coarse, fine, ratio)


# ------------------------------------------------------------------
# what the answer rests on
# ------------------------------------------------------------------

def test_the_independence_assumption_is_load_bearing():
    """Only the SUPPORT of the readout error is proved. Centering across
    observers annihilates any common-mode component outright, so the
    correlation structure -- an assumption -- moves the answer from ~2 to
    ~0. Independence is the conservative choice among these, which is why
    it is the one shipped, but it is a choice."""

    outcome = check_error_model_sensitivity()
    amps = outcome["amp_by_model"]
    assert amps["common_mode"] < 1e-6          # centering kills it
    assert 1.5 < amps["independent"] < 3.0
    assert amps["common_mode"] < amps["half_common"] < amps["independent"]
    assert outcome["independent_is_the_conservative_choice"]


def test_sigma_min_is_the_operative_quantity_not_a_formality():
    """amp * sigma_min stays inside one order of magnitude across four
    decades of sigma_min. The exact-model margin IS the error budget."""

    products = []
    for d, R, n in ((2, 8, 34), (2, 16, 48), (3, 5, 48), (3, 8, 48)):
        X, P = scene(n, R, seed=4321 + n + R, d=d)
        theta = flatten(X, P)
        reading = margin(theta, n, R, d)
        assert reading["margin_is_readable"], (d, R, n)
        products.append(
            amplification(theta, n, R, d, 1e-4)["amp"] * reading["sigma_min"]
        )
    assert max(products) / min(products) < 10


def test_crowding_observers_destroys_the_margin():
    """More observers lower the rigidity threshold and wreck the error
    budget. The exact model sees only the first half of that."""

    sigmas = {}
    for R in (8, 16, 24):
        X, P = scene(48, R, seed=6100 + R, d=2)
        reading = margin(flatten(X, P), 48, R, 2)
        assert reading["rigid"], R
        sigmas[R] = reading["sigma_min"]
    assert sigmas[8] > sigmas[16] > sigmas[24]
    assert sigmas[8] / sigmas[24] > 100


def test_the_proved_threshold_is_the_worst_place_to_operate():
    """R = d + 2 is exactly the corner G4c pins, and exactly the corner
    with the smallest margin."""

    outcome = check_minimum_observer_count_is_not_the_operating_point()
    assert outcome["passed"]
    assert outcome["threshold_is_always_worse"]
    by_d = {r["d"]: r["amp_penalty_for_sitting_at_the_threshold"]
            for r in outcome["rows"]}
    assert by_d[3] > 2.0          # 3+1D pays the most for minimality
    assert all(v > 1.0 for v in by_d.values())


def test_amplification_is_scale_free():
    """J is homogeneous of degree zero, so amp must not move when the
    whole scene is rescaled. This is the invariant that made the first
    version of this analysis wrong -- delta was compared to sigma_min
    without dividing by the scene size."""

    X, P = scene(34, 8, seed=4363, d=2)
    theta = flatten(X, P)
    big = theta * 7.5
    assert abs(diameter(big, 2) / diameter(theta, 2) - 7.5) < 1e-9
    small_reading = amplification(theta, 34, 8, 2, 1e-3)["amp"]
    big_reading = amplification(big, 34, 8, 2, 1e-3)["amp"]
    assert abs(small_reading - big_reading) / small_reading < 1e-6
