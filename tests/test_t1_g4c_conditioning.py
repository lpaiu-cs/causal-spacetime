"""Regressions for A1: the finite-resolution cost of the rigidity results.

Descriptive measurements, not preregistered and not gating anything. They
are pinned because the Section 8.6 metric claim now travels with operating
conditions, and a silent drift in those conditions would change what the
manuscript is entitled to say.
"""

from __future__ import annotations

import sys
from pathlib import Path

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "theory"
sys.path.insert(0, str(EXPERIMENT_DIR))

from t1_g4b_unlabeled_2plus1d import flatten, scene  # noqa: E402
from t1_g4c_conditioning import (  # noqa: E402
    amplification,
    check_frozen_instrument,
    check_minimum_observer_count_is_not_the_operating_point,
    diameter,
    smallest_nonzero_singular_value,
)


def test_the_frozen_instrument_survives_its_own_resolution():
    """The headline. At the configuration the pipeline actually uses, a
    delta/2 readout bound costs about two of itself in position, and the
    response is linear in delta -- so the 2+1D metric result is usable at
    the instrument's resolution rather than only in the exact model.
    """

    outcome = check_frozen_instrument()
    assert outcome["passed"]
    assert outcome["response_is_linear_in_delta"]
    assert 1.5 < outcome["amp_at_the_instruments_own_resolution"] < 3.0
    at_96 = next(r for r in outcome["rows"] if r["ticks"] == 96)
    assert at_96["position_error_over_L"] < 0.05


def test_position_error_falls_linearly_with_tick_count():
    """Halving delta halves the error. No floor, no threshold -- the
    exact-model regime is the one the instrument is actually in."""

    rows = {r["ticks"]: r for r in check_frozen_instrument()["rows"]}
    for coarse, fine in ((96, 192), (192, 384), (384, 768)):
        ratio = (rows[coarse]["position_error_over_L"]
                 / rows[fine]["position_error_over_L"])
        assert 1.8 < ratio < 2.2, (coarse, fine, ratio)


def test_sigma_min_is_the_operative_quantity_not_a_formality():
    """amp * sigma_min stays inside one order of magnitude across four
    decades of sigma_min. The exact-model margin IS the error budget."""

    products = []
    for d, R, n in ((2, 8, 34), (2, 16, 48), (3, 5, 48), (3, 8, 48)):
        X, P = scene(n, R, seed=4321 + n + R, d=d)
        theta = flatten(X, P)
        sigma = smallest_nonzero_singular_value(theta, n, R, d)
        products.append(amplification(theta, n, R, d, 1e-4) * sigma)
    assert max(products) / min(products) < 10


def test_crowding_observers_destroys_the_margin():
    """More observers lower the rigidity threshold and wreck the error
    budget. The exact model sees only the first half of that."""

    sigmas = {}
    for R in (8, 16, 24):
        X, P = scene(48, R, seed=6100 + R, d=2)
        sigmas[R] = smallest_nonzero_singular_value(flatten(X, P), 48, R, 2)
    assert sigmas[8] > sigmas[16] > sigmas[24]
    assert sigmas[8] / sigmas[24] > 100


def test_the_proved_threshold_is_the_worst_place_to_operate():
    """R = d + 2 is exactly the corner G4c pins, and exactly the corner
    with the smallest margin. The penalty is real in every dimension
    tested and worst in the ones that matter most."""

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
    a_small = amplification(theta, 34, 8, 2, 1e-3)
    a_big = amplification(big, 34, 8, 2, 1e-3)
    assert abs(a_small - a_big) / a_small < 1e-6
