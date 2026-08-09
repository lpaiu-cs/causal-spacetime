"""Regressions for S1: the Schwarzschild causality-cost measurement.

Pins the flight-time solver against every regime with an independent
answer -- the flat limit, the exact radial tortoise difference, the
weak-field Shapiro delay IN SCHWARZSCHILD COORDINATES (the isotropic
form differs at the same order and was the first draft's wrong
check), symmetry, and monotonicity -- plus the undecided band, the
patch-safety diagnostics, and the benchmark's deterministic counts.
Wall-clock SECONDS are host-dependent and are never pinned.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

import s1_schwarzschild_cost as s1  # noqa: E402


def test_the_flat_limit_is_the_euclidean_chord():
    """M = 0 must reproduce straight-line light: T equals the chord
    for both families (no-turn r1 != r2, one-turn r1 = r2)."""

    for r1, r2, dpsi in ((12.0, 17.0, 1.3), (15.0, 15.0, 1.0),
                         (10.0, 20.0, 0.4)):
        chord = math.sqrt(r1 * r1 + r2 * r2
                          - 2.0 * r1 * r2 * math.cos(dpsi))
        t, err = s1.flight_time(r1, r2, dpsi, m=0.0)
        assert abs(t - chord) <= max(err, 1e-6), (r1, r2, dpsi)


def test_the_radial_flight_is_the_exact_tortoise_difference():
    t, err = s1.flight_time(10.0, 20.0, 0.0)
    want = s1.tortoise(20.0, 1.0) - s1.tortoise(10.0, 1.0)
    assert err == 0.0
    assert t == pytest.approx(want, rel=1e-15)
    # and the closed form itself: r + 2M ln(r/2M - 1)
    assert s1.tortoise(10.0, 1.0) == pytest.approx(
        10.0 + 2.0 * math.log(4.0), rel=1e-15)


def test_the_weak_field_matches_shapiro_in_schwarzschild_coordinates():
    """First order in M, along the straight chord with impact
    distance p and foot segments s_i:

        T = d + 2M ln[(r1+r2+d)/(r1+r2-d)] - M (s1/r1 + s2/r2)

    The trailing term is the difference between Schwarzschild and
    isotropic coordinates at THIS order; the bare 2M-log (isotropic)
    form disagrees by ~0.2 M here and was the wrong first check.
    Residual must be O(M^2/p)."""

    r1 = r2 = 1e4
    dpsi = 0.2
    d = math.sqrt(r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * math.cos(dpsi))
    p = r1 * r2 * math.sin(dpsi) / d
    s_1 = math.sqrt(r1 * r1 - p * p)
    s_2 = math.sqrt(r2 * r2 - p * p)
    want = (d + 2.0 * math.log((r1 + r2 + d) / (r1 + r2 - d))
            - (s_1 / r1 + s_2 / r2))
    iso = d + 2.0 * math.log((r1 + r2 + d) / (r1 + r2 - d))
    t, _ = s1.flight_time(r1, r2, dpsi)
    assert abs(t - want) < 1e-3          # O(M^2/p) head-room
    assert abs(t - iso) > 0.1            # the wrong form IS wrong


def test_symmetry_and_monotonicity_on_the_patch():
    """T(x1, x2) = T(x2, x1) exactly, and T grows with the angular
    separation at fixed radii."""

    assert s1.flight_time(11.0, 19.0, 1.7)[0] == pytest.approx(
        s1.flight_time(19.0, 11.0, 1.7)[0], rel=1e-12)
    times = [s1.flight_time(12.0, 17.0, x)[0]
             for x in (0.2, 0.6, 1.0, 1.4, 1.8)]
    assert all(a < b for a, b in zip(times, times[1:]))
    # gravity delays: GR flight beats the flat chord everywhere tried
    for r1, r2, dpsi in ((12.0, 17.0, 1.3), (10.0, 10.0, 2.0)):
        chord = math.sqrt(r1 * r1 + r2 * r2
                          - 2.0 * r1 * r2 * math.cos(dpsi))
        assert s1.flight_time(r1, r2, dpsi)[0] > chord


def test_the_patch_stays_far_from_the_photon_sphere():
    """The frozen shell/cap keeps every direct perihelion above 5 M
    (photon sphere at 3 M): the worst case is the inner-radius pair
    at the maximal cap separation."""

    det: dict = {}
    t, _ = s1.flight_time(10.0, 10.0, s1.PSI_MAX, details=det)
    assert det["family"] == "one-turn"
    assert det["r_perihelion"] > 5.0
    det2: dict = {}
    s1.flight_time(12.0, 17.0, 0.3, details=det2)
    assert det2["family"] in ("no-turn", "one-turn")


def test_the_undecided_band_is_never_a_silent_verdict():
    """dt inside the error band returns None; clearly-before returns
    False; clearly-related returns True (P1's interval discipline)."""

    p = np.array([0.0, 12.0, 0.5, 0.3])
    q = np.array([0.0, 17.0, 0.7, 1.1])
    cosang = (math.sin(p[2]) * math.sin(q[2]) * math.cos(p[3] - q[3])
              + math.cos(p[2]) * math.cos(q[2]))
    dpsi = math.acos(cosang)
    t_min, err = s1.flight_time(p[1], q[1], dpsi)
    assert err > 0.0
    q_rel = q.copy()
    q_rel[0] = t_min * 2.0
    assert s1.causal_relation(p, q_rel) is True
    q_not = q.copy()
    q_not[0] = t_min * 0.5
    assert s1.causal_relation(p, q_not) is False
    q_amb = q.copy()
    q_amb[0] = t_min
    assert s1.causal_relation(p, q_amb) is None
    assert s1.causal_relation(q_rel, p) is False  # dt < 0


def test_the_sample_and_bench_are_deterministic():
    """Events respect the frozen patch (shell, cap, t-extent; any
    two cap directions are within PSI_MAX), and the benchmark's
    COUNTS reproduce from the seed -- seconds are host-dependent and
    unpinned."""

    rng = np.random.default_rng(s1.BENCH_SEED)
    ev = s1.sample_events(500, rng)
    assert (ev[:, 0] >= 0).all() and (ev[:, 0] <= s1.T_EXTENT).all()
    assert (ev[:, 1] >= s1.R_MIN).all() and (ev[:, 1] <= s1.R_MAX).all()
    assert (ev[:, 2] <= s1.CAP_HALF_ANGLE + 1e-12).all()
    b1 = s1.bench(n=60, seed=777_001, tol=1e-6, k_pairs=40)
    b2 = s1.bench(n=60, seed=777_001, tol=1e-6, k_pairs=40)
    assert b1["related"] == b2["related"]
    assert b1["undecided"] == b2["undecided"]
    assert b1["pairs_sampled"] == 40
