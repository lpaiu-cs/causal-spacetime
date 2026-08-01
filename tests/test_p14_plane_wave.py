"""Regressions for the P14 probe geometry (design §8 P1).

Every test names the case it is meant to catch. Two are load-bearing:

`test_closed_form_cost_matches_independent_minimisation` -- the causal
predicate rests entirely on a closed-form transverse action, and a
closed form nobody checked against an independent computation is the
defect class P12 §10.1 refused. The oracle here is a direct discrete
minimisation that never touches the closed form.

`test_error_bound_survives_the_cancellation_that_broke_v1` -- the two
legs of the cost have opposite signs and cancel, so the FIRST version of
the bound (scaled to the cancelled total) understated the real error by
`3.2e11` and decided pairs with the wrong sign. That case is pinned.

Nothing here runs a probe. These are structural checks on the geometry.
"""

from __future__ import annotations

import math
import sys
from decimal import Decimal, getcontext
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p14_plane_wave import (  # noqa: E402
    PlaneWaveGeometry,
    Slab,
    _exact_cost,
    arms,
    causal_relation,
    conjugate_du,
    cost_defocusing,
    cost_error_bound,
    cost_focusing,
    slab_within_conjugate_bound,
    sprinkle,
    transverse_cost,
)

# a slab comfortably inside the conjugate bound, reused below
W = 1.0
SLAB = slab_within_conjugate_bound(W, 0.6, dv=2.0, dx=1.5, dy=1.5)
CURVED, FLAT = arms(SLAB, W)


# --------------------------------------------------------------------
# the independent oracle: minimise the discretised action directly
# --------------------------------------------------------------------

def _min_action_discrete(s, ap, aq, w2, n=4000):
    """min over paths of  int_0^s (a'^2 + w2 a^2) du / 2,  a(0)=ap, a(s)=aq.

    Discretise, leave the interior nodes free, solve the resulting
    quadratic minimisation by one linear solve. Uses no closed form, no
    `sinh`, no `sin` -- it knows only the integrand.
    """

    h = s / n
    m = n - 1
    main = np.full(m, 2.0 / h + w2 * h / 2.0)
    off = np.full(m - 1, -1.0 / h + w2 * h / 4.0)
    a_mat = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    rhs = np.zeros(m)
    rhs[0] = ap * (1.0 / h - w2 * h / 4.0)
    rhs[-1] = aq * (1.0 / h - w2 * h / 4.0)
    a = np.concatenate(([ap], np.linalg.solve(a_mat, rhs), [aq]))
    da = np.diff(a) / h
    mid = 0.5 * (a[:-1] + a[1:])
    return float(0.5 * np.sum(h * (da ** 2 + w2 * mid ** 2)))


@pytest.mark.parametrize(
    ("s", "xp", "xq", "yp", "yq", "w"),
    [
        (1.0, 0.0, 0.0, 0.0, 0.0, 0.7),
        (1.0, 0.3, -0.2, 0.1, 0.4, 0.7),
        (0.5, -0.4, 0.4, -0.3, 0.2, 1.5),
        (2.0, 0.1, 0.1, 0.2, -0.2, 0.5),
        (1.3, 0.6, 0.0, 0.0, 0.5, 1.0),
    ],
)
def test_closed_form_cost_matches_independent_minimisation(
        s, xp, xq, yp, yq, w):
    """The case: the whole causal predicate is one closed-form action.
    If it is wrong, every relation is wrong TOGETHER and no downstream
    statistic can reveal it. Checked against a minimisation sharing no
    code with it.
    """

    closed = transverse_cost(s, xp, xq, yp, yq, w)
    direct = (_min_action_discrete(s, xp, xq, w * w)
              + _min_action_discrete(s, yp, yq, -w * w))
    assert closed == pytest.approx(direct, rel=2e-6, abs=1e-12)


def test_flat_arm_is_the_straight_line_and_not_a_numerical_limit():
    """The case: `w = 0` must give `(Dx^2 + Dy^2)/(2 s)` exactly, not a
    `sinh(0)/0` that happens to survive. The flat arm is the reference
    every paired statistic is measured against.
    """

    s, xp, xq, yp, yq = 1.7, 0.3, -0.4, 0.25, 0.6
    expected = ((xq - xp) ** 2 + (yq - yp) ** 2) / (2.0 * s)
    assert transverse_cost(s, xp, xq, yp, yq, 0.0) == expected


def test_series_branch_joins_the_closed_form_without_a_step(monkeypatch):
    """The case: both directions are `0/0` as `w s -> 0`. A step at the
    crossover would make `D` depend on which side of an arbitrary
    constant a pair fell.

    Both branches are evaluated at the SAME `w`, by moving the
    threshold. An earlier version compared `w` just below the cut with
    `w` just above and demanded 1e-12; the gap it found was the real
    `a^2` term responding to a 0.2% change in `w`. A continuity test
    that varies the physics is not a continuity test.
    """

    import p14_plane_wave as p14

    s, xp, xq, yp, yq = 1.0, 0.4, -0.3, 0.2, 0.5
    for w in (1e-5, 1e-4, 1e-3):
        monkeypatch.setattr(p14, "_SMALL_WS", 0.0)
        closed = p14.transverse_cost(s, xp, xq, yp, yq, w)
        monkeypatch.setattr(p14, "_SMALL_WS", 1.0)
        series = p14.transverse_cost(s, xp, xq, yp, yq, w)
        assert math.isfinite(closed) and math.isfinite(series)
        assert closed == pytest.approx(series, rel=max(1e-13, (w * s) ** 4))


def test_curvature_costs_more_in_one_direction_and_less_in_the_other():
    """The case: the profile is a QUADRUPOLE (design §4.5) -- `x`
    defocuses, `y` focuses. If both moved the same way the construction
    would be a scalar deformation and `D`'s gained/lost split would have
    nothing to separate.
    """

    s, d, w = 1.0, 0.5, 1.2
    flat_leg = d * d / (2.0 * s)
    assert cost_defocusing(s, 0.0, d, w) > flat_leg
    assert cost_focusing(s, 0.0, d, w) < flat_leg


# --------------------------------------------------------------------
# the error bound: R6.1
# --------------------------------------------------------------------

def _cancelling_case():
    """`S_x + S_y ~ 0` with both legs O(1). Solved in exact arithmetic
    so the fixture does not inherit the imprecision it exists to probe.
    """

    s, w, y = 1.0, 1.0, 1.0
    lo, hi = 0.5, 1.5
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _exact_cost(s, 0.0, mid, y, y, w) > 0:
            hi = mid
        else:
            lo = mid
    return s, w, y, 0.5 * (lo + hi)


def test_error_bound_survives_the_cancellation_that_broke_v1():
    """The case, review R6.1, measured rather than argued. `S_x > 0` and
    `S_y < 0` cancel: at O(1) coordinates the legs are `+0.5463` and
    `-0.5463` while the total is `-1.1e-16`. v1 budgeted the error as
    `rel_tol * (|Dv| + |C|)` -- against the magnitude AFTER cancellation
    -- giving `1.1e-28` against a true error of `3.6e-17`, understating
    it by `3.2e11` and deciding such pairs with the wrong sign.

    The bound must dominate the actual error. Nothing weaker is a bound.
    """

    s, w, y, x = _cancelling_case()
    leg_x = cost_defocusing(s, 0.0, x, w)
    leg_y = cost_focusing(s, y, y, w)
    total = transverse_cost(s, 0.0, x, y, y, w)

    # the configuration really does cancel
    assert abs(leg_x) > 0.5 and abs(leg_y) > 0.5
    assert abs(total) < 1e-14

    getcontext().prec = 60
    actual = abs(Decimal(repr(total)) - _exact_cost(s, 0.0, x, y, y, w))
    bound = cost_error_bound(s, 0.0, x, y, y, w)
    assert Decimal(repr(bound)) >= actual, (
        f"bound {bound:.3e} does not cover the actual error {actual:.3e}")

    # and the v1 model would NOT have covered it -- the regression
    v1_model = 1e-12 * abs(total)
    assert Decimal(repr(v1_model)) < actual


def test_error_bound_survives_the_series_cancellation_that_broke_v2():
    """The case, review R7.1, and it is R6.1's second face. In the
    series branch the two legs' `a^2` terms have OPPOSITE signs and
    cancel, while their `a^4` terms have the SAME sign and ADD. v2
    scaled its truncation budget by `a^4 * (leg_x + leg_y)` -- a
    quantity that shares the `a^2` cancellation -- so it vanished
    exactly where the real error did not.

    Measured at `s = 1`, `w = 9.9e-5`, all transverse endpoints `1`:
    the series returns `0.0`, the true cost is `-8.005e-18`, and v2
    reported `3.58e-23` -- understating by `2.2e5` and returning a
    definite `False` for a pair whose exact margin is positive.
    """

    s, w, val = 1.0, 9.9e-5, 1.0
    approx = transverse_cost(s, val, val, val, val, w)
    getcontext().prec = 60
    exact = _exact_cost(s, val, val, val, val, w)
    actual = abs(Decimal(repr(approx)) - exact)
    bound = cost_error_bound(s, val, val, val, val, w)

    # the configuration really is in the regime the bound must cover
    assert approx == 0.0 and abs(float(exact)) > 1e-18
    assert Decimal(repr(bound)) >= actual, (
        f"bound {bound:.3e} does not cover the actual error "
        f"{float(actual):.3e}")

    # the first omitted term, derived, reproduces the whole error
    a, p_sum, q_cross = w * s, 2.0 * val * val, 2.0 * val * val
    per_leg = a ** 4 * (-p_sum / 24.0
                        + 7.0 * (p_sum - q_cross) / 360.0) / (2.0 * s)
    assert abs(2.0 * per_leg - float(exact)) < 1e-24

    # and the predicate now escalates instead of deciding wrongly
    geom = PlaneWaveGeometry(Slab(du=2.0, dv=2.0, dx=2.0, dy=2.0), w)
    rel = causal_relation(geom, np.array([0.0, 0.0, val, val]),
                          np.array([s, 4e-18, val, val]))
    assert rel.escalated and rel.related is True


def test_a_pair_double_cannot_decide_escalates_instead_of_guessing():
    """The case, design §5.1: the response to an undecidable pair is
    precision, never a tolerance. A pair placed inside the double-
    precision error of the cone must come back either decided BY
    ESCALATION or `ambiguous`, and never decided by the double path.
    """

    s, w, y, x = _cancelling_case()
    getcontext().prec = 60
    exact = _exact_cost(s, 0.0, x, y, y, w)
    p = np.array([0.0, 0.0, 0.0, y])
    geom = PlaneWaveGeometry(Slab(du=2.0, dv=2.0, dx=2.0, dy=2.0), w)

    # sit just inside the true cone, by less than the double error
    q = np.array([s, float(-exact - Decimal("1e-17")), x, y])
    rel = causal_relation(geom, p, q)
    assert rel.escalated, "a pair this close must not be decided in double"
    assert rel.related is True

    q_out = np.array([s, float(-exact + Decimal("1e-17")), x, y])
    rel_out = causal_relation(geom, p, q_out)
    assert rel_out.escalated
    assert rel_out.related is False


def test_a_pair_exactly_on_the_cone_stays_ambiguous_after_escalation():
    """The case: escalation resolves precision, not degeneracy. A pair
    constructed exactly on the cone has no sign, and the predicate must
    say so rather than round it into a verdict.
    """

    w = 0.0                      # flat arm: the cone is exactly representable
    geom = PlaneWaveGeometry(Slab(du=2.0, dv=4.0, dx=2.0, dy=2.0), w)
    p = np.array([0.0, 0.0, 0.0, 0.0])
    cost = transverse_cost(1.0, 0.0, 0.5, 0.0, 0.5, w)
    q = np.array([1.0, -cost, 0.5, 0.5])
    rel = causal_relation(geom, p, q)
    assert rel.ambiguous and rel.related is None


# --------------------------------------------------------------------
# the causal predicate
# --------------------------------------------------------------------

def _minkowski_related(p, q):
    """Standard flat causality in (t, z, x, y), as an outside check."""

    du, dv = q[0] - p[0], q[1] - p[1]
    dt, dz = (du - dv) / math.sqrt(2.0), (du + dv) / math.sqrt(2.0)
    dx, dy = q[2] - p[2], q[3] - p[3]
    return dt > 0 and dt * dt >= dz * dz + dx * dx + dy * dy


@pytest.mark.parametrize("seed", range(6))
def test_flat_arm_predicate_agrees_with_minkowski_coordinates(seed):
    """The case: the flat arm must be ordinary Minkowski causality,
    written here in null coordinates with a cost function -- one
    algebraic slip from a subtly different cone that would bias `D`
    everywhere while looking self-consistent.
    """

    geom = PlaneWaveGeometry(Slab(du=4.0, dv=4.0, dx=4.0, dy=4.0), 0.0)
    rng = np.random.default_rng(seed)
    checked = 0
    for _ in range(400):
        p, q = rng.uniform(-1.0, 1.0, 4), rng.uniform(-1.0, 1.0, 4)
        if q[0] < p[0]:
            p, q = q, p
        if q[0] == p[0]:
            continue
        rel = causal_relation(geom, p, q)
        assert not rel.ambiguous
        assert rel.related == _minkowski_related(p, q)
        checked += 1
    assert checked > 350


def test_the_null_generator_at_zero_u_separation_is_causal():
    """The case, review R6.3: at `Du = 0` the interval is `dx^2 + dy^2`,
    so a transversally coincident pair with `Dv < 0` lies on the null
    generator of `∂v` and IS future-related. v1 rejected every `Du = 0`
    pair, which is measure zero in a sprinkling but wrong in the
    predicate's contract -- and fixed anchors can sit exactly here.
    """

    geom = PlaneWaveGeometry(Slab(du=1.0, dv=2.0, dx=1.0, dy=1.0), 0.9)
    p = np.array([0.0, 0.0, 0.0, 0.0])
    assert causal_relation(geom, p, np.array([0.0, -1.0, 0.0, 0.0])).related
    # past-directed along the same generator is not
    assert not causal_relation(
        geom, p, np.array([0.0, +1.0, 0.0, 0.0])).related
    # and transversally separated at Du = 0 is spacelike
    assert not causal_relation(
        geom, p, np.array([0.0, -1.0, 0.3, 0.0])).related


def test_past_directed_pairs_are_unrelated():
    """The case: `Du < 0` must be rejected before any cost is formed --
    `transverse_cost` divides by `s`.
    """

    geom = PlaneWaveGeometry(Slab(du=2.0, dv=2.0, dx=2.0, dy=2.0), 0.9)
    rel = causal_relation(geom, np.array([1.0, 0.0, 0.0, 0.0]),
                          np.array([0.5, -5.0, 0.0, 0.0]))
    assert rel.related is False and not rel.ambiguous


# --------------------------------------------------------------------
# the geometry object: the conjugate bound, enforced by construction
# --------------------------------------------------------------------

def test_the_conjugate_bound_is_enforced_by_construction_not_by_assert():
    """The case, review R6.2: `Slab.validate` existed and `sprinkle` did
    not call it, and the per-pair guard was an `assert` that vanishes
    under `python -O`. A bound the caller may skip is not a property of
    the box, and if it is not a property of the box then `D`'s
    denominator is not `C(N,2)`.
    """

    w = 1.5
    assert conjugate_du(w) == pytest.approx(math.pi / w)
    PlaneWaveGeometry(Slab(du=conjugate_du(w) * 0.99, dv=1, dx=1, dy=1), w)
    with pytest.raises(ValueError, match="conjugate point"):
        PlaneWaveGeometry(Slab(du=conjugate_du(w), dv=1, dx=1, dy=1), w)
    with pytest.raises(ValueError, match="conjugate point"):
        PlaneWaveGeometry(Slab(du=conjugate_du(w) * 1.01, dv=1, dx=1, dy=1), w)
    # a raise, not an assert: it must survive -O
    with pytest.raises(ValueError):
        PlaneWaveGeometry(Slab(du=1.0, dv=1, dx=1, dy=1), -1.0)


def test_the_safety_fraction_is_a_reported_choice_not_a_hidden_constant():
    """The case, review R6.2: v1 hard-coded `CONJUGATE_SAFETY = 0.8`, an
    operating-point choice that appears nowhere in the approved design.
    The core enforces only the exact bound; how much margin to leave is
    an argument the caller must pass and can report.
    """

    w = 1.2
    for fraction in (0.1, 0.5, 0.9):
        slab = slab_within_conjugate_bound(
            w, fraction, dv=1.0, dx=1.0, dy=1.0)
        assert slab.du == pytest.approx(fraction * conjugate_du(w))
        PlaneWaveGeometry(slab, w)
    for bad in (0.0, 1.0, 1.5, -0.2):
        with pytest.raises(ValueError, match="fraction"):
            slab_within_conjugate_bound(w, bad, dv=1.0, dx=1.0, dy=1.0)


# --------------------------------------------------------------------
# P0: the sprinkling
# --------------------------------------------------------------------

def test_uniform_coordinate_sprinkling_has_the_covariant_intensity():
    """The case, design §4 property 1: `sqrt(-g) = 1`, so a uniform
    COORDINATE sprinkling is the covariant Poisson process exactly. If
    the count did not track `rho * coordinate volume`, the shared-points
    construction would be unfounded and with it every paired statistic.
    """

    rho = 40.0
    rng = np.random.default_rng(20260801)
    counts = [len(sprinkle(CURVED, rho, rng)) for _ in range(300)]
    expected = rho * SLAB.coordinate_volume
    se = math.sqrt(expected / len(counts))
    assert abs(float(np.mean(counts)) - expected) < 4.0 * se
    assert float(np.var(counts)) == pytest.approx(expected, rel=0.25)


def test_the_two_arms_read_the_same_points_and_disagree():
    """The case, design §4.1: the control is not a second sprinkling, it
    is the SAME points read with `A = 0`.

    Review R6.4: this test used to count with `if relation.related`,
    which folds `ambiguous` back into False -- the exact pattern §5.1
    forbids, left in the tests as if it were correct. The three
    categories are counted separately, and the totals are asserted to
    partition every ordered pair.
    """

    pts = sprinkle(CURVED, 60.0, np.random.default_rng(7))
    n = len(pts)
    assert n > 30

    tallies = {}
    for geom, name in ((CURVED, "curved"), (FLAT, "flat")):
        related = unrelated = ambiguous = 0
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                rel = causal_relation(geom, pts[i], pts[j])
                if rel.related is True:
                    related += 1
                elif rel.related is False:
                    unrelated += 1
                else:
                    ambiguous += 1
        assert related + unrelated + ambiguous == n * (n - 1)
        tallies[name] = (related, unrelated, ambiguous)

    assert tallies["curved"][0] != tallies["flat"][0], (
        "the curved and flat readings of the same points must differ, "
        "or the probe has nothing to measure")
