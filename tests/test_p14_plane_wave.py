"""Regressions for the P14 probe geometry (design §8 P1).

Every test names the case it is meant to catch, per this repository's
convention. The load-bearing one is
`test_closed_form_cost_matches_independent_minimisation`: the causal
predicate rests entirely on a closed-form transverse action, and a
closed form nobody checked against an independent computation is the
defect class P12 §10.1 refused ("the kind of unverified constant that
killed Stage C v1"). Here the independent computation is a direct
discrete minimisation that never touches the closed form.

Nothing here runs a probe. These are structural checks on the geometry.
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

from p14_plane_wave import (  # noqa: E402
    CONJUGATE_SAFETY,
    Slab,
    causal_relation,
    cost_defocusing,
    cost_focusing,
    max_slab_du,
    sprinkle,
    transverse_cost,
)

# --------------------------------------------------------------------
# the independent oracle: minimise the discretised action directly
# --------------------------------------------------------------------

def _min_action_discrete(s, ap, aq, w2, n=4000):
    """min over paths of  int_0^s (a'^2 + w2 * a^2) du / 2,  a(0)=ap, a(s)=aq.

    `w2 = +w^2` is the defocusing direction (`H` term positive), `-w^2`
    the focusing one. Discretise with `n` intervals, leave the interior
    nodes free, and solve the resulting quadratic minimisation exactly
    by one linear solve. Uses no closed form, no `sinh`, no `sin` -- it
    knows only the integrand.
    """

    h = s / n
    # unknowns a_1..a_{n-1}; action = sum h/2 * [((a_{k+1}-a_k)/h)^2
    #                                            + w2 * ((a_k+a_{k+1})/2)^2]
    m = n - 1
    main = np.full(m, 2.0 / h + w2 * h / 2.0)
    off = np.full(m - 1, -1.0 / h + w2 * h / 4.0)
    a_mat = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    rhs = np.zeros(m)
    rhs[0] = ap * (1.0 / h - w2 * h / 4.0)
    rhs[-1] = aq * (1.0 / h - w2 * h / 4.0)
    interior = np.linalg.solve(a_mat, rhs)
    a = np.concatenate(([ap], interior, [aq]))
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
    If that formula is wrong, every relation is wrong and no statistic
    downstream can reveal it -- they would all be consistently wrong
    together. So it is checked against a minimisation that shares no
    code with it: discretise the integrand and solve for the minimum.
    """

    closed = transverse_cost(s, xp, xq, yp, yq, w)
    direct = (_min_action_discrete(s, xp, xq, w * w)
              + _min_action_discrete(s, yp, yq, -w * w))
    # the discretisation is O(h^2); n = 4000 leaves ~1e-7 relative
    assert closed == pytest.approx(direct, rel=2e-6, abs=1e-12)


def test_flat_arm_is_the_straight_line_and_not_a_numerical_limit():
    """The case: `w = 0` must give `(Dx^2 + Dy^2) / (2 s)` exactly, not
    a `sinh(0)/0` that happens to survive. The flat arm is the reference
    every paired statistic is measured against, so a limit taken
    numerically there would put noise into the baseline of everything.
    """

    s, xp, xq, yp, yq = 1.7, 0.3, -0.4, 0.25, 0.6
    expected = ((xq - xp) ** 2 + (yq - yp) ** 2) / (2.0 * s)
    assert transverse_cost(s, xp, xq, yp, yq, 0.0) == expected


def test_series_branch_joins_the_closed_form_without_a_step(monkeypatch):
    """The case: both directions are `0/0` as `w s -> 0`, so each has a
    series branch. A discontinuity at the crossover would make `D`
    depend on which side of an arbitrary constant a pair fell -- the
    same defect class as P12's C31, caught before it can bite.

    Both branches are evaluated at the SAME `w`, by moving the
    threshold rather than by moving `w`. The first version of this test
    compared `w` just below the cut against `w` just above it and
    demanded 1e-12 agreement -- but those are two different physical
    costs, and the 6e-12 gap it found was the real `a^2` term responding
    to a 0.2% change in `w`. A continuity test that varies the physics
    is not a continuity test.
    """

    import p14_plane_wave as P14

    s, xp, xq, yp, yq = 1.0, 0.4, -0.3, 0.2, 0.5
    for w in (1e-5, 1e-4, 1e-3):
        monkeypatch.setattr(P14, "_SMALL_WS", 0.0)      # force closed form
        closed = P14.transverse_cost(s, xp, xq, yp, yq, w)
        monkeypatch.setattr(P14, "_SMALL_WS", 1.0)      # force the series
        series = P14.transverse_cost(s, xp, xq, yp, yq, w)
        assert math.isfinite(closed) and math.isfinite(series)
        # they agree to the series' own truncation, O((w s)^4)
        assert closed == pytest.approx(series, rel=max(1e-13, (w * s) ** 4))


def test_curvature_costs_more_in_one_direction_and_less_in_the_other():
    """The case: the profile is a QUADRUPOLE (design §4.5) -- `x`
    defocuses and `y` focuses. If both directions moved the same way the
    construction would be a scalar deformation, and `D`'s gained/lost
    split (§5) would have nothing to separate. This pins the sign
    structure the whole probe design rests on.
    """

    s, d, w = 1.0, 0.5, 1.2
    flat_leg = d * d / (2.0 * s)
    assert cost_defocusing(s, 0.0, d, w) > flat_leg   # harder in x
    assert cost_focusing(s, 0.0, d, w) < flat_leg     # easier in y


# --------------------------------------------------------------------
# the causal predicate
# --------------------------------------------------------------------

def _minkowski_related(p, q):
    """Standard flat causality, in (t, z, x, y), as an outside check."""

    du, dv = q[0] - p[0], q[1] - p[1]
    dt, dz = (du - dv) / math.sqrt(2.0), (du + dv) / math.sqrt(2.0)
    dx, dy = q[2] - p[2], q[3] - p[3]
    return dt > 0 and dt * dt >= dz * dz + dx * dx + dy * dy


@pytest.mark.parametrize("seed", range(6))
def test_flat_arm_predicate_agrees_with_minkowski_coordinates(seed):
    """The case: the flat arm must be ordinary Minkowski causality. It
    is written in null coordinates with a cost function, which is one
    algebraic slip away from a subtly different cone -- and the flat arm
    is the reference, so such a slip would bias `D` everywhere at once
    while looking self-consistent.
    """

    rng = np.random.default_rng(seed)
    for _ in range(400):
        p = rng.uniform(-1.0, 1.0, 4)
        q = rng.uniform(-1.0, 1.0, 4)
        if q[0] <= p[0]:
            p, q = q, p
        if q[0] == p[0]:
            continue
        rel = causal_relation(p, q, 0.0)
        if rel.ambiguous:
            continue
        assert rel.related == _minkowski_related(p, q)


def test_undecided_pairs_are_withheld_rather_than_called_unrelated():
    """The case, design §5.1: a pair exactly on the cone is where `D`'s
    flips live. Reporting it as `related=False` would be a silent
    misclassification biased in one direction; the design requires it to
    be `ambiguous` so it widens `[D_lower, D_upper]` instead.
    """

    w, du = 0.9, 1.0
    p = np.array([0.0, 0.0, 0.1, -0.2])
    cost = transverse_cost(du, 0.1, 0.35, -0.2, 0.3, w)
    q = np.array([du, -cost, 0.35, 0.3])       # exactly on the cone
    rel = causal_relation(p, q, w)
    assert rel.ambiguous
    assert rel.related is None

    # and a pair comfortably inside or outside is decided
    inside = np.array([du, -cost - 0.05, 0.35, 0.3])
    outside = np.array([du, -cost + 0.05, 0.35, 0.3])
    assert causal_relation(p, inside, w).related is True
    assert causal_relation(p, outside, w).related is False


def test_past_directed_and_simultaneous_pairs_are_unrelated():
    """The case: `Du <= 0` must be rejected before any cost is formed.
    `transverse_cost` divides by `s`, so a zero or negative separation
    reaching it would be a division by zero or a sign-flipped action
    that still returns a number.
    """

    p = np.array([1.0, 0.0, 0.0, 0.0])
    for q in (np.array([1.0, -5.0, 0.0, 0.0]),
              np.array([0.5, -5.0, 0.0, 0.0])):
        rel = causal_relation(p, q, 0.9)
        assert rel.related is False and not rel.ambiguous


def test_a_pair_past_the_conjugate_point_aborts_rather_than_returning():
    """The case, design §4.3: `S_y` has poles at `w s = k pi`, and past
    the first one the closed form still RETURNS A NUMBER -- a finite,
    plausible, meaningless one. The slab is supposed to make this
    unreachable; the assert is what turns a silent wrong answer into a
    stop if a caller ever builds a slab that does not.
    """

    w = 1.0
    p = np.array([0.0, 0.0, 0.0, 0.0])
    q = np.array([math.pi + 0.1, -1.0, 0.0, 0.0])
    with pytest.raises(AssertionError, match="conjugate"):
        causal_relation(p, q, w)


# --------------------------------------------------------------------
# the slab, and P0's sprinkling
# --------------------------------------------------------------------

def test_the_slab_refuses_a_geometry_outside_the_conjugate_bound():
    """The case: the conjugate bound is a property of the BOX (design
    §4.3) -- that is what lets `D` keep every pair and a fixed `C(N,2)`.
    If a slab could be built past it, the bound would silently become a
    per-pair filter and the denominator would move.
    """

    w = 1.5
    assert max_slab_du(w) == pytest.approx(CONJUGATE_SAFETY * math.pi / w)
    Slab(du=max_slab_du(w) * 0.5, dv=1.0, dx=1.0, dy=1.0).validate(w)
    with pytest.raises(ValueError, match="conjugate bound"):
        Slab(du=max_slab_du(w) * 1.01, dv=1.0, dx=1.0, dy=1.0).validate(w)


def test_uniform_coordinate_sprinkling_has_the_covariant_intensity():
    """The case, design §4 property 1: `sqrt(-g) = 1`, so a uniform
    COORDINATE sprinkling is the covariant Poisson process exactly --
    and both arms can share one point set. If the count did not track
    `rho * coordinate volume` the shared-points construction would be
    unfounded, and with it every paired statistic.
    """

    slab = Slab(du=1.0, dv=2.0, dx=1.5, dy=1.5)
    rho = 40.0
    rng = np.random.default_rng(20260801)
    counts = [len(sprinkle(slab, rho, rng)) for _ in range(300)]
    expected = rho * slab.coordinate_volume
    mean = float(np.mean(counts))
    # Poisson: se of the mean is sqrt(expected/300)
    se = math.sqrt(expected / len(counts))
    assert abs(mean - expected) < 4.0 * se
    # and the variance tracks the mean, as a Poisson count must
    assert float(np.var(counts)) == pytest.approx(expected, rel=0.25)


def test_the_two_arms_see_the_same_points():
    """The case, design §4.1: the control is not a second sprinkling, it
    is the SAME points read with `A = 0`. This is what removes the twin
    calibration and the cross-arm gate P12 needed, so it deserves a test
    rather than a sentence: the sprinkling must not depend on `w` at
    all.
    """

    slab = Slab(du=0.8, dv=1.0, dx=1.0, dy=1.0)
    pts = sprinkle(slab, 60.0, np.random.default_rng(7))
    # the same array is handed to both arms; only the relation differs
    n_curved = sum(
        1 for i in range(len(pts)) for j in range(len(pts))
        if i != j and causal_relation(pts[i], pts[j], 1.0).related)
    n_flat = sum(
        1 for i in range(len(pts)) for j in range(len(pts))
        if i != j and causal_relation(pts[i], pts[j], 0.0).related)
    assert n_curved != n_flat, (
        "the curved and flat readings of the same points must differ, "
        "or the probe has nothing to measure")
