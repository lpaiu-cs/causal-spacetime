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
    GuardInsets,
    PlaneWaveGeometry,
    Slab,
    _defocusing_inset,
    _exact_cost,
    _focusing_excursion_factor,
    arms,
    causal_relation,
    class_c_eligible,
    conjugate_du,
    cost_defocusing,
    cost_error_bound,
    cost_focusing,
    guard_insets,
    guard_margins,
    max_negative_focusing,
    slab_within_conjugate_bound,
    sprinkle,
    trajectory_excursion,
    transverse_cost,
    worst_focusing_corner,
)

# a slab comfortably inside the conjugate bound, reused below
W = 1.0
SLAB = slab_within_conjugate_bound(W, 0.6, dv=2.0, dx=1.5, dy=1.5)
CURVED, FLAT = arms(SLAB, W)


# --------------------------------------------------------------------
# the independent oracle: minimise the discretised action directly
# --------------------------------------------------------------------

def _solve_tridiagonal(sub, diag, sup, rhs):
    """Thomas algorithm. `O(n)` time and memory, no dense matrix.

    Review R8.2: the first version of the oracle built the Hessian with
    `np.diag`, which materialises a dense `(n-1) x (n-1)` array -- 128 MB
    at `n = 4000`, several of them per call counting temporaries -- and
    then ran a dense LU at `O(n^3)`. The Hessian is tridiagonal, so none
    of that was necessary, and an oracle that can exhaust a runner's
    memory is not one anybody will keep running.
    """

    n = len(diag)
    c = np.empty(n - 1)
    d = np.empty(n)
    c[0] = sup[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for i in range(1, n):
        denom = diag[i] - sub[i - 1] * c[i - 1]
        if i < n - 1:
            c[i] = sup[i] / denom
        d[i] = (rhs[i] - sub[i - 1] * d[i - 1]) / denom
    x = np.empty(n)
    x[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


def _min_action_discrete(s, ap, aq, w2, n=20000):
    """min over paths of  int_0^s (a'^2 + w2 a^2) du / 2,  a(0)=ap, a(s)=aq.

    Discretise, leave the interior nodes free, solve the resulting
    quadratic minimisation exactly. Uses no closed form, no `sinh`, no
    `sin` -- it knows only the integrand.

    `n` is five times the dense version's and still far cheaper, because
    the solve is now tridiagonal: the discretisation is `O(h^2)`, so the
    extra nodes buy about `25x` in oracle accuracy for `O(n)` work.
    """

    h = s / n
    m = n - 1
    diag = np.full(m, 2.0 / h + w2 * h / 2.0)
    off = np.full(m - 1, -1.0 / h + w2 * h / 4.0)
    rhs = np.zeros(m)
    edge = 1.0 / h - w2 * h / 4.0
    rhs[0] = ap * edge
    rhs[-1] = aq * edge
    interior = _solve_tridiagonal(off, diag, off, rhs)
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
    If it is wrong, every relation is wrong TOGETHER and no downstream
    statistic can reveal it. Checked against a minimisation sharing no
    code with it.
    """

    closed = transverse_cost(s, xp, xq, yp, yq, w)
    direct = (_min_action_discrete(s, xp, xq, w * w)
              + _min_action_discrete(s, yp, yq, -w * w))
    # measured worst agreement over these cases is 2.4e-10, which is the
    # oracle's `O(h^2)` and not the closed form's; the tolerance is set
    # from that measurement rather than left at the dense version's 2e-6
    assert closed == pytest.approx(direct, rel=1.2e-9, abs=1e-12)


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


def test_error_bound_survives_proximity_to_the_conjugate_point():
    """The case, review R8.1, and the third face of the same defect. Near
    a conjugate point `sin(a)` is nearly zero, so the ROUNDING OF THE
    ARGUMENT `a = w s` is amplified twice over: once through
    `d(1/sin a) = |cos a| da / sin^2 a`. The bound had only the operand
    scale, which carries a single `1/sin a`, so it missed the second
    factor entirely.

    Measured at `w = 0.9`, `s = pi/0.9 - 1e-4` (`sin a = 9e-5`),
    endpoints `(0.3, -0.2, 0.4, 0.1)`: the true error is `1.80e-9` and
    the old bound reported `4.44e-12`, understating it 405-fold and
    deciding pairs on the wrong side of the cone without escalating.

    A slab this close to the bound is legal -- `PlaneWaveGeometry` only
    forbids reaching it -- so the bound has to cover it rather than the
    geometry being trusted to stay away.
    """

    w = 0.9
    s = math.pi / w - 1e-4
    xp, xq, yp, yq = 0.3, -0.2, 0.4, 0.1

    assert 0.0 < math.sin(w * s) < 1e-4      # really is near the pole
    approx = transverse_cost(s, xp, xq, yp, yq, w)
    getcontext().prec = 60
    actual = abs(Decimal(repr(approx)) - _exact_cost(s, xp, xq, yp, yq, w))
    bound = cost_error_bound(s, xp, xq, yp, yq, w)

    assert float(actual) > 1e-10, "the case must actually be hard"
    assert Decimal(repr(bound)) >= actual, (
        f"bound {bound:.3e} does not cover the actual error "
        f"{float(actual):.3e}")

    # and the bound is not vacuous far from the pole: it stays small
    # where the argument sensitivity is not amplified
    tame = cost_error_bound(1.0, xp, xq, yp, yq, 0.5)
    assert tame < 1e-14


# --------------------------------------------------------------------
# and the systematic version, because four rounds of counterexamples is
# enough evidence that waiting for the fifth is the wrong method
# --------------------------------------------------------------------

def _regime_sampler(rng, n):
    """Configurations spanning every regime that has broken the bound.

    Reviews R6.1 through R9.1 each found the bound in a place the
    previous fix did not reach: legs cancelling, the series' `a^2` terms
    cancelling between legs, the argument amplified near a conjugate
    point, and the separations rounded before exact arithmetic saw them.
    Each was repaired where the counterexample pointed. This draws from
    all of those regions at once instead.
    """

    for _ in range(n):
        kind = rng.integers(0, 6)
        if kind == 5:
            # deliberately near-cancelling: S_x > 0 always, and S_y < 0
            # whenever the y endpoints coincide, so solve for the x
            # separation that almost annihilates the sum. Perfect
            # cancellation is measure zero, which is why the uniform
            # regimes above hit R6.1's case only twice in 1400 draws --
            # the defect that started this series needs to be aimed at.
            s = float(rng.uniform(0.2, 2.0))
            w = float(rng.uniform(0.2, 1.4))
            if s >= math.pi / w:
                continue
            y = float(rng.normal(0.0, 1.0)) or 1.0
            leg_y = cost_focusing(s, y, y, w)
            coth = math.cosh(w * s) / math.sinh(w * s)
            target = -2.0 * leg_y / (w * coth)
            if target <= 0.0:
                continue
            xq = math.sqrt(target) * float(rng.uniform(0.999999, 1.000001))
            yield s, 0.0, xq, y, y, w
            continue
        if kind == 0:                                    # flat
            w = 0.0
            s = float(rng.uniform(1e-3, 5.0))
        elif kind == 1:                                  # series branch
            s = float(rng.uniform(0.1, 3.0))
            w = float(rng.uniform(1e-9, 0.9e-4)) / s
        elif kind == 2:                                  # mid range
            s = float(rng.uniform(0.1, 3.0))
            w = float(rng.uniform(0.05, 2.0))
        elif kind == 3:                                  # near conjugate
            w = float(rng.uniform(0.3, 2.0))
            s = (math.pi / w) * float(rng.uniform(0.90, 0.999999))
        else:                                            # near-coincident
            s = float(rng.uniform(0.1, 3.0))
            w = float(rng.uniform(0.0, 2.0))
        scale = 10.0 ** float(rng.uniform(-4.0, 1.0))
        xp = float(rng.normal(0.0, scale))
        yp = float(rng.normal(0.0, scale))
        if kind == 4:                                    # endpoints close
            xq = xp + float(rng.normal(0.0, scale * 1e-8))
            yq = yp + float(rng.normal(0.0, scale * 1e-8))
        else:
            xq = float(rng.normal(0.0, scale))
            yq = float(rng.normal(0.0, scale))
        if w > 0.0 and s >= math.pi / w:
            continue
        yield s, xp, xq, yp, yq, w


def test_the_error_bound_dominates_over_every_regime_that_has_broken_it():
    """The case: each of R6.1, R7.1, R8.1 and R9.1 found this bound
    wrong somewhere the previous repair did not reach, and each was
    fixed at the point the counterexample happened to expose. Four
    rounds of that is evidence the method is wrong, not just the code:
    a bound has to be checked against the thing it bounds, everywhere,
    rather than patched wherever the last reviewer looked.

    So: draw configurations from all four regimes and assert the bound
    dominates the measured error in every one. This is what should have
    been written after R6.1.

    NEGATIVE CONTROL, run before trusting this. Each historical bound
    model was reinstated and swept, and this test catches all three --
    a property test that passes for every version of the code it
    polices would be worth nothing:

        v1 (the R6.1 model)   271 violations / 1423, worst 7.1e4 too small
        v2 (the R7.1 model)    56 violations / 1423, worst 1.5e3 too small
        v3 (the R8.1 model)    17 violations / 1423, worst 9.7e2 too small
        current                 0 violations / 1423

    v1 registered only twice before the deliberately near-cancelling
    regime was added, which is itself the lesson: uniform sampling does
    not find a measure-zero configuration, so the regime that broke the
    code has to be aimed at rather than hoped for.
    """

    getcontext().prec = 60
    rng = np.random.default_rng(20260802)
    checked = failures = 0
    worst_ratio = 0.0
    for s, xp, xq, yp, yq, w in _regime_sampler(rng, 1500):
        approx = transverse_cost(s, xp, xq, yp, yq, w)
        if not math.isfinite(approx):
            continue
        exact = _exact_cost(s, xp, xq, yp, yq, w)
        bound = cost_error_bound(s, xp, xq, yp, yq, w)
        actual = abs(Decimal(approx) - exact)
        checked += 1
        if actual > 0:
            worst_ratio = max(worst_ratio, float(actual) / max(bound, 1e-300))
        if Decimal(bound) < actual:
            failures += 1
            if failures == 1:
                pytest.fail(
                    f"bound violated at s={s!r} w={w!r} "
                    f"x=({xp!r},{xq!r}) y=({yp!r},{yq!r}): "
                    f"bound {bound:.3e} < actual {float(actual):.3e}")
    assert checked > 1000, f"sampler produced only {checked} usable cases"
    assert failures == 0
    # and the bound is not vacuous: it should not be astronomically
    # larger than the error everywhere, or it would escalate every pair
    assert worst_ratio > 1e-6, (
        f"worst error/bound ratio {worst_ratio:.2e} suggests the bound "
        "is so loose it carries no information")


def test_the_verdict_agrees_with_exact_arithmetic_on_pairs_at_the_cone():
    """The case: the bound existing is not the claim -- the claim is that
    a returned verdict is right. Pairs are placed at the cone and nudged
    by amounts that straddle the double path's resolution, so most of
    them are decided by escalation, and every non-ambiguous verdict is
    checked against exact arithmetic that rebuilds the separations
    itself (R9.1).
    """

    getcontext().prec = 60
    rng = np.random.default_rng(31415)
    decided = escalated = 0
    for s, xp, xq, yp, yq, w in _regime_sampler(rng, 260):
        cost = transverse_cost(s, xp, xq, yp, yq, w)
        if not math.isfinite(cost):
            continue
        # the slab must hold the pair and stay inside the conjugate
        # bound; the sampler guarantees s < pi/w, so split the gap
        du = s * 1.5 if w == 0.0 else s + 0.5 * (math.pi / w - s)
        geom = PlaneWaveGeometry(
            Slab(du=du, dv=1.0, dx=1.0, dy=1.0), w)
        for nudge in (-1e-16, -1e-13, 0.0, 1e-13, 1e-16):
            p = np.array([0.0, 0.0, xp, yp])
            q = np.array([s, -cost + nudge, xq, yq])
            rel = causal_relation(geom, p, q)
            if rel.ambiguous:
                continue
            du_exact = Decimal(float(q[0])) - Decimal(float(p[0]))
            dv_exact = Decimal(float(q[1])) - Decimal(float(p[1]))
            truth = -(dv_exact
                      + _exact_cost(du_exact, xp, xq, yp, yq, w)) > 0
            assert rel.related == truth, (
                f"verdict {rel.related} != exact {truth} at s={s!r} "
                f"w={w!r} nudge={nudge!r}")
            decided += 1
            escalated += int(rel.escalated)
    assert decided > 400, f"only {decided} decided pairs exercised"
    assert escalated > 20, (
        f"only {escalated} escalations: the fixture is not reaching the "
        "regime where the double path has to hand over")


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


def test_escalation_rebuilds_the_separations_instead_of_inheriting_them():
    """The case, review R9.1, and the worst of the series: the double
    path escalated correctly and the ESCALATION returned the wrong
    answer, because it took `du` and `dv` as doubles that had already
    rounded in `q[i] - p[i]`.

    Flat arm, inside `Slab(du=1, dv=201, dx=10, dy=1)`. The `u`
    subtraction loses `3.1e-19`, and the `1/du` in the cost turns that
    into `1.2e-15` against a true margin of `+8.1e-15` -- so the verdict
    came back `False` with an `escalated=True` certificate on it. An
    escalation that begins by rounding its inputs is worse than no
    escalation: it certifies the wrong answer.
    """

    p = np.array([1.955619748358592e-07, 199.803718598693, 0.0, 0.0])
    q = np.array([0.053271600135298, 0.0, 4.613854078473958, 0.0])
    geom = PlaneWaveGeometry(Slab(du=1.0, dv=201.0, dx=10.0, dy=1.0), 0.0)

    getcontext().prec = 60
    du_exact = Decimal(float(q[0])) - Decimal(float(p[0]))
    dv_exact = Decimal(float(q[1])) - Decimal(float(p[1]))
    dx = Decimal(float(q[2])) - Decimal(float(p[2]))
    exact_margin = -(dv_exact + dx * dx / (2 * du_exact))

    # the case must really be one the double path cannot settle
    assert abs(float(exact_margin)) < 1e-13
    assert exact_margin > 0

    rel = causal_relation(geom, p, q)
    assert rel.escalated
    assert rel.related is True
    assert rel.margin_lower == pytest.approx(float(exact_margin), rel=1e-9)

    # the subtraction really does lose digits -- this is not a fixture
    # that would pass for an unrelated reason
    lost = abs(Decimal(float(q[0]) - float(p[0])) - du_exact)
    assert lost > 0


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


# --------------------------------------------------------------------
# P1: the Class C guard insets
# --------------------------------------------------------------------

def test_the_v_inset_matches_the_constrained_maximum_of_minus_the_cost():
    """The case, review R10.1: a diamond can rise above its own past
    endpoint because `C` goes negative, and the guard has to bound that.
    The closed form claims `max(-S_y) = w Y^2 tan(w s / 2)` at the box
    corner. Checked against a direct grid search, and separately that
    the defocusing direction never lifts `v` at all.

    The rectangle case of R12 has three regimes rather than one -- the
    interior optimum `y = y_p / cos a` is feasible only while
    `y_p <= Y cos a` -- so it is checked over random rectangles too,
    where getting a regime boundary wrong would show as the grid
    beating the formula.
    """

    for w in (0.3, 0.8, 1.4):
        for half_y in (0.4, 1.0):
            for frac in (0.2, 0.5, 0.9):
                s = frac * math.pi / w
                predicted = w * half_y ** 2 * math.tan(w * s / 2.0)
                grid = np.linspace(-half_y, half_y, 201)
                best = max(-cost_focusing(s, float(yp), float(y), w)
                           for yp in grid for y in grid)
                assert best == pytest.approx(predicted, rel=1e-9)
                # the symmetric case must be the formula's own corner
                assert max_negative_focusing(
                    w * s, w, half_y, half_y) == pytest.approx(
                        predicted, rel=1e-12)
                # and x can never contribute a rise
                worst_x = max(-cost_defocusing(s, float(xp), float(x), w)
                              for xp in grid for x in grid)
                assert worst_x <= 0.0

    rng = np.random.default_rng(20260802)
    worst_over = 0.0
    saw_interior = saw_corner = 0
    for _ in range(120):
        a = float(rng.uniform(0.05, 0.95)) * math.pi
        w = float(rng.uniform(0.3, 1.5))
        box = float(rng.uniform(0.5, 3.0))
        anchor = box * float(rng.uniform(0.05, 1.0))
        predicted = max_negative_focusing(a, w, anchor, box)
        if math.cos(a) > 0.0 and anchor <= box * math.cos(a):
            saw_interior += 1
        else:
            saw_corner += 1
        gp = np.linspace(-anchor, anchor, 121)
        gq = np.linspace(-box, box, 121)
        best = max(-cost_focusing(a / w, float(yp), float(y), w)
                   for yp in gp for y in gq)
        worst_over = max(worst_over,
                         (best - predicted) / max(abs(predicted), 1e-12))
    assert worst_over < 1e-12, (
        f"a grid beat the rectangle formula by {worst_over:.3e} relative")
    assert saw_interior > 10 and saw_corner > 10, (
        f"only one regime exercised: interior={saw_interior} "
        f"corner={saw_corner}")


def test_the_transverse_inset_is_the_alexandrov_radius_when_flat():
    """The case: the transverse guard is a chain of inequalities, so its
    one checkable anchor is the flat limit, where the answer is known
    exactly -- a diamond of proper time `tau` has transverse radius
    `tau/2`, and `tau^2 = 2 du dv`, giving `sqrt(du dv / 2)`. If the
    bound did not reduce to that, the inequalities would be wrong.
    """

    for du, dv in ((1.0, 1.0), (0.5, 2.0), (2.0, 0.3)):
        geom = PlaneWaveGeometry(Slab(du=du, dv=dv, dx=9.0, dy=9.0), 0.0)
        insets = guard_insets(geom)
        assert insets.v == 0.0
        assert insets.x == pytest.approx(math.sqrt(du * dv / 2.0))
        assert insets.y == pytest.approx(math.sqrt(du * dv / 2.0))


def test_the_focusing_direction_needs_its_own_inset():
    """The case, review R11.1: the first guard derived ONE transverse
    bound from `coth(z) >= 1/z` -- a fact about the DEFOCUSING direction
    -- and applied it to both. In the focusing direction the same slot
    holds `cot(z)`, where the inequality is false and the value is
    negative past `z = pi/2`, so nothing had been proved about `y` at
    all. A random containment sweep passing is not a proof, and this
    guard is about to be the ground-truth label for the order-invariant
    eligibility candidate.

    The excursion factor rises from `1.02 Y` at `w du = 0.4` to
    `5.88 Y` at `2.8`, so the geometry really does focus.

    What is NOT pinned here is which inset is larger. I asserted a
    fixed ordering twice and had it backwards both times; with the
    coupled solve of R12.1 it is regime-dependent, because `x` is
    computed from a budget that itself contains the `v` excursion that
    the `y` solve determines. There is no ordering to pin. What can be
    pinned is that the two are derived apart rather than one copied
    onto the other.
    """

    for a in (0.4, 1.0, 2.0, 2.8):
        assert _focusing_excursion_factor(a) > 1.0
    assert _focusing_excursion_factor(2.8) > 5.0
    assert (_focusing_excursion_factor(2.8)
            > _focusing_excursion_factor(0.4))

    for w, du in ((0.6, 0.5), (1.1, 0.9), (0.4, 0.3)):
        geom = PlaneWaveGeometry(Slab(du=du, dv=0.4, dx=3.0, dy=3.0), w)
        insets = guard_insets(geom)
        assert insets.x != insets.y, (
            f"w={w} du={du}: the two insets are equal, so they were not "
            "derived apart")


def test_the_asymmetric_v_bound_is_what_makes_the_guard_usable():
    """The case, review R12.1: the first guard bounded `-S_y` over the
    full box in BOTH arguments and solved the three components apart.
    That is why a table of `lim_x` and `lim_y` could report a trend
    while every row was already empty in `v` -- the analysis never
    consumed the component it had just added.

    The rectangle is the honest domain: the anchor is an eligible
    endpoint, hence already inset, while the intermediate point ranges
    over the whole box. This test pins that the distinction is
    load-bearing rather than cosmetic -- the symmetric bound empties
    all three geometries below, and the asymmetric one does not.
    """

    def admits(w, du):
        slab = Slab(du=du, dv=0.4, dx=3.0, dy=3.0)
        insets = guard_insets(PlaneWaveGeometry(slab, w))
        symmetric = w * (slab.dy / 2.0) ** 2 * math.tan(w * du / 2.0)
        eligible = (insets.x < 1.5 and insets.y < 1.5
                    and 2.0 * insets.v < slab.dv)
        return eligible, symmetric <= slab.dv / 2.0, insets.v, symmetric

    for w, du, sym_ratio in ((0.6, 0.5, 3.0), (1.1, 0.9, 20.0)):
        eligible, sym_ok, v_inset, symmetric = admits(w, du)
        assert not sym_ok, (
            f"w={w} du={du}: the symmetric bound admits this geometry, so "
            "it is not a case where the asymmetry matters")
        assert eligible, (
            f"w={w} du={du}: the asymmetric bound should have rescued this")
        assert v_inset <= symmetric / sym_ratio, (
            f"w={w} du={du}: asymmetric v={v_inset} saves less than "
            f"{sym_ratio}x on the symmetric {symmetric}")

    # It is not unconditional rescue: the saving is a bounded factor,
    # so a geometry far enough past the bound stays empty either way.
    assert not admits(0.9, 2.0)[0] and not admits(0.9, 2.0)[1]


def test_where_class_c_eligibility_is_non_empty_over_the_dimensionless_box():
    """The case, review R12.1 and R12.2: §4.6.3 reported a table of two
    margins at one box size, indexed by `a = w du` alone, and concluded
    from it that eligibility empties near `a ~ 1` because the
    defocusing condition closes first. Both halves were wrong. The
    table omitted `v`; and `a` does not index the problem, because
    scaling `w -> lambda w` with `du -> du/lambda` holds `a` fixed
    while changing `w dv` and `w dx/2`, which move the margins.

    Fixing `w = 1` (lengths in units of `1/w`) leaves four
    dimensionless quantities, `a`, `w dv`, `w dx/2`, `w dy/2`, and this
    pins the corrected picture over them.

    **Counted per condition, not per cell (R13.2).** The first version
    of this census recorded only the WORST margin in each empty cell,
    and I read `x` never being worst as `x` never binding -- concluding
    that tightening it would buy nothing, when it was in fact the
    ceiling on `a`. A condition can be violated in most of the empty
    cells and be the worst one in none of them. The counts below
    overlap on purpose.

    **And `X` and `Y` vary independently (R13.4).** Tying them makes
    the sweep a three-dimensional slice through a four-dimensional
    problem, which cannot support a claim about which direction binds.
    Untied, the box that admits is not square.
    """

    def margins(a, w_dv, w_x, w_y):
        slab = Slab(du=a, dv=w_dv, dx=2 * w_x, dy=2 * w_y)
        return guard_margins(PlaneWaveGeometry(slab, 1.0))

    def admits(a, w_dv, w_x, w_y):
        slab = Slab(du=a, dv=w_dv, dx=2 * w_x, dy=2 * w_y)
        insets = guard_insets(PlaneWaveGeometry(slab, 1.0))
        return (insets.x < w_x and insets.y < w_y
                and 2.0 * insets.v < slab.dv)

    # The slice §4.6.3 reported. The withdrawn claim said `lim_x` goes
    # negative at `a = 1.0` and eligibility empties there; it does not,
    # and it does not.
    eligible = {a: admits(a, 0.2, 3.0, 3.0) for a in (0.3, 0.6, 1.0, 1.6)}
    assert eligible == {0.3: True, 0.6: True, 1.0: True, 1.6: False}, (
        f"the corrected slice changed: {eligible}")
    # and `a = 1.6` fails on `v`, with both transverse margins positive
    m_x, m_y, m_v = margins(1.6, 0.2, 3.0, 3.0)
    assert m_x > 0.0 and m_y > 0.0 and m_v <= 0.0

    axes = (0.2, 0.4, 0.7, 1.0, 1.4, 2.0)
    dvs = (0.2, 1.0, 4.0, 16.0)
    halves = (0.3, 1.0, 3.0)

    # R14.3: all three sole-blocker counts, not just `x`. Pinning one
    # of them let a wrong value for another sit in the design table --
    # `v_only` was recorded as 0 and is 5, and "sole blocker in 0" is
    # exactly the reading that would justify dropping the `v` guard.
    def census(cells):
        counts = {"cells": 0, "ok": 0, "x": 0, "y": 0, "v": 0,
                  "x_only": 0, "y_only": 0, "v_only": 0}
        for a, w_dv, w_x, w_y in cells:
            counts["cells"] += 1
            if admits(a, w_dv, w_x, w_y):
                counts["ok"] += 1
                continue
            bad = [k for k, m in zip("xyv", margins(a, w_dv, w_x, w_y),
                                     strict=True)
                   if m <= 0.0]
            for key in bad:
                counts[key] += 1
            if len(bad) == 1:
                counts[bad[0] + "_only"] += 1
        return counts

    square = census([(a, w_dv, h, h)
                     for a in axes for w_dv in dvs for h in halves])
    assert square == {"cells": 72, "ok": 34, "x": 34, "y": 33, "v": 4,
                      "x_only": 1, "y_only": 0, "v_only": 4}, square

    rect = census([(a, w_dv, w_x, w_y) for a in axes for w_dv in dvs
                   for w_x in halves for w_y in halves])
    assert rect == {"cells": 216, "ok": 60, "x": 102, "y": 111, "v": 12,
                    "x_only": 34, "y_only": 48, "v_only": 5}, rect

    # the anisotropy is real: at this cell a wide `y` extent admits and
    # a wide `x` extent does not, which a square sweep cannot see
    assert admits(1.0, 1.0, 1.0, 3.0)
    assert not admits(1.0, 1.0, 3.0, 1.0)


def test_class_c_eligibility_cannot_be_taken_from_the_flat_arm():
    """The case, review R14.1: §4.6 says the eligibility rule is the
    same in both arms and comes from the curved one, because two arms
    with different eligible sets are not the same points and the
    paired comparison §4.1 exists for is gone. That was a sentence in
    a docstring; the function derived its insets from whatever geometry
    it was handed. A paired analysis passing each arm its own geometry
    got the flat arm's insets in the flat arm -- flatness collapses the
    guard to the Alexandrov radius, so they are much smaller -- and
    pairs would have entered one arm's denominator and not the other's.
    """

    slab = Slab(du=1.0, dv=0.4, dx=3.0, dy=3.0)
    curved, flat = arms(slab, 1.0)

    # the two guards really do differ, so this is not a moot contract
    assert guard_insets(flat).y < guard_insets(curved).y / 2.0

    # a point eligible under the flat guard but not the curved one
    y_flat = 0.5 * (1.5 - guard_insets(flat).y)
    p = np.array([0.0, 0.3, 0.0, y_flat])
    q = np.array([slab.du, 0.1, 0.0, y_flat])
    assert abs(y_flat) < 1.5 - guard_insets(flat).y
    assert abs(y_flat) > 1.5 - guard_insets(curved).y

    # both arms refuse it, because both read the curved guard
    assert not class_c_eligible(curved, p, q)
    assert not class_c_eligible(flat, p, q)

    # and a flat geometry with no curved partner is refused outright
    # rather than silently answering from its own insets
    orphan = PlaneWaveGeometry(slab, 0.0)
    with pytest.raises(ValueError, match="flat geometry"):
        class_c_eligible(orphan, p, q)

    # the link is checked, not trusted: a partner on a different slab
    # would mean the two arms are not the same points
    other = PlaneWaveGeometry(Slab(du=1.0, dv=0.5, dx=3.0, dy=3.0), 1.0)
    with pytest.raises(ValueError, match="same points"):
        PlaneWaveGeometry(slab, 0.0, guard_from=other)


def test_the_defocusing_chain_was_the_ceiling_on_a():
    """The case, review R13.2: the previous round concluded that the
    `x` inset never binds and would not be tightened. That was read off
    a census recording only the worst condition per cell, and it was
    backwards. With the old bound `cross = 2 w tanh(a/2) X^2` and the
    budget `dv + cross`,

        x_inset / X >= sqrt(du * cross / 2) / X = sqrt(a tanh(a/2)),

    with `dv` and the box size dropping out entirely. That exceeds 1 at
    `a = 1.5434`, so past there NO square box admitted a pair, whatever
    else was done -- while the effect being measured grows as `(wT)^4`.
    A ceiling on `a` is the most expensive thing this design can have,
    and it was hiding behind a statistic that never showed `x` as worst.

    Tightened the same two ways the rest of the guard already was --
    the anchor is an eligible endpoint, so the bilinear term carries
    `X_inner` and not `X`; and `tanh(a_p/2) + tanh(a_q/2)` over a free
    split peaks at `2 tanh(a/4)` by concavity, not `2 tanh(a/2)` -- the
    ratio drops below 1 everywhere tested and `a = 2` becomes
    reachable.
    """

    def old_ratio(a):
        return math.sqrt(a * math.tanh(a / 2.0))

    assert old_ratio(1.5) < 1.0 < old_ratio(1.6)
    assert old_ratio(2.0) == pytest.approx(1.2342, abs=5e-5)
    lo, hi = 1.0, 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid * math.tanh(mid / 2.0) < 1.0:
            lo = mid
        else:
            hi = mid
    assert lo == pytest.approx(1.5434, abs=5e-5)

    # the tightened inset, in the box-only limit where `dv` is
    # negligible, stays inside the box well past that
    for a, ceiling in ((1.6, 0.54), (2.0, 0.61), (2.5, 0.68), (3.0, 0.73)):
        slab = Slab(du=a, dv=1e-9, dx=2.0, dy=2.0)
        x_inset, _ = _defocusing_inset(slab, 1.0, 1.0)
        assert x_inset < ceiling, f"a={a}: x_inset/X = {x_inset}"
        assert x_inset < old_ratio(a)

    # and it is not merely a smaller number: `a = 2` now admits pairs,
    # which under the old bound was impossible at every box size
    admitting = []
    for w_dv in (0.1, 0.5, 2.0, 8.0):
        for w_x in (0.5, 1.0, 3.0, 9.0):
            for w_y in (0.5, 1.0, 3.0, 9.0):
                slab = Slab(du=2.0, dv=w_dv, dx=2 * w_x, dy=2 * w_y)
                ins = guard_insets(PlaneWaveGeometry(slab, 1.0))
                if (ins.x < w_x and ins.y < w_y
                        and 2.0 * ins.v < slab.dv):
                    admitting.append((w_dv, w_x, w_y))
    assert admitting, "a = 2 still admits nothing"

    # R14.4: and NO upper endpoint in `a` is claimed. A previous
    # version of the design read a miss on a coarse grid as a bound at
    # `2.4`; these three boxes admit, and the two larger ones need a
    # `w dv` two orders below anything that grid contained. Failing to
    # find a box is not evidence that none exists -- the only ceiling
    # established anywhere is the conjugate point, `a < pi`.
    for a, w_dv, w_x, w_y in ((2.4, 0.2, 0.6, 0.4),
                              (2.6, 0.005, 0.0835, 0.0519),
                              (2.8, 0.00875, 0.1112, 0.0519)):
        slab = Slab(du=a, dv=w_dv, dx=2 * w_x, dy=2 * w_y)
        ins = guard_insets(PlaneWaveGeometry(slab, 1.0))
        assert ins.x < w_x and ins.y < w_y and 2.0 * ins.v < slab.dv, (
            f"a={a} at ({w_dv}, {w_x}, {w_y}) no longer admits: {ins}")

    # the fixed point is a fixed point: `x_inset` reproduces itself
    for a, dv, half in ((0.7, 1.0, 2.0), (1.5, 0.3, 1.0), (2.5, 4.0, 5.0)):
        slab = Slab(du=a, dv=dv, dx=2 * half, dy=2 * half)
        x_inset, cross = _defocusing_inset(slab, 1.0, half)
        assert x_inset < half, f"a={a} dv={dv} half={half}: x closed"
        assert cross == pytest.approx(
            2.0 * math.tanh(a / 4.0) * (half - x_inset) * half)
        assert x_inset == pytest.approx(
            math.sqrt(slab.du * (slab.dv + cross) / 2.0), rel=1e-9)


def test_the_trajectory_excursion_is_exact_and_maximised_at_a_corner():
    """The case: the focusing inset rests on two structural claims --
    that `max_u |y(u)|` is computed exactly rather than sampled, and
    that its maximum over the endpoint box sits at a corner because the
    function is convex there. If the corner claim were wrong the inset
    would be too small everywhere, silently.
    """

    for a in (0.4, 1.0, 2.0, 2.8):
        grid = np.linspace(-1.0, 1.0, 121)
        best_interior = max(
            trajectory_excursion(float(p), float(q), a)
            for p in grid for q in grid)
        assert _focusing_excursion_factor(a) >= best_interior * (1 - 1e-12)

    # and the exact maximum really is the trajectory's, not an endpoint
    # value, wherever an interior peak falls inside the segment
    a, yp, yq = 2.8, 1.0, -1.0
    dense = max(abs(yp * math.cos(t) + ((yq - yp * math.cos(a))
                                        / math.sin(a)) * math.sin(t))
                for t in np.linspace(0.0, a, 200001))
    assert trajectory_excursion(yp, yq, a) == pytest.approx(dense, rel=1e-9)


def test_no_diamond_point_escapes_the_box_when_the_guard_admits_the_pair():
    """The case: the guard's whole purpose is that a Class C count is
    not truncated. That is a claim about points, so it is tested on
    points -- sample candidates from a region THREE TIMES the box, keep
    the ones that are causally between an eligible pair, and require
    every one of them to be inside the box.

    Sampling from the box itself would be circular: everything drawn
    would be inside by construction and the test would pass however
    wrong the inset was.

    NEGATIVE CONTROL, and it says something the guard's derivation did
    not. Shrinking each component and re-running this sweep:

        real insets        1724 interior points,    0 escapes
        transverse x0.9    1558,                    0
        transverse x0.7    1230,                    0
        transverse x0.5    1018,                    1  [caught]
        v inset x0.5       1748,                    0
        v inset x0.0       1783,                    0

    Two of those rows are executed below rather than quoted; the rest
    are a record.

    So this test polices the transverse component only loosely -- the
    bound is a chain of inequalities and is genuinely conservative --
    and **does not police the `v` component at all** in these
    geometries, where the excursion the guard reserves for it is
    `0.003` and `0.008` against `dv = 0.35`. That is now a statement
    about the geometries and not about the guard: the `v` component is
    computed over the honest rectangle since R12.1, and the earlier
    version of this docstring called the full-box bound "not worth
    iterating" one round before it turned out to be the difference
    between an empty eligible set and a usable one.

    What that leaves: the `v` formula itself is verified directly by
    `test_the_v_inset_matches_the_constrained_maximum_of_minus_the_cost`
    against a grid search, and this test verifies the two together do
    contain the diamond. Neither establishes that the `v` component is
    NECESSARY at the geometries the probe will use. §8 P1 should report
    whether it ever binds; if it never does, that is a simplification
    the design can take, and if it does, this control will start
    catching it.
    """

    interior, escapes = _containment_sweep(guard_insets)
    assert interior > 100, (
        f"only {interior} diamond points found; the sweep is not "
        "exercising containment")
    assert escapes == 0, f"{escapes}/{interior} diamond points escaped"

    def shrunk(fx=1.0, fv=1.0):
        def insets_fn(geom):
            real = guard_insets(geom)
            return GuardInsets(x=real.x * fx, y=real.y * fx, v=real.v * fv)
        return insets_fn

    caught_interior, caught = _containment_sweep(shrunk(fx=0.5))
    assert caught > 0, (
        f"halving the transverse inset leaked nothing in "
        f"{caught_interior} points; this sweep polices nothing")

    # and the row that says the sweep is blind to `v` here -- executed,
    # so that if a geometry change ever makes `v` bind, this assertion
    # fails and the docstring above stops being true silently
    blind_interior, blind = _containment_sweep(shrunk(fv=0.0))
    assert blind_interior > 100 and blind == 0, (
        f"removing the v inset now leaks ({blind}/{blind_interior}); the "
        "claim that these geometries do not exercise it is stale")


def _containment_sweep(insets_fn):
    """Body of the containment test, shared with its negative control.

    Taking the inset function as an argument is what lets the control
    run the SAME sweep against a weakened guard (R11.2) rather than
    quoting numbers in a docstring that nothing executes. Returns
    `(interior_points, escapes)` summed over the geometries.

    Review R12.3: eligibility is read off the insets THIS call was
    given, not off `class_c_eligible`, which always consults the real
    guard. The first version did the latter, so every weakened-`v` run
    placed its endpoints outside the real guard's window and was thrown
    out with zero eligible pairs -- the control reported "no escapes"
    for a guard it had never actually run. When the insets are the real
    ones the two predicates are asserted to agree, so reading them
    locally does not quietly diverge from the pipeline.
    """

    rng = np.random.default_rng(4242)
    total_interior = total_escapes = 0
    for w in (0.0, 0.6, 1.1):
        # a box the guard actually binds in: the flat inset alone is
        # sqrt(du dv / 2) = 0.247 against a half-extent of 0.45, so an
        # eligible pair sits in a genuinely small interior region
        slab = Slab(du=0.35, dv=0.35, dx=0.9, dy=0.9)
        # via arms(), so the geometry carries its own guard source and
        # Class C is legal on it even at w = 0 (R14.1)
        geom, _ = arms(slab, w)
        insets = insets_fn(geom)
        lim_x = slab.dx / 2.0 - insets.x
        lim_y = slab.dy / 2.0 - insets.y
        assert min(lim_x, lim_y) > 0.02, (
            f"w={w}: the geometry admits nothing to test")
        assert slab.dv - 2.0 * insets.v > 0.02, (
            f"w={w}: the v window is closed, nothing to test")
        real = insets_fn(geom) == guard_insets(geom)

        def admits(point, insets=insets, lim_x=lim_x, lim_y=lim_y,
                   slab=slab):
            return (abs(float(point[2])) <= lim_x
                    and abs(float(point[3])) <= lim_y
                    and insets.v <= float(point[1]) <= slab.dv - insets.v)

        # the FATTEST admissible diamonds, since those are the ones that
        # can reach a face: `p` at the top of the eligible `v` range and
        # `q` at the bottom, spanning almost the whole slab in `u`
        pairs = attempts = 0
        while pairs < 6:
            attempts += 1
            # R11.3: an explicit cap, so an eligibility or predicate
            # regression fails loudly instead of hanging CI
            assert attempts < 20000, (
                f"w={w}: {pairs} eligible causal pairs in {attempts} "
                "attempts -- eligibility or the predicate has regressed")
            p = np.array([rng.uniform(0.0, slab.du * 0.05),
                          slab.dv - insets.v,
                          rng.uniform(-lim_x, lim_x),
                          rng.uniform(-lim_y, lim_y)])
            q = np.array([rng.uniform(slab.du * 0.95, slab.du),
                          insets.v,
                          rng.uniform(-lim_x, lim_x),
                          rng.uniform(-lim_y, lim_y)])
            ok = admits(p) and admits(q)
            if real:
                assert ok == class_c_eligible(geom, p, q), (
                    "the local eligibility read has drifted from the "
                    "pipeline predicate")
            if not ok:
                continue
            if causal_relation(geom, p, q).related is not True:
                continue
            pairs += 1
            # candidates are drawn 30% BEYOND each transverse face and
            # past both `v` faces, so an escape is detectable; drawing
            # from the box itself would be circular, and drawing from a
            # much larger region finds nothing (two earlier attempts
            # returned 0 and 11 points)
            span_x, span_y = slab.dx * 0.65, slab.dy * 0.65
            for _ in range(6000):
                cand = np.array([rng.uniform(float(p[0]), float(q[0])),
                                 rng.uniform(float(q[1]) - 0.15,
                                             float(p[1]) + 0.15),
                                 rng.uniform(-span_x, span_x),
                                 rng.uniform(-span_y, span_y)])
                if causal_relation(geom, p, cand).related is not True:
                    continue
                if causal_relation(geom, cand, q).related is not True:
                    continue
                total_interior += 1
                inside = (0.0 <= cand[1] <= slab.dv
                          and abs(cand[2]) <= slab.dx / 2
                          and abs(cand[3]) <= slab.dy / 2)
                if not inside:
                    total_escapes += 1
    return total_interior, total_escapes


def test_a_weakened_focusing_guard_leaks_on_a_constructed_pair():
    """The case, review R11.2: the negative control existed only as
    numbers in a docstring, which nothing runs, so a weakened guard
    would still have passed the committed test.

    Executed here, and CONSTRUCTED rather than sampled -- random sweeps
    could not catch a shrunken guard at all, because at the small
    `w du` where eligibility is non-empty the focusing factor is barely
    above 1 and the geodesic hardly bulges. The counterexample instead
    uses the fact that makes it a counterexample: **every point of the
    geodesic joining a causally related pair lies in their diamond**, so
    if that geodesic leaves the box, the diamond does.

    The configuration is the maximising corner `(+Y, +Y)` -- same sign,
    where the trajectory bulges outward to `Y / cos(a/2)`. The opposite
    corner is the intuitive guess and reaches only `Y`, which is why an
    earlier version of this control found nothing.

    Review R12.3: it also has to be a control on the PIPELINE, not just
    on the arithmetic. The first version ran at `dy = 6`, where the
    guard admits no pair at all, so `class_c_eligible` would have
    rejected the constructed pair and a leak there could never have
    reached a Class C count. Same `a = 1.2`, wider box, and the pair is
    now put through the eligibility predicate itself.
    """

    w, du = 1.0, 1.2
    slab = Slab(du=du, dv=6.0, dx=16.0, dy=16.0)
    geom = PlaneWaveGeometry(slab, w)
    real = guard_insets(geom)
    half, a = slab.dy / 2.0, w * du
    sign_p, sign_q, alpha = worst_focusing_corner(a)
    assert alpha > 1.05, "the geometry must actually focus"
    assert (real.x < slab.dx / 2.0 and real.y < half
            and 2.0 * real.v < slab.dv), (
        f"the guard admits nothing here ({real}), so a leak could never "
        "reach a Class C count and this would not be a pipeline control")

    # rather than pick a magic shrink factor, find the largest one that
    # leaks. The claim being tested is that SOME shrinkage does -- i.e.
    # the guard is not arbitrarily loose in `y` -- and reporting the
    # threshold says how much slack it actually has.
    def leaks(f):
        lim = half - real.y * f
        return trajectory_excursion(sign_p * lim, sign_q * lim, a) > half

    assert not leaks(1.0), "the real guard must not leak"
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if leaks(mid):
            lo = mid
        else:
            hi = mid
    assert lo > 0.0, (
        "no shrinkage of the y inset makes the constructed geodesic "
        "leave the box, so this control polices nothing")
    lim = half - real.y * lo

    # and the pair really is causally related, so the peak really is a
    # diamond point rather than a curve nobody can traverse. The `v`
    # coordinates are placed inside the guard's own admissible window,
    # so the pair is a pair the pipeline would actually count.
    cost = transverse_cost(du, 0.0, 0.0, sign_p * lim, sign_q * lim, w)
    q_v = real.v
    p = np.array([0.0, q_v + cost + 1e-3, 0.0, sign_p * lim])
    q = np.array([du, q_v, 0.0, sign_q * lim])
    assert float(p[1]) <= slab.dv - real.v, (
        "the constructed pair does not fit the admissible v window")
    assert causal_relation(geom, p, q).related is True

    # THE POINT OF THE CONTROL: under the weakened guard this pair is
    # eligible, so its diamond would be counted while part of it lies
    # outside the box -- and under the real guard it is not.
    assert abs(float(p[3])) <= half - real.y * lo
    assert abs(float(p[3])) > half - real.y, (
        "the constructed pair survives the real guard too, so the "
        "control does not separate the two")
    assert not class_c_eligible(geom, p, q)

    # the real guard keeps the same construction inside
    lim_real = half - real.y
    assert trajectory_excursion(sign_p * lim_real, sign_q * lim_real,
                                a) <= half + 1e-9


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
