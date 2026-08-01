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
