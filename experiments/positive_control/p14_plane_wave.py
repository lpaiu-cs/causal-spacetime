"""P14 probe: does pure-Weyl curvature move the causal order, and how far?

Implements the EXPLORATORY probe of `docs/prereg/p14_weyl_curvature.md`
v0.8 (P0-P4), approved 2026-08-01 for exploration only. **No gate, no
threshold, no verdict, and no confirmation seed window appears in this
file, and none may be added to it** -- those belong to a freeze that
follows this probe's numbers and a separate review (design §8.2).

THE CONSTRUCTION (design §4). A vacuum plane wave in Brinkmann form

    ds^2 = 2 du dv + H(u,x,y) du^2 + dx^2 + dy^2,   H = A (x^2 - y^2)

with constant `A = w^2`. Three exact facts carry the design and are
checked by `docs/prereg/p14_checks/p14_brinkmann_vacuum_check.py`:
`sqrt(-g) = 1` for every profile, `Ricci = 0`, `Riemann != 0`. So the
sprinkling measure is Minkowski's EXACTLY -- a uniform-coordinate
Poisson sprinkling is the covariant one -- while the light cones are
curved. The two arms are therefore the SAME POINTS read with `A` and
with `0`, and any difference between them is a functional of the
relation change alone.

CAUSALITY, in closed form (design §4.2). `v` is cyclic so `u` is affine;
with `u` as parameter the transverse equations are linear,

    x'' = +w^2 x   (defocusing),      y'' = -w^2 y   (focusing).

For a causal curve `2 v' + H + x'^2 + y'^2 <= 0`, so `q` is in the
causal future of `p` iff `Du >= 0` and `Dv <= -C`, where `C` is the
minimal transverse action, the classical action of an inverted
oscillator in `x` and an ordinary one in `y`:

    C = S_x + S_y
    S_x = (w / (2 sinh(w s))) [ (x_p^2 + x_q^2) cosh(w s) - 2 x_p x_q ]
    S_y = (w / (2 sin (w s))) [ (y_p^2 + y_q^2) cos (w s) - 2 y_p y_q ]

`S_y` has poles at `w s = k pi`: the CONJUGATE POINTS. The slab is
frozen strictly inside the first one, so `sin(w s)` never vanishes on
admissible pairs and no pair needs a conjugate-point filter -- which is
what lets `D`'s denominator be `C(N,2)` unconditionally.

WHY THE ERROR BOUND IS BUILT THE WAY IT IS (review R6.1). `S_x` and
`S_y` have OPPOSITE SIGNS and cancel: the first version of this module
budgeted its rounding error as `rel_tol * (|Dv| + |C|)`, i.e. against
the magnitude AFTER cancellation. Measured, at O(1) coordinates:
`S_x = +0.5463`, `S_y = -0.5463`, `C = -1.1e-16` with a true error of
`3.6e-17` -- and a reported bound of `1.1e-28`, understating it by
`3.2e11` and confidently deciding pairs whose relation is the opposite.
So the bound below tracks the terms BEFORE they cancel, adds the series
truncation where that branch is used, and escalates to exact
arithmetic rather than guessing when the double result is inside it.

TWO ADMISSIBILITY RULES, NOT ONE (design §4.6). Class R statistics read
the predicate between two points and take EVERY pair, denominator
`C(N,2)`, no containment condition. Class C statistics count sprinkled
points inside a causally defined region and use a guard inset with their
own reported denominator.

AND THE CLAIM BOUNDARY (design §4.6.2). `D` consumes two orders plus the
point correspondence, so it is a PAIRED-ORDER SENSITIVITY DIAGNOSTIC and
not a property of any single causal set.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from decimal import Decimal, getcontext

import numpy as np

# --------------------------------------------------------------------
# exact arithmetic, for the pairs double precision cannot decide
# --------------------------------------------------------------------

#: Digits used when a pair escalates. At this precision the residual
#: error is ~1e-58 relative, i.e. far below anything the double path
#: can leave undecided, so an escalated pair is decided unless it sits
#: on the cone to within that -- which only an exact construction does.
ESCALATED_PRECISION = 60


def _series(x: Decimal, *, odd: bool, alternating: bool) -> Decimal:
    """`sinh`/`cosh`/`sin`/`cos` by Taylor series, summed to convergence.

    `odd` selects the odd-power family (`sin`, `sinh`), `alternating`
    the trigonometric one. Written out because `decimal` has `exp` and
    `ln` but no trigonometry, and the escalation path must not fall back
    on the same double-precision routines it exists to check.
    """

    term = x if odd else Decimal(1)
    total = term
    k = 1
    while True:
        if odd:
            term = term * x * x / (Decimal(2 * k) * Decimal(2 * k + 1))
        else:
            term = term * x * x / (Decimal(2 * k - 1) * Decimal(2 * k))
        contrib = -term if (alternating and k % 2 == 1) else term
        updated = total + contrib
        if updated == total:
            return updated
        total = updated
        k += 1


def _dec(value) -> Decimal:
    """A float's EXACT value as a `Decimal`.

    `Decimal(x)` is exact; `Decimal(repr(x))` is the shortest string
    that round-trips, which is a different number (R9.1). The
    escalation path exists to be exact about the stored coordinates, so
    it must not begin by rounding them.
    """

    return value if isinstance(value, Decimal) else Decimal(float(value))


def _exact_cost(s, xp, xq, yp, yq, w) -> Decimal:
    """`transverse_cost` in `decimal` at `ESCALATED_PRECISION` digits.

    `s` may be a `Decimal`, and in the escalation path it always is: the
    `u`-separation has to be formed as a Decimal DIFFERENCE of the two
    stored coordinates, never handed in as a double that already
    rounded (R9.1).
    """

    getcontext().prec = ESCALATED_PRECISION
    ds = _dec(s)
    dxp, dxq, dyp, dyq, dw = (_dec(xp), _dec(xq), _dec(yp), _dec(yq),
                              _dec(w))
    if dw == 0:
        return ((dxq - dxp) ** 2 + (dyq - dyp) ** 2) / (2 * ds)
    ws = dw * ds
    s_x = (dw / (2 * _series(ws, odd=True, alternating=False))) * (
        (dxp * dxp + dxq * dxq)
        * _series(ws, odd=False, alternating=False) - 2 * dxp * dxq)
    s_y = (dw / (2 * _series(ws, odd=True, alternating=True))) * (
        (dyp * dyp + dyq * dyq)
        * _series(ws, odd=False, alternating=True) - 2 * dyp * dyq)
    return s_x + s_y


# --------------------------------------------------------------------
# Geometry: the transverse cost
# --------------------------------------------------------------------

#: Below this `|w s|` the closed forms are evaluated by their series
#: instead: both legs are `0/0` as `w s -> 0`.
_SMALL_WS = 1e-4

#: Unit roundoff, and a generous operation count for one leg. The
#: budget is deliberately loose -- over-reporting ambiguity costs an
#: escalation, under-reporting it decides a pair wrongly (R6.1).
_EPS = 2.220446049250313e-16
_OPS_PER_LEG = 16.0


def _leg_cost_flat(s: float, dp: float, dq: float) -> float:
    return (dq - dp) ** 2 / (2.0 * s)


def cost_defocusing(s: float, xp: float, xq: float, w: float) -> float:
    """`x` direction, `x'' = +w^2 x`. Regular for every `s > 0`."""

    ws = w * s
    if abs(ws) < _SMALL_WS:
        return _leg_cost_flat(s, xp, xq) + ws * ws * (
            xp * xp + xq * xq + xp * xq) / (6.0 * s)
    return (w / (2.0 * math.sinh(ws))) * (
        (xp * xp + xq * xq) * math.cosh(ws) - 2.0 * xp * xq)


def cost_focusing(s: float, yp: float, yq: float, w: float) -> float:
    """`y` direction, `y'' = -w^2 y`. Singular at `w s = k pi`."""

    ws = w * s
    if abs(ws) < _SMALL_WS:
        return _leg_cost_flat(s, yp, yq) - ws * ws * (
            yp * yp + yq * yq + yp * yq) / (6.0 * s)
    return (w / (2.0 * math.sin(ws))) * (
        (yp * yp + yq * yq) * math.cos(ws) - 2.0 * yp * yq)


def transverse_cost(s: float, xp: float, xq: float,
                    yp: float, yq: float, w: float) -> float:
    """Minimal transverse action between two events at `u`-separation `s`."""

    if w == 0.0:
        return ((xq - xp) ** 2 + (yq - yp) ** 2) / (2.0 * s)
    return cost_defocusing(s, xp, xq, w) + cost_focusing(s, yp, yq, w)


def _leg_error_bound(s: float, ap: float, aq: float, w: float, *,
                     focusing: bool) -> float:
    """Roundoff and truncation for ONE leg, from the UNCANCELLED operands.

    Two separate cancellations happen here and each needs its own scale
    (R6.1, R7.1), which is why neither the leg's own value nor the
    two-leg total may be used:

    - inside a leg, the closed form's numerator
      `P cosh(a) - Q` with `P = ap^2 + aq^2`, `Q = 2 ap aq` cancels when
      the endpoints are close, so the operand scale is `P cosh + |Q|`;
    - between the legs, the `a^2` terms of the two series have OPPOSITE
      signs and cancel, while the `a^4` terms have the SAME sign and
      ADD. Scaling the truncation by anything that shares the `a^2`
      cancellation therefore vanishes exactly where the real error does
      not. The first omitted term is used explicitly instead:

          a^4 ( -P/24 + 7 (P - Q) / 360 ) / (2 s)

      identical in form for both legs, which is why they add. Measured
      at `s = 1`, `w = 9.9e-5`, all endpoints `1`: the series returns
      `0.0`, the true cost is `-8.005e-18`, and this term reproduces it
      to `1e-24`.
    """

    a = w * s
    p_sum = ap * ap + aq * aq
    q_cross = 2.0 * ap * aq

    if abs(a) < _SMALL_WS:
        flat = (aq - ap) ** 2 / (2.0 * s)
        second = a * a * abs(2.0 * p_sum + q_cross) / (12.0 * s)
        roundoff = _OPS_PER_LEG * _EPS * (flat + second)
        omitted = abs(a ** 4 * (-p_sum / 24.0
                                + 7.0 * (p_sum - q_cross) / 360.0)) / (2.0 * s)
        # the series converges fast at these `a`; twice the first
        # omitted term dominates the tail with room to spare
        return roundoff + 2.0 * omitted

    if focusing:
        denominator, trig = abs(math.sin(a)), abs(math.cos(a))
    else:
        denominator, trig = abs(math.sinh(a)), math.cosh(a)
    operands = (w / (2.0 * denominator)) * (p_sum * trig + abs(q_cross))
    roundoff = _OPS_PER_LEG * _EPS * operands

    # R8.1: the ARGUMENT `a = w s` carries its own rounding, and near a
    # conjugate point the nearly-zero `sin(a)` amplifies it. Propagating
    # `da` through the two places `a` enters:
    #
    #   d(1/sin a) -> |cos a| da / sin^2 a   =>  operands * |cot a| * da
    #   d(cos a)   -> |sin a| da             =>  (w/2) P da, the sin cancels
    #
    # so the sensitivity is `1/sin^2`, not `1/sin`. Measured at `w = 0.9`,
    # `s = pi/0.9 - 1e-4`: the true error is `1.80e-9` where the
    # roundoff term alone reports `4.44e-12`, understating it 405-fold
    # and deciding pairs on the wrong side without escalating.
    # `a` is wrong by two roundings of equal size, not one: `s` is
    # itself a subtraction of stored coordinates (R9.1) and `w * s` is a
    # product. Both give `eps |a|`; the trigonometric evaluation adds
    # about `eps`.
    argument_error = _EPS * (2.0 * abs(a) + 1.0)
    sensitivity = operands * (trig / denominator) + 0.5 * w * p_sum
    return roundoff + sensitivity * argument_error


def cost_error_bound(s: float, xp: float, xq: float,
                     yp: float, yq: float, w: float) -> float:
    """A bound on the double-precision error in `transverse_cost`.

    Built per leg from operands that have not yet cancelled, then added
    in absolute value so the legs' own cancellation cannot shrink it.
    `|S_x + S_y|` is never a scale for the error in forming it: R6.1
    measured that model `3.2e11` too small, and R7.1 measured its
    repaired version `2.2e5` too small in the series regime.
    """

    if w == 0.0:
        return _OPS_PER_LEG * _EPS * transverse_cost(s, xp, xq, yp, yq, w)
    return (_leg_error_bound(s, xp, xq, w, focusing=False)
            + _leg_error_bound(s, yp, yq, w, focusing=True))


# --------------------------------------------------------------------
# The geometry object: the conjugate bound is enforced by construction
# --------------------------------------------------------------------

@dataclass(frozen=True)
class Slab:
    """The sprinkled coordinate region: `u` in [0, du], and boxes in
    `v`, `x`, `y`."""

    du: float
    dv: float
    dx: float
    dy: float

    @property
    def coordinate_volume(self) -> float:
        return self.du * self.dv * self.dx * self.dy


@dataclass(frozen=True)
class PlaneWaveGeometry:
    """A slab paired with a curvature, validated on construction.

    This type exists so the conjugate bound is a property of the BOX and
    not a per-pair test (design §4.3): if a geometry cannot be built
    outside the bound, no pair inside one can straddle a conjugate
    point, and `D`'s denominator can be `C(N,2)` unconditionally.

    R6.2: sprinkling and relations take this object rather than a raw
    `(slab, w)`, because validation that a caller may skip is not
    enforcement -- and the check raises rather than asserting, so it
    survives `python -O`.

    The bound enforced here is the EXACT one, `du < pi / w`. How much
    margin to leave under it is an operating-point choice the probe
    makes and reports (see `slab_within_conjugate_bound`), not a
    constant this module gets to pick.
    """

    slab: Slab
    w: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.w) or self.w < 0.0:
            raise ValueError(f"w must be finite and non-negative, got {self.w}")
        if self.slab.du <= 0.0:
            raise ValueError(f"slab u-extent must be positive, got {self.slab.du}")
        if self.w > 0.0 and self.slab.du >= conjugate_du(self.w):
            raise ValueError(
                f"slab u-extent {self.slab.du} reaches the first conjugate "
                f"point {conjugate_du(self.w)} at w={self.w}: every pair in "
                "the slab must be strictly inside it")


def conjugate_du(w: float) -> float:
    """First conjugate `u`-separation, `pi / w`.

    Exact: the focusing Jacobi propagator is `sin(w s) / w` and its
    first zero is at `s = pi / w` (design §4.3, Harte & Drivas 2012).
    """

    return math.inf if w == 0.0 else math.pi / w


def slab_within_conjugate_bound(w: float, fraction: float, *,
                                dv: float, dx: float, dy: float) -> Slab:
    """A slab occupying `fraction` of the first conjugate separation.

    `fraction` is an OPERATING-POINT CHOICE the probe makes and reports,
    which is why it is an argument with no default rather than a module
    constant (R6.2). The design fixes only `du < pi / w`.
    """

    if not 0.0 < fraction < 1.0:
        raise ValueError(f"fraction must be in (0, 1), got {fraction}")
    return Slab(du=fraction * conjugate_du(w), dv=dv, dx=dx, dy=dy)


# --------------------------------------------------------------------
# The causal predicate
# --------------------------------------------------------------------

@dataclass(frozen=True)
class Relation:
    """The causal verdict for one ordered pair in one arm.

    `related` is None only when `[margin_lower, margin_upper]` fails to
    exclude zero even after escalation to exact arithmetic -- design
    §5.1's `ambiguous`. It is NEVER a silent False: an undecided pair
    stays in the denominator and widens `[D_lower, D_upper]`.
    """

    related: bool | None
    margin_lower: float
    margin_upper: float
    escalated: bool

    @property
    def ambiguous(self) -> bool:
        return self.related is None


def _verdict(margin: float, error: float, escalated: bool) -> Relation:
    lower, upper = margin - error, margin + error
    if lower > 0.0:
        return Relation(True, lower, upper, escalated)
    if upper < 0.0:
        return Relation(False, lower, upper, escalated)
    return Relation(None, lower, upper, escalated)


def causal_relation(geometry: PlaneWaveGeometry,
                    p: np.ndarray, q: np.ndarray) -> Relation:
    """Is `q` in the causal future of `p`? Returns a verified verdict.

    Events are `(u, v, x, y)`. The margin is `-(Dv + C)`: positive means
    related. A verdict is returned only when the error interval around
    the margin excludes zero; otherwise the pair is recomputed in exact
    arithmetic, and only a pair still undecided there is `ambiguous`.
    """

    w = geometry.w
    du = float(q[0]) - float(p[0])
    dv = float(q[1]) - float(p[1])
    dx = float(q[2]) - float(p[2])
    dy = float(q[3]) - float(p[3])

    if du < 0.0:
        return Relation(False, -math.inf, -math.inf, False)
    if du == 0.0:
        # R6.3: at zero u-separation `ds^2 = dx^2 + dy^2`, so the pair is
        # spacelike UNLESS it is transversally coincident, and then it
        # lies on the null generator of `∂v` -- future-directed for
        # decreasing `v`. Measure zero in a sprinkling, but the predicate
        # is defined on every pair and fixed anchors can sit here.
        if dx == 0.0 and dy == 0.0 and dv < 0.0:
            return Relation(True, -dv, -dv, False)
        return Relation(False, -math.inf, -math.inf, False)

    if w > 0.0 and du >= conjugate_du(w):
        raise ValueError(
            f"pair straddles a conjugate point: du={du} reaches "
            f"{conjugate_du(w)} at w={w}. PlaneWaveGeometry forbids such "
            "a slab, so this pair did not come from one.")

    xp, xq = float(p[2]), float(q[2])
    yp, yq = float(p[3]), float(q[3])
    cost = transverse_cost(du, xp, xq, yp, yq, w)
    margin = -(dv + cost)
    error = (cost_error_bound(du, xp, xq, yp, yq, w)
             + _OPS_PER_LEG * _EPS * abs(dv))

    verdict = _verdict(margin, error, escalated=False)
    if not verdict.ambiguous:
        return verdict

    # Escalate rather than guess (design §5.1: precision, not tolerance).
    #
    # R9.1: the separations are rebuilt HERE, as Decimal differences of
    # the stored coordinates. Passing the double `du` and `dv` down was
    # the defect: `q[0] - p[0]` has already rounded, and the cost
    # depends on `du` strongly enough for that to outweigh the margin.
    # Measured, flat arm, `p = (1.9556e-07, 199.8037, 0, 0)` and
    # `q = (0.0532716, 0, 4.6139, 0)`: the subtraction loses `3.1e-19`,
    # which the `1/du` in the cost turns into `1.2e-15` -- and the
    # escalated verdict came back `False` with a margin of `-1.7e-15`
    # where the exact one is positive. An escalation that begins by
    # rounding its inputs certifies the wrong answer more confidently
    # than the double path would have.
    getcontext().prec = ESCALATED_PRECISION
    du_exact = _dec(q[0]) - _dec(p[0])
    dv_exact = _dec(q[1]) - _dec(p[1])
    exact_margin = -(dv_exact
                     + _exact_cost(du_exact, xp, xq, yp, yq, w))
    exact_error = Decimal(10) ** (-(ESCALATED_PRECISION - 10))
    scale = abs(exact_margin) + abs(dv_exact) + _dec(abs(cost)) + 1
    return _verdict(float(exact_margin), float(exact_error * scale),
                    escalated=True)


# --------------------------------------------------------------------
# P0: the sprinkling
# --------------------------------------------------------------------

def sprinkle(geometry: PlaneWaveGeometry, rho: float,
             rng: np.random.Generator) -> np.ndarray:
    """Poisson sprinkling of the slab at density `rho`, as `(N, 4)`.

    Uniform in `(u, v, x, y)` with `N ~ Poisson(rho * coordinate
    volume)`. Because `sqrt(-g) = 1` this IS the covariant Poisson
    process of the true 4-volume, with nothing to calibrate -- and the
    SAME array serves both arms (design §4.1), which is why this takes a
    geometry only to be told the slab and to guarantee it was validated.
    """

    slab = geometry.slab
    n = int(rng.poisson(rho * slab.coordinate_volume))
    pts = np.empty((n, 4), dtype=float)
    pts[:, 0] = rng.uniform(0.0, slab.du, n)
    pts[:, 1] = rng.uniform(0.0, slab.dv, n)
    pts[:, 2] = rng.uniform(-slab.dx / 2.0, slab.dx / 2.0, n)
    pts[:, 3] = rng.uniform(-slab.dy / 2.0, slab.dy / 2.0, n)
    return pts


# --------------------------------------------------------------------
# P1: the Class C guard insets (design §4.6), in closed form
# --------------------------------------------------------------------

@dataclass(frozen=True)
class GuardInsets:
    """How far inside the box a Class C pair's endpoints must sit.

    Design §4.6: statistics that COUNT sprinkled points inside a
    causally defined region are biased if the region leaves the box, so
    their eligible pairs are decided by insetting the endpoints rather
    than by inspecting each diamond. Class R statistics use neither of
    these and keep every pair.
    """

    transverse: float
    v: float


def guard_insets(geometry: PlaneWaveGeometry) -> GuardInsets:
    """Both components, from the Jacobi propagators (design §4.6).

    The two costs factor usefully as

        S_x(s; a, b) = (w/2) coth(w s) (b-a)^2 + w tanh(w s/2) a b
        S_y(s; a, b) = (w/2) cot (w s) (b-a)^2 - w tan (w s/2) a b

    **The `v` inset** exists because `C` can be negative, so a diamond
    can rise ABOVE its own past endpoint (R10.1). Only the focusing
    direction does it: `max(-S_x) < 0` for every configuration, while
    maximising `-S_y` over the box `|y_p|, |y| <= Y` gives, at the
    corner `y_p = y = Y`,

        max(-S_y) = w Y^2 tan(w s / 2),

    increasing in `s`, so the bound is that at `s = du`. Note where it
    diverges: at `w du = pi`, the conjugate point -- NOT at `pi/2`,
    which is where the UNCONSTRAINED interior optimum
    `-(w y_p^2/2) tan(w s)` blows up. That optimum needs
    `|y| = |y_p / cos(w s)|`, which leaves the box first, and the
    constrained maximum is the one that governs. The v0.8 text said
    `pi/2`; it was describing an unattainable configuration.

    **The transverse inset** comes from the membership condition
    `C_p + C_q <= Dv`. Bounding the two negative focusing terms by the
    above and the cross terms by `w tanh(w du/2) X^2` each, and using
    `coth(z) >= 1/z` so `S_x >= (b-a)^2 / (2 s)` less that cross term,

        (|x - x_p| + |x - x_q|)^2 <= 2 du K,
        K = Dv + 2 w Y^2 tan(w du/2) + 2 w tanh(w du/2) X^2

    after minimising the left side over the split of `du` between the
    two legs. Hence an excursion of at most `sqrt(du K / 2)` beyond the
    endpoints. It is deliberately conservative -- a guard that is too
    tight silently truncates, a guard that is too loose only costs
    eligible pairs -- and the flat limit is exact: `K = dv` gives
    `sqrt(du dv / 2)`, the Alexandrov diamond's transverse radius.

    **[TO VERIFY at probe time]** what fraction of pairs this leaves
    eligible. The flat limit alone is `sqrt(du dv / 2)`, so the box has
    to be far wider transversally than the diamonds it holds, and §5's
    Class C statistics may end up living on a thin interior population.
    That is a measurement §8 P1 owes, not an assumption.
    """

    slab, w = geometry.slab, geometry.w
    half_x, half_y = slab.dx / 2.0, slab.dy / 2.0
    if w == 0.0:
        v_inset = 0.0
        budget = slab.dv
    else:
        a = w * slab.du
        v_inset = w * half_y * half_y * math.tan(a / 2.0)
        budget = (slab.dv + 2.0 * v_inset
                  + 2.0 * w * math.tanh(a / 2.0) * half_x * half_x)
    return GuardInsets(transverse=math.sqrt(slab.du * budget / 2.0),
                       v=v_inset)


def class_c_eligible(geometry: PlaneWaveGeometry,
                     p: np.ndarray, q: np.ndarray) -> bool:
    """Is this pair eligible for a Class C (counting) statistic?

    Both endpoints must sit inside the inset sub-box. The rule is a
    deterministic function of the frozen geometry and the coordinates,
    never of the realized data (design §4.6) -- and it is the SAME rule
    in both arms, taken from the curved one, so the two arms keep a
    common eligible set (§4.1).

    Class R statistics do not call this: they take every pair.
    """

    insets = guard_insets(geometry)
    slab = geometry.slab
    x_lim = slab.dx / 2.0 - insets.transverse
    y_lim = slab.dy / 2.0 - insets.transverse
    if x_lim <= 0.0 or y_lim <= 0.0:
        return False
    for point in (p, q):
        if abs(float(point[2])) > x_lim or abs(float(point[3])) > y_lim:
            return False
        if not (insets.v <= float(point[1]) <= slab.dv - insets.v):
            return False
    return True


def arms(slab: Slab, w: float) -> tuple[PlaneWaveGeometry, PlaneWaveGeometry]:
    """The curved and flat readings of one slab, both validated.

    They share the slab by construction, which is what makes the two
    arms the same points (design §4.1).
    """

    return PlaneWaveGeometry(slab, w), PlaneWaveGeometry(slab, 0.0)
