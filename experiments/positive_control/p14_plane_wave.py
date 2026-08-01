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

#: Relative room demanded of the guard's feasibility conditions, so
#: that the bisected inset errs in a FIXED direction -- too large,
#: costing eligible pairs -- rather than merely to within roundoff of
#: the true root (R12.4). Thousands of ulps, still far below the
#: sprinkling scale.
_GUARD_SLACK = 1.0 + 1e-12

#: Scan resolution for the guard's free parameter. The eligible volume
#: is not guaranteed unimodal in it, so the bracket is found by scan
#: before it is refined.
_GUARD_SCAN = 256


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

    **`x` and `y` need separate insets (R11.1).** The first version
    derived one transverse bound from `coth(z) >= 1/z` and applied it to
    both, but that inequality is a fact about the DEFOCUSING direction.
    The focusing one has `cot(z)` in the same slot, where it is false
    and in fact negative past `z = pi/2`. A random containment sweep
    passing is not a proof, and the guard is about to be used as the
    ground-truth label for the order-invariant eligibility candidate, so
    a certificate is what it needs.
    """

    x: float
    y: float
    v: float


def trajectory_excursion(yp: float, yq: float, a: float) -> float:
    """`max |y(u)|` on the focusing geodesic between two endpoints.

    With `u` as parameter the focusing trajectory is
    `y(u) = y_p cos(w u) + C sin(w u)`, `C = (y_q - y_p cos a)/sin a`,
    i.e. `R sin(w u + phi)` -- so the maximum over the segment is either
    at an endpoint or at an interior peak, and both are exact. Nothing
    is sampled.
    """

    if a <= 0.0:
        return max(abs(yp), abs(yq))
    c = (yq - yp * math.cos(a)) / math.sin(a)
    r = math.hypot(yp, c)
    if r == 0.0:
        return 0.0
    phi = math.atan2(yp, c)
    best = max(abs(math.sin(phi)), abs(math.sin(phi + a)))
    first = math.ceil((phi - math.pi / 2.0) / math.pi)
    for k in (first, first + 1):
        peak = math.pi / 2.0 + k * math.pi
        if phi <= peak <= phi + a:
            best = 1.0
            break
    return r * best


def worst_focusing_corner(a: float) -> tuple[float, float, float]:
    """The unit-box corner with the largest focusing excursion, and it.

    `y(u)` is affine in `(y_p, y_q)`, `|.|` is convex and a maximum of
    convex functions is convex, so `trajectory_excursion` is convex on
    the box and its maximum is attained at a CORNER. Four evaluations,
    no grid, and homogeneous of degree one so the result scales.

    The maximising corner is the SAME-SIGN one, `(+Y, +Y)`, where the
    trajectory bulges outward and reaches `Y / cos(a/2)` -- not the
    opposite-sign corner, which is the intuitive guess and is wrong.
    Returned rather than just its value because the constructed
    counterexample in the tests needs the configuration, not the number.
    """

    best = (1.0, 1.0, trajectory_excursion(1.0, 1.0, a))
    for sy in (-1.0, 1.0):
        for sq in (-1.0, 1.0):
            value = trajectory_excursion(sy, sq, a)
            if value > best[2]:
                best = (sy, sq, value)
    return best


def _focusing_excursion_factor(a: float) -> float:
    """`max |y(u)| / Y` over endpoints in `[-Y, Y]^2`, exactly."""

    return worst_focusing_corner(a)[2]


def max_negative_focusing(a: float, w: float,
                          anchor_limit: float, box_limit: float) -> float:
    """`max(-S_y)` over `|y_p| <= anchor_limit`, `|y| <= box_limit`.

    Maximising `2 y_p y - cos(a)(y_p^2 + y^2)` over the RECTANGLE, which
    is what the geometry actually offers: the anchor is an eligible
    endpoint and is already inset, while the intermediate point ranges
    over the whole box. R12 asked for this; the symmetric version, with
    both limits at `dy/2`, is between 1.8x and 12x too large.

    Three regimes, all exact and verified against a grid to 3e-15:
    `cos a <= 0` puts the maximum at the corner; otherwise the interior
    optimum `y = y_p / cos a` is feasible while `y_p <= box cos a` and
    gives `y_p^2 sin^2(a)/cos(a)`; past that it is the corner again.
    """

    cos_a, sin_a = math.cos(a), math.sin(a)
    if cos_a > 0.0 and anchor_limit <= box_limit * cos_a:
        bracket = anchor_limit * anchor_limit * sin_a * sin_a / cos_a
    else:
        bracket = (2.0 * anchor_limit * box_limit
                   - cos_a * (anchor_limit ** 2 + box_limit ** 2))
    return w * bracket / (2.0 * sin_a)


def guard_insets(geometry: PlaneWaveGeometry) -> GuardInsets:
    """Both components, from the Jacobi propagators (design §4.6).

    The two costs factor usefully as

        S_x(s; a, b) = (w/2) coth(w s) (b-a)^2 + w tanh(w s/2) a b
        S_y(s; a, b) = (w/2) cot (w s) (b-a)^2 - w tan (w s/2) a b

    **The `v` inset** exists because `C` can be negative, so a diamond
    can rise ABOVE its own past endpoint (R10.1). Only the focusing
    direction does it: `max(-S_x) < 0` for every configuration, while
    `max(-S_y)` over the RECTANGLE `|y_p| <= Y_inner`, `|y| <= Y` is
    `max_negative_focusing`. The rectangle, not the square, is what the
    geometry offers: the anchor is an eligible endpoint and so already
    inset, while the intermediate point ranges over the whole box
    (R12). Its symmetric special case is the corner value
    `w Y^2 tan(w s / 2)`, increasing in `s`, so the bound is that at
    `s = du`. Note where it diverges: at `w du = pi`, the conjugate
    point -- NOT at `pi/2`, which is where the UNCONSTRAINED interior
    optimum `-(w y_p^2/2) tan(w s)` blows up. That optimum needs
    `|y| = |y_p / cos(w s)|`, which leaves the box first, and the
    constrained maximum is the one that governs. The v0.8 text said
    `pi/2`; it was describing an unattainable configuration.

    **The transverse inset** comes from the membership condition
    `C_p + C_q <= Dv`. Bounding the two negative focusing terms by the
    above and the cross terms by `w tanh(w du/2) X^2` each, and using
    `coth(z) >= 1/z` so `S_x >= (b-a)^2 / (2 s)` less that cross term,

        (|x - x_p| + |x - x_q|)^2 <= 2 du K,
        K = Dv + 2 max(-S_y) + 2 w tanh(w du/2) X^2

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
        flat = math.sqrt(slab.du * slab.dv / 2.0)
        return GuardInsets(x=flat, y=flat, v=0.0)

    a = w * slab.du

    # The focusing direction (R11.1). The reachable `y` at split
    # `a_p + a_q = a` is bounded by the trajectory plus a spread:
    #
    #   f(y) >= S_y(du; y_p, y_q) + Q (y - y*)^2,
    #   Q = (w/2) sin(a) / (sin a_p sin a_q) >= w cot(a/2) > 0
    #
    # the identity holding even where the individual cotangents are
    # negative, and `min_y f = S_y(du; y_p, y_q)` by composition of the
    # classical action. So the spread is at most
    # `sqrt(D tan(a/2) / w)` with `D` the budget above the direct
    # action, and the trajectory term is `alpha(a) * Y_inner` with
    # `alpha` the corner maximum.
    #
    # `Y_inner = half_y - y_inset` appears on both sides, so the inset
    # is NOT a closed form (R12.4). Nor is it the root of anything:
    # R12.1 asked for the three components to be solved TOGETHER, and
    # once they are, the guard turns out to have a free parameter.
    alpha = _focusing_excursion_factor(a)
    tan_half = math.tan(a / 2.0)
    cross = 2.0 * w * math.tanh(a / 2.0) * half_x * half_x

    # The coupling runs one way. Admitting anchors further out in `y`
    # raises the asymmetric `max(-S_y)` -- the anchor is an eligible
    # endpoint and so already inset, while the intermediate point
    # ranges over the whole box -- and that excursion appears twice:
    # once as the `v` window `dv - 2 v` the endpoints must fit in, and
    # again inside the budget that sets the defocusing inset. So `y`
    # room is BOUGHT with `x` room and `v` room, and there is no
    # distinguished point on that trade-off curve.
    #
    # Maximising `y_inner` -- the first thing I tried -- spends the
    # whole `v` window and leaves a region of zero volume, which
    # satisfies every inequality and admits no pairs at all. The
    # criterion has to be the quantity the statistic actually consumes,
    # so it is the eligible coordinate volume itself, which is also
    # what §8 P1 owes as a measurement.
    def v_excursion(y_inner: float) -> float:
        return max_negative_focusing(a, w, y_inner, half_y)

    def room(y_inner: float) -> tuple[float, float, float, float]:
        """`(x_lim, y_room, v_window, v_inset)` at this anchor limit.

        Every margin carries `_GUARD_SLACK` of relative room, so the
        insets built from it can only err LARGE -- costing eligible
        pairs, never admitting a pair whose diamond leaves the box.
        Margins may be non-positive; the caller decides.
        """

        v_here = v_excursion(y_inner) * _GUARD_SLACK
        budget = slab.dv + cross + 2.0 * v_here
        reach = alpha * y_inner + math.sqrt(budget * tan_half / w)
        return (half_x - math.sqrt(slab.du * budget / 2.0) * _GUARD_SLACK,
                half_y - reach * _GUARD_SLACK,
                slab.dv - 2.0 * v_here,
                v_here)

    def volume(y_inner: float) -> float:
        x_lim, y_room, v_window, _ = room(y_inner)
        if x_lim <= 0.0 or y_room <= 0.0 or v_window <= 0.0:
            return 0.0
        return 2.0 * x_lim * 2.0 * y_inner * v_window * slab.du

    # Coarse scan then local refinement: the objective vanishes at both
    # ends -- no `y` room at zero, no `v` window at the far end -- and
    # is not guaranteed unimodal, so a bracket is found by scan rather
    # than assumed.
    grid = [half_y * i / _GUARD_SCAN for i in range(_GUARD_SCAN + 1)]
    best = max(range(len(grid)), key=lambda i: volume(grid[i]))
    if volume(grid[best]) <= 0.0:
        # Nothing is eligible at any anchor limit. Block on all three,
        # not just the one that happens to bind: reporting a positive
        # `y` limit here would let pairs through the very condition
        # that emptied the set.
        return GuardInsets(x=half_x, y=half_y, v=slab.dv)

    lo = grid[max(best - 1, 0)]
    hi = grid[min(best + 1, _GUARD_SCAN)]
    for _ in range(80):  # ternary search, bracket below an ulp
        m1, m2 = lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0
        if volume(m1) < volume(m2):
            lo = m1
        else:
            hi = m2
    y_inner = 0.5 * (lo + hi)
    if volume(y_inner) <= 0.0:  # refinement wandered off the feasible set
        y_inner = grid[best]
    x_lim, _, _, v_inset = room(y_inner)
    return GuardInsets(x=half_x - x_lim, y=half_y - y_inner, v=v_inset)


def guard_margins(geometry: PlaneWaveGeometry) -> tuple[float, float, float]:
    """The three interior margins, as fractions of the box's own scale.

    `guard_insets` blocks on all three components once the eligible set
    is empty, which is right for the predicate and useless for asking
    WHICH condition emptied it. This reports that separately: the
    margins at the guard's operating point, or -- when the set is
    empty -- at the anchor limit that comes closest to feasible.

    R12.2: they are reported as fractions so that the answer depends on
    the four dimensionless ratios `w du`, `w dv`, `w dx/2`, `w dy/2`
    and not on the scale, which is what makes a sweep over them a
    complete statement rather than one slice of it.
    """

    slab, w = geometry.slab, geometry.w
    half_x, half_y = slab.dx / 2.0, slab.dy / 2.0
    insets = guard_insets(geometry)
    if insets.x < half_x:  # non-empty: read the operating point off it
        return ((half_x - insets.x) / half_x,
                (half_y - insets.y) / half_y,
                (slab.dv - 2.0 * insets.v) / slab.dv)

    if w == 0.0:
        flat = math.sqrt(slab.du * slab.dv / 2.0)
        return ((half_x - flat) / half_x, (half_y - flat) / half_y, 1.0)

    a = w * slab.du
    alpha = _focusing_excursion_factor(a)
    tan_half = math.tan(a / 2.0)
    cross = 2.0 * w * math.tanh(a / 2.0) * half_x * half_x

    def margins(y_inner: float) -> tuple[float, float, float]:
        v_here = max_negative_focusing(a, w, y_inner, half_y)
        budget = slab.dv + cross + 2.0 * v_here
        reach = alpha * y_inner + math.sqrt(budget * tan_half / w)
        return ((half_x - math.sqrt(slab.du * budget / 2.0)) / half_x,
                (half_y - reach) / half_y,
                (slab.dv - 2.0 * v_here) / slab.dv)

    grid = [half_y * i / _GUARD_SCAN for i in range(_GUARD_SCAN + 1)]
    return max((margins(y) for y in grid), key=min)


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
    x_lim = slab.dx / 2.0 - insets.x
    y_lim = slab.dy / 2.0 - insets.y
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
