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

import functools
import math
from dataclasses import dataclass, field
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


#: Held by `arms()` alone, so that a linked geometry cannot be
#: hand-built and "came from arms()" is a property a check can
#: read rather than a convention a caller must follow (R14.7).
_ARMS_KEY = object()


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

    `guard_from` carries the arm whose guard governs eligibility, and
    exists because a docstring saying "the same rule in both arms" was
    not a mechanism (R14.1). `class_c_eligible` derived its insets from
    whatever geometry it was handed, so a paired analysis passing each
    arm its own geometry got the flat arm's much smaller insets -- at
    `Slab(1, 0.4, 3, 3)` the `y` inset is `1.207` curved against
    `0.447` flat -- and the two arms would have had different eligible
    sets, which is the one thing §4.1 exists to prevent.

    **The first version of the link checked the wrong things (R14.7.)**
    It validated the shared slab and forbade chaining, neither of which
    constrains `guard_from.w`, so a `w = 0.9` arm could be pointed at a
    `w = 0.2` guard and would apply it -- admitting `|y| <= 1.0419`
    where the real guard allows `0.3429`, and leaking diamonds out of
    the box in 71% of causal pairs. And a hand-built flat-to-flat link
    passed every check and compared equal to `arms()`' own output, so
    "came from `arms()`" was not a property any check could read.

    **A LINKED ARM IS ALWAYS THE FLAT READING (R17.1.)** The second
    version admitted `w in (0, guard_from.w)`, which also lets a curved
    arm link to itself -- a third shape `arms()` never emits, since the
    curved arm it returns is unlinked. Nothing needed that shape and
    `__reduce__` was written assuming it away, so a self-linked curved
    geometry came back from a round trip as the flat arm, silently
    swapping its guard: the R14.1 defect arriving through
    serialisation. Enforcing `w == 0` is what the two `arms()` shapes
    already satisfy -- a flat arm reading its curved partner, and, at
    `w = 0`, a flat arm reading a flat source -- and it makes that
    premise true rather than assumed.

    A link can only be made with the key `arms()` holds, so "came from
    `arms()`" and "may be used for Class C" are one condition.
    """

    slab: Slab
    w: float
    guard_from: PlaneWaveGeometry | None = None
    link_key: object | None = field(
        default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.w) or self.w < 0.0:
            raise ValueError(f"w must be finite and non-negative, got {self.w}")
        if self.guard_from is not None:
            if self.link_key is not _ARMS_KEY:
                raise ValueError(
                    "a linked geometry is built by arms() and nowhere else, "
                    "so that 'came from arms()' is a readable property")
            if self.guard_from.slab != self.slab:
                raise ValueError(
                    "guard_from must share this arm's slab, else the two "
                    "arms are not the same points")
            if self.guard_from.guard_from is not None:
                raise ValueError(
                    "guard_from must be a guard source itself, so that "
                    "one arm's eligibility cannot be read off another's")
            if self.w != 0.0:
                raise ValueError(
                    f"only a flat arm carries a link, and this one has "
                    f"w={self.w}: an arm reads its own guard unless it is "
                    "the flat reading of a curved partner")
        if self.slab.du <= 0.0:
            raise ValueError(f"slab u-extent must be positive, got {self.slab.du}")
        if self.w > 0.0 and self.slab.du >= conjugate_du(self.w):
            raise ValueError(
                f"slab u-extent {self.slab.du} reaches the first conjugate "
                f"point {conjugate_du(self.w)} at w={self.w}: every pair in "
                "the slab must be strictly inside it")

    def __reduce__(self):
        """Rebuild through the checks, so the link survives a round trip.

        R15.14: `pickle` and `deepcopy` restore a dataclass by setting
        `__dict__` directly, which skips `__post_init__`. A linked
        geometry came back with its `guard_from` intact and its key
        gone -- comparing and hashing equal to the original, accepted
        by `class_c_eligible`, but rejected by `dataclasses.replace`.
        The capability was lost silently, which is the one way it must
        not be lost. §8 P3/P4 make independent sprinklings the
        replication unit and `multiprocessing` pickles its arguments,
        so this lands the moment that is parallelised.
        """

        if self.guard_from is None:
            return (PlaneWaveGeometry, (self.slab, self.w))
        return (_arm_of, (self.slab, self.guard_from.w, 1))


def _arm_of(slab: Slab, w: float, index: int) -> PlaneWaveGeometry:
    """Reconstructor for `PlaneWaveGeometry.__reduce__` (R15.14).

    **A selector on `arms()`, not a linker (R16.2.)** The first version
    called the constructor with `_ARMS_KEY` itself, which made it a
    general-purpose linking function that happened to be private by
    naming convention, and rebuilt the flat-to-flat forgery R14.7
    exists to refuse -- admitting a pair both real arms reject. A
    `__reduce__` reconstructor has to reproduce an object that already
    passed the checks, never to hand out the capability to make one.

    The fixed index `1` is safe because `PlaneWaveGeometry` now
    enforces what this assumed: a linked arm has `w == 0`, so it is
    `arms()`' second output and never its first (R17.1). Before that
    was enforced, a self-linked curved arm reduced to `(slab, w, 1)`
    and came back as the flat arm with the flat guard.
    """

    return arms(slab, w)[index]


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
# P1: the Class C guard insets (design §4.6)
#
# NOT in closed form, and §8 P1 requires reporting it as the frozen
# deterministic construction it is. This banner said "in closed form"
# from `9c63ff7`, whose subject made the same claim, and survived the
# `a1f4a2a` audit because that audit's scope was the design document
# (R14.6). Every other "closed form" in this module is about the §4.2
# transverse action, which genuinely has one.
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
    both limits at `dy/2`, is 3.06x and 20.04x too large on the two
    geometries §4.6 and `test_the_asymmetric_v_bound_...` both use.

    **`0 <= anchor_limit <= box_limit` is a PRECONDITION, and checked
    (R14.11).** Zero is admissible and does occur: `guard_insets`
    passes it from the scan grid's first term and from the span
    bisection's lower bracket, so the strict form this docstring
    first stated would have made the module raise on its own first
    call (R15.3). The three branches below cover the rectangle only on
    that side of the diagonal. Beyond it, at `anchor > box / cos a`,
    the maximum in `y_p` is interior at `y_p = box / cos a` and the
    corner form is not merely loose but wrong: at
    `(a, w, anchor, box) = (0.6, 1, 5, 2)` it returns `-3.484` where
    the true maximum is `+1.368` -- and no maximum of `-S_y` over a
    rectangle containing the origin can be negative at all. Nothing
    reaches it today (`v_excursion` passes `anchor <= half_y = box`,
    `direct_bound` passes `anchor == box`, so it would need
    `cos a > 1`), and the error is in the unsafe direction, so the
    precondition is enforced rather than documented: the natural next
    tightening is to inset the anchor harder, which walks toward this
    edge, not away from it.

    Three regimes on that domain, all exact and verified against a
    grid to 3e-15: `cos a <= 0` puts the maximum at the corner;
    otherwise the interior optimum `y = y_p / cos a` is feasible while
    `y_p <= box cos a` and gives `y_p^2 sin^2(a)/cos(a)`; past that it
    is the corner again.
    """

    if not 0.0 <= anchor_limit <= box_limit:
        raise ValueError(
            f"anchor_limit {anchor_limit} must lie in [0, {box_limit}]: "
            "the anchor is an eligible endpoint and so is inset at least "
            "as far as the box, and the branches below assume it")

    cos_a, sin_a = math.cos(a), math.sin(a)
    if cos_a > 0.0 and anchor_limit <= box_limit * cos_a:
        bracket = anchor_limit * anchor_limit * sin_a * sin_a / cos_a
    else:
        bracket = (2.0 * anchor_limit * box_limit
                   - cos_a * (anchor_limit ** 2 + box_limit ** 2))
    return w * bracket / (2.0 * sin_a)


def _defocusing_inset(slab: Slab, w: float,
                      half_x: float) -> tuple[float, float]:
    """`(x_inset, cross)` for the defocusing direction (R13.2).

    The membership budget for the two `x` legs is `Dv + 2V + cross`,
    and for an ELIGIBLE pair `Dv <= dv - 2V`, so `V` cancels and the
    budget is `dv + cross` (R13.1). What remains is the cross term,
    which the first three versions bounded as `2 w tanh(a/2) X^2`.
    Both factors were loose, and the looseness was not harmless: with
    that bound

        x_inset / X >= sqrt(a tanh(a/2)),

    which exceeds 1 at `a = 1.5434` **regardless of `dv` or box size**,
    so every square box closed on `x` alone past that point. I read a
    census that never showed `x` as the WORST condition as showing `x`
    did not bind, and concluded that tightening it would buy nothing.
    It is the ceiling on `a`, and the effect being measured grows as
    `(w T)^4`.

    Two corrections, both already made elsewhere in this guard:

    **The anchor is inset.** The bilinear terms sum to
    `w x_m [tanh(a_p/2) x_p + tanh(a_q/2) x_q]`, where `x_p, x_q` are
    eligible ENDPOINTS -- bounded by `X_inner` -- while only the
    intermediate `x_m` ranges over the box. This is the same asymmetry
    that R12.1 found in `-S_y`, in the term that turned out to set the
    ceiling.

    **The split is not free.** `tanh` is concave, so
    `tanh(a_p/2) + tanh(a_q/2)` over `a_p + a_q = a` is largest at the
    midpoint, giving `2 tanh(a/4)` and not `2 tanh(a/2)`.

    So `cross = 2 w tanh(a/4) X_inner X` with `X_inner = X - x_inset`,
    which puts `x_inset` on both sides. The map is decreasing in
    `x_inset` while the identity is increasing, so the fixed point is
    unique and bisection finds it; the root is approached from the
    conservative side. In the `dv -> 0`, `X = 1` limit it is the root
    of `r^2 = a tanh(a/4) (1 - r)`, giving `r = 0.6046` at `a = 2` and
    `0.7245` at `a = 3` -- so the ceiling is gone. (An earlier draft of
    this note quoted `0.689` and `0.778`, which solve the same equation
    with `tanh(a/2)`: they are the intermediate result from correcting
    only the first of the two factors, and describe no version of this
    code. R14.2.)
    """

    if w == 0.0:
        return math.sqrt(slab.du * slab.dv / 2.0), 0.0

    a = w * slab.du
    tanh_quarter = math.tanh(a / 4.0)

    def cross_at(x_inset: float) -> float:
        return (2.0 * w * tanh_quarter
                * max(half_x - x_inset, 0.0) * half_x)

    def implied(x_inset: float) -> float:
        return math.sqrt(slab.du * (slab.dv + cross_at(x_inset)) / 2.0)

    if implied(half_x) >= half_x:
        # Even a zero-width interior does not satisfy it. Report the
        # inset the condition actually demands -- `sqrt(du dv / 2)`,
        # the flat Alexandrov radius, since `cross` vanishes at zero
        # width -- rather than clamping to `half_x` (R14.8). Clamping
        # made `guard_margins`' `x` component come out at exactly
        # `0.0` for every closed cell, so a geometry one ulp from
        # admitting and one missing by 89x read identically, and the
        # census's `<= 0.0` was the only thing distinguishing them
        # from a cell where `x` was satisfied.
        #
        # `cross` is reported at FULL box width for the same reason.
        # It is vacuous here -- there are no eligible `x` endpoints to
        # cross with -- and the zero it would otherwise take is the
        # most optimistic value in a chain that feeds the `y`
        # diagnostic monotonically.
        return implied(half_x) * _GUARD_SLACK, cross_at(0.0)

    lo, hi = 0.0, half_x
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if implied(mid) > mid:
            lo = mid
        else:
            hi = mid
    x_inset = hi * _GUARD_SLACK  # the side that over-insets
    return x_inset, cross_at(x_inset)


def _blocked(half_x: float, half_y: float,
             x_inset: float, v_inset: float) -> GuardInsets:
    """Insets for a geometry that admits nothing.

    `y = half_y` alone already refuses every pair -- `class_c_eligible`
    returns on `y_lim <= 0` before it looks at a coordinate -- so the
    other two components do not also have to be jammed, and jamming
    them cost information. The first version returned `v = slab.dv`,
    which made the `v` window the rejection that actually fired for
    every blocked geometry and left the extent branch unreachable and
    therefore untested (R15.4). Reporting the `x` inset the condition
    DEMANDS, which may exceed `half_x`, both keeps the extent branch
    live and says by how much the geometry missed.

    R16.6: the `x` component is the RAW demanded inset, not
    `max(x_inset, half_x)`. Clamping reported exactly `half_x` in 559
    of 974 blocked geometries -- every one where `x` was not the
    binding condition -- which is R14.8's saturated zero returning in
    the cells it was not looking at. Nothing rests on the clamp:
    `y_lim = slab.dy/2 - half_y` is a bitwise zero in every blocked
    case, computed from the same `slab.dy` on both sides, so the extent
    branch fires on `y` whatever `x` says.
    """

    return GuardInsets(x=x_inset, y=half_y, v=v_inset)


#: `guard_insets` is a constant of the geometry -- it reads only the
#: frozen `slab` and `w` -- and `class_c_eligible` calls it once per
#: pair, inside what will be an `O(N^2)` loop. One solve is 383 us and
#: dominates the predicate entirely; over a 7140-pair loop caching it
#: is 456x with identical eligible counts (R14.18). Sized to hold every
#: geometry a sweep touches.
#:
#: R16.1: this decorator has to sit against `guard_insets` itself.
#: Inserting `_blocked` between the two put it on a four-float
#: constructor instead -- 3 hits against 9 misses -- and reverted
#: R14.18 entirely, at 285.5 us per call against 0.211 us cached.
@functools.lru_cache(maxsize=4096)
def guard_insets(geometry: PlaneWaveGeometry) -> GuardInsets:
    """All three components, from the Jacobi propagators (§4.6).

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

    **The defocusing inset** comes from the membership condition
    `C_p + C_q <= Dv`, using `coth(z) >= 1/z` so that
    `S_x >= (b-a)^2 / (2 s)` less its cross term, and minimising the
    left side over the split of `du` between the two legs:

        (|x - x_p| + |x - x_q|)^2 <= 2 du K,
        K = dv + cross,   cross = 2 w tanh(a/4) X_inner X

    so the excursion is at most `sqrt(du K / 2)`. See
    `_defocusing_inset`, which solves it -- `X_inner` is on both sides.

    **`K` carries no `max(-S_y)` term**, which is the whole of R13.1:
    an eligible pair sits inside `[V, dv - V]`, so `Dv <= dv - 2V`,
    and the `2V` the focusing legs return cancels it exactly. This
    docstring stated `K = Dv + 2 max(-S_y) + 2 w tanh(w du/2) X^2`
    two rounds after that cancellation and one after `tanh(a/2) X^2`
    became `tanh(a/4) X_inner X` -- on the sweep slab it overstated
    `cross` by 4.5x before the missing term, and a reader reconciling
    it against the design document would have concluded the code was
    wrong (R14.12).

    It is deliberately conservative -- a guard that is too tight
    silently truncates, a guard that is too loose only costs eligible
    pairs -- and the flat limit is exact: `K = dv` gives
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

    # **`V` IS NOT PAID TWICE (R13.1).** The first joint version wrote
    # the budget as `dv + cross + 2V` in both conditions, on the reading
    # that the diamond may use the whole slab in `v` AND rise by `V` at
    # each end. It may not: an ELIGIBLE pair already sits inside the
    # window `[V, dv - V]`, so the drop it has to spend is
    #
    #   Dv = v_p - v_q <= dv - 2 V,
    #
    # and the two focusing legs return at most `2V` of it. In the
    # defocusing condition the two cancel exactly,
    #
    #   (x-x_p)^2/(2 s_p) + (x-x_q)^2/(2 s_q)
    #       <= Dv + 2V + cross <= dv + cross,
    #
    # so **the `x` inset does not depend on `y_inner` at all** and the
    # guard is far less coupled than R12 concluded. Charging `V` to the
    # window and again to the budget shrank the eligible volume and the
    # census with it.
    #
    # The `y` spread does not get the same cancellation, because what
    # it must subtract is the direct action between the two ENDPOINTS,
    # both of which are inset -- an inner-inner bound `W`, not the
    # inner-outer `V` that governs an intermediate point:
    #
    #   Q (y - y*)^2 <= Dv - (S_x,p + S_x,q) - S_y(du; y_p, y_q)
    #                <= (dv - 2V) + cross + W.
    def v_excursion(y_inner: float) -> float:
        """`max(-S_y)` for one leg: anchor inset, intermediate free."""

        return max_negative_focusing(a, w, y_inner, half_y)

    def direct_bound(y_inner: float) -> float:
        """`max(-S_y)` between two ENDPOINTS, so both are inset."""

        return max_negative_focusing(a, w, y_inner, y_inner)

    x_inset, cross = _defocusing_inset(slab, w, half_x)
    x_lim = half_x - x_inset

    def room(y_inner: float) -> tuple[float, float, float]:
        """`(y_room, v_window, v_inset)` at this anchor limit.

        Every margin carries `_GUARD_SLACK` of relative room, so the
        insets built from it can only err LARGE -- costing eligible
        pairs, never admitting a pair whose diamond leaves the box.
        Margins may be non-positive; the caller decides.
        """

        v_here = v_excursion(y_inner) * _GUARD_SLACK
        v_window = slab.dv - 2.0 * v_here
        budget = max(v_window, 0.0) + cross + direct_bound(y_inner)
        reach = alpha * y_inner + math.sqrt(budget * tan_half / w)
        return half_y - reach * _GUARD_SLACK, v_window, v_here

    def volume(y_inner: float) -> float:
        y_room, v_window, _ = room(y_inner)
        if x_lim <= 0.0 or y_room <= 0.0 or v_window <= 0.0:
            return 0.0
        return 2.0 * x_lim * 2.0 * y_inner * v_window * slab.du

    # **The scan is bracketed by the `v` window, not by the box.** The
    # first version scanned `[0, half_y]` uniformly, and the feasible
    # set can be a sliver far narrower than one step of that: at
    # `Slab(1.59, 0.2, 6, 6)`, `w = 1`, it is `[0, 0.00452]` against a
    # step of `0.0117`, so every sample missed it and the geometry was
    # reported as admitting nothing when it admits a thin region. The
    # `v` window closes monotonically in `y_inner` -- `max(-S_y)` grows
    # with its domain -- so its upper end is found exactly by bisection
    # and the scan gets the interval that can contain the optimum.
    #
    # R13.3: within that bracket the objective is still not known to be
    # unimodal, so every interior local maximum is refined and the best
    # kept. A peak is missed only if it is narrower than one step of
    # the BRACKET, which is the guarantee claimed -- not "the maximum".
    span = half_y
    if slab.dv - 2.0 * v_excursion(span) * _GUARD_SLACK <= 0.0:
        lo_v, hi_v = 0.0, span
        for _ in range(80):
            mid = 0.5 * (lo_v + hi_v)
            if slab.dv - 2.0 * v_excursion(mid) * _GUARD_SLACK > 0.0:
                lo_v = mid
            else:
                hi_v = mid
        span = lo_v
    if span <= 0.0:
        return _blocked(half_x, half_y, x_inset,
                        v_excursion(0.0) * _GUARD_SLACK)

    grid = [span * i / _GUARD_SCAN for i in range(_GUARD_SCAN + 1)]
    vols = [volume(y) for y in grid]
    if max(vols) <= 0.0:
        return _blocked(half_x, half_y, x_inset,
                        v_excursion(0.0) * _GUARD_SLACK)

    # The global argmax is always an interior peak already -- `vols[0]`
    # is identically zero and the far end is either zero or a boundary
    # maximum -- so it is added only when the endpoints could carry it,
    # rather than unconditionally re-refining a bracket already in the
    # list (R14.9).
    peaks = {i for i in range(1, _GUARD_SCAN)
             if vols[i] >= vols[i - 1] and vols[i] >= vols[i + 1]
             and vols[i] > 0.0}
    peaks.add(max(range(len(vols)), key=lambda i: vols[i]))
    y_inner = 0.0
    best_vol = 0.0
    for peak in sorted(peaks):
        lo = grid[max(peak - 1, 0)]
        hi = grid[min(peak + 1, _GUARD_SCAN)]
        for _ in range(80):  # ternary search, bracket below an ulp
            m1, m2 = lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0
            if volume(m1) < volume(m2):
                lo = m1
            else:
                hi = m2
        for candidate in (0.5 * (lo + hi), grid[peak]):
            here = volume(candidate)
            if here > best_vol:
                y_inner, best_vol = candidate, here

    _, _, v_inset = room(y_inner)
    return GuardInsets(x=x_inset, y=half_y - y_inner, v=v_inset)


def guard_margins(geometry: PlaneWaveGeometry) -> tuple[float, float, float]:
    """The three interior margins, as fractions of the box's own scale.

    `guard_insets` jams `y` once the eligible set is empty -- that
    alone refuses every pair -- which is right for the predicate and
    useless for asking WHICH condition emptied it. (It used to jam all
    three; R15.4 stopped, because jamming them cost exactly the
    information this function exists to recover. The sentence here
    still described the old sentinel, which is R14.6's shape again:
    the site that changed was updated and the sibling claim elsewhere
    was not.) This reports separately: the
    margins at the guard's operating point, or -- when the set is
    empty -- at the anchor limit that comes closest to feasible.

    R12.2: they are reported as fractions so that the answer depends on
    the four dimensionless ratios `w du`, `w dv`, `w dx/2`, `w dy/2`
    and not on the scale, which is what makes a sweep over them a
    complete statement rather than one slice of it.

    **Read all three (R13.2).** A census that records only the WORST
    condition per empty cell answers a different question from the one
    it looks like it answers: I read `x` never being worst as `x` never
    binding, and concluded that tightening it would buy nothing. It
    binds in most of the empty cells. Counting has to be per condition,
    not per cell, and every caller of this gets all three margins for
    that reason.

    **All three mean CONTAINMENT SLACK, in both branches (R14.10.)**
    The first version reported the `y` component as the eligible width
    `y_inner / half_y` when the set was non-empty and as the slack
    `(half_y - reach) / half_y` when it was not, so the number jumped
    by 205x across the ulp of `a` where a cell closes and the two
    branches could not be compared. Slack is the one of the two that
    answers "which condition emptied this cell", so it is what both
    branches report, and `_GUARD_SLACK` is applied here exactly as the
    predicate applies it.

    **And the empty branch no longer reports from a zero-volume
    anchor.** Every margin is non-increasing in `y_inner`, so
    `max(key=min)` always landed on `y_inner = 0` -- where the eligible
    width is zero by construction and all three margins can be
    positive at once, giving an empty cell with no blocker. The scan
    runs over anchors that would carry positive volume if the margins
    allowed it, so the reported margins are the ones that stopped it.
    """

    # R15.1: redirect the whole GEOMETRY, once, exactly as
    # `class_c_eligible` does. The first version unpacked `slab, w`
    # from the source and then called `guard_insets` on the original,
    # so the anchor handed to the curved closure was the flat arm's --
    # `1.0528` against `0.3429` at `Slab(1, 0.4, 3, 3)`, w = 0.9 -- and
    # in 146 of 2400 sampled geometries the two arms took different
    # branches and named different blockers for one shared guard. That
    # is R13.2's mis-attribution arriving through the diagnostic.
    if geometry.guard_from is not None:
        geometry = geometry.guard_from
    slab, w = geometry.slab, geometry.w
    half_x, half_y = slab.dx / 2.0, slab.dy / 2.0

    if w == 0.0:
        # R15.2: no `_GUARD_SLACK` here. `guard_insets`' flat branch
        # applies none, and a diagnostic that is pessimistic relative
        # to the predicate disagrees with it inside a 1e-12 band.
        flat = math.sqrt(slab.du * slab.dv / 2.0)
        return ((half_x - flat) / half_x, (half_y - flat) / half_y, 1.0)

    a = w * slab.du
    alpha = _focusing_excursion_factor(a)
    tan_half = math.tan(a / 2.0)
    x_inset, cross = _defocusing_inset(slab, w, half_x)
    x_margin = (half_x - x_inset) / half_x

    def margins(y_inner: float) -> tuple[float, float, float]:
        v_here = max_negative_focusing(a, w, y_inner, half_y) * _GUARD_SLACK
        v_window = slab.dv - 2.0 * v_here
        budget = (max(v_window, 0.0) + cross
                  + max_negative_focusing(a, w, y_inner, y_inner))
        reach = alpha * y_inner + math.sqrt(budget * tan_half / w)
        return (x_margin,
                (half_y - reach * _GUARD_SLACK) / half_y,
                v_window / slab.dv)

    insets = guard_insets(geometry)
    if insets.y < half_y:  # non-empty: read the operating point off it
        return margins(half_y - insets.y)

    # Empty. The anchor that comes closest to feasible, over anchors
    # that carry width -- `y_inner = 0` carries none, so its margins
    # describe a configuration that admits nothing whatever they say.
    grid = [half_y * i / _GUARD_SCAN for i in range(1, _GUARD_SCAN + 1)]
    return max((margins(y) for y in grid), key=min)


def class_c_eligible(geometry: PlaneWaveGeometry,
                     p: np.ndarray, q: np.ndarray) -> bool:
    """Is this pair eligible for a Class C (counting) statistic?

    Both endpoints must sit inside the inset sub-box. The rule is a
    deterministic function of the frozen geometry and the coordinates,
    never of the realized data (design §4.6) -- and it is the SAME rule
    in both arms, taken from the curved one, so the two arms keep a
    common eligible set (§4.1).

    **That last sentence used to be only a sentence (R14.1.)** This
    derived its insets from whatever geometry it was handed, so a
    paired analysis passing each arm its own geometry got the flat
    arm's insets in the flat arm -- much smaller, since flatness
    reduces the guard to the Alexandrov radius -- and the two arms
    ended up with different eligible sets and different denominators.
    A flat geometry is refused here unless `arms()` built it and linked
    it to its curved partner, which is the only construction where the
    shared set is guaranteed.

    Class R statistics do not call this: they take every pair.
    """

    if geometry.guard_from is not None:
        geometry = geometry.guard_from
    elif geometry.w == 0.0:
        raise ValueError(
            "Class C eligibility cannot be derived from a flat geometry: "
            "the guard must come from the curved arm so both arms share "
            "one eligible set (design §4.1). Build the pair with arms().")

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
    arms the same points (design §4.1). And the flat arm is linked to
    the curved one, so Class C eligibility is the SAME set in both --
    a property of the objects rather than of the caller's discipline
    (R14.1).

    At `w = 0` both arms are the same geometry and nothing can
    disagree, but they are still linked, so that "came from `arms()`"
    and "may be used for Class C" are the same condition and there is
    no second rule to remember.
    """

    source = PlaneWaveGeometry(slab, w)
    flat = PlaneWaveGeometry(slab, 0.0, guard_from=source,
                             link_key=_ARMS_KEY)
    if w == 0.0:
        return flat, PlaneWaveGeometry(slab, 0.0, guard_from=source,
                                       link_key=_ARMS_KEY)
    return source, flat
