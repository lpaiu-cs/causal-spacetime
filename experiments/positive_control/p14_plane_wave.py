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

What that does NOT give (design §4.4, review R1.2): the volume of a
causal DIAMOND is not flat, because its boundary moves with the cones.
At equal proper time `V_A/V_0 = 1 + (wT)^4/252 + O((wT)^8)`. Interval
cardinality is therefore a signal, not a control.

CAUSALITY, in closed form (design §4.2). `v` is cyclic so `u` is affine;
with `u` as parameter the transverse equations are linear,

    x'' = +w^2 x   (defocusing),      y'' = -w^2 y   (focusing).

For a causal curve `2 v' + H + x'^2 + y'^2 <= 0`, so `q` is in the
causal future of `p` iff `Du > 0` and

    Dv <= -C(Du; x_p,x_q, y_p,y_q),

where `C` is the minimal transverse action, attained on the geodesic and
available as the classical action of a (inverted) harmonic oscillator:

    C = S_x + S_y
    S_x = (w / (2 sinh(w s))) [ (x_p^2 + x_q^2) cosh(w s) - 2 x_p x_q ]
    S_y = (w / (2 sin (w s))) [ (y_p^2 + y_q^2) cos (w s) - 2 y_p y_q ]

with `s = Du`. The flat arm is the `w -> 0` limit, `C = (Dx^2+Dy^2)/(2s)`,
used directly rather than by taking a numerical limit.

`S_y` has poles at `w s = k pi`: those are the CONJUGATE POINTS, and the
slab is frozen strictly inside the first one (design §4.3), so `sin(w s)`
never vanishes on admissible pairs and no pair needs a conjugate-point
filter -- which is what lets `D`'s denominator be `C(N,2)` unconditionally.

TWO ADMISSIBILITY RULES, NOT ONE (design §4.6, review R3.1). Class R
statistics read the predicate between two points and take EVERY pair,
denominator `C(N,2)`, no containment condition. Class C statistics count
sprinkled points inside a causally defined region and use a guard inset
with their own reported denominator. Collapsing them was a defect.

AND THE CLAIM BOUNDARY (design §4.6.2, review R4.1). `D` consumes two
orders plus the point correspondence, so it is a PAIRED-ORDER
SENSITIVITY DIAGNOSTIC and not a property of any single causal set.
Nothing computed here licenses a single-poset claim except the global
relation fraction over all elements.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# --------------------------------------------------------------------
# Geometry: the transverse cost, and the causal predicate
# --------------------------------------------------------------------

#: Below this `|w s|` the closed forms are evaluated by their series
#: instead. Both `S_x` and `S_y` are `0/0` at `w s -> 0` -- the flat
#: limit is finite but the expression is not -- and the cancellation
#: costs digits long before it costs the answer. Chosen by measuring
#: both forms against the flat limit (see tests).
_SMALL_WS = 1e-4


def _leg_cost_flat(s: float, dp: float, dq: float) -> float:
    """One transverse direction's cost in the flat arm.

    `(dq - dp)^2 / (2 s)`, the straight line. Written as its own
    function because the flat arm is used directly and not as a
    numerical `w -> 0` limit of the curved one.
    """

    return (dq - dp) ** 2 / (2.0 * s)


def cost_defocusing(s: float, xp: float, xq: float, w: float) -> float:
    """Transverse cost in the `x` (defocusing) direction, `x'' = +w^2 x`.

    The classical action of the inverted oscillator. Regular for every
    `s > 0`: `sinh` has no zero, so this direction contributes no
    conjugate points.
    """

    ws = w * s
    if abs(ws) < _SMALL_WS:
        # sinh/cosh series to the order the flat limit needs
        return _leg_cost_flat(s, xp, xq) + ws * ws * (
            xp * xp + xq * xq + xp * xq) / (6.0 * s)
    return (w / (2.0 * math.sinh(ws))) * (
        (xp * xp + xq * xq) * math.cosh(ws) - 2.0 * xp * xq)


def cost_focusing(s: float, yp: float, yq: float, w: float) -> float:
    """Transverse cost in the `y` (focusing) direction, `y'' = -w^2 y`.

    The classical action of the ordinary oscillator. **Singular at
    `w s = k pi`** -- the conjugate points. The caller is responsible
    for staying inside the first one; `causal_relation` asserts it
    rather than trusting it.
    """

    ws = w * s
    if abs(ws) < _SMALL_WS:
        return _leg_cost_flat(s, yp, yq) - ws * ws * (
            yp * yp + yq * yq + yp * yq) / (6.0 * s)
    return (w / (2.0 * math.sin(ws))) * (
        (yp * yp + yq * yq) * math.cos(ws) - 2.0 * yp * yq)


def transverse_cost(s: float, xp: float, xq: float,
                    yp: float, yq: float, w: float) -> float:
    """Minimal transverse action between two events at `u`-separation `s`.

    `q` is causally after `p` iff `Dv <= -transverse_cost(...)`. In the
    flat arm (`w = 0`) this is `(Dx^2 + Dy^2) / (2 s)`.
    """

    if w == 0.0:
        return ((xq - xp) ** 2 + (yq - yp) ** 2) / (2.0 * s)
    return (cost_defocusing(s, xp, xq, w)
            + cost_focusing(s, yp, yq, w))


#: Fraction of the first conjugate `u`-separation the slab may occupy.
#: The bound is EXACT -- the focusing Jacobi propagator is `sin(w s)/w`,
#: first zero at `s = pi / w` (design §4.3, Harte & Drivas 2012) -- and
#: this is the safety margin under it, not the bound itself. It is a
#: property of the BOX: every pair in a slab this thin satisfies it by
#: construction, so no pair is ever excluded on conjugate grounds.
CONJUGATE_SAFETY = 0.8


def max_slab_du(w: float) -> float:
    """Largest `u`-extent a slab may have at curvature `w`."""

    if w == 0.0:
        return math.inf
    return CONJUGATE_SAFETY * math.pi / w


@dataclass(frozen=True)
class Relation:
    """The causal verdict for one ordered pair in one arm.

    `related` is None when the sign of `margin` is not resolved at the
    achieved precision -- design §5.1's `ambiguous`. It is NEVER a
    silent False: an undecided pair stays in the denominator and widens
    `[D_lower, D_upper]` instead of quietly reading as "unrelated".
    """

    related: bool | None
    margin: float
    error: float

    @property
    def ambiguous(self) -> bool:
        return self.related is None


def causal_relation(p: np.ndarray, q: np.ndarray, w: float,
                    rel_tol: float = 1e-12) -> Relation:
    """Is `q` in the causal future of `p`? Returns a verified verdict.

    Events are `(u, v, x, y)`. The margin is `-(Dv + C)`: positive means
    related, negative means not, and the verdict is withheld when the
    margin is inside the accumulated floating-point error rather than
    being guessed (design §5.1).

    The error estimate is deliberately crude and CONSERVATIVE -- it
    scales the terms that actually cancel, so it over-reports rather
    than under-reports ambiguity. A too-large ambiguous fraction is
    visible in the reported bracket and is answered with precision, not
    with a wider tolerance.
    """

    du = float(q[0] - p[0])
    if du <= 0.0:
        return Relation(False, -math.inf, 0.0)

    if w != 0.0:
        assert du < math.pi / w, (
            f"pair straddles a conjugate point: w*du = {w * du:.6f} "
            f">= pi. The slab must be frozen thinner (see max_slab_du).")

    cost = transverse_cost(du, float(p[2]), float(q[2]),
                           float(p[3]), float(q[3]), w)
    dv = float(q[1] - p[1])
    margin = -(dv + cost)

    scale = abs(dv) + abs(cost)
    error = rel_tol * max(scale, 1e-300)
    if abs(margin) <= error:
        return Relation(None, margin, error)
    return Relation(margin > 0.0, margin, error)


# --------------------------------------------------------------------
# P0: the sprinkling. Uniform in coordinates IS the covariant Poisson
# process, exactly, because sqrt(-g) = 1 for every profile.
# --------------------------------------------------------------------

@dataclass(frozen=True)
class Slab:
    """The sprinkled region: `u` in [0, du], `v`, `x`, `y` in boxes.

    `du` must respect `max_slab_du(w)`; `Slab.validate` asserts it so a
    geometry that violates the conjugate bound cannot be constructed
    silently.
    """

    du: float
    dv: float
    dx: float
    dy: float

    @property
    def coordinate_volume(self) -> float:
        return self.du * self.dv * self.dx * self.dy

    def validate(self, w: float) -> None:
        if w != 0.0 and self.du >= max_slab_du(w):
            raise ValueError(
                f"slab u-extent {self.du} is not inside the conjugate "
                f"bound {max_slab_du(w)} at w={w}")


def sprinkle(slab: Slab, rho: float, rng: np.random.Generator
             ) -> np.ndarray:
    """Poisson sprinkling of the slab at density `rho`, as `(N, 4)`.

    Uniform in `(u, v, x, y)` with `N ~ Poisson(rho * coordinate
    volume)`. Because `sqrt(-g) = 1` this IS the covariant Poisson
    process of the true 4-volume, with nothing to calibrate -- and the
    same points serve both arms (design §4.1).
    """

    n = int(rng.poisson(rho * slab.coordinate_volume))
    pts = np.empty((n, 4), dtype=float)
    pts[:, 0] = rng.uniform(0.0, slab.du, n)
    pts[:, 1] = rng.uniform(0.0, slab.dv, n)
    pts[:, 2] = rng.uniform(-slab.dx / 2.0, slab.dx / 2.0, n)
    pts[:, 3] = rng.uniform(-slab.dy / 2.0, slab.dy / 2.0, n)
    return pts
