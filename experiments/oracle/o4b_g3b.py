"""G3b's inside probe: an interval, not a ray.

The two outside probes move one way and the margin grows without
bound, so `o4b_g3a.place_row` finds their first satisfying placement by
bracketing exponentially and bisecting. The inside probe is not like
that, in two ways, and both of them break that routine.

FIRST, THE SUCCESS REGION IS AN INTERVAL. One intermediate time `t_x`
carries two legs: `p -> x` needs `t_x` late enough, `x -> q` needs it
early enough. Raising `t_x` widens one margin and narrows the other.
The satisfying set is therefore the INTERSECTION of a ray up and a ray
down, and an exponential bracket that overshoots lands outside it and
concludes there is nothing there.

SECOND, AND LESS OBVIOUS: THE FREE VARIABLE IS `t_x`, NOT `dt`.
`place_row` searches in `dt`. Here each leg's `dt` is a DIFFERENCE
formed from `t_x` -- `t_x - t_p` and `t_q - t_x` -- and each of those
subtractions rounds. A `dt1` that satisfies leg one need not be
reachable by any representable `t_x`, and the `t_x` that comes nearest
produces a `dt2` that was never the one examined. Searching in `dt`
and mapping back would be reasoning in the wrong coordinate. So both
endpoints are located in `t_x` bit space, where the thing being chosen
actually lives.

WHAT MONOTONICITY SURVIVES. Each leg separately is still monotone in
`t_x`: `t_x - t_p` is non-decreasing in `t_x` and `t_q - t_x` is
non-increasing, and `t_min` and `err` do not move with `t_x` at all --
they are fixed by the geometry. So each ENDPOINT can be bisected; it
is only their conjunction that cannot. That is the whole difference,
and it is why this module bisects twice and intersects rather than
searching once.

NOTHING IS ASSUMED FROM THE ENDPOINTS. Having found `[lo, hi]` and
picked a time inside it, both conditions are recomputed from the
formed differences at that time. The endpoints are located by a
predicate; the verdict is taken from a recomputation. An earlier stage
of this program computed a realized margin and then reported the
answer instead, and the two are not the same claim.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from o4b_g3a import _bits, _from_bits  # noqa: E402


def leg_margin(dt: float, t_min: float, err: float) -> float:
    """The realized decision margin of one leg, in the predicate's own
    terms: how far past the error band the formed `dt` actually sits.

    `>= eta` is the condition; `> 0` alone is what the frozen G3 asked
    for, and that is the margin the abort found sitting at zero."""

    return (dt - t_min) - err


def _bisect_endpoint(seed: float, ok, upward: bool,
                     steps: int = 4096) -> dict:
    """The first representable `t_x` on `seed`'s satisfying side.

    `ok(t_x)` must be monotone in the float order along `upward` --
    true for each leg separately, false for their conjunction, which is
    the reason this is called twice rather than once."""

    comparisons = 1
    if not math.isfinite(seed) or seed < 0.0:
        return {"t_x": seed, "reached": False,
                "comparisons": comparisons,
                "why": "seed time is not a usable float"}
    if ok(seed):
        return {"t_x": seed, "reached": True,
                "comparisons": comparisons, "ulp_distance": 0}

    base = _bits(seed)
    sign = 1 if upward else -1

    def at(k: int) -> float:
        n = base + sign * k
        return _from_bits(n) if n >= 0 else math.nan

    span = 1
    while span < 2 ** 62:
        candidate = at(span)
        comparisons += 1
        if not math.isfinite(candidate) or candidate < 0.0:
            return {"t_x": candidate, "reached": False,
                    "comparisons": comparisons,
                    "why": ("ran out of representable times before "
                            "the margin was reached")}
        if ok(candidate):
            break
        span *= 2
    else:                                          # pragma: no cover
        return {"t_x": at(span), "reached": False,
                "comparisons": comparisons,
                "why": "bracket exhausted the exponent range"}

    low, high = span // 2, span           # low fails, high satisfies
    while high - low > 1 and comparisons < steps:
        mid = (low + high) // 2
        comparisons += 1
        if ok(at(mid)):
            high = mid
        else:
            low = mid
    return {"t_x": at(high), "reached": True,
            "comparisons": comparisons, "ulp_distance": high}


def inside_range(t_p: float, t_q: float,
                 t_min_1: float, err_1: float,
                 t_min_2: float, err_2: float,
                 eta: float) -> dict:
    """The representable times `t_x` at which BOTH legs clear `eta`.

    Returns the interval and whether it is non-empty. Emptiness is a
    fact about the cluster, not a failure: it is exactly the case
    `W_robust > 0` was meant to exclude in advance, and a cluster that
    reaches here empty is one whose eligibility predicate and whose
    construction disagree -- which the caller must record, not paper
    over."""

    def leg_1(t_x: float) -> bool:
        return leg_margin(t_x - t_p, t_min_1, err_1) >= eta

    def leg_2(t_x: float) -> bool:
        return leg_margin(t_q - t_x, t_min_2, err_2) >= eta

    low = _bisect_endpoint(t_p + t_min_1 + err_1 + eta, leg_1,
                           upward=True)
    high = _bisect_endpoint(t_q - t_min_2 - err_2 - eta, leg_2,
                            upward=False)
    comparisons = low["comparisons"] + high["comparisons"]
    if not (low["reached"] and high["reached"]):
        return {"non_empty": False, "lo": low, "hi": high,
                "comparisons": comparisons,
                "why": ("one leg has no representable satisfying "
                        "time, so their intersection is empty")}
    if low["t_x"] > high["t_x"]:
        return {"non_empty": False, "lo": low, "hi": high,
                "comparisons": comparisons,
                "why": ("the two legs' ranges do not overlap: "
                        "widening one margin narrows the other past "
                        "eta, and no single time satisfies both")}
    return {"non_empty": True, "lo": low, "hi": high,
            "comparisons": comparisons}


def place_inside(t_p: float, t_q: float,
                 t_min_1: float, err_1: float,
                 t_min_2: float, err_2: float,
                 eta: float) -> dict:
    """One inside-probe time, with BOTH conditions re-verified at it.

    The chosen time is the bit-space midpoint of the interval, not an
    endpoint: an endpoint is by construction the first value that
    passes, so it sits one representable step from failing, and a
    probe placed there tests the search as much as the wrapper."""

    span = inside_range(t_p, t_q, t_min_1, err_1, t_min_2, err_2, eta)
    if not span["non_empty"]:
        return {"reached": False, "why": span["why"],
                "search_comparisons": span["comparisons"],
                "range": span}

    lo, hi = span["lo"]["t_x"], span["hi"]["t_x"]
    t_x = _from_bits((_bits(lo) + _bits(hi)) // 2)

    # recomputed from the FORMED differences at the chosen time --
    # not carried over from the endpoints that located it
    dt_1, dt_2 = t_x - t_p, t_q - t_x
    m_1 = leg_margin(dt_1, t_min_1, err_1)
    m_2 = leg_margin(dt_2, t_min_2, err_2)
    reached = m_1 >= eta and m_2 >= eta
    return {
        "reached": reached,
        "t_x": t_x, "lo": lo, "hi": hi,
        "dt_1": dt_1, "dt_2": dt_2,
        "realized_margin_1": m_1, "realized_margin_2": m_2,
        "realized_margin": min(m_1, m_2),
        "search_comparisons": span["comparisons"],
        "ulp_width": _bits(hi) - _bits(lo),
        "why": None if reached else (
            "the interval endpoints satisfied their own legs but the "
            "midpoint does not satisfy both -- the recomputation, not "
            "the search, is what decides"),
        "range": span,
    }
