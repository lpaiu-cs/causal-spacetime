"""Frozen constants and arithmetic for the G3 redesign.

The design is in `docs/prereg/p14_o4_g3_redesign.md`; this module is
its numeric half, kept apart from both the frozen campaign runner and
the replay diagnostic so that a constant has ONE home and the document
can be pinned against it.

There is no runner here. G3a and G3b are not implemented yet -- their
sizes and thresholds are frozen only after the census -- but every
value that must NOT be chosen after seeing the census lives here now,
which is the whole point of the split:

  * `ETA` is a float-separation margin (design section 4).
  * `MAX_NUDGES` bounds the placement search (design section 4.3).

Both are determined by rounding arguments alone. Neither may be
revisited in light of a frequency the census reports.
"""

from __future__ import annotations

import math

#: The separation margin, in the same units as the flight times.
#: Chosen so that it dominates the rounding in forming a probe time
#: (all magnitudes <= dt = 8.5, `ulp(8.5) = 1.78e-15`, three or four
#: operations) while sitting six decades below the `err` scale, so it
#: cannot drive eligibility. It is NOT chosen from the frequency of
#: `err2 >= 1e-6`, and it is not a bare addend: the runner must CHECK
#: that the realized margin reaches it (design section 4.2).
ETA = 1e-12

#: The exploratory grid the census scans. Eligibility is strict:
#: `L - err1 - err2 - 2*eta > 0`.
ETA_GRID = (0.0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9,
            1e-8, 1e-7, 1e-6)

#: How many `nextafter` steps the probe placement may take before the
#: construction is declared unavailable, frozen HERE rather than after
#: the census -- the cap decides whether a given cluster yields a valid
#: probe or an unavailable one, so choosing it later would let the
#: census pick which clusters count.
#:
#: The bound comes from rounding alone. A candidate is placed at a
#: nominal margin of exactly `ETA`, and the arithmetic that forms it
#: (`t_x` in two or three operations, then `dt`, then the difference)
#: rounds by at most `0.5 * ulp(8.5) = 8.9e-16` each, so the realized
#: margin can fall short by at most about `3.6e-15`. One step moves
#: `dt` by one ulp of the probe time, worst case `ulp(1.0) = 2.2e-16`,
#: so about 17 steps suffice. 64 leaves a factor of ~3.7.
MAX_NUDGES = 64


#: G3b availability, frozen as a DOUBLE structure (prereg re-opening
#: section 2). `N_AVAIL` is a fixed denominator: the first this many
#: `L_S1 > 0` candidates of the G1 stream, on which the eligible and
#: fully-testable rates are reported. `K_G3B` is a separate completion
#: target: keep scanning the SAME stream until this many fully-testable
#: clusters exist.
#:
#: The two cannot be one number. `K / scanned` from a stopping-time
#: sample is not a fixed-denominator availability estimate -- the
#: denominator is itself chosen by the outcome. So the prefix reports
#: availability and the scan merely procures the contract sample.
N_AVAIL = 100_000
K_G3B = 100_000

#: One-sided Clopper-Pearson level for the G3b mismatch CHARACTERISATION.
#: Not an accuracy gate: promoting "mismatch rate <= x" to a verdict
#: sentence would need coverage, sidedness and power frozen separately,
#: and none of those are frozen here.
ALPHA_G3B = 0.05

#: The scan has no cap of its own. It ends where the frozen G1 sample
#: ends, so no new arbitrary number enters the protocol; the value is
#: read from the campaign's own sizing module rather than copied.
def scan_cap() -> int:
    """The end of the frozen G1 sample -- the scan's only bound."""

    import o4_sizing as sz
    return sz.N_G1


def w_robust(length: float, err1: float, err2: float,
             eta: float = ETA) -> float:
    """`L - err1 - err2 - 2*eta`, in the predicate's coordinates."""

    return length - err1 - err2 - 2.0 * eta


def is_eligible(length: float, err1: float, err2: float,
                eta: float = ETA) -> bool:
    """Strictly positive robust window. Strict, not `>=`: a zero window
    admits no probe with a margin to spare."""

    return w_robust(length, err1, err2, eta) > 0.0


def lower_probe_time(lo: float, err1: float,
                     eta: float = ETA) -> float:
    """Probe (iii): strictly below the window's lower edge.

    Split out because eligibility needs it on its own -- whether this
    time is non-negative decides between a lower-boundary solver probe
    and the predicate's negative-`dt` short circuit, which are
    different checks and are counted apart (design section 3.4)."""

    return lo - err1 - eta


def probe_times(lo: float, hi: float, err1: float, err2: float,
                eta: float = ETA) -> dict:
    """The three nominal probe times of design section 3.3.

    Nominal: each is where the algebra puts the probe. Whether the
    REALIZED margin reaches `eta` is a separate check the runner makes
    against the predicate itself, and is not decided here."""

    return {
        "inside": 0.5 * ((lo + err1 + eta) + (hi - err2 - eta)),
        "outside_above": hi + err2 + eta,
        "outside_below": lower_probe_time(lo, err1, eta),
    }


# ------------------------------------- G3a: the deterministic case table

#: The solver's branch structure, which the G3a table must cover. Read
#: off `s1_schwarzschild_cost.flight_time`: `dpsi < 1e-12` is radial and
#: reports `err` exactly 0; `|dpsi - a_eq| <= tol` is the equal-perihelion
#: arc; above `a_eq` is one-turn and below it no-turn. Equal-radius pairs
#: have a zero arc and are therefore always one-turn.
FAMILIES = ("radial", "no-turn", "equal-perihelion", "one-turn")

#: The radial predicate's threshold, quoted so the boundary rows can
#: straddle it. Strict `<`, so the threshold itself is NOT radial.
RADIAL_THRESHOLD = 1e-12

#: Geometries: radius direction, equal radii, and the shell / angular
#: patch edges. Radii are the frozen anchors and sampling-box edges;
#: `None` is filled in from `o4_sizing` so no boundary is hand-copied.
G3A_GEOMETRIES = (
    ("increasing-interior", 13.0, 17.0),
    ("decreasing-interior", 17.0, 13.0),
    ("equal-radius", 15.0, 15.0),
    ("anchor-pair", None, None),           # (R_IN, R_OUT)
    ("shell-inner-edge", "R_LO", 17.0),
    ("shell-outer-edge", 13.0, "R_HI"),
    ("shell-both-edges", "R_LO", "R_HI"),
)

#: Cases are specified in THETA, the coordinate, not in `dpsi`.
#:
#: The wrapper never receives `dpsi`; it recovers it as
#: `acos(clamp(cos theta))`, and that composition does not reach every
#: float. Its smallest nonzero value is `acos(nextafter(1, -inf))` =
#: 1.4901e-08, so NOTHING lies between 0 and that -- and the solver's
#: own radial threshold, 1e-12, sits inside the gap (review R1). A case
#: frozen at an unreachable `dpsi` would pre-compute `(t_min, err)` on
#: one branch while the wrapper took another, and G3a would report
#: `INVALID` against itself.
#:
#: So every case names a `theta` the wrapper can actually be handed,
#: and the branch edges are straddled by ADJACENT binary64 thetas found
#: by bisecting on the family the solver reports.
G3A_ANGLE_SPECS = (
    ("radial", "theta", 0.0),
    ("radial-reachable-edge-below", "radial-edge", "below"),
    ("radial-reachable-edge-above", "radial-edge", "above"),
    ("no-turn-mid", "band-fraction", 0.5),
    ("equal-perihelion-inside", "band-centre", None),
    ("equal-perihelion-lower-edge-below", "band-edge", "lower-below"),
    ("equal-perihelion-lower-edge-above", "band-edge", "lower-above"),
    ("equal-perihelion-upper-edge-below", "band-edge", "upper-below"),
    ("equal-perihelion-upper-edge-above", "band-edge", "upper-above"),
    ("one-turn-beyond", "band-multiple", 1.5),
    ("patch-max-angle", "psi-max", None),
)

#: Which wrapper rows each case runs. D is not here: it is a single
#: negative-`dt` case that must not consult the solver at all, and is
#: checked once rather than per geometry.
G3A_ROWS = ("A", "B", "C")


def family_at(r1: float, r2: float, dpsi: float, m: float,
              tol: float) -> str:
    """The family the solver REPORTS at these arguments.

    Asked rather than re-derived: a second copy of the branch algebra
    here could drift from the one that decides, and then the table
    would be straddling an edge that is not the edge."""

    import s1_schwarzschild_cost as s1
    details: dict = {}
    s1.flight_time(r1, r2, dpsi, m, tol, details)
    return details["family"]


def equal_perihelion_band(r1: float, r2: float, m: float, tol: float,
                          steps: int = 200) -> tuple[float, float]:
    """The closed interval of angles the solver calls
    `equal-perihelion`, located by bisection on its own label.

    Returns `(lower_edge, upper_edge)`. Raises if the geometry has no
    such band at all -- equal-radius pairs have a zero arc and are
    one-turn everywhere, which is a fact about the geometry and must
    surface rather than be papered over."""

    def fam(x: float) -> str:
        return family_at(r1, r2, x, m, tol)

    lo, hi = 1e-9, 3.0
    if fam(lo) != "no-turn" or fam(hi) != "one-turn":
        raise ValueError(
            f"({r1}, {r2}) has no no-turn/one-turn transition to "
            f"bracket: the solver reports {fam(lo)!r} at {lo} and "
            f"{fam(hi)!r} at {hi}")
    inside = None
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        found = fam(mid)
        if found == "equal-perihelion":
            inside = mid
            break
        if found == "no-turn":
            lo = mid
        else:
            hi = mid
    if inside is None:
        raise ValueError(
            f"({r1}, {r2}): bisection never landed in the "
            f"equal-perihelion band between {lo!r} and {hi!r}")

    def edge(low: float, high: float, want_low: bool) -> float:
        """Bisect the transition between the band and its neighbour."""

        for _ in range(steps):
            mid = 0.5 * (low + high)
            in_band = fam(mid) == "equal-perihelion"
            if in_band == want_low:
                high = mid
            else:
                low = mid
            if high - low <= 0.0 or math.nextafter(low, high) >= high:
                break
        return high if want_low else low

    return edge(lo, inside, True), edge(inside, hi, False)


def wrapper_dpsi(theta: float) -> float:
    """The angle the wrapper recovers from a `theta` coordinate.

    Written out rather than imported so this module states the
    composition it is reasoning about: for the G3 geometry the two
    events differ only in polar angle, so `cosang` is `cos(theta)`
    exactly and the recovery is `acos(clamp(cos theta))`."""

    return math.acos(max(-1.0, min(1.0, math.cos(theta))))


def straddle_theta(predicate, low: float, high: float,
                   steps: int = 200) -> tuple[float, float]:
    """Adjacent binary64 thetas whose RECOVERED angles fall on opposite
    sides of `predicate`.

    Adjacent in theta, because theta is what the wrapper is handed.
    Straddling in `dpsi` space would name floats the wrapper cannot
    produce."""

    if predicate(wrapper_dpsi(low)) == predicate(wrapper_dpsi(high)):
        raise ValueError(
            f"theta bracket [{low!r}, {high!r}] does not straddle the "
            f"predicate: both sides agree")
    for _ in range(steps):
        mid = 0.5 * (low + high)
        if not low < mid < high:
            break
        if predicate(wrapper_dpsi(mid)) == predicate(
                wrapper_dpsi(low)):
            low = mid
        else:
            high = mid
        if math.nextafter(low, high) >= high:
            break
    return low, high


def resolve_theta(spec: tuple, r1: float, r2: float, m: float,
                  tol: float, psi_max: float) -> float:
    """One case specification, resolved to a binary64 `theta`.

    Deterministic: the same table, the same solver, the same numbers.
    Freezing the SPEC freezes the case."""

    _, kind, value = spec
    if kind == "theta":
        return float(value)
    if kind == "psi-max":
        return psi_max
    if kind == "radial-edge":
        # the reachability edge: where cos(theta) stops rounding to 1
        below, above = straddle_theta(lambda d: d == 0.0, 0.0, 1e-6)
        return below if value == "below" else above

    lower, upper = equal_perihelion_band(r1, r2, m, tol)
    centre = 0.5 * (lower + upper)
    if kind == "band-fraction":
        return lower * float(value)
    if kind == "band-multiple":
        return upper * float(value)
    if kind == "band-centre":
        return centre
    if kind == "band-edge":
        which, side = value.split("-")
        bracket = ((lower * 0.5, centre) if which == "lower"
                   else (centre, upper * 1.5))
        below, above = straddle_theta(
            lambda d: family_at(r1, r2, d, m, tol) == "equal-perihelion",
            *bracket)
        return below if side == "below" else above
    raise ValueError(f"unknown case spec {spec!r}")

