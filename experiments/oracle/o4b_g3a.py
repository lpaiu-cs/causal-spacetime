"""G3a: the wrapper semantics contract, run BEFORE any seed is drawn.

The O4 campaign met the wrapper contract for the first time after
spending twelve hours of sample, and died there with no verdict. G3a
is the repair: a deterministic, branch-stratified case table, executed
on the exact-freeze checkout, whose failure ends the stage `INVALID`
while every campaign scalar is still unspent.

Each case fixes the answer BY CONSTRUCTION rather than by geometric
intuition. The predicate's own `(t_min, err)` is computed first -- at
the angle the wrapper will recover, never at the `theta` handed in --
and `dt` is then built from it:

    A   dt = t_min + err + eta   ->  True     (decided, above)
    B   dt = t_min - err - eta   ->  False    (decided, below)
    C   dt = t_min               ->  None     (legitimately undecided)
    D   dt < 0                   ->  False    (short circuit, no solver)

Row C passing is the whole repair. The frozen G3 aborted on any
`None`; a tri-state predicate must be allowed to be undecided exactly
where its own contract says it is.

Row B runs only where `dt >= 0`. Where `t_min - err - eta` goes
negative the predicate short-circuits and B would pass without ever
exercising the lower error-band path -- a silent duplicate of D. Those
cases are recorded `construction-unavailable` and counted, not skipped
quietly.

Row D is checked once, and by proving the solver is NOT CALLED: a
returned `False` cannot distinguish the short circuit from a solver
that ran and agreed.
"""

from __future__ import annotations

import math
import struct
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parents[0] / "positive_control"))

import o4_g3_redesign as g3  # noqa: E402
import o4_sizing as base  # noqa: E402
import s1_schwarzschild_cost as s1  # noqa: E402

#: Row D's fixed case: a pair whose time separation is negative. The
#: radii and angle are immaterial -- the predicate must not look.
ROW_D_CASE = {"dt": -1.0, "r1": 13.0, "r2": 17.0, "theta": 0.3}


def geometry_radii(name: str, r1, r2,
                   box: tuple[float, float] | None = None,
                   ) -> tuple[float, float]:
    """Resolve a table geometry's radii, filling the frozen anchors and
    sampling-box edges from `o4_sizing` rather than transcribing them.

    `box` (review PR #90 R2): a rung's OWN sampling-box edges
    (r_lo, r_hi); None keeps the frozen M = 1 edges. The anchors are
    absolute and rung-independent either way."""

    if r1 is None:
        return base.R_IN, base.R_OUT
    lo, hi = box if box is not None else (base.R_LO, base.R_HI)
    edges = {"R_LO": lo, "R_HI": hi}
    return edges.get(r1, r1), edges.get(r2, r2)


def _events(r1: float, r2: float, theta: float,
            dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Two events separated by `dt` in time and `theta` in polar angle.

    The wrapper recovers its own angle from these coordinates; the
    caller must have used `wrapper_dpsi(theta)` to build `dt`."""

    return (np.array([0.0, r1, 0.0, 0.0]),
            np.array([dt, r2, theta, 0.0]))


def recovered_dpsi(p: np.ndarray, q: np.ndarray) -> float:
    """The angle recovered FROM THE COORDINATES, by the expression
    `causal_relation` uses.

    This is the half of the bit-identity argument that determinism
    cannot supply: it establishes that the wrapper's `flight_time` call
    carries the same angle the pre-computation used. Comparing
    `wrapper_dpsi(theta)` against itself would establish nothing, which
    is exactly what an earlier draft of this function did."""

    cosang = (math.sin(p[2]) * math.sin(q[2]) * math.cos(p[3] - q[3])
              + math.cos(p[2]) * math.cos(q[2]))
    return math.acos(max(-1.0, min(1.0, cosang)))


def _bits(x: float) -> int:
    """A non-negative double's bit pattern, which orders exactly as the
    value does -- so "the next representable value" is "+1 here"."""

    return struct.unpack("<Q", struct.pack("<d", x))[0]


def _from_bits(n: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", n))[0]


def place_row(t_min: float, err: float, eta: float,
              above: bool) -> dict:
    """Place row A or row B at the FIRST representable time whose
    realized margin reaches `eta`.

    The nominal time is `t_min +/- (err + eta)`, but forming it rounds
    and the margin recomputed from the formed value can land short --
    it does, on more than half the table.

    THERE IS NO STEP CAP. An earlier draft capped the search at 64 and
    recorded anything longer as unavailable; one legitimate case needs
    8,042 single-ulp steps, and discarding it would have been one more
    non-adaptive margin in a stage whose whole failure was
    non-adaptive margins.

    Nor is the search performed one ulp at a time. Walking 8,042 steps
    is harmless once, but G3b places three probes at 100,000 clusters
    and a linear walk there is a new way to stall for hours. The margin
    is monotone in the float ORDER, so the search brackets
    exponentially in bit space and then bisects, which returns the same
    first-satisfying value in a few dozen comparisons instead of
    thousands of moves.

    Termination is on representability, never on a budget:

      * a non-finite `dt` -- nothing further exists to move to;
      * for the lower probe, `dt` reaching zero -- below it the
        predicate short-circuits and the probe stops testing the
        solver at all, which is row D's job.

    `ulp_distance` is the DIAGNOSTIC distance from the nominal
    placement, in representable steps -- not a count of work done, so
    it must not be read as cost. `search_comparisons` is the work.

    MONOTONICITY IS A PRECONDITION, NOT A PROPERTY OF ALL PROBES. This
    routine may be used only where moving one way increases the margin
    without bound: G3a's rows A and B, and G3b's two outside probes.
    G3b's INSIDE probe is different -- widening one leg's margin
    narrows the other's, so its success region is an INTERVAL, and an
    exponential jump can step clean over it. That probe must intersect
    the two legs' representable ranges directly and re-verify both
    conditions at whatever it picks; it must NOT call this function."""

    start = t_min + err + eta if above else t_min - err - eta

    def margin(x: float) -> float:
        return abs(x - t_min) - err

    def reached(x: float) -> bool:
        return math.isfinite(x) and margin(x) >= eta

    if not math.isfinite(start):
        return {"dt": start, "ulp_distance": 0,
                "search_comparisons": 1, "reached": False,
                "why": "non-finite probe time",
                "realized_margin": float("nan")}
    if not above and start < 0.0:
        return {"dt": start, "ulp_distance": 0,
                "search_comparisons": 1, "reached": False,
                "why": ("dt crossed zero, where the predicate "
                        "short-circuits and row B would duplicate "
                        "row D"),
                "realized_margin": margin(start)}
    if reached(start):
        return {"dt": start, "ulp_distance": 0,
                "search_comparisons": 1, "reached": True,
                "realized_margin": margin(start)}

    base = _bits(start)
    sign = 1 if above else -1
    comparisons = 1                      # the `start` check above

    def at(k: int) -> float:
        return _from_bits(base + sign * k)

    # exponential bracket: double the offset until the margin is met.
    # Valid ONLY because the margin is monotone in this direction --
    # see `place_row`'s note on the inside probe, where it is not.
    hi = 1
    while True:
        comparisons += 1
        if not above and base - hi <= 0:          # dt would reach zero
            return {"dt": at(base), "ulp_distance": base,
                    "search_comparisons": comparisons,
                    "reached": False,
                    "why": ("dt reached zero before the margin, where "
                            "the predicate short-circuits"),
                    "realized_margin": margin(_from_bits(0))}
        candidate = at(hi)
        if not math.isfinite(candidate):
            return {"dt": candidate, "ulp_distance": hi,
                    "search_comparisons": comparisons,
                    "reached": False,
                    "why": "non-finite probe time",
                    "realized_margin": float("nan")}
        if reached(candidate):
            break
        hi *= 2

    # bisect for the FIRST offset that reaches it
    lo = hi // 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        comparisons += 1
        if reached(at(mid)):
            hi = mid
        else:
            lo = mid
    return {"dt": at(hi), "ulp_distance": hi,
            "search_comparisons": comparisons, "reached": True,
            "realized_margin": margin(at(hi)),
            "search": "exponential bracket then bisection in float order"}


def check_case(name: str, r1: float, r2: float, theta: float,
               tol: float, eta: float = g3.ETA,
               m: float = s1.M) -> dict:
    """One case of the table: rows A, B and C at one (geometry, angle).

    Every quantity the decision turns on is recorded, so a failure is
    a coordinate rather than a bare disagreement. A row passes only if
    the tri-state is right AND -- for A and B -- the realized margin
    reaches `eta`; getting the right answer with too little margin is
    not the contract."""

    dpsi = g3.wrapper_dpsi(theta)
    t_min, err = s1.flight_time(r1, r2, dpsi, m, tol)
    recovered = recovered_dpsi(*_events(r1, r2, theta, 0.0))
    rows = []
    for row, above, want in (("A", True, "true"),
                             ("B", False, "false")):
        placed = place_row(t_min, err, eta, above)
        dt = placed["dt"]
        if not placed["reached"]:
            rows.append({"row": row,
                         "outcome": "construction-unavailable",
                         "want": want, **placed})
            continue
        p, q = _events(r1, r2, theta, dt)
        got = s1.causal_relation(p, q, m, tol)
        label = ("undecided" if got is None
                 else ("true" if got else "false"))
        rows.append({
            "row": row, "dt": dt, "want": want, "got": label,
            "outcome": "pass" if label == want else "FAIL",
            "abs_dt_minus_t_min": abs(dt - t_min),
            "ulp_distance": placed["ulp_distance"],
            "search_comparisons": placed["search_comparisons"],
            "realized_margin": placed["realized_margin"],
            "margin_reaches_eta": placed["realized_margin"] >= eta,
        })

    # row C wants the undecided answer, so it carries no margin at all
    p, q = _events(r1, r2, theta, t_min)
    got = s1.causal_relation(p, q, m, tol)
    label = ("undecided" if got is None
             else ("true" if got else "false"))
    rows.append({
        "row": "C", "dt": t_min, "want": "undecided", "got": label,
        "outcome": "pass" if label == "undecided" else "FAIL",
        "abs_dt_minus_t_min": abs(t_min - t_min),
        "why_no_margin": ("row C asks for the undecided answer, so a "
                          "separation margin would defeat it"),
    })
    return {
        "case": name, "r1": r1, "r2": r2, "theta": theta,
        "dpsi": dpsi,
        "recovery_is_bit_identical": recovered == dpsi,
        "family": _family_of(r1, r2, dpsi, tol, m),
        "t_min": t_min, "err": err, "eta": eta,
        "rows": rows,
    }


def _family_of(r1: float, r2: float, dpsi: float, tol: float,
               m: float = s1.M) -> str:
    details: dict = {}
    s1.flight_time(r1, r2, dpsi, m, tol, details)
    return details["family"]


def check_row_d(tol: float, m: float = s1.M) -> dict:
    """Row D: the predicate must answer without consulting the solver.

    Proved by making the solver fail. A returned `False` alone cannot
    tell a short circuit from a solver that ran and agreed."""

    case = ROW_D_CASE
    p, q = _events(case["r1"], case["r2"], case["theta"], case["dt"])
    called = {"n": 0}
    original = s1.flight_time

    def forbidden(*args, **kwargs):
        called["n"] += 1
        raise AssertionError(
            "row D: causal_relation consulted the solver on a "
            "negative dt, so the short circuit is not what answered")

    s1.flight_time = forbidden
    try:
        got = s1.causal_relation(p, q, m, tol)
    finally:
        s1.flight_time = original
    return {"case": "negative-dt-short-circuit", **case,
            "got": "false" if got is False else repr(got),
            "solver_calls": called["n"],
            "outcome": ("pass" if got is False and called["n"] == 0
                        else "FAIL")}


def run_g3a(tol: float = s1.DEFAULT_TOL,
            eta: float = g3.ETA, m: float = s1.M,
            geometry: dict | None = None) -> dict:
    """The whole table. Returns the record; the caller decides what to
    do with a failure, because that decision is the freeze's, not the
    checker's."""

    box = ((geometry["r_lo"], geometry["r_hi"])
           if geometry is not None else None)
    cap = (geometry["psi_max"] if geometry is not None
           else base.PSI_MAX)
    cases, unreachable, build_failures = [], [], []
    for name, r1_spec, r2_spec in g3.G3A_GEOMETRIES:
        r1, r2 = geometry_radii(name, r1_spec, r2_spec, box)
        for spec in g3.G3A_ANGLE_SPECS:
            label = f"{name}/{spec[0]}"
            try:
                theta = g3.resolve_theta(spec, r1, r2, m, tol,
                                         cap)
            except g3.ExpectedUnreachable as exc:
                # the geometry says it has no such branch. Still
                # checked against the frozen list: the raiser says
                # WHY, the list says WHERE, and both have to agree.
                unreachable.append({
                    "case": label, "why": str(exc),
                    "expected": label in g3.EXPECTED_UNREACHABLE_LABELS,
                    "note": g3.WHY_EXPECTED_UNREACHABLE})
                continue
            except ValueError as exc:
                # the construction failed. Not an omission at all, and
                # not excused by the label being on the list -- a
                # bracket that went wrong on `equal-radius` would
                # otherwise be waved through as the expected absence.
                build_failures.append({
                    "case": label, "why": str(exc),
                    "on_the_expected_list": (
                        label in g3.EXPECTED_UNREACHABLE_LABELS),
                    "note": ("the case could not be BUILT; the freeze "
                             "accepts a branch the geometry does not "
                             "have, not a search that failed")})
                continue
            cases.append(check_case(label, r1, r2, theta, tol, eta,
                                    m))

    row_d = check_row_d(tol, m)
    failures = [
        {"case": c["case"], "row": r["row"], "want": r["want"],
         "got": r["got"], "family": c["family"], "r1": c["r1"],
         "r2": c["r2"], "theta": c["theta"], "dpsi": c["dpsi"],
         "t_min": c["t_min"], "err": c["err"],
         "realized_margin": r["realized_margin"]}
        for c in cases for r in c["rows"] if r["outcome"] == "FAIL"]
    failures += [{"case": row_d["case"], "row": "D",
                  "got": row_d["got"],
                  "solver_calls": row_d["solver_calls"]}
                 ] if row_d["outcome"] == "FAIL" else []
    recovery_failures = [c["case"] for c in cases
                         if not c["recovery_is_bit_identical"]]
    margin_failures = [
        {"case": c["case"], "row": r["row"], "family": c["family"],
         "realized_margin": r["realized_margin"], "eta": c["eta"]}
        for c in cases for r in c["rows"]
        if r["outcome"] == "pass" and "margin_reaches_eta" in r
        and not r["margin_reaches_eta"]]
    families = sorted({c["family"] for c in cases})
    # Two cases in one geometry can name different thetas that the
    # wrapper recovers to the SAME dpsi -- the recovery lattice is
    # coarse near zero. Both are legitimately named rows and both run,
    # but they exercise one solver state, so they are recorded here
    # rather than counted as independent coverage (review R8).
    by_angle: dict[tuple, list] = {}
    for c in cases:
        by_angle.setdefault((c["r1"], c["r2"], c["dpsi"]),
                            []).append(c["case"])
    # The `radial` / `radial-reachable-edge-below` pairs collapse BY
    # CONSTRUCTION: the case below the reachability edge is defined as
    # the theta whose recovery is still exactly 0, so agreeing with
    # `theta = 0` is what it is there to demonstrate. The rest are
    # incidental -- two independently specified angles that the coarse
    # lattice happens to merge.
    deliberate = {"radial", "radial-reachable-edge-below"}
    duplicates = [
        {"r1": r1, "r2": r2, "dpsi": dpsi, "cases": sorted(names),
         "family": next(c["family"] for c in cases
                        if c["case"] == names[0]),
         "by_construction": {n.split("/", 1)[1]
                             for n in names} == deliberate}
        for (r1, r2, dpsi), names in by_angle.items() if len(names) > 1]
    # Every (geometry, spec) pair lands in exactly one of `cases` or
    # `unreachable`, so these two lists together say whether the table
    # that ran is the table that was frozen.
    unexpected_unreachable = [u for u in unreachable
                              if not u["expected"]]
    built = {c["case"] for c in cases}
    unexpectedly_reachable = sorted(
        g3.EXPECTED_UNREACHABLE_LABELS & built)
    return {
        "cases": len(cases),
        "unreachable": unreachable,
        "unexpected_unreachable": unexpected_unreachable,
        "unexpectedly_reachable": unexpectedly_reachable,
        "case_build_failures": build_failures,
        "duplicate_recovered_angles": duplicates,
        "why_duplicates_are_recorded": (
            "distinct thetas that the wrapper recovers to one dpsi. "
            "Each is a valid named row and each runs, but together "
            "they exercise a single solver state, so they must not be "
            "counted as independent family or state coverage. Those "
            "marked by_construction are the radial reachability-edge "
            "pairs, where agreeing is the demonstration; the rest are "
            "the coarse lattice merging two independent specs"),
        "distinct_solver_states": len(by_angle),
        "families_covered": families,
        "covers_every_family": set(families) == set(g3.FAMILIES),
        "construction_unavailable": sum(
            1 for c in cases for r in c["rows"]
            if r["outcome"] == "construction-unavailable"),
        "recovery_failures": recovery_failures,
        "margin_failures": margin_failures,
        "ulp_distance_total": sum(r.get("ulp_distance", 0)
                                  for c in cases for r in c["rows"]),
        "search_comparisons_total": sum(
            r.get("search_comparisons", 0)
            for c in cases for r in c["rows"]),
        "distance_is_not_cost": ("ulp_distance measures how far the "
                                 "nominal placement sat from a "
                                 "satisfying one; search_comparisons "
                                 "is the work actually done"),
        "row_d": row_d,
        "failures": failures,
        "passed": (not failures and not recovery_failures
                   and not margin_failures
                   and not unexpected_unreachable
                   and not unexpectedly_reachable
                   and not build_failures
                   and set(families) == set(g3.FAMILIES)),
        "detail": cases,
    }


def solver_determinism(tol: float = s1.DEFAULT_TOL,
                       repeats: int = 3,
                       m: float = s1.M) -> dict:
    """The premise the whole separation argument rests on: identical
    arguments give bit-identical results."""

    args = (13.0, 17.0, g3.wrapper_dpsi(0.3), m, tol)
    first = s1.flight_time(*args)
    same = all(s1.flight_time(*args) == first for _ in range(repeats))
    return {"repeats": repeats, "bit_identical": same,
            "probe": {"r1": args[0], "r2": args[1], "dpsi": args[2]}}


def run_preflight(tol: float = s1.DEFAULT_TOL,
                  eta: float = g3.ETA, m: float = s1.M,
                  geometry: dict | None = None) -> dict:
    """THE G3a verdict. One place, composing every condition.

    An earlier draft ran the determinism probe beside the table and
    printed it, so a solver that was not deterministic would have been
    reported and then ignored -- and the realized-margin contract was
    computed and likewise not composed. Anything the freeze requires
    has to be in the `passed` that gates the stage, or it is decoration.

    The stage proceeds only on this function's `passed`, and sizing
    measures this function, so the number and the verdict come from the
    same execution."""

    table = run_g3a(tol, eta, m, geometry)
    determinism = solver_determinism(tol, m=m)
    conditions = {
        "tri_state_rows": not table["failures"],
        "recovered_dpsi_bit_identical": not table["recovery_failures"],
        "families_covered": table["covers_every_family"],
        "ab_realized_margin": not table["margin_failures"],
        "only_expected_unreachable": (
            not table["unexpected_unreachable"]
            and not table["unexpectedly_reachable"]),
        "every_case_built": not table["case_build_failures"],
        "row_d_no_solver_call": table["row_d"]["outcome"] == "pass",
        "solver_determinism": determinism["bit_identical"],
    }
    return {
        "passed": all(conditions.values()),
        "conditions": conditions,
        "failed_conditions": sorted(k for k, ok in conditions.items()
                                    if not ok),
        "on_failure": ("INVALID, before any fresh seed is touched; "
                       "G3b is not run"),
        "table": table,
        "determinism": determinism,
    }


def main() -> None:
    preflight = run_preflight()
    result = preflight["table"]
    print(f"G3a: {result['cases']} cases, "
          f"{len(result['unreachable'])} unreachable "
          f"({len(result['unexpected_unreachable'])} unexpected, "
          f"{len(result['unexpectedly_reachable'])} expected but "
          f"built), {len(result['case_build_failures'])} build "
          f"failures, families {result['families_covered']}")
    print(f"  construction-unavailable rows: "
          f"{result['construction_unavailable']}")
    print(f"  row D: {result['row_d']['outcome']} "
          f"(solver calls {result['row_d']['solver_calls']})")
    print(f"  margin failures: {len(result['margin_failures'])}")
    for name, ok in preflight["conditions"].items():
        print(f"    {'ok ' if ok else 'FAIL'} {name}")
    print(f"  {'PASS' if preflight['passed'] else 'INVALID'}")
    if not preflight["passed"]:
        print(f"  failed: {preflight['failed_conditions']}")


if __name__ == "__main__":
    main()
