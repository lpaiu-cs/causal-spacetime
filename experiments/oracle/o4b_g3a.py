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


def geometry_radii(name: str, r1, r2) -> tuple[float, float]:
    """Resolve a table geometry's radii, filling the frozen anchors and
    sampling-box edges from `o4_sizing` rather than transcribing them."""

    if r1 is None:
        return base.R_IN, base.R_OUT
    edges = {"R_LO": base.R_LO, "R_HI": base.R_HI}
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


def place_row(t_min: float, err: float, eta: float,
              above: bool) -> dict:
    """Place row A or row B so the REALIZED margin reaches `eta`.

    The nominal time is `t_min +/- (err + eta)`, but forming it rounds,
    and the margin recomputed from the formed value can land short --
    it does, on more than half the table, by a few ulps. So the
    placement is nudged outward one ulp at a time with `nextafter`,
    within a PRE-FROZEN SEARCH BUDGET of `MAX_NUDGES` steps.

    The budget is not a guarantee that any shortfall closes. One step
    gains `ulp(dt)`, and where `t_min` nearly cancels `err` the probe
    time collapses toward zero and that gain collapses with it. Not
    reaching the margin inside the budget is recorded as an
    availability outcome, `construction-unavailable`, and is never
    reclassified as a mismatch or an instrument failure.

    The budget counts MOVES: the initial placement is checked first,
    then at most `MAX_NUDGES` outward steps. A failed final check does
    not move again.

    Checking the margin against the pre-computed `(t_min, err)` is
    checking it against the wrapper's own: the events carry the same
    `theta`, so the wrapper recovers the same `dpsi`, and `flight_time`
    is deterministic. Those two facts are themselves checked."""

    direction = math.inf if above else -math.inf
    dt = t_min + err + eta if above else t_min - err - eta
    for nudges in range(g3.MAX_NUDGES + 1):
        margin = abs(dt - t_min) - err
        if margin >= eta:
            return {"dt": dt, "nudges": nudges, "reached": True,
                    "realized_margin": margin}
        if nudges == g3.MAX_NUDGES:
            break                 # the budget is spent; do NOT move
        dt = math.nextafter(dt, direction)
    return {"dt": dt, "nudges": g3.MAX_NUDGES, "reached": False,
            "realized_margin": abs(dt - t_min) - err}


def check_case(name: str, r1: float, r2: float, theta: float,
               tol: float, eta: float = g3.ETA) -> dict:
    """One case of the table: rows A, B and C at one (geometry, angle).

    Every quantity the decision turns on is recorded, so a failure is
    a coordinate rather than a bare disagreement. A row passes only if
    the tri-state is right AND -- for A and B -- the realized margin
    reaches `eta`; getting the right answer with too little margin is
    not the contract."""

    dpsi = g3.wrapper_dpsi(theta)
    t_min, err = s1.flight_time(r1, r2, dpsi, s1.M, tol)
    recovered = recovered_dpsi(*_events(r1, r2, theta, 0.0))
    rows = []
    for row, above, want in (("A", True, "true"),
                             ("B", False, "false")):
        placed = place_row(t_min, err, eta, above)
        dt = placed["dt"]
        if not above and dt < 0.0:
            rows.append({
                "row": row, "outcome": "construction-unavailable",
                "why": ("t_min - err - eta is negative, so the "
                        "predicate would short-circuit and this row "
                        "would duplicate row D"),
                "dt": dt, **placed})
            continue
        if not placed["reached"]:
            rows.append({
                "row": row, "outcome": "construction-unavailable",
                "why": (f"the realized margin stayed below eta after "
                        f"{g3.MAX_NUDGES} nudges"),
                "want": want, **placed})
            continue
        p, q = _events(r1, r2, theta, dt)
        got = s1.causal_relation(p, q, s1.M, tol)
        label = ("undecided" if got is None
                 else ("true" if got else "false"))
        rows.append({
            "row": row, "dt": dt, "want": want, "got": label,
            "outcome": "pass" if label == want else "FAIL",
            "abs_dt_minus_t_min": abs(dt - t_min),
            "nudges": placed["nudges"],
            "realized_margin": placed["realized_margin"],
            "margin_reaches_eta": placed["realized_margin"] >= eta,
        })

    # row C wants the undecided answer, so it carries no margin at all
    p, q = _events(r1, r2, theta, t_min)
    got = s1.causal_relation(p, q, s1.M, tol)
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
        "family": _family_of(r1, r2, dpsi, tol),
        "t_min": t_min, "err": err, "eta": eta,
        "rows": rows,
    }


def _family_of(r1: float, r2: float, dpsi: float, tol: float) -> str:
    details: dict = {}
    s1.flight_time(r1, r2, dpsi, s1.M, tol, details)
    return details["family"]


def check_row_d(tol: float) -> dict:
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
        got = s1.causal_relation(p, q, s1.M, tol)
    finally:
        s1.flight_time = original
    return {"case": "negative-dt-short-circuit", **case,
            "got": "false" if got is False else repr(got),
            "solver_calls": called["n"],
            "outcome": ("pass" if got is False and called["n"] == 0
                        else "FAIL")}


def run_g3a(tol: float = s1.DEFAULT_TOL,
            eta: float = g3.ETA) -> dict:
    """The whole table. Returns the record; the caller decides what to
    do with a failure, because that decision is the freeze's, not the
    checker's."""

    cases, unreachable = [], []
    for name, r1_spec, r2_spec in g3.G3A_GEOMETRIES:
        r1, r2 = geometry_radii(name, r1_spec, r2_spec)
        for spec in g3.G3A_ANGLE_SPECS:
            label = f"{name}/{spec[0]}"
            try:
                theta = g3.resolve_theta(spec, r1, r2, s1.M, tol,
                                         base.PSI_MAX)
            except ValueError as exc:
                unreachable.append({
                    "case": label, "why": str(exc),
                    "note": ("the geometry has no such branch -- "
                             "equal radii have a zero arc, so no-turn "
                             "and equal-perihelion do not exist")})
                continue
            cases.append(check_case(label, r1, r2, theta, tol, eta))

    row_d = check_row_d(tol)
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
    return {
        "cases": len(cases),
        "unreachable": unreachable,
        "families_covered": families,
        "covers_every_family": set(families) == set(g3.FAMILIES),
        "construction_unavailable": sum(
            1 for c in cases for r in c["rows"]
            if r["outcome"] == "construction-unavailable"),
        "recovery_failures": recovery_failures,
        "margin_failures": margin_failures,
        "nudges_used": sum(r.get("nudges", 0) for c in cases
                           for r in c["rows"]),
        "row_d": row_d,
        "failures": failures,
        "passed": (not failures and not recovery_failures
                   and not margin_failures
                   and set(families) == set(g3.FAMILIES)),
        "detail": cases,
    }


def solver_determinism(tol: float = s1.DEFAULT_TOL,
                       repeats: int = 3) -> dict:
    """The premise the whole separation argument rests on: identical
    arguments give bit-identical results."""

    args = (13.0, 17.0, g3.wrapper_dpsi(0.3), s1.M, tol)
    first = s1.flight_time(*args)
    same = all(s1.flight_time(*args) == first for _ in range(repeats))
    return {"repeats": repeats, "bit_identical": same,
            "probe": {"r1": args[0], "r2": args[1], "dpsi": args[2]}}


def run_preflight(tol: float = s1.DEFAULT_TOL,
                  eta: float = g3.ETA) -> dict:
    """THE G3a verdict. One place, composing every condition.

    An earlier draft ran the determinism probe beside the table and
    printed it, so a solver that was not deterministic would have been
    reported and then ignored -- and the realized-margin contract was
    computed and likewise not composed. Anything the freeze requires
    has to be in the `passed` that gates the stage, or it is decoration.

    The stage proceeds only on this function's `passed`, and sizing
    measures this function, so the number and the verdict come from the
    same execution."""

    table = run_g3a(tol, eta)
    determinism = solver_determinism(tol)
    conditions = {
        "tri_state_rows": not table["failures"],
        "recovered_dpsi_bit_identical": not table["recovery_failures"],
        "families_covered": table["covers_every_family"],
        "ab_realized_margin": not table["margin_failures"],
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
          f"{len(result['unreachable'])} unreachable, "
          f"families {result['families_covered']}")
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
