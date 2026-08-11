"""O2 price ladder: the [TO SIZE] measurement of the volume oracle.

Measures what tightening the certified interval COSTS, on the
NEIGHBOR configuration (12, 18, 8.0)M -- deliberately NOT the frozen
(12, 18, 8.5)M. The oracle is deterministic, but the freeze
discipline's point is that nothing about the frozen run is tuned
after seeing its number, so the price is measured next door and the
frozen volume stays unobserved until the PR-O3 execution.

Design notes (both revisions were forced by what the first attempts
actually did, and are recorded so the artifact is readable without
this history):

1. The ladder is ONE adaptive pass that records the cost each time
   the running interval first crosses a target -- not one restarted
   run per rung. The restart version discarded every earlier rung's
   refinement, cost ~3x for the same numbers, and emitted nothing
   until the end, so a slow rung was indistinguishable from a hung
   one. This version streams its curve and publishes the full
   (calls, cells, ratio, wall) trace, so the convergence exponent is
   auditable rather than asserted.

2. The convergence is FLOOR-LIMITED, so a pure power-law
   extrapolation is the wrong model: every cell inside the diamond
   carries its center flight-time's uncertainty however small the
   cell becomes, and only a larger `n_sub` shrinks that. The
   analysis therefore fits `ratio = floor + C * calls^slope` and
   reports, per quadrature setting, both the asymptote and the cost
   to reach the frozen 0.01 -- which is the configuration guidance
   PR-O3 actually needs.

Targets not reached inside the cost caps are reported as
EXTRAPOLATIONS, in fields that cannot be mistaken for measured
crossings. The Monte Carlo cross-check (seed 40000271, a spent
diagnostic stream in the probe ledger) is DIAGNOSTIC only:
disjointness would raise an investigation flag, never a verdict.

Run:      python experiments/oracle/o2_price_ladder.py
Reanalyse: python experiments/oracle/o2_price_ladder.py --reanalyze
  (recomputes the derived fields from the artifact's own stored
  curve, without repeating the multi-hour measurement)
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from volume_oracle import (  # noqa: E402
    OracleConfig,
    assemble,
    mc_diagnostic,
    quadrature_floor_diagnostic,
)

TARGETS = (0.10, 0.05, 0.03, 0.02, 0.01)
FROZEN_TARGET = 0.01
NEIGHBOR = (12.0, 18.0, 8.0)   # NOT the frozen (12, 18, 8.5)
MC_SEED = 40_000_271           # observed diagnostic stream (ledger)
MC_N = 20_000
N_SUB = 16
N_SUB_LADDER = (16, 32, 64)
MAX_CALLS = 120_000
MAX_WALL_S = 9_000.0           # 150 min global cap

_ARTIFACT = (Path(__file__).resolve().parents[2]
             / "docs" / "prereg" / "p14_oracle_price.json")


# fit windows: fraction of the trace (from the end) used for the
# log-log fit. The spread of the extrapolation across windows is the
# honest error bar of the model -- a single window makes the number
# look sharper than the model is (review ruling on PR #65: the
# window-sensitivity range moved the frozen-target calls by tens of
# thousands, so the range is reported, never one value alone)
FIT_WINDOWS = (0.50, 0.33, 0.67)


def _fit_refine(curve: list[dict], floor: float, window: float,
                ) -> tuple[float, float] | None:
    """Least-squares slope of log(ratio - floor) vs log(calls) over
    the trailing `window` fraction of the trace. Fitting the
    REFINABLE part, not the total, is what keeps the extrapolation
    consistent with the floor diagnostic."""

    pts = [(math.log(s["calls"]), math.log(s["ratio"] - floor))
           for s in curve
           if s["calls"] > 0 and s["ratio"] - floor > 0.0]
    pts = pts[int(len(pts) * (1.0 - window)):]
    if len(pts) < 3:
        return None
    n = len(pts)
    mx = sum(x for x, _ in pts) / n
    my = sum(y for _, y in pts) / n
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    if sxx == 0.0:
        return None
    slope = sum((x - mx) * (y - my) for x, y in pts) / sxx
    return slope, my - slope * mx


def _calls_range(target: float, fit_floor: float, proj_floor: float,
                 curve: list[dict]) -> dict | None:
    """The target extrapolation across every fit window: {primary,
    min, max} calls, or None when no window can fit or the target
    sits at/below the projection floor.

    Two DISTINCT floors (review R1 on this PR): the refinable
    component E(calls) = ratio - floor is fitted with the floor of
    the quadrature the curve was actually MEASURED at (`fit_floor`,
    always the N_SUB floor) -- refitting the same curve with a
    candidate's smaller floor absorbs the floor difference into the
    power law and flattens the slope artificially. The candidate's
    floor (`proj_floor`) enters only when the fixed fit is projected
    onto the target: calls = E^-1(target - proj_floor)."""

    if target <= proj_floor:
        return None
    vals = []
    for w in FIT_WINDOWS:
        fit = _fit_refine(curve, fit_floor, w)
        if fit is None:
            continue
        v = _calls_for(target, proj_floor, fit)
        if v is not None:
            vals.append(v)
    if not vals:
        return None
    primary_fit = _fit_refine(curve, fit_floor, FIT_WINDOWS[0])
    primary = (_calls_for(target, proj_floor, primary_fit)
               if primary_fit else vals[0])
    return {"primary": primary, "min": min(vals), "max": max(vals)}


def _calls_for(target: float, floor: float,
               fit: tuple[float, float]) -> float | None:
    """Calls at which floor + C calls^slope first reaches `target`,
    or None when the target lies at or below the floor (unreachable
    by cell refinement at that quadrature)."""

    if target <= floor:
        return None
    slope, intercept = fit
    return math.exp((math.log(target - floor) - intercept) / slope)


def _summarize(curve: list[dict], crossings: dict, final: dict,
               floors: list[dict], rate: float) -> dict:
    """Derived analysis -- a pure function of the stored trace, so
    `--reanalyze` reproduces it without repeating the run."""

    measured = [{"target_ratio": t, **crossings[t]}
                for t in TARGETS if t in crossings]
    base = next((f for f in floors if f["n_sub"] == N_SUB), None)
    floor = base["floor_ratio"] if base else 0.0
    fit = _fit_refine(curve, floor, FIT_WINDOWS[0])

    extrapolated = []
    for t in TARGETS:
        if t in crossings or fit is None:
            continue
        rng = _calls_range(t, floor, floor, curve)
        row = {"target_ratio": t}
        if rng is None:
            row["estimated_reachable_at_this_n_sub"] = False
            row["reason"] = (
                f"target {t:g} is at or below the estimated "
                f"quadrature floor {floor:.5f} at n_sub={N_SUB}: "
                f"cell refinement alone is not expected to reach "
                f"it, a larger n_sub would be required")
        else:
            row["estimated_reachable_at_this_n_sub"] = True
            row["calls_extrapolated"] = rng["primary"]
            row["calls_extrapolated_range"] = [rng["min"],
                                               rng["max"]]
            row["wall_s_extrapolated"] = (rng["primary"]
                                          / max(rate, 1e-9))
        extrapolated.append(row)

    plan = []
    for f in floors:
        # fit with the MEASURED floor (N_SUB), project with the
        # candidate's floor -- see _calls_range
        rng = _calls_range(FROZEN_TARGET, floor, f["floor_ratio"],
                           curve)
        # a call costs roughly linearly in n_sub (the quadrature is
        # the inner loop), so scale the measured rate accordingly --
        # a MODEL assumption, like everything else in this table
        scaled_rate = rate * N_SUB / f["n_sub"]
        plan.append({
            "n_sub": f["n_sub"],
            "floor_ratio": f["floor_ratio"],
            "frozen_target_estimated_reachable": rng is not None,
            "calls_extrapolated": (rng["primary"] if rng else None),
            "calls_extrapolated_range": ([rng["min"], rng["max"]]
                                         if rng else None),
            "wall_s_extrapolated": (rng["primary"]
                                    / max(scaled_rate, 1e-9)
                                    if rng else None),
        })

    return {
        "measured_crossings": measured,
        "extrapolated_targets": extrapolated,
        "convergence_fit": (
            {"model": "ratio = floor + exp(intercept) * calls^slope",
             "log_log_slope": fit[0], "log_intercept": fit[1],
             "floor_ratio_used": floor,
             "fit_windows": list(FIT_WINDOWS),
             "note": ("fitted on the trailing part of the trace's "
                      "REFINABLE component (ratio - floor); the "
                      "window-to-window spread is reported as "
                      "calls_extrapolated_range because the single "
                      "number is sharper than the model deserves. "
                      "An extrapolation aid, not a certified claim "
                      "and not an execution cap.")}
            if fit else None),
        "quadrature_floor": {
            "rows": floors,
            "binding_lever_estimate": (
                "quadrature" if floor > 0.5 * final["ratio"]
                else "cell-refinement"),
            "note": ("DIAGNOSTIC ESTIMATE. floor_ratio is the ratio "
                     "cell refinement alone is modelled to converge "
                     "to at that n_sub, because every interior cell "
                     "keeps its center flight-time's uncertainty. "
                     "The n_sub=32/64 rows are NOT full-integrator "
                     "runs: they combine flight-time widths at a "
                     "few probe points with a linear cost model. "
                     "Estimated standalone (support integral on a "
                     "fixed grid with the fast non-certified "
                     "solver, widths from the certified solver); "
                     "never part of a certified interval."),
        },
        "frozen_target_plan": {
            "target_ratio": FROZEN_TARGET,
            "rows": plan,
            "note": ("configuration GUIDANCE for PR-O3, not frozen "
                     "numbers: every row is a model-based "
                     "extrapolation (floor model + fitted exponent "
                     "+ linear cost scaling), and the review ruling "
                     "on PR #65 forbids freezing any of them as an "
                     "execution cap before the mode-width "
                     "diagnostic has run"),
            "reading": ("n_sub=16 is the current RECOMMENDED PLAN. "
                        "'A larger n_sub is counterproductive' is a "
                        "MODEL-BASED EXPECTATION -- it buys a lower "
                        "floor that the model says never binds, "
                        "while every call costs proportionally "
                        "more -- not a full-integrator measurement. "
                        "The refinement exponent degrades over the "
                        "run; the allowed statement is that cell "
                        "refinement is the LIKELY binding lever and "
                        "the first-order mode a CANDIDATE carrier "
                        "of the remaining width (O(h) per cell "
                        "against the tangent mode's O(h^2)). The "
                        "mode-width diagnostic on the neighbor "
                        "configuration (o2_mode_width.py, "
                        "p14_oracle_mode_width.json) is what "
                        "settles it, and it must run BEFORE any "
                        "O3 budget or algorithm decision."),
        },
    }


def _measure() -> dict:
    cfg = OracleConfig(*NEIGHBOR, m=1.0, target_ratio=FROZEN_TARGET,
                       n_sub=N_SUB, max_calls=MAX_CALLS,
                       max_wall_s=MAX_WALL_S, max_depth=16,
                       init_rho=12, init_psi=12)
    last = [0.0]

    def show(s: dict) -> None:
        if s["wall_s"] - last[0] < 30.0:
            return
        last[0] = s["wall_s"]
        print(f"  calls={s['calls']:6d} cells={s['cells']:6d} "
              f"ratio={s['ratio']:.5f} "
              f"V=[{s['v_lo']:.4f}, {s['v_hi']:.4f}] "
              f"t={s['wall_s']:.0f}s", flush=True)

    print(f"O2 price ladder on neighbor anchors {NEIGHBOR} "
          f"(frozen 8.5 stays unobserved)", flush=True)
    res = assemble(cfg, targets=list(TARGETS), progress=show)
    v = res["v"]
    print(f"assembly {res['status']}: V=[{float(v.lo):.5f}, "
          f"{float(v.hi):.5f}] ratio={res['ratio']:.5f} "
          f"calls={res['calls']} cells={res['cells']} "
          f"wall={res['wall_s']:.0f}s", flush=True)
    mc = mc_diagnostic(OracleConfig(*NEIGHBOR, m=1.0), MC_N, MC_SEED)
    overlap = not (mc["ci95"][1] < float(v.lo)
                   or mc["ci95"][0] > float(v.hi))
    print(f"MC diagnostic: {mc['estimate']:.4f} +- {mc['se']:.4f} "
          f"overlap={overlap}", flush=True)
    return {
        # OUTWARD endpoints: serialized intervals must still contain
        # the true value (PR #67 review R1)
        "final": {"status": res["status"], "v_lo": v.lo_float(),
                  "v_hi": v.hi_float(), "ratio": res["ratio"],
                  "calls": res["calls"], "cells": res["cells"],
                  "modes": res["modes"],
                  "raw_width_by_mode": res["raw_width_by_mode"],
                  "raw_total_before_intersection":
                      res["raw_total_before_intersection"],
                  "certified_total_after_intersection":
                      res["certified_total_after_intersection"],
                  "intersection_active":
                      res["intersection_active"],
                  "wall_s": res["wall_s"]},
        "crossings": res["crossings"],
        "curve": res["curve"],
        "mc_diagnostic": {**mc, "overlap_with_certified": overlap,
                          "role": ("diagnostic only -- disjointness "
                                   "raises an investigation flag, "
                                   "never a verdict")},
    }


def main() -> None:
    ap = argparse.ArgumentParser(allow_abbrev=False)
    ap.add_argument("--reanalyze", action="store_true",
                    help=("recompute derived fields from the stored "
                          "curve without repeating the measurement"))
    args = ap.parse_args()

    if args.reanalyze:
        prior = json.loads(_ARTIFACT.read_text(encoding="utf-8"))
        run = {k: prior[k] for k in
               ("final", "curve", "mc_diagnostic")}
        run["crossings"] = {r["target_ratio"]: {
            k: r[k] for k in
            ("calls", "cells", "ratio", "v_lo", "v_hi", "wall_s")}
            for r in prior["measured_crossings"]}
        print("reanalyzing the stored curve "
              f"({len(run['curve'])} samples)", flush=True)
    else:
        run = _measure()

    cfg = OracleConfig(*NEIGHBOR, m=1.0)
    v_mid = 0.5 * (run["final"]["v_lo"] + run["final"]["v_hi"])
    print("measuring the quadrature floor "
          f"at n_sub in {N_SUB_LADDER}", flush=True)
    floors = quadrature_floor_diagnostic(cfg, v_mid,
                                         n_subs=N_SUB_LADDER)
    for f in floors:
        print(f"  n_sub={f['n_sub']:3d}: mean flight-time width "
              f"{f['mean_flight_time_width']:.3e} -> floor ratio "
              f"{f['floor_ratio']:.5f}", flush=True)

    rate = run["final"]["calls"] / max(run["final"]["wall_s"], 1e-9)
    derived = _summarize(run["curve"], run["crossings"],
                         run["final"], floors, rate)
    for row in derived["measured_crossings"]:
        print(f"  MET     {row['target_ratio']:g}: {row['calls']} "
              f"calls, {row['cells']} cells, {row['wall_s']:.0f}s",
              flush=True)
    for row in derived["extrapolated_targets"]:
        if row["estimated_reachable_at_this_n_sub"]:
            lo, hi = row["calls_extrapolated_range"]
            print(f"  NOT MET {row['target_ratio']:g}: extrapolated "
                  f"{row['calls_extrapolated']:.0f} calls "
                  f"(window range {lo:.0f}-{hi:.0f})", flush=True)
        else:
            print(f"  NOT MET {row['target_ratio']:g}: below the "
                  f"estimated quadrature floor at n_sub={N_SUB}",
                  flush=True)
    print("  binding lever (estimate): "
          f"{derived['quadrature_floor']['binding_lever_estimate']}",
          flush=True)

    art = {
        "scope": ("[TO SIZE] price ladder on the NEIGHBOR anchors "
                  "(12, 18, 8.0)M; the frozen (12, 18, 8.5)M volume "
                  "is deliberately unobserved until the PR-O3 "
                  "execution. ONE adaptive pass records each "
                  "crossing; targets past the cost caps are "
                  "EXTRAPOLATED, never reported as measured."),
        "neighbor_anchors": {"r_in": NEIGHBOR[0],
                             "r_out": NEIGHBOR[1],
                             "dt": NEIGHBOR[2], "m": 1.0},
        "config": {"n_sub": N_SUB, "k_micro": cfg.k_micro,
                   "d_switch": cfg.d_switch,
                   "max_calls": MAX_CALLS, "max_wall_s": MAX_WALL_S,
                   "max_depth": 16, "init_grid": [12, 12]},
        "host": {"machine": platform.machine(),
                 "python": platform.python_version()},
        "final": run["final"],
        **derived,
        "curve": run["curve"],
        "mc_diagnostic": run["mc_diagnostic"],
    }
    with open(_ARTIFACT, "w", encoding="utf-8", newline="\n") as f:
        json.dump(art, f, ensure_ascii=False, indent=1)
        f.write("\n")
    print(f"artifact: {_ARTIFACT}", flush=True)


if __name__ == "__main__":
    main()
