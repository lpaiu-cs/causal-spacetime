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


def _fit_refine(curve: list[dict], floor: float,
                ) -> tuple[float, float] | None:
    """Least-squares slope of log(ratio - floor) vs log(calls) over
    the latter half of the trace. Fitting the REFINABLE part, not
    the total, is what makes the extrapolation consistent with the
    floor diagnostic."""

    pts = [(math.log(s["calls"]), math.log(s["ratio"] - floor))
           for s in curve
           if s["calls"] > 0 and s["ratio"] - floor > 0.0]
    pts = pts[len(pts) // 2:]
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
    fit = _fit_refine(curve, floor)

    extrapolated = []
    for t in TARGETS:
        if t in crossings or fit is None:
            continue
        calls = _calls_for(t, floor, fit)
        row = {"target_ratio": t}
        if calls is None:
            row["reachable_at_this_n_sub"] = False
            row["reason"] = (
                f"target {t:g} is at or below the quadrature floor "
                f"{floor:.5f} at n_sub={N_SUB}: no amount of cell "
                f"refinement reaches it, a larger n_sub is required")
        else:
            row["reachable_at_this_n_sub"] = True
            row["calls_extrapolated"] = calls
            row["wall_s_extrapolated"] = calls / max(rate, 1e-9)
        extrapolated.append(row)

    plan = []
    for f in floors:
        calls = (_calls_for(FROZEN_TARGET, f["floor_ratio"], fit)
                 if fit else None)
        # a call costs roughly linearly in n_sub (the quadrature is
        # the inner loop), so scale the measured rate accordingly
        scaled_rate = rate * N_SUB / f["n_sub"]
        plan.append({
            "n_sub": f["n_sub"],
            "floor_ratio": f["floor_ratio"],
            "frozen_target_reachable": calls is not None,
            "calls_extrapolated": calls,
            "wall_s_extrapolated": (calls / max(scaled_rate, 1e-9)
                                    if calls else None),
        })

    return {
        "measured_crossings": measured,
        "extrapolated_targets": extrapolated,
        "convergence_fit": (
            {"model": "ratio = floor + exp(intercept) * calls^slope",
             "log_log_slope": fit[0], "log_intercept": fit[1],
             "floor_ratio_used": floor,
             "note": ("fitted on the latter half of the trace's "
                      "REFINABLE part (ratio - floor); an "
                      "extrapolation aid, not a certified claim")}
            if fit else None),
        "quadrature_floor": {
            "rows": floors,
            "binding_lever": ("quadrature" if
                              floor > 0.5 * final["ratio"]
                              else "cell-refinement"),
            "note": ("DIAGNOSTIC. floor_ratio is the ratio that cell "
                     "refinement alone converges to at that n_sub, "
                     "because every interior cell keeps its center "
                     "flight-time's uncertainty. Estimated "
                     "standalone (support integral on a fixed grid "
                     "with the fast non-certified solver, widths "
                     "from the certified solver); never part of a "
                     "certified interval."),
        },
        "frozen_target_plan": {
            "target_ratio": FROZEN_TARGET,
            "rows": plan,
            "note": ("what PR-O3 must configure to reach the frozen "
                     "target; wall times scale the measured rate by "
                     "N_SUB/n_sub and are EXTRAPOLATIONS"),
            "reading": ("A larger n_sub does NOT pay here: it buys "
                        "a lower floor that was never binding, "
                        "while every call costs proportionally "
                        "more, so the extrapolated wall clock gets "
                        "worse with n_sub. The cost is cell "
                        "refinement, whose exponent degrades over "
                        "the run. NOT YET MEASURED, and the next "
                        "diagnostic PR-O3 should run: which cell "
                        "MODE carries the remaining width -- the "
                        "anchor neighbourhoods use the first-order "
                        "mode, whose per-cell width is O(h) against "
                        "the tangent mode's O(h^2), so they are the "
                        "suspect. `assemble` now reports "
                        "width_by_mode for exactly this."),
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
        "final": {"status": res["status"], "v_lo": float(v.lo),
                  "v_hi": float(v.hi), "ratio": res["ratio"],
                  "calls": res["calls"], "cells": res["cells"],
                  "modes": res["modes"],
                  "width_by_mode": res["width_by_mode"],
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
        if row["reachable_at_this_n_sub"]:
            print(f"  NOT MET {row['target_ratio']:g}: extrapolated "
                  f"{row['calls_extrapolated']:.0f} calls",
                  flush=True)
        else:
            print(f"  NOT MET {row['target_ratio']:g}: below the "
                  f"quadrature floor at n_sub={N_SUB}", flush=True)
    print("  binding lever: "
          f"{derived['quadrature_floor']['binding_lever']}",
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
