"""O2b mode-width diagnostic: WHICH cell mode carries the width.

The PR #65 review ruling: before any O3 budget or algorithm decision,
the bottleneck hypothesis ("the first-order anchor-neighbourhood
cells, O(h) per cell, carry the remaining width against the tangent
mode's O(h^2)") must be MEASURED, on the NEIGHBOR configuration
(12, 18, 8.0)M only -- running the frozen (12, 18, 8.5)M first and
then choosing budgets would collapse the freeze boundary.

This runner executes one adaptive assembly on the neighbor anchors
and records, at every checkpoint and at the end, the decomposition
the ruling asked for:

- raw_width_by_mode        (per-mode split of the raw directed sum)
- raw_total_before_intersection
- certified_total_after_intersection
- intersection_active      (has the L6d nesting ever clipped)
- per-mode cell counts

The split decomposes the RAW sum: when the intersection has clipped,
the certified total is tighter than the sum of its parts and each
mode's share is an upper bound on its part of the certified width --
the artifact stores both totals so that gap is visible.

Everything here is DIAGNOSTIC: it informs the O3 configuration, it
never enters a certified interval, and the verdict wording it emits
is capped at what a single neighbor-configuration run supports.

Run:  python experiments/oracle/o2_mode_width.py
"""

from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from volume_oracle import OracleConfig, assemble  # noqa: E402

NEIGHBOR = (12.0, 18.0, 8.0)   # NOT the frozen (12, 18, 8.5)
N_SUB = 16
MAX_CALLS = 60_000
MAX_WALL_S = 4_500.0           # 75 min: deep enough to see the
#                                late-run regime where the exponent
#                                degrades, without repeating the
#                                full 150-min ladder

_ARTIFACT = (Path(__file__).resolve().parents[2]
             / "docs" / "prereg" / "p14_oracle_mode_width.json")


def _shares(sample: dict) -> dict:
    total = sum(sample["raw_width_by_mode"].values())
    if total <= 0.0:
        return {}
    return {mode: w / total
            for mode, w in sample["raw_width_by_mode"].items()}


def main() -> None:
    cfg = OracleConfig(*NEIGHBOR, m=1.0, target_ratio=1e-9,
                       n_sub=N_SUB, max_calls=MAX_CALLS,
                       max_wall_s=MAX_WALL_S, max_depth=16,
                       init_rho=12, init_psi=12)
    last = [0.0]

    def show(s: dict) -> None:
        if s["wall_s"] - last[0] < 60.0:
            return
        last[0] = s["wall_s"]
        sh = _shares(s)
        print(f"  calls={s['calls']:6d} ratio={s['ratio']:.5f} "
              f"first-order={sh.get('first-order', 0.0):.3f} "
              f"tangent={sh.get('tangent', 0.0):.3f} "
              f"t={s['wall_s']:.0f}s", flush=True)

    print(f"O2b mode-width diagnostic on neighbor anchors "
          f"{NEIGHBOR} (frozen 8.5 stays unobserved)", flush=True)
    res = assemble(cfg, progress=show)
    v = res["v"]
    fo_share = _shares({"raw_width_by_mode":
                        res["raw_width_by_mode"]})
    print(f"final ratio={res['ratio']:.5f} "
          f"V=[{float(v.lo):.5f}, {float(v.hi):.5f}] "
          f"calls={res['calls']}", flush=True)
    print(f"raw width split: {res['raw_width_by_mode']}",
          flush=True)
    print(f"shares: {fo_share}", flush=True)
    print(f"intersection_active={res['intersection_active']} "
          f"raw={res['raw_total_before_intersection']:.5f} "
          f"certified={res['certified_total_after_intersection']:.5f}",
          flush=True)

    fo = fo_share.get("first-order", 0.0)
    if fo >= 0.5:
        verdict = (
            f"MEASURED on the neighbor configuration: the "
            f"first-order cells carry {fo:.0%} of the raw remaining "
            f"width at ratio {res['ratio']:.4f}. Improving the "
            f"anchor-neighbourhood treatment is the indicated next "
            f"lever for O3.")
    elif fo >= 0.2:
        verdict = (
            f"MEASURED on the neighbor configuration: the "
            f"first-order cells carry {fo:.0%} of the raw remaining "
            f"width at ratio {res['ratio']:.4f} -- a material but "
            f"not dominant share; an anchor-neighbourhood "
            f"improvement helps but does not remove the cost of "
            f"tangent-cell refinement.")
    else:
        verdict = (
            f"MEASURED on the neighbor configuration: the "
            f"first-order cells carry only {fo:.0%} of the raw "
            f"remaining width at ratio {res['ratio']:.4f}. The "
            f"candidate hypothesis is REFUTED at this depth: the "
            f"width lives in the tangent cells, and the O3 budget "
            f"should assume plain cell-refinement cost.")
    print(f"verdict: {verdict}", flush=True)

    art = {
        "scope": ("mode-width bottleneck diagnostic on the NEIGHBOR "
                  "anchors (12, 18, 8.0)M, per the PR #65 review "
                  "ruling: measured BEFORE any O3 budget or "
                  "algorithm decision, and the frozen (12, 18, "
                  "8.5)M volume stays unobserved. DIAGNOSTIC only "
                  "-- never part of a certified interval."),
        "neighbor_anchors": {"r_in": NEIGHBOR[0],
                             "r_out": NEIGHBOR[1],
                             "dt": NEIGHBOR[2], "m": 1.0},
        "config": {"n_sub": N_SUB, "k_micro": cfg.k_micro,
                   "d_switch": cfg.d_switch,
                   "max_calls": MAX_CALLS, "max_wall_s": MAX_WALL_S,
                   "max_depth": 16, "init_grid": [12, 12]},
        "host": {"machine": platform.machine(),
                 "python": platform.python_version()},
        "final": {
            "status": res["status"],
            # OUTWARD endpoints (PR #67 review R1)
            "v_lo": v.lo_float(), "v_hi": v.hi_float(),
            "ratio": res["ratio"], "calls": res["calls"],
            "cells": res["cells"],
            "cell_counts_by_mode": res["modes"],
            "raw_width_by_mode": res["raw_width_by_mode"],
            "raw_width_shares": fo_share,
            "raw_total_before_intersection":
                res["raw_total_before_intersection"],
            "certified_total_after_intersection":
                res["certified_total_after_intersection"],
            "intersection_active": res["intersection_active"],
            "wall_s": res["wall_s"],
        },
        "decomposition_note": (
            "raw_width_by_mode decomposes the RAW directed sum of "
            "live cells (x 2 pi). When intersection_active, the "
            "certified total is tighter than the sum of its parts "
            "and each mode's share is an UPPER BOUND on its part of "
            "the certified width; both totals are stored so the gap "
            "is visible."),
        "verdict": verdict,
        "curve": res["curve"],
    }
    with open(_ARTIFACT, "w", encoding="utf-8", newline="\n") as f:
        json.dump(art, f, ensure_ascii=False, indent=1)
        f.write("\n")
    print(f"artifact: {_ARTIFACT}", flush=True)


if __name__ == "__main__":
    main()
