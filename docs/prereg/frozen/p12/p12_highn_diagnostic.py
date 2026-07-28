"""Quarantined high-n diagnostic (design-check space 2000000+).

Adjudicates whether Stage A-12's non-monotone top rung is a curvature
plateau or n=12 fluctuation, adds an N=4800 rung, and reconciles the
flat-arm level using the RUNNER's own code path.
"""
import sys

import numpy as np

sys.path.insert(0, "experiments/positive_control")
from p12_curved import run_sample  # noqa: E402

NS = (600, 1200, 2400, 4800)
NSAMP = 200
BASE = 2_000_000

rows = {}
for k, n in enumerate(NS):
    ys, curved, flat, ms = [], [], [], []
    skips = 0
    for i in range(NSAMP):
        rec, ok = run_sample(n, BASE + 100_000 * k + 200 * i)
        if not ok:
            skips += 1
            continue
        ys.append(rec["y"])
        curved.append(rec["median_relerr_curved"])
        flat.append(rec["median_relerr_vs_flat_template"])
        ms.append(rec["mean_m_open"])
    ys = np.array(ys)
    rows[n] = dict(mean_y=ys.mean(), se=ys.std(ddof=1)/np.sqrt(ys.size),
                   curved=float(np.median(curved)), flat=float(np.median(flat)),
                   m=float(np.mean(ms)), n_ok=ys.size, skips=skips)
    r = rows[n]
    print(f"n={n:5d}: mean y {r['mean_y']:+.4f} ± {r['se']:.4f} | curved "
          f"{r['curved']:.4f} | flat {r['flat']:.4f} | m {r['m']:.1f} | "
          f"ok {r['n_ok']} skips {r['skips']}", flush=True)

logn = np.log10(np.array(NS, dtype=float))
means = np.array([rows[n]["mean_y"] for n in NS])
print(f"\nslope over all four rungs: {np.polyfit(logn, means, 1)[0]:+.4f}")
print(f"slope over 600-2400 only : {np.polyfit(logn[:3], means[:3], 1)[0]:+.4f}")
print(f"Delta 600->2400: {means[2]-means[0]:+.4f} | "
      f"600->4800: {means[3]-means[0]:+.4f}")
mono = all(means[i] > means[i+1] for i in range(3))
print("monotone decreasing across all four rungs:", mono)
