"""Isolated diagnostic: what makes the FLAT twin's contrast nonzero?

In exact flat space the chain statistic is boost-invariant, so given
the interval count m it should not know the box's aspect ratio or
absolute size. If conditioning on m removes the rung difference, the
artifact is an m-distribution effect; if it survives, it is not.

Design-check space 7500000+ (never an experimental window).
"""
import sys

import numpy as np

sys.path.insert(0, "experiments/positive_control")
import p13_tau_ell as P  # noqa: E402

from causal_spacetime_lab.estimators import (  # noqa: E402
    estimate_tau_from_longest_chain_1p1,
)

NSAMP = 120
rows = {}
for tau_c in P.RUNGS:
    eta_lo, xhalf, _ = P.PATCH[tau_c]
    centre = P.TWIN_BAND_CENTRE[tau_c]
    band = (centre * 0.9, centre * 1.1)
    rho_twin = P.TWIN_RHO[tau_c]
    recs = []
    rng = np.random.default_rng(7_500_000 + int(tau_c * 100))
    for _ in range(NSAMP):
        eta, x = P.sprinkle(eta_lo, xhalf, rho_twin, rng, flat=True)
        u, v = eta + x, eta - x
        pool = P.eligible_pairs(eta, x, eta_lo, xhalf, band, flat=True)
        if pool.shape[0] == 0:
            continue
        pairs, ok, _ = P.draw_disjoint(pool, u, v, rng)
        if not ok:
            continue
        rho_hat = eta.size / ((abs(eta_lo) - 1.0) * 2.0 * xhalf)
        for i, j in pairs:
            inside = (u > u[i]) & (u < u[j]) & (v > v[i]) & (v < v[j])
            idxs = np.concatenate(([i], np.flatnonzero(inside), [j]))
            su, sv = u[idxs] / 2.0, v[idxs] / 2.0
            from causal_spacetime_lab.chains import longest_chain_length
            L = longest_chain_length(
                (su[:, None] < su[None, :]) & (sv[:, None] < sv[None, :]),
                start=0, end=idxs.size - 1, event_times=su)
            th = float(estimate_tau_from_longest_chain_1p1(L, rho=rho_hat))
            tt = float(np.sqrt((u[j] - u[i]) * (v[j] - v[i])))
            m_open = int(idxs.size - 2)
            aspect = (u[j] - u[i]) / (v[j] - v[i])
            recs.append((m_open, abs(th - tt) / tt, L, aspect, rho_hat))
    arr = np.array(recs)
    rows[tau_c] = arr
    m, e, L, asp = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    print(f"tau/ell={tau_c}: pairs {arr.shape[0]:5d} | m {m.mean():6.1f}"
          f" +/- {m.std():5.1f} | relerr {np.median(e):.4f} | L {L.mean():5.1f}"
          f" | aspect {np.median(asp):8.2f} [{np.percentile(asp,5):7.2f},"
          f" {np.percentile(asp,95):8.2f}]", flush=True)

print("\nCONDITIONED on m in [70, 82] -- flat space should be blind to"
      " everything else given m:")
for tau_c in P.RUNGS:
    a = rows[tau_c]
    sel = (a[:, 0] >= 70) & (a[:, 0] <= 82)
    e = a[sel, 1]
    print(f"  tau/ell={tau_c}: n {sel.sum():5d} | median relerr"
          f" {np.median(e):.4f} | mean L {a[sel,2].mean():.2f}")

print("\nmean of log10(relerr) per pair, all m vs conditioned:")
for tau_c in P.RUNGS:
    a = rows[tau_c]
    sel = (a[:, 0] >= 70) & (a[:, 0] <= 82)
    print(f"  tau/ell={tau_c}: all {np.mean(np.log10(a[:,1])):+.4f}"
          f" | cond {np.mean(np.log10(a[sel,1])):+.4f}")
