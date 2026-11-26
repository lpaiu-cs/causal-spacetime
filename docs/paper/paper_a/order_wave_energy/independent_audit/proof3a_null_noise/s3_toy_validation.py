# ruff: noqa: E501, UP031
"""Validation of the referee's 1+1 analogue (toy_common.py) and of the
probe-sampling term.

V1-V4  the counting code equals brute force; the functional is exactly zero
       on population moments and the reconstructed phase coordinate is the
       true u there.
V5     the tensor-grid probe quadrature agrees with exact cell integration.
V6     finite random probes: E|S_n - Sbar|^2 = E[int|wH|^2 - |Sbar|^2]/n_P,
       i.e. the probe term of (5) is 1/sqrt(n n_P) and nothing larger.

Deterministic; prints PASS/FAIL; exits nonzero on any FAIL.
"""

from __future__ import annotations

import math
import sys

import numpy as np
import toy_common as tc

FAILS = []


def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL") + " " + name + (" :: " + detail if detail else ""))
    if not ok:
        FAILS.append(name)


def brute_moments(u, v, xu, xv):
    """Raw moments at arbitrary probe points by direct O(n) counting per probe."""
    n = u.size
    past = (u[None, :] < xu[:, None]) & (v[None, :] < xv[:, None])
    fut = (u[None, :] > xu[:, None]) & (v[None, :] > xv[:, None])
    anc = (u > tc.UA) & (v > tc.VA)
    comp = (u[:, None] < u[None, :]) & (v[:, None] < v[None, :])
    dom = comp.sum(axis=0)  # dom[k] = number of references preceding k
    p = past.sum(axis=1) / n
    f = fut.sum(axis=1) / n
    b = (past & anc[None, :]).sum(axis=1) / n
    a = (past.astype(float) @ dom.astype(float)) / (0.5 * n * (n - 1))
    pa = np.count_nonzero((u < tc.UA) & (v < tc.VA)) / n
    fa = np.count_nonzero(anc) / n
    return p, f, a, b, pa, fa, dom


def v1_to_v4():
    rng = np.random.default_rng(101)
    u, v = rng.random(700), rng.random(700)
    brute = ((u[None, :] < u[:, None]) & (v[None, :] < v[:, None])).sum(axis=1)
    check(
        "V1 dominance counts equal brute force",
        np.array_equal(brute, tc.dominance_counts(u, v)),
    )

    ug, vg, hu, hv = tc.grid(400, 60)
    amp, uh = tc.calibrate(*tc.population_moments(ug, vg))
    check(
        "V2 population identity: w*H = 0 and u_hat = u on population moments",
        np.abs(amp).max() < 1e-14 and np.abs(uh - ug[:, None]).max() < 1e-12,
        f"max|wH|={np.abs(amp).max():.1e} max|u_hat-u|={np.abs(uh - ug[:, None]).max():.1e}",
    )

    n = 300
    u, v = rng.random(n), rng.random(n)
    ug, vg, hu, hv = tc.grid(50, 12)
    emp = tc.empirical_moments(u, v, ug, vg)
    uu, vv = np.meshgrid(ug, vg, indexing="ij")
    ref = brute_moments(u, v, uu.ravel(), vv.ravel())
    gap = max(float(np.abs(emp[k].ravel() - ref[k]).max()) for k in range(4))
    gap = max(gap, abs(emp[4] - ref[4]), abs(emp[5] - ref[5]))
    check(
        "V3 grid moments equal brute-force counts (incl. U-statistic)",
        gap < 1e-12,
        f"max gap {gap:.1e}",
    )

    omegas = np.array([1.0, 16.0, 64.0])
    exact, second = tc.exact_cell_scores(u, v, omegas)
    us, vs = np.sort(u), np.sort(v)
    eu = np.clip(np.r_[0.0, us, 1.0], tc.U_LO, tc.U_HI)
    ev = np.clip(np.r_[0.0, vs, 1.0], tc.V_LO, tc.V_HI)
    mid_u = 0.5 * (np.r_[0.0, us] + np.r_[us, 1.0])
    mid_v = 0.5 * (np.r_[0.0, vs] + np.r_[vs, 1.0])
    uu, vv = np.meshgrid(mid_u, mid_v, indexing="ij")
    p, f, a, b, pa, fa, _ = brute_moments(u, v, uu.ravel(), vv.ravel())
    amp, uh = tc.calibrate(p, f, a, b, pa, fa)
    area = (np.diff(eu)[:, None] * np.diff(ev)[None, :]).ravel()
    brute_score = np.array([np.sum(area * amp * np.exp(-1j * w * uh)) for w in omegas])
    gap = float(np.abs(exact - brute_score).max() / np.abs(brute_score).max())
    check(
        "V4 exact cell integration equals brute-force cell evaluation",
        gap < 1e-10,
        f"rel gap {gap:.1e}",
    )


def v5_grid_against_exact():
    rng = np.random.default_rng(202)
    omegas = np.array([1.0, 16.0, 64.0, 128.0])
    ug, vg, hu, hv = tc.grid(4200, 40)
    worst = 0.0
    for _ in range(4):
        u, v = rng.random(2000), rng.random(2000)
        exact, _ = tc.exact_cell_scores(u, v, omegas)
        amp, uh = tc.calibrate(*tc.empirical_moments(u, v, ug, vg))
        approx = tc.fourier_sums(amp, uh, omegas, hu * hv)
        worst = max(worst, float(np.max(np.abs(exact - approx) / np.abs(exact))))
    check(
        "V5 grid quadrature (4200x40) vs exact cells, n=2000, omega<=128",
        worst < 0.06,
        f"worst relative difference {worst:.4f} (tolerance 0.06)",
    )


def random_probe_score(rng, u, v, dom, n_probe, omegas):
    n = u.size
    x, y = rng.random(n_probe), rng.random(n_probe)
    inside = (x > tc.U_LO) & (x < tc.U_HI) & (y > tc.V_LO) & (y < tc.V_HI)
    xu, xv = x[inside], y[inside]
    past = (u[None, :] < xu[:, None]) & (v[None, :] < xv[:, None])
    fut = (u[None, :] > xu[:, None]) & (v[None, :] > xv[:, None])
    anc = (u > tc.UA) & (v > tc.VA)
    p = past.sum(axis=1) / n
    f = fut.sum(axis=1) / n
    b = (past & anc[None, :]).sum(axis=1) / n
    a = (past.astype(float) @ dom.astype(float)) / (0.5 * n * (n - 1))
    pa = np.count_nonzero((u < tc.UA) & (v < tc.VA)) / n
    fa = np.count_nonzero(anc) / n
    amp, uh = tc.calibrate(p, f, a, b, pa, fa)
    return (amp[None, :] * np.exp(-1j * omegas[:, None] * uh[None, :])).sum(
        axis=1
    ) / n_probe


def v6_probe_term():
    omegas = np.array([8.0, 64.0])
    ug, vg, hu, hv = tc.grid(4200, 40)
    rows = []
    for n, reps in ((1000, 300), (2000, 300), (4000, 200)):
        rng = np.random.default_rng(303 + n)
        diff2, pred = [], []
        for _ in range(reps):
            u, v = rng.random(n), rng.random(n)
            amp, uh = tc.calibrate(*tc.empirical_moments(u, v, ug, vg))
            sbar = tc.fourier_sums(amp, uh, omegas, hu * hv)
            second = hu * hv * float(np.sum(amp * amp))
            dom = tc.dominance_counts(u, v)
            sn = random_probe_score(rng, u, v, dom, n, omegas)
            diff2.append(np.abs(sn - sbar) ** 2)
            pred.append((second - np.abs(sbar) ** 2) / n)
        diff2, pred = np.array(diff2), np.array(pred)
        mean = diff2.mean(axis=0)
        se = diff2.std(axis=0, ddof=1) / math.sqrt(reps)
        ratio = mean / pred.mean(axis=0)
        rows.append((n, mean, se))
        check(
            "V6 n=n_P=%d: E|S_n-Sbar|^2 matches E[int|wH|^2-|Sbar|^2]/n_P" % n,
            bool(np.all(np.abs(ratio - 1) < 0.3)),
            f"ratio at omega=8,64: {ratio[0]:.3f} {ratio[1]:.3f} (rel SE {se[0] / mean[0]:.2f} {se[1] / mean[1]:.2f}); n*n_P*E|S_n-Sbar|^2 = {n * n * mean[0]:.3e} {n * n * mean[1]:.3e}",
        )
    # scaling in n with n_P = n: E|S_n - Sbar|^2 ~ n^-2
    logn = np.log([r[0] for r in rows])
    for k, omega in enumerate(omegas):
        logm = np.log([r[1][k] for r in rows])
        sig = np.array([r[2][k] / r[1][k] for r in rows])
        wts = 1 / sig**2
        xbar = np.sum(wts * logn) / np.sum(wts)
        slope = np.sum(wts * (logn - xbar) * logm) / np.sum(wts * (logn - xbar) ** 2)
        slope_se = math.sqrt(1 / np.sum(wts * (logn - xbar) ** 2))
        check(
            f"V6 probe term exponent at omega={omega:g}: E|S_n-Sbar|^2 ~ n^(-2) (n_P=n)",
            abs(slope + 2) < max(3 * slope_se, 0.15) and slope_se < 0.15,
            f"fitted exponent {slope:.3f} +- {slope_se:.3f} (L2 norm exponent {slope / 2:.3f} +- {slope_se / 2:.3f}, claim -1)",
        )


if __name__ == "__main__":
    v1_to_v4()
    v5_grid_against_exact()
    v6_probe_term()
    print("SUMMARY failures=%d" % len(FAILS))
    sys.exit(1 if FAILS else 0)
