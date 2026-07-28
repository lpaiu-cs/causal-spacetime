"""P13 design check: choose per-rung constants so that m is held FIXED
while tau/ell varies, and six disjoint boxes still pack.

Design-check seed space 7000000+ (never an experimental window).
"""
import numpy as np

ELL = 1.0
K = 6
CAP = 200
M_TARGET = 76.0          # P12's top-rung interval count
BAND_REL = 0.10          # +/- 10% band around each rung's tau/ell


def tau_of(e1, x1, e2, x2):
    z = (e1**2 + e2**2 - (x1 - x2)**2) / (2 * e1 * e2)
    return ELL * np.arccosh(np.clip(z, 1.0, None))


def sprinkle(rho, eta_lo, xhalf, rng):
    """Poisson with intensity rho*Omega^2, Omega = ell/|eta|, over
    eta in [eta_lo, -1], |x| <= xhalf."""
    omega2_max = 1.0
    vol_coord = (abs(eta_lo) - 1.0) * 2 * xhalf
    n_prop = rng.poisson(rho * vol_coord * omega2_max)
    eta = rng.uniform(eta_lo, -1.0, n_prop)
    x = rng.uniform(-xhalf, xhalf, n_prop)
    keep = rng.random(n_prop) < (ELL**2 / eta**2) / omega2_max
    return eta[keep], x[keep]


def probe(tau_c, eta_lo, xhalf, rho, trials, seed0):
    rng = np.random.default_rng(seed0)
    lo, hi = tau_c * (1 - BAND_REL), tau_c * (1 + BAND_REL)
    done, ms, ns = 0, [], []
    for _ in range(trials):
        eta, x = sprinkle(rho, eta_lo, xhalf, rng)
        n = eta.size
        ns.append(n)
        if n < 10:
            continue
        u, v = eta + x, eta - x
        du = u[None, :] - u[:, None]
        dv = v[None, :] - v[:, None]
        rel = (du > 0) & (dv > 0)
        tau = tau_of(eta[:, None], x[:, None], eta[None, :], x[None, :])
        keep = rel & (tau >= lo) & (tau <= hi)
        a, b = np.nonzero(keep)
        if a.size == 0:
            continue
        e1, x1 = (u[b] + v[a]) / 2, (u[b] - v[a]) / 2
        e2, x2 = (u[a] + v[b]) / 2, (u[a] - v[b]) / 2
        ok = ((e1 >= eta_lo) & (e1 <= -1) & (np.abs(x1) <= xhalf)
              & (e2 >= eta_lo) & (e2 <= -1) & (np.abs(x2) <= xhalf))
        a, b = a[ok], b[ok]
        if a.size == 0:
            continue
        boxes, rej, chosen = [], 0, []
        for idx in rng.permutation(a.size):
            i, j = a[idx], b[idx]
            box = (u[i], u[j], v[i], v[j])
            if all((box[1] < o[0] or o[1] < box[0]
                    or box[3] < o[2] or o[3] < box[2]) for o in boxes):
                boxes.append(box)
                chosen.append((i, j))
                if len(chosen) == K:
                    break
            else:
                rej += 1
                if rej >= CAP:
                    break
        if len(chosen) < K:
            continue
        done += 1
        for i, j in chosen:
            inside = ((u > u[i]) & (u < u[j]) & (v > v[i]) & (v < v[j]))
            ms.append(int(inside.sum()))
    return done / trials, (float(np.mean(ms)) if ms else float('nan')), \
        float(np.mean(ns))


# candidate constants per rung: (tau/ell, eta_lo, xhalf, rho)
CANDIDATES = [
    (0.30, -1.8, 1.5, 2.0 * M_TARGET / 0.30**2),
    (0.60, -2.6, 2.0, 2.0 * M_TARGET / 0.60**2),
    (1.00, -4.0, 3.0, 2.0 * M_TARGET / 1.00**2),
    (1.50, -7.0, 4.5, 2.0 * M_TARGET / 1.50**2),
]
for tau_c, eta_lo, xhalf, rho in CANDIDATES:
    comp, m, n = probe(tau_c, eta_lo, xhalf, rho, 30, 7_000_000)
    print(f"tau/ell={tau_c:.2f} eta_lo={eta_lo:5.1f} X={xhalf:.1f} "
          f"rho={rho:7.1f} -> complete {comp:5.0%} | mean m {m:6.1f} "
          f"| mean N {n:7.0f}", flush=True)
