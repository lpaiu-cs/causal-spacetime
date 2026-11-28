"""Referee check (e), part 2: the smoothing lemma with NON-uniform densities.

The proof's density is a product of k one-dimensional densities with
rho_i <= C/h and |rho_i'| <= C/h on an interval of length h.  On the unit cube
(s = t/h) this is rho_i <= C, |rho_i'| <= C h.  Here the harsher case
rho_i(s) = 1 + kappa_i (2 s - 1), |kappa_i| up to 0.9 (derivative O(1), not
O(h)) is tested.

Referee's explicit Lipschitz constant (referee report, finding on Section 4):
    rbar   = prod_i sup rho_i,   rbar'_j = sup|rho_j'| prod_{l != j} sup rho_l
    mixed entries    <= rbar
    repeated entries <= 2 rbar + rbar'_j        (cube faces + density derivative)
    L_rho = (m+1) rbar + max_j rbar'_j

Tests
  F1  k=2, SymPy closed form with linear densities: p is C^1 across tau=0,
      |p''| <= L_rho on a grid of (kappa_1, kappa_2, tau), p(1)=p'(1)=0.
  F2  k=3 evaluator (nested Gauss-Legendre between all breakpoints, exact for
      these piecewise polynomials) agrees with the exact rational Lasserre
      volume when kappa = 0.
  F3  k=3, chain / star / vee / triangle, all sign patterns, random kappa:
      second differences |p(tau+d)+p(tau-d)-2p(tau)| <= L_rho |d|^2 for
      |d| from 1e-1 to 1e-4 (no blow-up), including near the zero set.
  F4  same configurations: |sqrt p(tau') - sqrt p(tau)| <= sqrt(L_rho/2) |tau'-tau|,
      including pairs with p = 0 on one side.

Deterministic (fixed seed).  Run: python check_e_lemma_density.py
"""

from __future__ import annotations

import functools
import itertools
from fractions import Fraction as Fr

import numpy as np
import sympy as sp
from audit_common import Reporter, pattern_prob

rep = Reporter("check (e2): smoothing lemma with non-uniform product densities")
rng = np.random.default_rng(20261124)
GLX, GLW = np.polynomial.legendre.leggauss(8)


def L_rho(kap, m):
    kap = np.abs(np.asarray(kap, dtype=float))
    sup = 1.0 + kap
    rbar = float(np.prod(sup))
    rprime = max(2.0 * kap[j] * rbar / sup[j] for j in range(len(kap)))
    return (m + 1) * rbar + rprime


# ---------------------------------------------------------------- F1
tau, t1, t2, k1, k2 = sp.symbols("tau t1 t2 k1 k2", real=True)
r1 = 1 + k1 * (2 * t1 - 1)
r2 = 1 + k2 * (2 * t2 - 1)
p_pos = sp.expand(
    sp.integrate(r1 * sp.integrate(r2, (t2, t1 + tau, 1)), (t1, 0, 1 - tau))
)  # 0<=tau<=1
p_neg = sp.expand(
    1 - sp.integrate(r1 * sp.integrate(r2, (t2, 0, t1 + tau)), (t1, -tau, 1))
)  # -1<=tau<=0
c1_ok = (
    sp.simplify(p_pos.subs(tau, 0) - p_neg.subs(tau, 0)) == 0
    and sp.simplify(sp.diff(p_pos, tau).subs(tau, 0) - sp.diff(p_neg, tau).subs(tau, 0))
    == 0
)
edge_ok = (
    sp.simplify(p_pos.subs(tau, 1)) == 0
    and sp.simplify(sp.diff(p_pos, tau).subs(tau, 1)) == 0
    and sp.simplify(p_neg.subs(tau, -1) - 1) == 0
    and sp.simplify(sp.diff(p_neg, tau).subs(tau, -1)) == 0
)
d2pos = sp.lambdify((tau, k1, k2), sp.diff(p_pos, tau, 2), "numpy")
d2neg = sp.lambdify((tau, k1, k2), sp.diff(p_neg, tau, 2), "numpy")
kk = np.linspace(-0.95, 0.95, 39)
tt = np.linspace(0.0, 1.0, 201)
worst = 0.0
for a in kk:
    for b in kk:
        lr = L_rho([a, b], 1)
        worst = max(
            worst,
            float(np.max(np.abs(d2pos(tt, a, b))) / lr),
            float(np.max(np.abs(d2neg(-tt, a, b))) / lr),
        )
rep.check(
    "F1 k=2 linear densities: C^1 at tau=0, flat at tau=+-1, |p''| <= L_rho",
    bool(c1_ok and edge_ok and worst <= 1.0),
    f"max |p''|/L_rho = {worst:.4f} over 39x39x201 grid",
)


# ---------------------------------------------------------------- k=3 evaluator
def cdf(kap, t):
    t = np.clip(t, 0.0, 1.0)
    return t + kap * (t * t - t)


def p3(rows, kap):
    """P(pattern) for t in [0,1]^3 with densities 1 + kap_i (2 t_i - 1).

    rows: dict (i,j) -> (dirn, theta), meaning dirn * (t_j - t_i - theta) >= 0.
    Nested Gauss-Legendre between every breakpoint: exact up to rounding.
    """
    r01, r12, r02 = rows.get((0, 1)), rows.get((1, 2)), rows.get((0, 2))
    bps = [0.0, 1.0]
    if r01:
        bps += [-r01[1], 1 - r01[1]]
    if r02:
        bps += [-r02[1], 1 - r02[1]]
    if r01 and r12:
        bps += [-r01[1] - r12[1], 1 - r01[1] - r12[1]]
    if r02 and r12:
        bps += [r12[1] - r02[1], 1 + r12[1] - r02[1]]
    bps = np.unique(np.clip(np.array(bps), 0.0, 1.0))
    left, right = bps[:-1], bps[1:]
    t0 = (left[:, None] + (GLX[None, :] + 1) * 0.5 * (right - left)[:, None]).ravel()
    w0 = (GLW[None, :] * 0.5 * (right - left)[:, None]).ravel()
    # t1 domain
    lo1 = np.zeros_like(t0)
    hi1 = np.ones_like(t0)
    if r01:
        if r01[0] > 0:
            lo1 = np.maximum(lo1, t0 + r01[1])
        else:
            hi1 = np.minimum(hi1, t0 + r01[1])
    hi1 = np.maximum(hi1, lo1)
    cands = [lo1, hi1]
    if r12:
        cands += [np.full_like(t0, -r12[1]), np.full_like(t0, 1 - r12[1])]
        if r02:
            cands.append(t0 + r02[1] - r12[1])
    B = np.sort(np.stack([np.clip(c, lo1, hi1) for c in cands], axis=1), axis=1)
    bl, br = B[:, :-1], B[:, 1:]
    t1n = bl[:, :, None] + (GLX[None, None, :] + 1) * 0.5 * (br - bl)[:, :, None]
    w1n = GLW[None, None, :] * 0.5 * (br - bl)[:, :, None]
    t0b = t0[:, None, None]
    lo2 = np.zeros_like(t1n)
    hi2 = np.ones_like(t1n)
    if r02:
        if r02[0] > 0:
            lo2 = np.maximum(lo2, t0b + r02[1])
        else:
            hi2 = np.minimum(hi2, t0b + r02[1])
    if r12:
        if r12[0] > 0:
            lo2 = np.maximum(lo2, t1n + r12[1])
        else:
            hi2 = np.minimum(hi2, t1n + r12[1])
    inner2 = np.where(hi2 > lo2, cdf(kap[2], hi2) - cdf(kap[2], lo2), 0.0)
    inner2 = np.maximum(inner2, 0.0)
    rho1 = 1 + kap[1] * (2 * t1n - 1)
    inner1 = (rho1 * inner2 * w1n).sum(axis=(1, 2))
    rho0 = 1 + kap[0] * (2 * t0 - 1)
    return float((rho0 * inner1 * w0).sum())


def pattern_value(pairs, signs, kap, tau_signed):
    """p3 for signed thresholds; arguments bound explicitly per test case."""
    return p3(make_rows(pairs, signs, tau_signed), kap)


def make_rows(pairs, signs, tau_signed):
    """sigma (t_j - t_i) >= tau  <=>  dirn = sigma, theta = sigma * tau."""
    return {
        pr: (sg, sg * th) for pr, sg, th in zip(pairs, signs, tau_signed, strict=False)
    }


CONFIGS = [
    ("chain", [(0, 1), (1, 2)]),
    ("star", [(0, 1), (0, 2)]),
    ("vee", [(0, 2), (1, 2)]),
    ("triangle", [(0, 1), (1, 2), (0, 2)]),
]

# ---------------------------------------------------------------- F2
worst = 0.0
for _name, pairs in CONFIGS:
    m = len(pairs)
    for sg in itertools.product((1, -1), repeat=m):
        for _ in range(12):
            tq = [Fr(int(rng.integers(-80, 81)), 64) for _ in range(m)]
            exact = float(
                pattern_prob(
                    3, [(i, j, s) for (i, j), s in zip(pairs, sg, strict=False)], tq
                )
            )
            num = p3(make_rows(pairs, sg, [float(x) for x in tq]), [0.0, 0.0, 0.0])
            worst = max(worst, abs(exact - num))
rep.check(
    "F2 k=3 Gauss-Legendre evaluator equals exact Lasserre volume at kappa=0",
    worst < 5e-14,
    f"max abs diff {worst:.2e}",
)

# ---------------------------------------------------------------- F3, F4
for name, pairs in CONFIGS:
    m = len(pairs)
    max_q = 0.0
    max_s = 0.0
    ok3 = ok4 = True
    nzero = 0
    ntest = 0
    for sg in itertools.product((1, -1), repeat=m):
        for _trial in range(5):
            kap = rng.uniform(-0.9, 0.9, 3)
            lr = L_rho(kap, m)
            f = functools.partial(pattern_value, pairs, sg, kap)
            points = [rng.uniform(-1.2, 1.2, m) for _ in range(3)]
            # near the zero set along the diagonal
            base = None
            for _ in range(30):
                cand = rng.uniform(-0.7, 0.7, m)
                if f(cand) > 1e-6:
                    base = cand
                    break
            if base is not None:
                lo, hi = 0.0, 4.0
                for _ in range(40):
                    mid = 0.5 * (lo + hi)
                    if f(base + mid) > 0.0:
                        lo = mid
                    else:
                        hi = mid
                for back in (1e-1, 1e-2, 1e-3):
                    if hi - back > 0:
                        points.append(base + hi - back)
                # straddling pair
                pin, pout = f(base + hi - 2e-3), f(base + hi + 2e-3)
                dist = 4e-3 * np.sqrt(m)
                if pout == 0.0 and pin > 0.0:
                    nzero += 1
                if abs(np.sqrt(pin) - np.sqrt(pout)) > np.sqrt(lr / 2) * dist + 1e-7:
                    ok4 = False
            for tv in points:
                p0 = f(tv)
                for step in (1e-1, 1e-2, 1e-3, 1e-4):
                    d = rng.normal(size=m)
                    d *= step / np.linalg.norm(d)
                    pp, pm = f(tv + d), f(tv - d)
                    quo = abs(pp + pm - 2 * p0) / step**2
                    tol = 4e-15 / step**2
                    max_q = max(max_q, (quo - tol) / lr)
                    if quo > lr + tol:
                        ok3 = False
                    slope = abs(np.sqrt(pp) - np.sqrt(p0)) / step
                    max_s = max(max_s, (slope - 1e-7 / step) / np.sqrt(lr / 2))
                    if slope > np.sqrt(lr / 2) + 1e-7 / step:
                        ok4 = False
                    ntest += 1
    rep.check(
        "F3 second differences bounded by L_rho |d|^2, k=3 " + name,
        ok3,
        f"max quotient/L_rho = {max_q:.4f} over {ntest} tests",
    )
    rep.check(
        "F4 sqrt(p) Lipschitz with sqrt(L_rho/2), k=3 " + name,
        ok4,
        f"max slope/sqrt(L_rho/2) = {max_s:.4f}; zero-set crossings {nzero}",
    )

rep.finish()
