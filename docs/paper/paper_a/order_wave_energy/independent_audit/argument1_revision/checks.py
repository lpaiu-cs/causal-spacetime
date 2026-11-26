"""Deterministic checks for quantized_order_bound_r2.md (argument 1 written out).

No sampling of orders. Seeded draws are used only to test inequalities.
Run: python checks.py
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np

failed = 0


def check(name, ok, detail=""):
    global failed
    print(("PASS " if ok else "FAIL ") + name + (" | " + detail if detail else ""))
    failed += 0 if ok else 1


def f_max(q):
    return 1.0 / math.cos(math.sqrt(q) / 2) ** 2


def c_e(q):
    return 2 * math.sqrt(q) * (4 * f_max(q) - 2)


# 1. The constants C_E and c0 = 2 sqrt 2 C_E (referee: 6.39, 18.1, 132, 13532).
for q, want in ((0.25, 6.39), (1.0, 18.1), (4.0, 132.4), (9.0, 13532)):
    c0 = 2 * math.sqrt(2) * c_e(q)
    check(f"c0 at q={q}", abs(c0 - want) / want < 2e-3, f"c0 = {c0:.4g}")

# 2. Lemma 2.1: e <= C_E lam min(l, lam) and e < b, by composite Gauss quadrature.
gx, gw = np.polynomial.legendre.leggauss(12)


def oscillatory(u1, u2, lam, q):
    """I = int_u1^u2 sec^2(sqrt q u) e^{iu/lam} du, panels no wider than lam/4."""
    panels = max(1, math.ceil((u2 - u1) / (lam / 4)))
    edges = np.linspace(u1, u2, panels + 1)
    mid, half = (edges[1:] + edges[:-1]) / 2, (edges[1:] - edges[:-1]) / 2
    u = (mid[:, None] + half[:, None] * gx[None, :]).ravel()
    w = (half[:, None] * gw[None, :]).ravel()
    return np.sum(w * np.exp(1j * u / lam) / np.cos(math.sqrt(q) * u) ** 2)


def b_exact(u1, u2, q):
    r = math.sqrt(q)
    return (math.tan(r * u2) - math.tan(r * u1)) / r


worst_ratio, worst_b = 0.0, 0.0
ends = np.linspace(-0.5, 0.5, 17)
for q in (0.25, 1.0, 4.0, 9.0):
    for lam in (1 / 8, 1 / 32, 1 / 128):
        s = 2 * lam * math.sqrt(q)
        t = s / math.sqrt(1 + s * s)
        for i, u1 in enumerate(ends):
            for u2 in list(ends[i + 1 :]) + [u1 + 1e-3, u1 + 1e-5]:
                if u2 > 0.5:
                    continue
                e = t * abs(oscillatory(u1, u2, lam, q))
                bound = c_e(q) * lam * min(u2 - u1, lam)
                worst_ratio = max(worst_ratio, e / bound)
                worst_b = max(worst_b, e / b_exact(u1, u2, q))
check(
    "Lemma 2.1: e <= C_E lam min(l, lam)",
    worst_ratio <= 1.0,
    f"largest e / bound = {worst_ratio:.3f}",
)
check("Lemma 2.1: e < b", worst_b < 1.0, f"largest e / b = {worst_b:.3f}")

# 3. Lemma 2.2 and Lemma 2.4 on a deterministic grid near the null cone.
#    Tested with the smallest constant the proof needs for that pair,
#    c = 2 sqrt 2 e / delta, which is at most c0 by Lemma 2.1.
violations = {"first": 0, "second": 0, "time": 0}
tested = 0
for q in (0.25, 1.0, 4.0):
    for lam in (1 / 8, 1 / 32):
        delta = lam * lam
        s = 2 * lam * math.sqrt(q)
        t = s / math.sqrt(1 + s * s)
        for u1, u2 in ((-0.5, 0.5), (-0.3, -0.29), (0.1, 0.1 + 1e-4), (-0.2, 0.3)):
            b = b_exact(u1, u2, q)
            e = t * abs(oscillatory(u1, u2, lam, q))
            c_pair = 2 * math.sqrt(2) * e / delta
            for psi in np.linspace(0, math.pi, 7):
                cs, sn = math.cos(psi), math.sin(psi)
                rot = np.array([[cs, -sn], [sn, cs]])
                emat = e * rot @ np.diag([1.0, -1.0]) @ rot.T
                inv = np.linalg.inv(b * np.eye(2) + emat)
                for ang in np.linspace(0, 2 * math.pi, 13)[:-1]:
                    for r in (1e-4, 0.03, 0.3, 0.9):
                        dy = r * np.array([math.cos(ang), math.sin(ang)])
                        dv0 = dy @ dy / (2 * b)
                        for shift in np.linspace(-3, 3, 25):
                            dv = dv0 + shift * e * max(dv0, 1e-12) / b
                            rel = 2 * dv >= dy @ inv @ dy
                            dt = (b + dv) / math.sqrt(2)
                            dz = math.sqrt(((b - dv) / math.sqrt(2)) ** 2 + dy @ dy)
                            g = dt - dz
                            tested += 1
                            if g >= c_pair * delta and not rel:
                                violations["first"] += 1
                            if rel and g < -c_pair * delta * (1 + 1e-12):
                                violations["second"] += 1
                            if rel and dt < 0:
                                violations["time"] += 1
check(
    "Lemma 2.2: both inclusions with c = 2 sqrt 2 e / delta",
    violations["first"] == 0 and violations["second"] == 0,
    f"{tested} pairs, violations {violations['first']}, {violations['second']}",
)

#    Lemma 2.4 needs Delta v of both signs (the grid above stays near the
#    null cone, where Delta v > 0, so the check there cannot fail).
#    Added after review (finding M9).
related = 0
for q in (0.25, 1.0, 9.0):
    lam = 1 / 8
    s = 2 * lam * math.sqrt(q)
    t = s / math.sqrt(1 + s * s)
    for u1, u2 in ((-0.5, 0.5), (0.1, 0.1 + 1e-4), (-0.45, -0.4)):
        b = b_exact(u1, u2, q)
        e = t * abs(oscillatory(u1, u2, lam, q))
        for psi in np.linspace(0, math.pi, 5):
            cs, sn = math.cos(psi), math.sin(psi)
            rot = np.array([[cs, -sn], [sn, cs]])
            inv = np.linalg.inv(b * np.eye(2) + e * rot @ np.diag([1.0, -1.0]) @ rot.T)
            for dv in np.linspace(-1.0, 1.0, 81):
                for r in (0.0, 1e-3, 0.1, 1.0):
                    for ang in (0.0, 1.0, 2.5):
                        dy = r * np.array([math.cos(ang), math.sin(ang)])
                        if 2 * dv >= dy @ inv @ dy:
                            related += 1
                            violations["time"] += (b + dv) / math.sqrt(2) < 0 or dv < 0
check(
    "Lemma 2.4: a wave relation has Delta v >= 0 and Delta T >= 0",
    violations["time"] == 0 and related > 0,
    f"{related} related pairs among Delta v in [-1, 1]",
)

# 4. Lemmas 4.1 and 4.2 on the cyclic triangle (uniform densities, m = 3,
#    Lambda = m + 1 = 4). At tau = -eps (1,1,1), p = 9/2 eps^2 - 4 eps^3 and,
#    by the cyclic symmetry, each partial derivative is -(3 eps - 4 eps^2).
ratios = []
for eps in (1e-6, 1e-3, 0.01, 0.1, 0.4):
    p = 4.5 * eps**2 - 4 * eps**3
    grad2 = 3 * (3 * eps - 4 * eps**2) ** 2
    ratios.append(grad2 / (2 * p))
check(
    "Lemma 4.2: |grad p|^2 / (2p) <= Lambda = 4 on the cyclic triangle",
    max(ratios) <= 4.0,
    f"largest ratio {max(ratios):.4f} (tends to 3)",
)
#    The frozen control value H^2/eps^2 = 13.44677 at eps = 0.001 (h = 1,
#    c0 delta = eps) against Lemma 4.3's constant 2^(m+1) m Lambda.
check(
    "Lemma 4.3 constant covers the frozen cyclic control",
    13.44677 <= 2**4 * 3 * 4,
    f"13.44677 against {2**4 * 3 * 4}",
)

#    Step (ii) for k = 2 with non-uniform densities: G(tau) = int rho1(s)
#    rho2(s + tau) 1{0 <= s + tau <= 1} ds must be Lipschitz with constant
#    2 R1 R2 + R2' R1 (integrating s2 first).
nodes, weights = np.polynomial.legendre.leggauss(64)


def g_slice(tau, beta1, beta2):
    lo, hi = max(0.0, -tau), min(1.0, 1.0 - tau)
    if hi <= lo:
        return 0.0
    s = lo + (hi - lo) * (nodes + 1) / 2
    w = (hi - lo) * weights / 2
    rho1 = 1 + beta1 * (s - 0.5)
    rho2 = 1 + beta2 * (s + tau - 0.5)
    return float(np.sum(w * rho1 * rho2))


worst = 0.0
for beta1, beta2 in ((0.0, 0.0), (0.8, -0.6), (-1.5, 1.9), (1.9, 1.9)):
    r1, r2 = 1 + abs(beta1) / 2, 1 + abs(beta2) / 2
    allowed = 2 * r1 * r2 + abs(beta2) * r1
    taus = np.linspace(-1.2, 1.2, 4801)
    vals = np.array([g_slice(x, beta1, beta2) for x in taus])
    slope = np.max(np.abs(np.diff(vals) / np.diff(taus)))
    worst = max(worst, slope / allowed)
check(
    "Lemma 4.1 (ii): slice Lipschitz constant, k = 2",
    worst <= 1.0,
    f"largest ratio {worst:.3f}",
)

# 5. Hellinger convention: TV <= H and H^2 subadditive over products.
rng = np.random.default_rng(7)
bad_tv, bad_prod = 0, 0
for _ in range(2000):
    n = int(rng.integers(2, 9))
    p, q = rng.dirichlet(np.ones(n)), rng.dirichlet(np.ones(n))
    if rng.random() < 0.3:
        p[0], q[1] = 0.0, 0.0
        p, q = p / p.sum(), q / q.sum()
    h2 = np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)
    bad_tv += 0.5 * np.sum(np.abs(p - q)) > math.sqrt(h2) + 1e-15
    p2, q2 = rng.dirichlet(np.ones(3)), rng.dirichlet(np.ones(3))
    h2b = np.sum((np.sqrt(p2) - np.sqrt(q2)) ** 2)
    joint = np.sum((np.sqrt(np.outer(p, p2)) - np.sqrt(np.outer(q, q2))) ** 2)
    bad_prod += joint > h2 + h2b + 1e-15
check("TV <= H with H^2 = sum (sqrt p - sqrt q)^2", bad_tv == 0, "2000 pairs")
check("H^2 of a product <= sum of H^2", bad_prod == 0, "2000 pairs")

# 6. Section 5 algebra, in logarithms (powers like delta^(2K) underflow).
worst_id = 0.0
for k in (1, 2, 3, 7, 14):
    for ln_n in (0.0, math.log(37.0), math.log(1e4), math.log(1e9)):
        for ln_d in (math.log(1e-2), math.log(1e-6), math.log(1e-12)):
            ln_h = 2 * (ln_d - k * ln_n) / (2 * k + 1)
            target = ((3 * k + 1) * ln_n + 2 * k * ln_d) / (2 * k + 1)
            for value in ((k + 1) * ln_n + k * ln_h, ln_n + ln_d - ln_h / 2):
                worst_id = max(worst_id, abs(value - target))
            worst_id = max(worst_id, abs(ln_n + ln_h - (ln_n + 2 * ln_d) / (2 * k + 1)))
check("choice of h: the three identities", worst_id < 1e-9, f"{worst_id:.1e} in logs")

#    If h < c1 delta then the third term of (8) is at least c1^(-(3K+1)/(2K)).
worst_low = math.inf
for k in (1, 2, 5, 14):
    for c1 in (1.0, 5.0, 14.2, 1e4):
        ln_c1 = math.log(c1)
        for ln_d in np.linspace(math.log(1e-12), 0.0, 25):
            for ln_n in np.linspace(0.0, math.log(1e30), 61):
                ln_h = 2 * (ln_d - k * ln_n) / (2 * k + 1)
                if ln_h < ln_c1 + ln_d:
                    term = ((3 * k + 1) * ln_n + 2 * k * ln_d) / (2 * k + 1)
                    worst_low = min(worst_low, term + (3 * k + 1) / (2 * k) * ln_c1)
check(
    "h below h_min forces a large third term",
    worst_low >= -1e-12,
    f"smallest log ratio {worst_low:.3f}",
)

#    Coverage of alpha < 4/3 and the alpha = 1.3 example, in exact rationals
#    (in floating point, alpha/(4 - 3 alpha) at an exact integer can round
#    down and make floor(x) + 1 equal x).
cover_ok = True
for j in range(1, 4000):
    alpha = Fraction(j, 3000)
    k = math.floor(alpha / (4 - 3 * alpha)) + 1
    cover_ok &= 4 * k - (3 * k + 1) * alpha > 0
a13 = Fraction(13, 10)
k13 = math.floor(a13 / (4 - 3 * a13)) + 1
margin = (4 * k13 - (3 * k13 + 1) * a13) / (2 * k13 + 1)
check("every alpha = j/3000 < 4/3 is covered", bool(cover_ok))
check(
    "alpha = 1.3 needs K = 14, margin <= 0.0035",
    k13 == 14 and margin <= Fraction(35, 10000),
    f"K={k13}, margin {float(margin):.5f}",
)

print(f"SUMMARY failed={failed}")
raise SystemExit(1 if failed else 0)
