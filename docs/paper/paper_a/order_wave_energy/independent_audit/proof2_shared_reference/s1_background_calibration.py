# ruff: noqa: E501, UP031
"""Referee check of Section 1 of local_null_challenge.md and of the calibration
V(p,f,c;q) defined in null_coordinate_challenge.md Section 2.

Symbolic (SymPy):
  S1  R_uu = 0 for the circular family with sinh(2 beta) = 2 lambda sqrt(q); 2q for the background.
  S2  C_2, C_3 closed forms of the background kernel sin(k l)/k; alpha = (1 - k l cot k l)/(2q).
  S3  transverse integrals behind Vol = pi eta^2 C_2/K and P_2 = K C_3/(3 C_2^2).
  S4  V(p,f,c;q) = 1 - t + t eta - (1-t) zeta equals v identically.
Numeric:
  N1  p0(0) = 1/10, d log p0/dz <= -1/50 on [0,1] (dense grid), alpha increasing.
  N2  F(population moments) = v to rounding at random interior events (q = 1, 0.25).
Monte Carlo (independent of the volume/chain hierarchy; exact pair relation only):
  M1  background past/future volume fractions of an off-axis event (q=1).
  M2  a_X = p_X^2 c_X with denominator choose(M,2), comparable = either order.
  M3  background pair probability at q=6 equals p0(6) (clearly different from 1/10).

Deterministic (fixed seeds).  Prints PASS/FAIL lines; exits nonzero on any FAIL.
Run:  python s1_background_calibration.py     (system Python: numpy + sympy + mpmath)
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_geometry import (  # noqa: E402
    A_TIP,
    Wave,
    alpha0,
    c_series,
    calibration_f,
    p0,
    t0,
)

FAILS = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" | " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


def symbolic():
    # ---- S1: Ricci component of the Rosen metric, evaluated numerically from exact derivatives
    u = sp.symbols("u", real=True)
    k, lam, phi, beta = sp.symbols("k lambda phi beta", positive=True)
    theta = (u / lam + phi) / 2
    rot = sp.Matrix([[sp.cos(theta), -sp.sin(theta)], [sp.sin(theta), sp.cos(theta)]])

    def ruu(gamma):
        inv = gamma.inv()
        first, second = gamma.diff(u), gamma.diff(u, 2)
        return (
            -sp.Rational(1, 2) * (inv * second).trace()
            + sp.Rational(1, 4) * (inv * first * inv * first).trace()
        )

    wave = (
        sp.cos(k * u) ** 2 * rot * sp.diag(sp.exp(2 * beta), sp.exp(-2 * beta)) * rot.T
    )
    flat = sp.cos(k * u) ** 2 * sp.eye(2)
    f_wave = sp.lambdify((u, k, lam, phi, beta), ruu(wave), "mpmath")
    f_flat = sp.lambdify((u, k), ruu(flat), "mpmath")
    mp.mp.dps = 30
    rng = np.random.default_rng(5)
    worst_w, worst_b = 0, 0
    for _ in range(40):
        uu, kk = rng.uniform(-0.5, 0.5), rng.uniform(0.2, 1.0)
        ll, pp = rng.uniform(0.001, 0.7), rng.uniform(0, 6.28)
        bb = mp.asinh(2 * ll * kk) / 2
        worst_w = max(
            worst_w, abs(f_wave(mp.mpf(uu), mp.mpf(kk), mp.mpf(ll), mp.mpf(pp), bb))
        )
        worst_b = max(worst_b, abs(f_flat(mp.mpf(uu), mp.mpf(kk)) - 2 * kk * kk))
    check(
        "S1 exact vacuum: R_uu=0 for sinh(2beta)=2 lambda sqrt(q); background R_uu=2q",
        worst_w < 1e-15 and worst_b < 1e-15,
        f"max|R_uu| wave {mp.nstr(worst_w, 3)}, max|R_uu-2q| background {mp.nstr(worst_b, 3)}",
    )

    # ---- S2: chain convolutions of K(l) = sin(k l)/k
    w, ell = sp.symbols("w l", positive=True)
    kern = lambda x: sp.sin(k * x) / k  # noqa: E731
    c2 = sp.integrate(kern(w) * kern(ell - w), (w, 0, ell))
    c2_claim = (sp.sin(k * ell) - k * ell * sp.cos(k * ell)) / (2 * k**3)
    ok2 = sp.simplify(c2 - c2_claim) == 0
    alpha = sp.simplify(c2_claim / kern(ell))
    ok_alpha = sp.simplify(alpha - (1 - k * ell * sp.cot(k * ell)) / (2 * k**2)) == 0
    c3 = sp.integrate(c2_claim.subs(ell, w) * kern(1 - w), (w, 0, 1))
    c3_claim = ((3 - k**2) * sp.sin(k) - 3 * k * sp.cos(k)) / (8 * k**5)
    ok3 = sp.simplify(c3 - c3_claim) == 0
    series_ok = True
    for z in (0.01, 0.3, 1.0, 4.0, 6.0):
        kk = math.sqrt(z)
        vals = [
            math.sin(kk) / kk,
            float(c2_claim.subs({k: kk, ell: 1})),
            float(c3_claim.subs({k: kk})),
        ]
        for r, val in zip((1, 2, 3), vals, strict=False):
            series_ok &= abs(c_series(r, z) - val) < 1e-11
    check(
        "S2 C_2, C_3, alpha closed forms and the 0F1 series agree",
        ok2 and ok3 and ok_alpha and series_ok,
        f"C2 {ok2}, C3 {ok3}, alpha {ok_alpha}, series {series_ok}",
    )

    # ---- S3: transverse integrals (after completing the square, w = Q^(1/2) x)
    eta, r = sp.symbols("eta r", positive=True)
    one = sp.integrate((eta - r**2 / 2) * 2 * sp.pi * r, (r, 0, sp.sqrt(2 * eta)))
    two = sp.integrate(
        (eta - r**2 / 2) ** 3 / 3 * 2 * sp.pi * r, (r, 0, sp.sqrt(2 * eta))
    )
    kk, cc2, cc3 = sp.symbols("K C2 C3", positive=True)
    # ordered 2-chain measure pi*(pi eta^4/6) C3/K; volume pi eta^2 C2/K; either order: factor 2
    p2 = sp.simplify(
        2 * (sp.pi * two.subs(eta, eta) * cc3 / kk) / (one * cc2 / kk) ** 2
    )
    check(
        "S3 Vol = pi eta^2 C2/K and P_2(comparable) = K C3/(3 C2^2)",
        sp.simplify(one - sp.pi * eta**2) == 0
        and sp.simplify(two - sp.pi * eta**4 / 6) == 0
        and sp.simplify(p2 - kk * cc3 / (3 * cc2**2)) == 0,
        f"P_2 = {p2}",
    )

    # ---- S4: the calibration identity
    m, big, v, rad = sp.symbols("m M v rho", positive=True)
    eta_b = v - rad**2 / (2 * m)
    zeta_b = 1 - v - rad**2 / (2 * (big - m))
    t = m / big
    check(
        "S4 1 - t + t eta - (1-t) zeta = v in the background (any transverse radius)",
        sp.simplify(1 - t + t * eta_b - (1 - t) * zeta_b - v) == 0,
    )
    # t(1-s) = 1 - t(s): the form F = 1 - t(s) + A(s) sqrt p - A(1-s) sqrt f of the proof
    s = sp.symbols("s", positive=True)

    def tt(x):
        return (sp.tan(k * (x - sp.Rational(1, 2))) + sp.tan(k / 2)) / (
            2 * sp.tan(k / 2)
        )  # noqa: E731

    check("S5 t(1-s) = 1 - t(s)", sp.simplify(tt(1 - s) + tt(s) - 1) == 0)


def numeric():
    check("N1a p0(0) = 1/10", abs(p0(0.0) - 0.1) < 1e-15, f"p0(0)={p0(0.0):.16f}")
    mp.mp.dps = 30

    def c(r, z):
        return mp.hyp0f1(r + mp.mpf(1) / 2, -z / 4) / mp.factorial(2 * r - 1)

    def logp0(z):
        return mp.log(c(1, z) * c(3, z) / (3 * c(2, z) ** 2))

    worst, worst_formula = -1e9, 0.0
    for j in range(0, 2001):
        z = mp.mpf(j) / 2000
        formula = -c(2, z) / c(1, z) - 3 * c(4, z) / c(3, z) + 4 * c(3, z) / c(2, z)
        if 0 < j < 2000:
            numeric_derivative = mp.diff(logp0, z)
            worst_formula = max(worst_formula, abs(numeric_derivative - formula))
        worst = max(worst, float(formula))
    check(
        "N1b d log p0/dz <= -1/50 on [0,1] and the derivative formula is right",
        worst < -1 / 50 and worst_formula < 1e-18,
        f"max over grid {worst:.8f} (author certifies -0.03555); formula error {mp.nstr(worst_formula, 3)}",
    )
    # informational: how far beyond q=1 the inverse calibration would remain monotone
    last = None
    for j in range(1, 960):
        z = mp.mpf(j) / 100
        formula = -c(2, z) / c(1, z) - 3 * c(4, z) / c(3, z) + 4 * c(3, z) / c(2, z)
        if formula >= 0:
            last = float(z)
            break
    print(
        "INFO p0 stays strictly decreasing on (0, 9.59]"
        if last is None
        else f"INFO p0 stops decreasing near z={last}"
    )
    lengths = np.linspace(0.01, 1.0, 400)
    inc = all(np.all(np.diff([alpha0(x, q) for x in lengths]) > 0) for q in (0.25, 1.0))
    check("N1c alpha(l) strictly increasing on (0,1] for q = 0.25, 1", inc)

    rng = np.random.default_rng(11)
    worst = 0.0
    for q in (1.0, 0.25):
        k = math.sqrt(q)
        big = 2 * math.tan(k / 2) / k
        gamma = alpha0(1.0, q)
        for _ in range(300):
            s = rng.uniform(0.05, 0.95)
            m = (math.tan(k * (s - 0.5)) + math.tan(k / 2)) / k
            v = rng.uniform(0.1, 0.9)
            rad2 = rng.uniform(0, 0.9) * 2 * min(v * m, (1 - v) * (big - m))
            eta, zeta = v - rad2 / (2 * m), 1 - v - rad2 / (2 * (big - m))
            p = alpha0(s, q) * eta**2 / gamma
            f = alpha0(1 - s, q) * zeta**2 / gamma
            a = p * p * p0(q * s * s)
            worst = max(worst, abs(calibration_f(p, f, a, q) - v))
            assert abs(t0(s, q) - m / big) < 1e-14
    check(
        "N2 F(p,f,a) = V(p,f,a/p^2;q) returns v at 600 random background events",
        worst < 1e-9,
        f"max |F - v| = {worst:.3e}",
    )


def monte_carlo():
    # M1/M2: background, q=1
    q = 1.0
    bg = Wave(q, 0.25, background=True)
    rng = np.random.default_rng(424242)
    n = 2000000
    pts = bg.sample(n, rng)
    event = np.array([0.15, 0.52, 0.20, -0.10])
    eta, zeta = bg.slacks(event)
    s = event[0] - A_TIP
    gamma = alpha0(1.0, q)
    p_formula = alpha0(s, q) * eta**2 / gamma
    f_formula = alpha0(1 - s, q) * zeta**2 / gamma
    anc, dec = bg.precedes(pts, event), bg.follows(pts, event)
    zp = (anc.mean() - p_formula) / math.sqrt(p_formula * (1 - p_formula) / n)
    zf = (dec.mean() - f_formula) / math.sqrt(f_formula * (1 - f_formula) / n)
    check(
        "M1 background p = alpha(s) eta^2/gamma, f = alpha(1-s) zeta^2/gamma (2e6 events)",
        abs(zp) < 4.5 and abs(zf) < 4.5,
        f"p formula {p_formula:.6f} MC {anc.mean():.6f} (z={zp:.2f}); "
        f"f formula {f_formula:.6f} MC {dec.mean():.6f} (z={zf:.2f})",
    )
    refs = pts[:150000]
    ancestors = refs[bg.precedes(refs, event)]
    count = bg.comparable_pairs(ancestors)
    a_hat = count / math.comb(len(refs), 2)
    a_formula = p_formula**2 * p0(q * s * s)
    wrong_denominator = count / math.comb(len(ancestors), 2)
    check(
        "M2 a_X = p_X^2 c_X with denominator choose(M,2) (150000 references)",
        abs(a_hat / a_formula - 1) < 0.08
        and abs(wrong_denominator / a_formula - 1) > 10,
        f"a formula {a_formula:.6e}, MC {a_hat:.6e} (ratio {a_hat / a_formula:.4f}); "
        f"ancestors {len(ancestors)}, c_hat {wrong_denominator:.4f} vs p0 {p0(q * s * s):.4f}",
    )
    # M3: q = 6 (outside the theorem's range, used only to test the pair formula where
    # it differs strongly from the flat value 1/10)
    q6 = 6.0
    bg6 = Wave(q6, 0.25, background=True)
    rng6 = np.random.default_rng(99)
    pts6 = bg6.sample(16000, rng6)
    frac = bg6.comparable_pairs(pts6) / math.comb(len(pts6), 2)
    check(
        "M3 background pair probability at q=6 matches p0(6), not 1/10",
        abs(frac - p0(q6)) < 0.006 and abs(frac - 0.1) > 0.02,
        f"p0(6) = {p0(q6):.5f}, MC {frac:.5f}",
    )


if __name__ == "__main__":
    symbolic()
    numeric()
    monte_carlo()
    print(f"SUMMARY s1_background_calibration: {len(FAILS)} FAIL")
    sys.exit(1 if FAILS else 0)
