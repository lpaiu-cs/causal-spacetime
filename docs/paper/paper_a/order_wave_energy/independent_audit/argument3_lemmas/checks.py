"""Deterministic checks for gram_fourier_lemmas.md.

No order is sampled. Seeded random draws are used only to place test points.
Every integral is a fixed composite Gauss-Legendre rule (16 nodes per panel,
panels at most lambda/4 wide, so the phase e^(iu/lambda) turns by at most a
quarter radian per panel).

Pass criteria (written before the full run):
  G1  Background. Lambda_0(s, t) by quadrature equals alpha(t - s) to relative
      1e-10 (q = 0.01, 0.25, 1; four intervals). The pair fraction
      Lambda_20/(3 Lambda_0^2) equals p0(q l^2) = C1 C3/(3 C2^2), with
      C_r(z) = 0F1(r + 1/2; -z/4)/(2r - 1)!, to 1e-10 (l = 0.2, 0.5, 0.8, 1).
      At q = 1e-8 it equals 1/10 to 1e-6.
  G2  Lemma K for the wave (q = 1; lam = 1/32, 1/128; two phases). The
      composition identity at 40 seeded points holds to 1e-12 (absolute). The
      slice length from the two slacks of Theorem 2 equals the closed form
      h(x) at the same points to 1e-12. The slice area by a 1601 x 1601 grid
      equals pi kappa^2 / sqrt(det Q) to relative 3e-3.
  G3  Corollary K2(ii). |Lambda_lam/c - Lambda_0| <= 3 f_max C_E^2 lam^2
      Phi(l) for four intervals, q = 0.25, 1, two phases, lam = 1/32 to
      1/512. |Lambda_2lam/c^2 - Lambda_20| <= C lam^4 (1 + log(1/lam)) with
      the ratio at every lam at most 3 times its value at lam = 1/32, for
      lam = 1/32, 1/64, 1/128, interval (-1/2, 0.2).
  G4  Lemma E (q = 0.25, 1; lam = 1/32, 1/64, 1/256 with lam <= lambda_1;
      phases 0, 0.37, 2.1; 17 values of u reaching s = 1e-7; four transverse
      directions). zeta preserved to 1e-12. det T equals
      (alpha_0^2 - |E|^2)^(-1/2) to relative 1e-12. |det T - 1| <= 4 e^2,
      ||T - I|| <= 4 e, |delta v| <= 12 f_max C_E lam^2 / M,
      |delta y| <= 4 e sqrt(2 f_max s).
  G5  Lemma F. The numerical (pi/Z_0) int a^2 W_0 |det T - 1| du is below
      the bound (7) (q = 0.25, 1; lam = 1/32, 1/128, 1/512 where lam <=
      lambda_1; two phases). Its ratio to lam^4 log(1/lam) is printed.
  G6  The interior expansion of Delta_A(X): max |exact - first order| /
      lam^4 <= 50 over the three frozen anchors (rotations 0, pi/8), probes
      u in {0.12, 0.2, 0.28}, v~ = 0.7, |z| = 0.15 at 6 angles, lam = 1/32 to
      1/256, two phases. The beta_L expansion: residual / lam^3 <= 20 on
      u in [0.1, 0.3].
  G7  Angular bookkeeping. <sum ell0_i Delta_i> over 8 angles, minus
      Re[beta_L D_A], minus S~: max / lam^4 <= 50 (same anchors, phases,
      lam; u on 9 points of [0.1, 0.3]).
  G8  The Fourier coefficient. E_num = int a^2 w_u <sum ell0 Delta> e^(-i
      omega u) du against E_pred = -(i k/c) lam^2 e^(i phi) int w_u D_A
      e^(-i delta u) du, delta in {0, 1/8}, q = 1, phase 0.37, rotation 0:
      the relative difference decreases from lam = 1/64 to 1/128 to 1/256
      and is below 0.05 at 1/256.
  G9  Exact rational re-certification on the outer box (Fractions; square
      roots and sec^2(1/2) enclosed outward). Pass if: a_3x >= 0.079625,
      a_3y >= 0.079374, |ell_1|, |ell_2| <= 1.009434, |a_2x a_2y| <=
      2.0063e-5, Im D_A >= 0.0106978, |D_A| <= 0.0424, grid coefficient
      >= 0.0101, anchor-pair slack >= 0.0289166, active-probe slack >=
      0.100197.
  G10 Lemma B on a grid: along 30 background null generators, 1{R <_0 X(u)}
      never decreases for 400 grid references (2001 u values; at least
      one switch-on must occur, so the test is not vacuous). Lemma J: E(sum Dbar)^4 =
      n E Dbar^4 + 3 n (n - 1)(E Dbar^2)^2 in exact arithmetic, n = 1..5.

Before the full run (debugging, recorded here): G1, G9 and G10 were run in
full and passed; the helpers behind G2 to G8 were run on single small cases
(lam = 1/32) and not compared with any criterion.

Run: python checks.py   (system Python 3.11 with numpy and mpmath)
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction as Fr

import mpmath as mp
import numpy as np

NODES, WEIGHTS = np.polynomial.legendre.leggauss(16)
RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, info: str) -> None:
    RESULTS.append((name, bool(ok), info))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {info}", flush=True)


# ---------------------------------------------------------------- geometry


class Geo:
    """Background quantities for fixed q."""

    def __init__(self, q: float) -> None:
        self.q = q
        self.k = math.sqrt(q)
        self.fmax = 1.0 / math.cos(self.k / 2) ** 2
        self.CE = 2 * self.k * (4 * self.fmax - 2)
        self.M = self.m(-0.5, 0.5)
        self.gamma = self.alpha(1.0)
        self.lam1 = 1.0 / (4 * self.CE)

    def a(self, u):
        return np.cos(self.k * u)

    def f(self, u):
        return 1.0 / np.cos(self.k * u) ** 2

    def m(self, s, t):
        return (np.tan(self.k * t) - np.tan(self.k * s)) / self.k

    def alpha(self, ell):
        x = self.k * ell
        return (1 - x / np.tan(x)) / (2 * self.q)


def panels(s: float, t: float, width: float):
    n = max(1, math.ceil((t - s) / width - 1e-12))
    edges = np.linspace(s, t, n + 1)
    lo, hi = edges[:-1, None], edges[1:, None]
    x = 0.5 * (hi - lo) * NODES[None, :] + 0.5 * (hi + lo)
    w = 0.5 * (hi - lo) * WEIGHTS[None, :]
    return edges, x.ravel(), w.ravel()


class Wave:
    """beta(s, t) through an accurate primitive of f e^(iu/lam)."""

    def __init__(self, geo: Geo, lam: float, phi: float) -> None:
        self.g = geo
        self.lam = lam
        self.phi = phi
        self.c = math.sqrt(1 + 4 * geo.q * lam**2)
        self.pref = (2 * lam * geo.k / self.c) * np.exp(1j * phi)
        self.width = min(lam / 4, 1 / 64)
        edges, x, w = panels(-0.5, 0.5, self.width)
        self.edges = edges
        vals = (w * self.integrand(x)).reshape(-1, 16).sum(axis=1)
        self.Gedge = np.concatenate([[0.0], np.cumsum(vals)])

    def integrand(self, x):
        return self.g.f(x) * np.exp(1j * x / self.lam)

    def G(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        j = np.clip(np.searchsorted(self.edges, x, side="right") - 1, 0, None)
        j = np.minimum(j, len(self.edges) - 2)
        a0 = self.edges[j]
        span = (x - a0)[:, None]
        nodes = a0[:, None] + span * (NODES[None, :] + 1) / 2
        part = (span * WEIGHTS[None, :] / 2 * self.integrand(nodes)).sum(axis=1)
        return self.Gedge[j] + part

    def beta(self, s, t):
        return self.pref * (self.G(t) - self.G(s))


def bmat(b: complex) -> np.ndarray:
    return np.array([[b.real, b.imag], [b.imag, -b.real]])


def sym_pow(A: np.ndarray, p: float) -> np.ndarray:
    w, V = np.linalg.eigh(A)
    return V @ np.diag(w**p) @ V.T


def frame(wv: Wave, u: float):
    """M_L, M_R, T and the scalars of Lemma E at u."""
    g = wv.g
    m = float(g.m(-0.5, u))
    bL = complex(wv.beta(-0.5, u)[0])
    bR = complex(wv.beta(u, 0.5)[0])
    ML = m * np.eye(2) - bmat(bL)
    MR = (g.M - m) * np.eye(2) - bmat(bR)
    W0 = m * (g.M - m) / g.M
    Q = np.linalg.inv(ML) + np.linalg.inv(MR)
    T = sym_pow(W0 * Q, -0.5)
    return dict(m=m, bL=bL, bR=bR, ML=ML, MR=MR, W0=W0, Q=Q, T=T)


def to_phys(wv: Wave, u: float, vt: float, z: np.ndarray):
    fr = frame(wv, u)
    y = fr["T"] @ z
    v = vt - z @ z / (2 * fr["m"]) + y @ np.linalg.solve(fr["ML"], y) / 2
    return v, y, fr


# ------------------------------------------------------------------ G1


def C_r(r: int, z: float):
    return mp.hyp0f1(r + mp.mpf(1) / 2, -mp.mpf(z) / 4) / mp.factorial(2 * r - 1)


def g1() -> None:
    worst = 0.0
    for q in (0.01, 0.25, 1.0):
        geo = Geo(q)
        for s, t in ((-0.5, 0.5), (-0.3, 0.2), (0.1, 0.45), (-0.5, -0.2)):
            _, x, w = panels(s, t, 1 / 64)
            lam0 = np.sum(w * geo.a(x) ** 2 * geo.m(s, x) * geo.m(x, t) / geo.m(s, t))
            worst = max(worst, abs(lam0 / geo.alpha(t - s) - 1))
    record("G1 Lambda_0 = alpha", worst < 1e-10, f"max rel err {worst:.2e}")

    worst = 0.0
    for q in (0.25, 1.0):
        geo = Geo(q)
        for ell in (0.2, 0.5, 0.8, 1.0):
            s, t = -0.5, -0.5 + ell
            _, x, w = panels(s, t, 1 / 64)
            kern = geo.a(x) ** 2 * geo.m(s, x) * geo.m(x, t) / geo.m(s, t)
            lam2 = np.sum(w * geo.alpha(t - x) * kern)
            c = lam2 / (3 * geo.alpha(ell) ** 2)
            z = q * ell**2
            p0 = float(C_r(1, z) * C_r(3, z) / (3 * C_r(2, z) ** 2))
            worst = max(worst, abs(c - p0))
    record("G1 pair fraction = p0", worst < 1e-10, f"max abs err {worst:.2e}")

    geo = Geo(1e-8)
    _, x, w = panels(-0.5, 0.5, 1 / 64)
    kern = geo.a(x) ** 2 * geo.m(-0.5, x) * geo.m(x, 0.5) / geo.m(-0.5, 0.5)
    lam2 = np.sum(w * geo.alpha(0.5 - x) * kern)
    c = lam2 / (3 * geo.alpha(1.0) ** 2)
    record("G1 flat limit 1/10", abs(c - 0.1) < 1e-6, f"c = {c:.9f}")


# ------------------------------------------------------------------ G2


def kappa(wv: Wave, P1, P2) -> float:
    (u1, v1, y1), (u2, v2, y2) = P1, P2
    Mt = float(wv.g.m(u1, u2)) * np.eye(2) - bmat(complex(wv.beta(u1, u2)[0]))
    dy = y2 - y1
    return v2 - v1 - dy @ np.linalg.solve(Mt, dy) / 2


def g2() -> None:
    rng = np.random.default_rng(20261126)
    geo = Geo(1.0)
    werr_comp = werr_len = werr_area = 0.0
    for lam in (1 / 32, 1 / 128):
        for phi in (0.37, 2.1):
            wv = Wave(geo, lam, phi)
            u1, u2 = -0.3, 0.25
            y1 = np.array([0.02, -0.05])
            y2 = np.array([0.1, 0.04])
            v1 = 0.2
            Mt = float(geo.m(u1, u2)) * np.eye(2) - bmat(complex(wv.beta(u1, u2)[0]))
            dy = y2 - y1
            v2 = v1 + dy @ np.linalg.solve(Mt, dy) / 2 + 0.3
            P1, P2 = (u1, v1, y1), (u2, v2, y2)
            kap = kappa(wv, P1, P2)
            for _ in range(40):
                u = rng.uniform(u1 + 0.01, u2 - 0.01)
                x = rng.normal(0, 0.1, 2)
                M1 = float(geo.m(u1, u)) * np.eye(2) - bmat(complex(wv.beta(u1, u)[0]))
                M2 = float(geo.m(u, u2)) * np.eye(2) - bmat(complex(wv.beta(u, u2)[0]))
                Qu = np.linalg.inv(M1) + np.linalg.inv(M2)
                d1 = np.linalg.solve(Qu, np.linalg.solve(M2, dy))
                xs = y1 + d1
                vv = rng.uniform(0, 1)
                lhs = kappa(wv, P1, (u, vv, x)) + kappa(wv, (u, vv, x), P2)
                rhs = kap - (x - xs) @ Qu @ (x - xs) / 2
                werr_comp = max(werr_comp, abs(lhs - rhs))
                lo = v1 + (x - y1) @ np.linalg.solve(M1, x - y1) / 2
                hi = v2 - (y2 - x) @ np.linalg.solve(M2, y2 - x) / 2
                werr_len = max(werr_len, abs(max(0.0, hi - lo) - max(0.0, rhs)))
            # slice area at one u by a grid, against the closed form
            u = 0.0
            M1 = float(geo.m(u1, u)) * np.eye(2) - bmat(complex(wv.beta(u1, u)[0]))
            M2 = float(geo.m(u, u2)) * np.eye(2) - bmat(complex(wv.beta(u, u2)[0]))
            Qu = np.linalg.inv(M1) + np.linalg.inv(M2)
            xs = y1 + np.linalg.solve(Qu, np.linalg.solve(M2, dy))
            R = math.sqrt(2 * kap / np.linalg.eigvalsh(Qu).min()) * 1.02
            g = np.linspace(-R, R, 1601)
            X1, X2 = np.meshgrid(g + xs[0], g + xs[1], indexing="ij")
            pts = np.stack([X1 - y1[0], X2 - y1[1]], axis=-1)
            pts2 = np.stack([y2[0] - X1, y2[1] - X2], axis=-1)
            q1 = np.einsum("...i,ij,...j", pts, np.linalg.inv(M1), pts) / 2
            q2 = np.einsum("...i,ij,...j", pts2, np.linalg.inv(M2), pts2) / 2
            hlen = np.maximum(0.0, (v2 - q2) - (v1 + q1))
            area = hlen.sum() * (g[1] - g[0]) ** 2
            exact = math.pi * kap**2 / math.sqrt(np.linalg.det(Qu))
            werr_area = max(werr_area, abs(area / exact - 1))
    record("G2 composition identity", werr_comp < 1e-12, f"max err {werr_comp:.1e}")
    record("G2 slice length", werr_len < 1e-12, f"max err {werr_len:.1e}")
    record("G2 slice area", werr_area < 3e-3, f"max rel err {werr_area:.2e}")


# ------------------------------------------------------------------ G3


def lam_wave_over_c(wv: Wave, s: float, t: float) -> float:
    g = wv.g
    _, x, w = panels(s, t, wv.width)
    Gs, Gt, Gx = wv.G(s)[0], wv.G(t)[0], wv.G(x)
    b1 = np.abs(wv.pref * (Gx - Gs)) ** 2
    b2 = np.abs(wv.pref * (Gt - Gx)) ** 2
    b = abs(wv.pref * (Gt - Gs)) ** 2
    m1, m2, mm = g.m(s, x), g.m(x, t), g.m(s, t)
    integ = g.a(x) ** 2 * np.sqrt((m1**2 - b1) * (m2**2 - b2) / (mm**2 - b))
    return float(np.sum(w * integ))


def lam2_wave_over_c2(wv: Wave, s: float, t: float) -> float:
    """Outer w over global nodes; inner over [w, t] = partial panel + rest."""
    g = wv.g
    edges, x, w = panels(s, t, wv.width)
    Gx = wv.G(x)
    Gs, Gt = wv.G(s)[0], wv.G(t)[0]
    pan = np.repeat(np.arange(len(edges) - 1), 16)
    total = 0.0
    for i in range(len(x)):
        wi, Gw = x[i], Gx[i]
        b_end = edges[pan[i] + 1]
        span = b_end - wi
        pn = wi + span * (NODES + 1) / 2
        pw = span * WEIGHTS / 2
        # G at partial nodes: G(w) + int_w^node
        sub = (pn - wi)[:, None]
        inner_nodes = wi + sub * (NODES[None, :] + 1) / 2
        Gp = Gw + (sub * WEIGHTS[None, :] / 2 * wv.integrand(inner_nodes)).sum(axis=1)
        mask = pan > pan[i]
        xn = np.concatenate([pn, x[mask]])
        wn = np.concatenate([pw, w[mask]])
        Gn = np.concatenate([Gp, Gx[mask]])
        b1 = np.abs(wv.pref * (Gn - Gw)) ** 2
        b2 = np.abs(wv.pref * (Gt - Gn)) ** 2
        b = abs(wv.pref * (Gt - Gw)) ** 2
        m1, m2, mm = g.m(wi, xn), g.m(xn, t), g.m(wi, t)
        integ = g.a(xn) ** 2 * np.sqrt((m1**2 - b1) * (m2**2 - b2) / (mm**2 - b))
        lam_wt = float(np.sum(wn * integ))
        o1 = abs(wv.pref * (Gw - Gs)) ** 2
        o2 = abs(wv.pref * (Gt - Gw)) ** 2
        ob = abs(wv.pref * (Gt - Gs)) ** 2
        M1, M2, MM = g.m(s, wi), g.m(wi, t), g.m(s, t)
        outer = g.a(wi) ** 2 * math.sqrt((M1**2 - o1) * (M2**2 - o2) / (MM**2 - ob))
        total += w[i] * lam_wt * outer
    return total


def lam2_background(geo: Geo, s: float, t: float) -> float:
    _, x, w = panels(s, t, 1 / 256)
    kern = geo.a(x) ** 2 * geo.m(s, x) * geo.m(x, t) / geo.m(s, t)
    return float(np.sum(w * geo.alpha(t - x) * kern))


def Phi(ell: float, lam: float) -> float:
    if ell <= lam:
        return ell**2 / 2
    return lam**2 / 2 + lam**2 * math.log(ell / lam)


def g3() -> None:
    ok = True
    ratios = []
    for q in (0.25, 1.0):
        geo = Geo(q)
        for phi in (0.37, 2.1):
            for lam in (1 / 32, 1 / 64, 1 / 128, 1 / 256, 1 / 512):
                if lam > geo.lam1:
                    continue
                wv = Wave(geo, lam, phi)
                for s, t in ((-0.5, 0.2), (-0.5, 0.5), (-0.25, 0.2), (0.2, 0.5)):
                    dev = abs(lam_wave_over_c(wv, s, t) - geo.alpha(t - s))
                    bound = 3 * geo.fmax * geo.CE**2 * lam**2 * Phi(t - s, lam)
                    ok &= dev <= bound
                    ratios.append(dev / (lam**4 * (1 + math.log(1 / lam))))
    record(
        "G3 Lambda bound",
        ok,
        f"max dev/(lam^4(1+log)) = {max(ratios):.3f}",
    )
    geo = Geo(1.0)
    rat = []
    for lam in (1 / 32, 1 / 64, 1 / 128):
        wv = Wave(geo, lam, 0.37)
        dev = abs(lam2_wave_over_c2(wv, -0.5, 0.2) - lam2_background(geo, -0.5, 0.2))
        rat.append(dev / (lam**4 * (1 + math.log(1 / lam))))
    ok2 = all(r <= 3 * rat[0] for r in rat)
    record(
        "G3 Lambda_2 rate",
        ok2,
        "dev/(lam^4(1+log)) = " + ", ".join(f"{r:.3f}" for r in rat),
    )


# ------------------------------------------------------------------ G4, G5


def g4() -> None:
    ok = {"zeta": True, "det": True, "dT": True, "TI": True, "dv": True, "dy": True}
    worst = {"zeta": 0.0, "det": 0.0, "dT": 0.0, "TI": 0.0, "dv": 0.0, "dy": 0.0}
    for q in (0.25, 1.0):
        geo = Geo(q)
        for lam in (1 / 32, 1 / 64, 1 / 256):
            if lam > geo.lam1:
                continue
            for phi in (0.0, 0.37, 2.1):
                wv = Wave(geo, lam, phi)
                svals = (1e-7, 1e-5, 1e-3, lam / 3, lam, 3 * lam, 0.05, 0.2)
                us = [0.0] + [-0.5 + sv for sv in svals] + [0.5 - sv for sv in svals]
                for u in us:
                    s = min(u + 0.5, 0.5 - u)
                    e = geo.CE * lam * min(1.0, lam / s)
                    fr = frame(wv, u)
                    m, W0, T = fr["m"], fr["W0"], fr["T"]
                    eL = abs(fr["bL"]) / m
                    eR = abs(fr["bR"]) / (geo.M - m)
                    sig = m / geo.M
                    a0 = (1 - sig) / (1 - eL**2) + sig / (1 - eR**2)
                    EL, ER = bmat(fr["bL"]) / m, bmat(fr["bR"]) / (geo.M - m)
                    E = (1 - sig) * EL / (1 - eL**2) + sig * ER / (1 - eR**2)
                    En = np.linalg.norm(E, 2)
                    detT = np.linalg.det(T)
                    closed = (a0**2 - En**2) ** -0.5
                    worst["det"] = max(worst["det"], abs(detT / closed - 1))
                    worst["dT"] = max(worst["dT"], abs(detT - 1) / (4 * e**2))
                    TI = np.linalg.norm(T - np.eye(2), 2)
                    worst["TI"] = max(worst["TI"], TI / (4 * e))
                    for ang in (0.3, 1.4, 2.6, 4.0):
                        r = 0.95 * math.sqrt(2 * W0) * math.sqrt(0.5)
                        z = r * np.array([math.cos(ang), math.sin(ang)])
                        vt = 0.25 + z @ z / (2 * m)
                        v, y, _ = to_phys(wv, u, vt, z)
                        zeta_l = 1 - v - y @ np.linalg.solve(fr["MR"], y) / 2
                        zeta_0 = 1 - vt - z @ z / (2 * (geo.M - m))
                        worst["zeta"] = max(worst["zeta"], abs(zeta_l - zeta_0))
                        bdv = 12 * geo.fmax * geo.CE * lam**2 / geo.M
                        worst["dv"] = max(worst["dv"], abs(v - vt) / bdv)
                        bdy = 4 * e * math.sqrt(2 * geo.fmax * s)
                        worst["dy"] = max(worst["dy"], np.linalg.norm(y - z) / bdy)
    ok["zeta"] = worst["zeta"] < 1e-12
    ok["det"] = worst["det"] < 1e-12
    for key in ("dT", "TI", "dv", "dy"):
        ok[key] = worst[key] <= 1.0
    record("G4 zeta preserved", ok["zeta"], f"max err {worst['zeta']:.1e}")
    record("G4 det T closed form", ok["det"], f"max rel err {worst['det']:.1e}")
    record("G4 |det T - 1| <= 4e^2", ok["dT"], f"max ratio {worst['dT']:.3f}")
    record("G4 ||T - I|| <= 4e", ok["TI"], f"max ratio {worst['TI']:.3f}")
    record("G4 |delta v| bound", ok["dv"], f"max ratio {worst['dv']:.4f}")
    record("G4 |delta y| bound", ok["dy"], f"max ratio {worst['dy']:.3f}")


def g5() -> None:
    ok = True
    info = []
    for q in (0.25, 1.0):
        geo = Geo(q)
        for lam in (1 / 32, 1 / 128, 1 / 512):
            if lam > geo.lam1:
                continue
            for phi in (0.37, 2.1):
                wv = Wave(geo, lam, phi)
                _, x, w = panels(-0.5, 0.5, wv.width)
                m = geo.m(-0.5, x)
                bL = wv.beta(-0.5, x)
                bR = wv.beta(x, 0.5)
                eL = np.abs(bL) / m
                eR = np.abs(bR) / (geo.M - m)
                sig = m / geo.M
                a0 = (1 - sig) / (1 - eL**2) + sig / (1 - eR**2)
                # |E| for E = (1-sig) EL/(1-eL^2) + sig ER/(1-eR^2), complex rep
                Ec = (1 - sig) * (bL / m) / (1 - eL**2) + sig * (bR / (geo.M - m)) / (
                    1 - eR**2
                )
                detT = (a0**2 - np.abs(Ec) ** 2) ** -0.5
                W0 = m * (geo.M - m) / geo.M
                Z0 = math.pi * geo.gamma
                tv = math.pi / Z0 * np.sum(w * geo.a(x) ** 2 * W0 * np.abs(detT - 1))
                bound = (
                    4
                    * geo.fmax
                    * geo.CE**2
                    / geo.gamma
                    * lam**4
                    * (1 + 2 * math.log(1 / (2 * lam)))
                )
                ok &= tv <= bound
                info.append(f"{tv / (lam**4 * math.log(1 / lam)):.2f}")
    record("G5 TV bound (7)", ok, "tv/(lam^4 log) = " + ", ".join(info))


# ------------------------------------------------------------------ G6-G8


def anchors(theta: float):
    d = 0.08
    R = np.array(
        [[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]]
    )
    zs = [R @ np.array(v) for v in ((d, 0.0), (0.0, d), (d, d))]
    us = (-0.35, -0.25, -0.15)
    vs = (0.15, 0.25, 0.35)
    GA = np.array([[zs[0] @ zs[0], zs[0] @ zs[1]], [zs[0] @ zs[1], zs[1] @ zs[1]]])
    l12 = -np.array([zs[0] @ zs[2], zs[1] @ zs[2]]) @ np.linalg.inv(GA)
    ell = np.array([l12[0], l12[1], 1.0])
    return us, vs, zs, ell


def delta_exact_and_first(wv: Wave, A, X):
    g = wv.g
    (uA, vA, zA), (uX, vX, zX) = A, X
    vpA, yA, frA = to_phys(wv, uA, vA, zA)
    vpX, yX, frX = to_phys(wv, uX, vX, zX)
    dm = float(g.m(uA, uX))
    BAX = bmat(complex(wv.beta(uA, uX)[0]))
    dy = yX - yA
    kap_l = vpX - vpA - dy @ np.linalg.solve(dm * np.eye(2) - BAX, dy) / 2
    dz = zX - zA
    kap_c = vX - vA - dz @ dz / (2 * dm)
    exact = dm * (kap_l - kap_c)
    tX, tA = frX["T"] - np.eye(2), frA["T"] - np.eye(2)
    first = (
        dm * ((vpX - vX) - (vpA - vA))
        - dz @ (tX @ zX - tA @ zA)
        - dz @ BAX @ dz / (2 * dm)
    )
    return exact, first


def g6() -> None:
    worst = 0.0
    for phi in (0.37, 2.1):
        for lam in (1 / 32, 1 / 64, 1 / 128, 1 / 256):
            wv = Wave(Geo(1.0), lam, phi)
            for theta in (0.0, math.pi / 8):
                us, vs, zs, _ = anchors(theta)
                for i in range(3):
                    for uX in (0.12, 0.2, 0.28):
                        for ang in np.linspace(0, 2 * math.pi, 7)[:-1]:
                            zX = 0.15 * np.array([math.cos(ang), math.sin(ang)])
                            ex, fi = delta_exact_and_first(
                                wv, (us[i], vs[i], zs[i]), (uX, 0.7, zX)
                            )
                            worst = max(worst, abs(ex - fi) / lam**4)
    record("G6 Delta first order", worst <= 50, f"max resid/lam^4 = {worst:.3f}")
    worst = 0.0
    geo = Geo(1.0)
    for phi in (0.37, 2.1):
        for lam in (1 / 32, 1 / 64, 1 / 128, 1 / 256):
            wv = Wave(geo, lam, phi)
            u = np.linspace(0.1, 0.3, 41)
            b = wv.beta(-0.5, u)
            lead = (
                (2 * geo.k / wv.c)
                * np.exp(1j * phi)
                * (-1j * lam**2 * geo.f(u) * np.exp(1j * u / lam))
            )
            cst = wv.pref * 1j * lam * geo.f(-0.5) * np.exp(-1j / (2 * lam))
            worst = max(worst, float(np.max(np.abs(b - lead - cst))) / lam**3)
    record("G6 beta_L expansion", worst <= 20, f"max resid/lam^3 = {worst:.3f}")


def D_A(geo: Geo, u, us, zs, ell):
    out = 0
    for i in range(3):
        zeta = zs[i][0] + 1j * zs[i][1]
        out = out - 0.5 * ell[i] * np.conj(zeta) ** 2 / (geo.m(us[i], u))
    return out


def angle_avg(wv: Wave, uX: float, theta: float, nang: int = 8) -> float:
    us, vs, zs, ell = anchors(theta)
    tot = 0.0
    for ang in np.arange(nang) * 2 * math.pi / nang:
        zX = 0.15 * np.array([math.cos(ang), math.sin(ang)])
        for i in range(3):
            ex, _ = delta_exact_and_first(wv, (us[i], vs[i], zs[i]), (uX, 0.7, zX))
            tot += ell[i] * ex
    return tot / nang


def S_tilde(wv: Wave, uX: float, theta: float) -> float:
    g = wv.g
    us, vs, zs, ell = anchors(theta)
    out = 0.0
    for i in range(3):
        vp, y, fr = to_phys(wv, us[i], vs[i], zs[i])
        BLi = bmat(fr["bL"])
        dmi = float(g.m(us[i], uX))
        tau = fr["T"] - np.eye(2)
        out += ell[i] * (
            zs[i] @ BLi @ zs[i] / (2 * dmi) - dmi * (vp - vs[i]) - zs[i] @ tau @ zs[i]
        )
    return out


def g7() -> None:
    worst = 0.0
    geo = Geo(1.0)
    for phi in (0.37, 2.1):
        for lam in (1 / 32, 1 / 64, 1 / 128, 1 / 256):
            wv = Wave(geo, lam, phi)
            for theta in (0.0, math.pi / 8):
                us, vs, zs, ell = anchors(theta)
                for uX in np.linspace(0.1, 0.3, 9):
                    avg = angle_avg(wv, uX, theta)
                    bL = complex(wv.beta(-0.5, uX)[0])
                    loc = (bL * D_A(geo, uX, us, zs, ell)).real
                    res = avg - loc - S_tilde(wv, uX, theta)
                    worst = max(worst, abs(res) / lam**4)
    record("G7 angular bookkeeping", worst <= 50, f"max resid/lam^4 = {worst:.3f}")


def w_u(u):
    x = (u - 0.2) / 0.1
    out = np.zeros_like(u)
    inside = np.abs(x) < 1
    out[inside] = np.exp(-1 / (1 - x[inside] ** 2))
    return out


def g8() -> None:
    geo = Geo(1.0)
    phi, theta = 0.37, 0.0
    us, vs, zs, ell = anchors(theta)
    rels = {0.0: [], 0.125: []}
    for lam in (1 / 64, 1 / 128, 1 / 256):
        wv = Wave(geo, lam, phi)
        _, x, w = panels(0.1, 0.3, lam / 4)
        avg = np.array([angle_avg(wv, float(u), theta) for u in x])
        for dl in (0.0, 0.125):
            om = 1 / lam + dl
            e_num = np.sum(w * geo.a(x) ** 2 * w_u(x) * avg * np.exp(-1j * om * x))
            da = D_A(geo, x, us, zs, ell)
            e_pred = (
                -(1j * geo.k / wv.c)
                * lam**2
                * np.exp(1j * phi)
                * np.sum(w * w_u(x) * da * np.exp(-1j * dl * x))
            )
            rels[dl].append(abs(e_num - e_pred) / abs(e_pred))
    ok = True
    for r in rels.values():
        ok &= r[0] > r[1] > r[2] and r[2] < 0.05
    info = "; ".join(
        f"delta={dl}: " + ", ".join(f"{v:.2e}" for v in r) for dl, r in rels.items()
    )
    record("G8 Fourier coefficient", ok, info)


# ------------------------------------------------------------------ G9


SCALE = 10**40


def sqrt_lo(x: Fr) -> Fr:
    n = math.isqrt(x.numerator * SCALE * SCALE // x.denominator)
    return Fr(n, SCALE)


def sqrt_hi(x: Fr) -> Fr:
    return sqrt_lo(x) + Fr(1, SCALE)


def g9() -> None:
    t = Fr(2, 100000)  # Gram tolerance 0.00002
    du = Fr(2, 1000)  # u and v tolerance 0.002
    G11 = (Fr(64, 10000) - t, Fr(64, 10000) + t)
    G22 = G11
    G33 = (Fr(128, 10000) - t, Fr(128, 10000) + t)
    G13 = G11
    G23 = G11
    G12abs = t
    a3x = G13[0] / sqrt_hi(G11[1])
    a3y = (G23[0] - G12abs * G13[1] / G11[0]) / sqrt_hi(G22[1])
    det_lo = G11[0] * G22[0] - G12abs**2
    ell2 = G23[1] / (G22[0] - G12abs**2 / G11[0]) + G13[1] * G12abs / det_lo
    ell1 = G13[1] / (G11[0] - G12abs**2 / G22[0]) + G23[1] * G12abs / det_lo
    a2xy = G12abs * sqrt_hi(G22[1] / G11[0])
    # sec^2(1/2) upper: cos(1/2) >= 1 - x^2/2 + x^4/24 - x^6/720 at x = 1/2
    xh = Fr(1, 2)
    cos_lo = 1 - xh**2 / 2 + xh**4 / 24 - xh**6 / 720
    sec2_hi = 1 / cos_lo**2
    u3_lo = Fr(-15, 100) - du
    u2_hi = Fr(-25, 100) + du
    u1_hi = Fr(-35, 100) + du
    u3_hi = Fr(-15, 100) + du
    up_lo, up_hi = Fr(1, 10), Fr(3, 10)
    gap3_hi = sec2_hi * (up_hi - u3_lo)  # m(u) - m_3 <= f_max (u - u_3)
    gap2_lo = up_lo - u2_hi  # m' >= 1
    gap1_lo = up_lo - u1_hi
    gap3_lo = up_lo - u3_hi
    im_da = a3x * a3y / gap3_hi - ell2 * a2xy / gap2_lo
    abs_da = (ell1 * G11[1] / gap1_lo + ell2 * G22[1] / gap2_lo + G33[1] / gap3_lo) / 2
    grid = im_da - abs_da * Fr(1, 10) / 8
    # anchor-pair slack: min over pairs of dv_lo - |z_i - z_j|^2_hi / (2 du_lo)
    pairs = {
        (1, 2): (G11[1] + G22[1] + 2 * G12abs, Fr(1, 10) - 2 * du),
        (2, 3): (G22[1] + G33[1] - 2 * G23[0], Fr(1, 10) - 2 * du),
        (1, 3): (G11[1] + G33[1] - 2 * G13[0], Fr(2, 10) - 2 * du),
    }
    pair_slack = min(dd - dz2 / (2 * dd) for dz2, dd in pairs.values())
    # active-probe slack, probes v in [0.65, 0.75], R^2 in [0.0025, 0.04]
    vX_lo, R2_hi = Fr(65, 100), Fr(4, 100)
    probe = []
    for vi, ui, Gii in (
        (Fr(15, 100), Fr(-35, 100), G11[1]),
        (Fr(25, 100), Fr(-25, 100), G22[1]),
        (Fr(35, 100), Fr(-15, 100), G33[1]),
    ):
        dz_hi = (sqrt_hi(R2_hi) + sqrt_hi(Gii)) ** 2
        probe.append((vX_lo - (vi + du)) - dz_hi / (2 * (up_lo - (ui + du))))
    probe_slack = min(probe)
    checks = [
        ("a_3x", a3x, Fr(79625, 10**6), True),
        ("a_3y", a3y, Fr(79374, 10**6), True),
        ("|ell_1|", ell1, Fr(1009434, 10**6), False),
        ("|ell_2|", ell2, Fr(1009434, 10**6), False),
        ("|a_2x a_2y|", a2xy, Fr(20063, 10**9), False),
        ("Im D_A", im_da, Fr(106978, 10**7), True),
        ("|D_A|", abs_da, Fr(424, 10**4), False),
        ("grid coefficient", grid, Fr(101, 10**4), True),
        ("anchor-pair slack", pair_slack, Fr(289166, 10**7), True),
        ("active-probe slack", probe_slack, Fr(100197, 10**6), True),
    ]
    for name, val, target, lower in checks:
        ok = val >= target if lower else val <= target
        rel = ">=" if lower else "<="
        record(f"G9 {name}", ok, f"{float(val):.10f} {rel} {float(target)}")


# ------------------------------------------------------------------ G10


def g10() -> None:
    geo = Geo(1.0)

    def before(R, u, v, y1, y2):
        # background relation R <=_0 X(u) in the (T, Z, y) chart, vector in u
        mR, mX = geo.m(-0.5, R[0]), geo.m(-0.5, u)
        dT = (mX + v - mR - R[1]) / math.sqrt(2)
        dZ = (mX - v - mR + R[1]) / math.sqrt(2)
        dy = math.hypot(y1 - R[2], y2 - R[3])
        return dT >= np.hypot(dZ, dy)

    grid = [
        (u, v, y1, 0.0)
        for u in np.linspace(-0.45, 0.45, 10)
        for v in np.linspace(0.05, 0.95, 8)
        for y1 in np.linspace(-0.2, 0.2, 5)
    ]
    violations = 0
    switches = 0
    rng = np.random.default_rng(7)
    ugrid = np.linspace(-0.5, 0.5, 2001)
    for _ in range(30):
        v, y1, y2 = (
            rng.uniform(0.1, 0.9),
            rng.uniform(-0.2, 0.2),
            rng.uniform(-0.2, 0.2),
        )
        for R in grid:
            ind = before(R, ugrid, v, y1, y2).astype(int)
            steps = np.diff(ind)
            violations += int(np.sum(steps < 0))
            switches += int(np.sum(steps > 0))
    record(
        "G10 monotone past indicators",
        violations == 0 and switches > 0,
        f"violations {violations}, switches on {switches}",
    )

    vals = [Fr(-1), Fr(0), Fr(2)]
    probs = [Fr(1, 6), Fr(1, 2), Fr(1, 3)]
    mean = sum(v * p for v, p in zip(vals, probs, strict=True))
    cen = [v - mean for v in vals]
    m2 = sum(c**2 * p for c, p in zip(cen, probs, strict=True))
    m4 = sum(c**4 * p for c, p in zip(cen, probs, strict=True))
    ok = True
    for n in range(1, 6):
        # exact E (sum of n iid centred)^4 by enumeration
        tot = Fr(0)
        for idx in np.ndindex(*([3] * n)):
            s = sum(cen[i] for i in idx)
            p = Fr(1)
            for i in idx:
                p *= probs[i]
            tot += s**4 * p
        ok &= tot == n * m4 + 3 * n * (n - 1) * m2**2
    record("G10 fourth-moment identity", ok, "n = 1..5 exact")


def main() -> int:
    for fn in (g1, g2, g3, g4, g5, g6, g7, g8, g9, g10):
        fn()
    n_pass = sum(ok for _, ok, _ in RESULTS)
    n_fail = len(RESULTS) - n_pass
    print(f"\n{n_pass} pass, {n_fail} fail")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
