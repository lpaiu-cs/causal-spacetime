# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written in this folder by the fresh-context AI referee of local_null_challenge_r2.md.
"""Referee check r4: Lemma 4.1 (the uniform bias bound) on the whole of K.

Lemma 4.1 imports [GFL] Corollary K2 and claims, uniformly on
K = {|u| <= 0.3, eta_0 >= 1/4, zeta_0 >= 2/5},

  |theta_lam(xi) - theta_0(X~)|_inf + |F(theta_lam(xi)) - v~(X~)|
      <= C_b lam^4 (1 + log(1/lam)).

The audit referee measured the F-part only at u in {-0.15, 0.05, 0.25}
(and once at -0.2). Here the exact interval volumes of [GFL] Lemma K are
computed by composite Gauss-Legendre quadrature, wave and background with the
same nodes. theta depends on the canonical point only through
(u, eta, zeta); the wave/background ratios depend on u only, and neither
depends on the phase (only |beta| enters). Grid: u in {-0.3, ..., 0.3}
(step 0.1), (eta, zeta) on 7 points of K's range, lam in {1/40, ..., 1/320};
then u = -0.3 (where the sup sits) down to lam = 1/2560.

  Q0  the background quadrature reproduces alpha(u + 1/2) and
      p0(q (u + 1/2)^2) (0F1 series, mpmath) to 1e-12 relative.
  B2  single logarithm: at u = -0.3 and the worst (eta, zeta), the slope of
      y = bias / lam^4 against L = log(1/lam) over the last two doublings
      agree within 25%, for theta and for F. (A lam^4 log^2 law would raise
      the slope by about 13% per doubling at these L; a lam^3 law would
      double it.) The implied constant C_b (the slope) is printed.

Revision note. The first version of this script (output r4_output_v1.txt,
kept) used lam in {1/32, ..., 1/256} and the criterion "the ratio
bias / (lam^4 (1 + log)) at 1/256 is at most twice its value at 1/64". It
FAILED for F at both q (factor 2.04): the ratio approaches its limit from
below, because bias / lam^4 = A + B log(1/lam) with A < 0 here. That
criterion did not test what it meant to. This version tests the log-slope
instead, over a wider range of lam.
Run: python r4_bias_on_K.py   (a few minutes)
"""

from __future__ import annotations

import math
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 40
XG, WG = np.polynomial.legendre.leggauss(10)
OK = []


def rec(name, ok, msg):
    OK.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


class Vol:
    """Interval volumes Lambda(s,t), Lambda_2(-1/2,u) on panel edges."""

    def __init__(self, q, lam, wave):
        self.k = math.sqrt(q)
        self.lam = lam
        self.wave = wave
        c = math.sqrt(1 + 4 * q * lam * lam)
        self.pref = (2 * lam * self.k / c) if wave else 0.0
        h = lam / 2
        self.n = int(round(1 / h))
        self.h = h
        self.e = -0.5 + h * np.arange(self.n + 1)
        lo = self.e[:-1, None]
        self.x = lo + h * (XG[None, :] + 1) / 2  # (n, 10)
        self.w = np.broadcast_to(h * WG / 2, self.x.shape)
        # G at edges and at nodes
        full = (self.w * self.g(self.x)).sum(axis=1)
        self.Ge = np.concatenate([[0], np.cumsum(full)])
        self.Gx = self.Ge[:-1, None] + self.partial(lo, self.x)

    def g(self, x):
        return np.exp(1j * x / self.lam) / np.cos(self.k * x) ** 2

    def partial(self, a, b):
        """int_a^b g, a and b broadcastable, b - a <= one panel."""
        a = np.asarray(a)[..., None]
        b = np.asarray(b)[..., None]
        nodes = a + (b - a) * (XG + 1) / 2
        return ((b - a) * WG / 2 * self.g(nodes)).sum(axis=-1)

    def tan_(self, x):
        return np.tan(self.k * x) / self.k

    def D(self, Ta, Ga, Tb, Gb):
        m = Tb - Ta
        return m * m - (self.pref**2) * np.abs(Gb - Ga) ** 2

    def a2(self, x):
        return np.cos(self.k * x) ** 2

    def lam1(self, i, j):
        """Lambda(e_i, e_j)/c over panels i..j-1."""
        x, w = self.x[i:j].ravel(), self.w[i:j].ravel()
        Gx = self.Gx[i:j].ravel()
        Ts, Tt, Tx = self.tan_(self.e[i]), self.tan_(self.e[j]), self.tan_(x)
        Gs, Gt = self.Ge[i], self.Ge[j]
        val = self.a2(x) * np.sqrt(
            self.D(Ts, Gs, Tx, Gx) * self.D(Tx, Gx, Tt, Gt) / self.D(Ts, Gs, Tt, Gt)
        )
        return float(np.sum(w * val))

    def lam2(self, j):
        """Lambda_2(-1/2, e_j)/c^2."""
        Tt, Gt = self.tan_(self.e[j]), self.Ge[j]
        Ts, Gs = self.tan_(self.e[0]), self.Ge[0]
        Dst = self.D(Ts, Gs, Tt, Gt)
        total = 0.0
        for p in range(j):
            wo = self.x[p]  # outer nodes in panel p
            Two, Gwo = self.tan_(wo), self.Gx[p]
            # inner partial panel [w, e_{p+1}]
            right = self.e[p + 1]
            xi = wo[:, None] + (right - wo)[:, None] * (XG[None, :] + 1) / 2
            wi = (right - wo)[:, None] * WG[None, :] / 2
            Gxi = self.Ge[p] + self.partial(self.e[p], xi)
            Txi = self.tan_(xi)
            vals = self.a2(xi) * np.sqrt(
                self.D(Two[:, None], Gwo[:, None], Txi, Gxi)
                * self.D(Txi, Gxi, Tt, Gt)
                / self.D(Two[:, None], Gwo[:, None], Tt, Gt)
            )
            inner = (wi * vals).sum(axis=1)
            if p + 1 < j:
                xf = self.x[p + 1 : j].ravel()
                wf = self.w[p + 1 : j].ravel()
                Gf = self.Gx[p + 1 : j].ravel()
                Tf = self.tan_(xf)
                vf = self.a2(xf)[None, :] * np.sqrt(
                    self.D(Two[:, None], Gwo[:, None], Tf[None, :], Gf[None, :])
                    * self.D(Tf[None, :], Gf[None, :], Tt, Gt)
                    / self.D(Two[:, None], Gwo[:, None], Tt, Gt)
                )
                inner = inner + (vf * wf[None, :]).sum(axis=1)
            outer = self.a2(wo) * np.sqrt(
                self.D(Ts, Gs, Two, Gwo) * self.D(Two, Gwo, Tt, Gt) / Dst
            )
            total += float(np.sum(self.w[p] * inner * outer))
        return total


def p0_mp(z):
    def C(r):
        return mp.hyp0f1(r + mp.mpf(1) / 2, -mp.mpf(z) / 4) / mp.factorial(2 * r - 1)

    return C(1) * C(3) / (3 * C(2) ** 2)


def alpha(k, ell):
    x = k * mp.mpf(ell)
    return (1 - x * mp.cot(x)) / (2 * k * k)


def F_value(q, p, f, c, z_start):
    """Calibration (2): V(p, f, c; q), with s = (p0^-1(c)/q)^(1/2)."""
    k = mp.sqrt(q)
    z = mp.findroot(lambda zz: p0_mp(zz) - c, mp.mpf(z_start))
    s = mp.sqrt(z / q)
    M = 2 * mp.tan(k / 2) / k

    def t_of(sv):
        return (mp.tan(k * (sv - mp.mpf(1) / 2)) + mp.tan(k / 2)) / k / M

    gam = alpha(k, 1)
    t = t_of(s)
    return 1 - t + t * mp.sqrt(p * gam / alpha(k, s)) - (1 - t) * mp.sqrt(f * gam / alpha(k, 1 - s))


def bias_at(q, wv, bg, tot_w, tot_b, u, grid):
    k = math.sqrt(q)
    kq = mp.sqrt(q)
    M = 2 * math.tan(k / 2) / k
    gam = float(alpha(kq, 1))
    n = wv.n
    j = int(round((u + 0.5) / wv.h))
    assert abs(j * wv.h - (u + 0.5)) < 1e-12
    s = u + 0.5
    Lw, Lb = wv.lam1(0, j), bg.lam1(0, j)
    Rw, Rb = wv.lam1(j, n), bg.lam1(j, n)
    L2w, L2b = wv.lam2(j), bg.lam2(j)
    c0q = L2b / (3 * Lb * Lb)
    c0 = float(p0_mp(q * s * s))
    qerr = max(abs(Lb - float(alpha(kq, s))) / float(alpha(kq, s)), abs(c0q - c0) / c0)
    Rp = (Lw / tot_w) / (Lb / tot_b)
    Rf = (Rw / tot_w) / (Rb / tot_b)
    Rc = (L2w / (3 * Lw * Lw)) / c0q
    out = []
    for eta, zeta in grid:
        p0v = eta**2 * float(alpha(kq, s)) / gam
        f0v = zeta**2 * float(alpha(kq, 1 - s)) / gam
        a0v = p0v**2 * c0
        pl, fl = p0v * Rp, f0v * Rf
        al = pl**2 * c0 * Rc
        th = max(abs(pl - p0v), abs(fl - f0v), abs(al - a0v))
        t = ((math.tan(k * u) + math.tan(k / 2)) / k) / M
        vt = 1 - t + t * eta - (1 - t) * zeta
        Fl = F_value(q, mp.mpf(pl), mp.mpf(fl), mp.mpf(c0) * mp.mpf(Rc), q * s * s)
        out.append((th, abs(float(Fl - mp.mpf(vt))), (eta, zeta)))
    return out, qerr


def main():
    grid = [(e, z) for e in (0.25, 0.4, 0.6) for z in (0.4, 0.5, 0.75) if e + z <= 1 + 1e-12]
    us_full = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]
    lams_full = (1 / 40, 1 / 80, 1 / 160, 1 / 320)
    lams_edge = (1 / 640, 1 / 1280, 1 / 2560)
    for q in (0.25, 1.0):
        qerr = 0.0
        edge = {}  # lam -> (theta bias, F bias) at u = -0.3, worst over grid
        for lam in lams_full + lams_edge:
            wv, bg = Vol(q, lam, True), Vol(q, lam, False)
            tot_w, tot_b = wv.lam1(0, wv.n), bg.lam1(0, bg.n)
            L = math.log(1 / lam)
            us = us_full if lam in lams_full else [-0.3]
            th_max = F_max = F_log = 0.0
            where = None
            for u in us:
                out, e = bias_at(q, wv, bg, tot_w, tot_b, u, grid)
                qerr = max(qerr, e)
                for th, Fb, ez in out:
                    th_max = max(th_max, th / (lam**4 * (1 + L)))
                    if Fb / (lam**4 * (1 + L)) > F_max:
                        F_max, where = Fb / (lam**4 * (1 + L)), (u,) + ez
                    F_log = max(F_log, Fb / (lam**4 * L))
                if u == -0.3:
                    edge[lam] = (max(o[0] for o in out) / lam**4, max(o[1] for o in out) / lam**4)
            scope = "u grid" if lam in lams_full else "u = -0.3 only"
            print(
                f"[INFO] q={q} lam=1/{round(1 / lam)} ({scope}): max |theta bias|/(lam^4(1+L)) = "
                f"{th_max:.3e}; max |F bias|/(lam^4(1+L)) = {F_max:.3e} at (u, eta, zeta) = "
                f"{where}; max |F bias|/(lam^4 L) = {F_log:.3e}",
                flush=True,
            )
        rec(f"Q0 q={q} background quadrature", qerr < 1e-12, f"max rel err {qerr:.1e}")
        lams = lams_full + lams_edge
        Ls = [math.log(1 / x) for x in lams]
        for idx, name in ((0, "theta"), (1, "F")):
            ys = [edge[x][idx] for x in lams]
            slopes = [(ys[i + 1] - ys[i]) / (Ls[i + 1] - Ls[i]) for i in range(len(ys) - 1)]
            ok = abs(slopes[-1] - slopes[-2]) <= 0.25 * abs(slopes[-2])
            rec(
                f"B2 q={q} {name} bias at u=-0.3 is lam^4 (A + B log(1/lam))",
                ok,
                "bias/lam^4: " + ", ".join(f"{y:.4g}" for y in ys)
                + "; slopes per unit log: " + ", ".join(f"{s:.4g}" for s in slopes),
            )
    n_fail = OK.count(False)
    print(f"\n{OK.count(True)} PASS, {n_fail} FAIL")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
