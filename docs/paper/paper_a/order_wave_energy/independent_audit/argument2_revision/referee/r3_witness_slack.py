# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written in this folder by the fresh-context AI referee of local_null_challenge_r2.md.
"""Referee check r3: Lemmas 3.1-3.3 on ADVERSARIAL witness pairs.

checks.py N8 samples W_lam at random plus 16 corners. Here the pairs are
chosen at the worst edges of W_lam: |z_1| in {rho_1, rho_2}, the direction
of z_1 at the edges +/- pi/6 of both admissible arcs, the largest gap
v~_2 - v~_1 = -g_2 lam^2, |w| = c' lam^(3/2) with 24 directions, eta in
{0.3, 0.4}, tau on five points of [pi/2, 3 pi/2], u_1 on nine points of
[-1/4, 1/4]. Code written for this review (no helper of checks.py is used);
the formulas are [GFL] Section 1 and Lemma E(i), [WCR] Theorem 2.

  S0  constants c_L, g_1, g_2, c' at q = 1 and 0.25 against the text.
  S1  the canonical map preserves the slacks: eta_lam(xi) = eta_0(X~) and
      zeta_lam(xi) = zeta_0(X~), relative error < 1e-9.
  S2  Lemma 3.3(i),(ii): at lam in {1/2048, 1/8192, 1/32768} every pair is
      related, both points are in K, and min kappa/lam^2 >= (3/64) c_L rho1^2.
      lam = 1/512 is information.
  S3  Lemma 3.2: the action A <= 2 |dy|^2 / (3 tau lam) on every pair, and
      |kappa - kappa_lead| / lam^3 is bounded (printed; must not grow by more
      than a factor 2 from 1/2048 to 1/32768), where
      kappa_lead = dv~ - (lam^2/2) z1^T [l] z1 - A.
Run: python r3_witness_slack.py
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

XG, WG = np.polynomial.legendre.leggauss(12)
OK = []
RHO1, RHO2 = 0.1, 0.2


def rec(name, ok, msg):
    OK.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def tm(w):
    return np.array([[w.real, w.imag], [w.imag, -w.real]])


class Model:
    def __init__(self, q, lam, phi):
        self.q, self.k, self.lam, self.phi = q, math.sqrt(q), lam, phi
        self.c = math.sqrt(1 + 4 * q * lam * lam)
        self.pref = 2 * lam * self.k / self.c * np.exp(1j * phi)
        self.M = self.m(-0.5, 0.5)
        h = lam / 8
        n = int(round(1 / h))
        self.edges = np.linspace(-0.5, 0.5, n + 1)
        lo, hi = self.edges[:-1, None], self.edges[1:, None]
        x = 0.5 * (hi - lo) * XG + 0.5 * (hi + lo)
        w = 0.5 * (hi - lo) * WG
        self.cum = np.concatenate([[0], np.cumsum((w * self.g(x)).sum(axis=1))])

    def f(self, u):
        return 1 / np.cos(self.k * u) ** 2

    def g(self, x):
        return self.f(x) * np.exp(1j * x / self.lam)

    def m(self, s, t):
        return (math.tan(self.k * t) - math.tan(self.k * s)) / self.k

    def seg(self, a, b):
        x = 0.5 * (b - a) * XG + 0.5 * (b + a)
        return complex(np.sum(0.5 * (b - a) * WG * self.g(x)))

    def G(self, x):
        j = min(int(np.searchsorted(self.edges, x, side="right")) - 1, len(self.edges) - 2)
        return self.cum[j] + self.seg(self.edges[j], x)

    def beta(self, s, t):
        if t - s < 4 * self.lam:
            return self.pref * self.seg(s, t)
        return self.pref * (self.G(t) - self.G(s))

    def frame(self, u):
        m = self.m(-0.5, u)
        ML = m * np.eye(2) - tm(self.beta(-0.5, u))
        MR = (self.M - m) * np.eye(2) - tm(self.beta(u, 0.5))
        W0 = m * (self.M - m) / self.M
        Q = np.linalg.inv(ML) + np.linalg.inv(MR)
        ev, V = np.linalg.eigh(W0 * Q)
        T = V @ np.diag(ev**-0.5) @ V.T
        return m, ML, MR, T

    def phys(self, u, vt, z):
        m, ML, MR, T = self.frame(u)
        y = T @ z
        v = vt - z @ z / (2 * m) + (T @ z) @ np.linalg.solve(ML, T @ z) / 2
        return v, y, (m, ML, MR, T)


def constants(q):
    k = math.sqrt(q)
    M = 2 * math.tan(k / 2) / k
    cL = 8 * math.sqrt(2) * k / (1.008 * M * M)
    return cL, cL * RHO1**2 / 16, cL * RHO1**2 / 8, math.sqrt(3 * math.pi * cL * RHO1**2 / 128)


def main():
    text = {1.0: (9.40, 0.00588, 0.0118, 0.0832), 0.25: (5.38, 0.00336, 0.00672, 0.0629)}
    for q, want in text.items():
        got = constants(q)
        ok = all(abs(a - b) <= 0.006 * abs(b) for a, b in zip(got, want))
        rec(f"S0 constants q={q}", ok, "c_L, g1, g2, c' = " + ", ".join(f"{v:.5g}" for v in got))

    us = np.linspace(-0.25, 0.25, 9)
    taus = np.linspace(math.pi / 2, 3 * math.pi / 2, 5)
    offs = (-math.pi / 6, math.pi / 6, math.pi - math.pi / 6, math.pi + math.pi / 6)
    wdirs = np.linspace(0, 2 * math.pi, 24, endpoint=False)
    for q in (0.25, 1.0):
        cL, g1, g2, cp = constants(q)
        floor = 3 / 64 * cL * RHO1**2
        lead_floor = 3 / 32 * cL * RHO1**2
        res_by_lam = {}
        for lam in (1 / 512, 1 / 2048, 1 / 8192, 1 / 32768):
            mod = Model(q, lam, 0.37)
            kmin, related, total, inK, slack_err = math.inf, 0, 0, True, 0.0
            act_ratio, resid = 0.0, 0.0
            for u1, tau in itertools.product(us, taus):
                u2 = u1 + tau * lam
                m1 = mod.m(-0.5, u1)
                ell = (
                    2j * mod.k * mod.f(u1) / (mod.c * m1 * (mod.M - m1))
                    * np.exp(1j * (u1 / lam + mod.phi)) * (np.exp(1j * tau) - 1)
                )
                psi_m = np.angle(ell) / 2 + math.pi / 2
                fr1, fr2 = mod.frame(u1), mod.frame(u2)
                b12 = mod.beta(u1, u2)
                M12 = mod.m(u1, u2) * np.eye(2) - tm(b12)
                M12i = np.linalg.inv(M12)
                for off, r1, eta in itertools.product(offs, (RHO1, RHO2), (0.3, 0.4)):
                    z1 = r1 * np.array([math.cos(psi_m + off), math.sin(psi_m + off)])
                    vt1 = eta + z1 @ z1 / (2 * m1)
                    for th in wdirs:
                        w = cp * lam**1.5 * np.array([math.cos(th), math.sin(th)])
                        z2 = z1 + w
                        vt2 = vt1 - g2 * lam * lam
                        out = []
                        for (u, vt, z, (m, ML, MR, T)) in ((u1, vt1, z1, fr1), (u2, vt2, z2, fr2)):
                            y = T @ z
                            v = vt - z @ z / (2 * m) + y @ np.linalg.solve(ML, y) / 2
                            eta_l = v - y @ np.linalg.solve(ML, y) / 2
                            zeta_l = 1 - v - y @ np.linalg.solve(MR, y) / 2
                            eta0 = vt - z @ z / (2 * m)
                            zeta0 = 1 - vt - z @ z / (2 * (mod.M - m))
                            slack_err = max(
                                slack_err,
                                abs(eta_l - eta0) / eta0, abs(zeta_l - zeta0) / zeta0,
                            )
                            inK = inK and abs(u) <= 0.3 and eta0 >= 0.25 and zeta0 >= 0.4
                            out.append((v, y))
                        (v1, y1), (v2, y2) = out
                        dy = y2 - y1
                        A = dy @ M12i @ dy / 2
                        kap = (v2 - v1) - A
                        related += kap > 0
                        total += 1
                        kmin = min(kmin, kap / lam**2)
                        act_ratio = max(act_ratio, A / (2 * (dy @ dy) / (3 * tau * lam)))
                        lead = (vt2 - vt1) - lam**2 / 2 * z1 @ tm(ell) @ z1 - A
                        resid = max(resid, abs(kap - lead) / lam**3)
            res_by_lam[lam] = resid
            msg = (
                f"lam=1/{round(1 / lam)}: related {related}/{total}, min kappa/lam^2 = "
                f"{kmin:.6f} (floor {floor:.6f}, leading-order floor {lead_floor:.6f}), "
                f"both in K: {inK}; slack err {slack_err:.1e}; max A/bound {act_ratio:.3f}; "
                f"max |kappa - lead|/lam^3 = {resid:.2f}"
            )
            if lam == 1 / 512:
                print(f"[INFO] S2 q={q} {msg}", flush=True)
                continue
            rec(
                f"S1+S2+S3a q={q}",
                related == total and inK and kmin >= floor and slack_err < 1e-9 and act_ratio <= 1,
                msg,
            )
        rec(
            f"S3b q={q} residual O(lam^3)",
            res_by_lam[1 / 32768] <= 2 * res_by_lam[1 / 2048],
            ", ".join(f"1/{round(1 / k)}: {v:.2f}" for k, v in res_by_lam.items()),
        )
    n_fail = OK.count(False)
    print(f"\n{OK.count(True)} PASS, {n_fail} FAIL")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
