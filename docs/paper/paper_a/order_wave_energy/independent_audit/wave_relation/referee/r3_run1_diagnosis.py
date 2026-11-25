# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of wave_causal_relation.md.
"""Referee check R3: audit of checks.py's W1 step rule and the run-1 diagnosis.

The finite-difference functions below are COPIED from checks.py (at 1ba2e6e)
so that the author's numerics can be audited; they are not used for any
independent claim. Questions:
  1. Does halving the run-1 step cut the error ~16x, and converge to 0 / 2q?
  2. Is the run-1 failure at q=9, lam=1, u=0.45 caused by the metric samples
     crossing the pole, or by truncation of the OUTER derivative? (The metric
     a^2 Gamma is entire in u; only gamma^-1 is singular, and it is evaluated
     only at Christoffel centres.) Test: apply the same outer 5-point stencil
     to EXACT Christoffels (closed form) and compare.
  3. What do d/100, d/200, d/400 give on the W1 grid? Is the d/200 pass a
     convergence plateau or a just-passing tuning?
"""

import math

import numpy as np


def tmat(psi):
    return np.array([[math.cos(psi), math.sin(psi)], [math.sin(psi), -math.cos(psi)]])


def rot(th):
    return np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])


class Family:  # copied
    def __init__(self, q, lam, phi, s=None, background=False):
        self.q, self.lam, self.phi = q, lam, phi
        self.rq = math.sqrt(q)
        self.s = 2 * lam * self.rq if s is None else s
        self.c = math.sqrt(1 + self.s**2)
        self.t = self.s / self.c
        self.beta = 0.5 * math.asinh(self.s)
        self.background = background

    def a(self, u):
        return math.cos(self.rq * u)

    def gamma(self, u):
        if self.background:
            return self.a(u) ** 2 * np.eye(2)
        r = rot(0.5 * (u / self.lam + self.phi))
        d = np.diag([math.exp(2 * self.beta), math.exp(-2 * self.beta)])
        return self.a(u) ** 2 * r @ d @ r.T

    def metric(self, p):
        g = np.zeros((4, 4))
        g[0, 1] = g[1, 0] = -1.0
        g[2:, 2:] = self.gamma(p[0])
        return g


def d1(fun, p, k, h):  # copied
    e = np.zeros(4)
    e[k] = h
    return (-fun(p + 2 * e) + 8 * fun(p + e) - 8 * fun(p - e) + fun(p - 2 * e)) / (12 * h)


def christoffel(fam, p, hs):  # copied
    g = fam.metric(p)
    gi = np.linalg.inv(g)
    dg = np.array([d1(fam.metric, p, k, hs[k]) for k in range(4)])
    t1 = np.einsum("bdc->dbc", dg)
    t2 = np.einsum("cdb->dbc", dg)
    t3 = dg
    return 0.5 * np.einsum("ad,dbc->abc", gi, t1 + t2 - t3)


def ricci_fd(fam, p, hs, chris=None):  # copied, with optional exact Christoffels
    chris = chris or (lambda pp: christoffel(fam, pp, hs))
    gam = chris(p)
    dgam = np.array([d1(chris, p, k, hs[k]) for k in range(4)])
    term1 = np.einsum("aabd->bd", dgam)
    tr = np.einsum("aab->b", gam)
    dtr = np.einsum("daab->bd", dgam)
    term3 = np.einsum("e,ebd->bd", tr, gam)
    term4 = np.einsum("ade,eab->bd", gam, gam)
    return term1 - dtr + term3 - term4


def exact_christoffel(fam):
    """Closed-form Christoffels of Lemma 4 (verified independently in r1)."""

    def chris(p):
        u = p[0]
        k = fam.rq
        if fam.background:
            G, Gp = np.eye(2), np.zeros((2, 2))
        else:
            psi = u / fam.lam + fam.phi
            G = fam.c * np.eye(2) + fam.s * tmat(psi)
            Gp = (fam.s / fam.lam) * np.array(
                [[-math.sin(psi), math.cos(psi)], [math.cos(psi), math.sin(psi)]]
            )
        A, Ap = math.cos(k * u) ** 2, -k * math.sin(2 * k * u)
        gam, gamp = A * G, Ap * G + A * Gp
        L = np.zeros((4, 4, 4))
        L[1, 2:, 2:] = 0.5 * gamp
        P = np.linalg.solve(gam, gamp)
        L[2:, 0, 2:] = 0.5 * P
        L[2:, 2:, 0] = 0.5 * P
        return L

    return chris


# ---- 1 and 2: the run-1 worst point
q, lam, phi, u = 9.0, 1.0, 0.0, 0.45
p = np.array([u, 0.3, -0.2, 0.7])
pole = math.pi / (2 * math.sqrt(q))
print(f"pole = {pole:.4f}, distance d = {pole - u:.4f}")
for phi in (0.0, 1.1):
    fam = Family(q, lam, phi)
    scale = 2 * q + fam.s**2 / (2 * lam**2) + math.sqrt(q) / lam
    prev = None
    for h in (lam / 40, lam / 80, lam / 160, lam / 320, (pole - u) / 200):
        hs = [h, 1e-2, 1e-2, 1e-2]
        r_fd = ricci_fd(fam, p, hs)
        r_ex = ricci_fd(fam, p, hs, chris=exact_christoffel(fam))
        err = np.max(np.abs(r_fd)) / scale
        err_ex = np.max(np.abs(r_ex)) / scale
        ratio = "" if prev is None else f" ratio {prev / err:6.2f}"
        print(
            f"phi={phi} h={h:.5f} reach(metric)={u + 4 * h:.3f}  |R|/scale FD-of-FD={err:.3e}"
            f"  outer-FD of EXACT Christoffels={err_ex:.3e}{ratio}"
        )
        prev = err

# ---- 3: the W1 grid under several step rules
print("\nW1 grid maxima (vacuum |R|/scale, control err/scale) under step rules:")
for label, rule in (
    ("run1 lam/40", lambda lam, d: lam / 40),
    ("lam/80 only", lambda lam, d: lam / 80),
    ("min(lam/80,d/100)", lambda lam, d: min(lam / 80, d / 100)),
    ("min(lam/80,d/200) [run2]", lambda lam, d: min(lam / 80, d / 200)),
    ("min(lam/80,d/400)", lambda lam, d: min(lam / 80, d / 400)),
    ("min(lam/160,d/400)", lambda lam, d: min(lam / 160, d / 400)),
):
    w1 = 0.0
    w1c = 0.0
    arg = None
    for q in (0.25, 1.0, 4.0, 9.0):
        for lam in (1.0, 1 / 8, 1 / 64):
            for phi in (0.0, 1.1):
                for u in (-0.4, 0.1, 0.45):
                    pole = math.pi / (2 * math.sqrt(q))
                    hs = [rule(lam, pole - abs(u)), 1e-2, 1e-2, 1e-2]
                    p = np.array([u, 0.3, -0.2, 0.7])
                    fam = Family(q, lam, phi)
                    scale = 2 * q + fam.s**2 / (2 * lam**2) + math.sqrt(q) / lam
                    val = float(np.max(np.abs(ricci_fd(fam, p, hs)))) / scale
                    if val > w1:
                        w1, arg = val, (q, lam, phi, u)
                    ctrl = Family(q, lam, phi, s=1.1 * fam.s)
                    sc = 2 * q + ctrl.s**2 / (2 * lam**2) + math.sqrt(q) / lam
                    rc = ricci_fd(ctrl, p, hs)
                    want = 2 * q - ctrl.s**2 / (2 * lam**2)
                    rest = rc.copy()
                    rest[0, 0] = 0
                    w1c = max(w1c, max(abs(rc[0, 0] - want), np.max(np.abs(rest))) / sc)
    print(f"  {label:28s} vacuum {w1:.2e} at {arg}   control {w1c:.2e}")

print("\nW1 background max rel err under step rules:")
for label, rule in (
    ("run1 1e-3", lambda d: 1e-3),
    ("min(1e-3,d/100)", lambda d: min(1e-3, d / 100)),
    ("min(1e-3,d/200) [run2]", lambda d: min(1e-3, d / 200)),
    ("min(1e-3,d/400)", lambda d: min(1e-3, d / 400)),
):
    bg = 0.0
    for q in (0.25, 1.0, 4.0, 9.0):
        fam = Family(q, 1.0, 0.0, background=True)
        for u in (-0.4, 0.1, 0.45):
            pole = math.pi / (2 * math.sqrt(q))
            rc = ricci_fd(fam, np.array([u, 0.1, 0.2, -0.3]), [rule(pole - abs(u)), 1e-2, 1e-2, 1e-2])
            rest = rc.copy()
            rest[0, 0] = 0
            bg = max(bg, abs(rc[0, 0] - 2 * q) / (2 * q), np.max(np.abs(rest)) / (2 * q))
    print(f"  {label:24s} {bg:.2e}")
