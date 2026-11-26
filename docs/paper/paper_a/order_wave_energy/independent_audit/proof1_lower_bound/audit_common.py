"""Shared helpers for the independent audit of quantized_order_bound.md.

Written by the referee from the definitions in full_order_resolution_probe.md
and coordinate_oracle_barrier.md.  Nothing here is imported from the author's
code.  ASCII only.

Conventions (y chart, x = sqrt(c) y):
  a(u) = cos(sqrt(q) u),  f = a^-2,  m(u) = integral_{-1/2}^u f
  s = sinh(2 beta) = 2 lambda sqrt(q),  c = cosh(2 beta) = sqrt(1 + s^2)
  Mtilde(u1,u2) = b I + E,  b = m(u2) - m(u1),
  E = -(s/c) [[Re G, Im G], [Im G, -Re G]],
  G = integral_{u1}^{u2} f(u) exp(i (u/lambda + phi)) du
  e = (s/c) |G| is the positive eigenvalue of the trace-free symmetric E.
  wave relation (u1 < u2):  2 (v2 - v1) >= dy^T (b I + E)^-1 dy
  background relation:      2 (v2 - v1) >= |dy|^2 / b
  Minkowski chart of the background: T = (m + v)/sqrt2, Z = (m - v)/sqrt2.
"""

from __future__ import annotations

import sys
from fractions import Fraction as Fr

import numpy as np

SQRT2 = float(np.sqrt(2.0))


class Reporter:
    """Collect PASS/FAIL lines; exit nonzero if any FAIL."""

    def __init__(self, title):
        self.fails = 0
        self.count = 0
        print("=== " + title + " ===")

    def check(self, name, ok, detail=""):
        self.count += 1
        tag = "PASS" if ok else "FAIL"
        if not ok:
            self.fails += 1
        line = tag + " " + name
        if detail:
            line += " | " + detail
        print(line)
        sys.stdout.flush()
        return ok

    def info(self, text):
        print("INFO " + text)
        sys.stdout.flush()

    def finish(self):
        verdict = "ALL PASS" if self.fails == 0 else "FAILURES"
        print(f"SUMMARY checks={self.count} fails={self.fails} -> {verdict}")
        sys.exit(1 if self.fails else 0)


class Family:
    """Exact circular vacuum family in the y chart, with an accurate primitive."""

    def __init__(self, q, lam, phi=0.37, slack=1.0, gl=12):
        self.q = float(q)
        self.lam = float(lam)
        self.phi = float(phi)
        self.h = float(slack)  # axial slack: tips at v=0 and v=h, y=0
        self.rq = float(np.sqrt(q))
        self.s = 2.0 * self.lam * self.rq
        self.c = float(np.sqrt(1.0 + self.s**2))
        self.tanh2b = self.s / self.c
        self.delta = self.lam**2
        self.fmax = 1.0 / np.cos(0.5 * self.rq) ** 2
        self.Mtot = 2.0 * np.tan(0.5 * self.rq) / self.rq
        # explicit constants re-derived by the referee
        #   |int f e^{iu/lam}| <= lam (f(s)+f(t)+TV_f) <= lam (4 fmax - 2)
        #   e <= C_E lam min(ell, lam),  C_E = 2 sqrt(q) (4 fmax - 2)
        #   sandwich constant c0 = 2 sqrt(2) C_E
        self.C_E = 2.0 * self.rq * (4.0 * self.fmax - 2.0)
        self.c0 = 2.0 * SQRT2 * self.C_E
        n = max(512, int(np.ceil(8.0 / self.lam)))
        self.n = n
        self.dx = 1.0 / n
        self.grid = -0.5 + self.dx * np.arange(n + 1)
        x, w = np.polynomial.legendre.leggauss(gl)
        self.glx, self.glw = x, w
        nodes = self.grid[:-1, None] + (x[None, :] + 1.0) * 0.5 * self.dx
        vals = self.f(nodes) * np.exp(1j * nodes / self.lam)
        cell = (vals * w[None, :]).sum(axis=1) * 0.5 * self.dx
        self.G = np.concatenate([[0.0 + 0.0j], np.cumsum(cell)])

    def a(self, u):
        return np.cos(self.rq * u)

    def f(self, u):
        return 1.0 / np.cos(self.rq * u) ** 2

    def m(self, u):
        """m(-1/2, u) in closed form."""
        return (np.tan(self.rq * u) + np.tan(0.5 * self.rq)) / self.rq

    def u_of_m(self, m):
        return np.arctan(self.rq * m - np.tan(0.5 * self.rq)) / self.rq

    def prim(self, u):
        """integral_{-1/2}^u f(u') exp(i u'/lambda) du' (phase phi excluded)."""
        u = np.asarray(u, dtype=float)
        k = np.clip(np.floor((u + 0.5) / self.dx).astype(np.int64), 0, self.n - 1)
        a = self.grid[k]
        half = 0.5 * (u - a)
        nodes = a[..., None] + (self.glx + 1.0) * half[..., None]
        vals = self.f(nodes) * np.exp(1j * nodes / self.lam)
        return self.G[k] + (vals * self.glw).sum(axis=-1) * half

    def pair_G(self, prim1, prim2):
        return np.exp(1j * self.phi) * (prim2 - prim1)

    def seg(self, u1, u2, panels=16):
        """integral_{u1}^{u2} f exp(i u/lambda) du, cancellation-free.

        Segments shorter than 4 lambda are integrated directly (composite
        Gauss-Legendre, phase advance <= 0.25 rad per panel); longer ones use
        the primitive.  Phase phi excluded.
        """
        u1 = np.asarray(u1, dtype=float)
        u2 = np.asarray(u2, dtype=float)
        out = self.prim(u2) - self.prim(u1)
        short = (u2 - u1) <= 4.0 * self.lam
        if np.any(short):
            a = u1[short]
            ell = u2[short] - a
            acc = np.zeros(a.shape, dtype=complex)
            hp = ell / panels
            for p in range(panels):
                left = a + p * hp
                nodes = left[:, None] + (self.glx[None, :] + 1.0) * 0.5 * hp[:, None]
                vals = self.f(nodes) * np.exp(1j * nodes / self.lam)
                acc += (vals * self.glw[None, :]).sum(axis=1) * 0.5 * hp
            out[short] = acc
        return out

    def wave_form(self, b, G, dy1, dy2):
        """dy^T (b I + E)^-1 dy by explicit 2x2 inversion (no closed form used)."""
        k = self.s / self.c
        e11 = b - k * G.real
        e22 = b + k * G.real
        e12 = -k * G.imag
        det = e11 * e22 - e12 * e12
        return (e22 * dy1 * dy1 - 2.0 * e12 * dy1 * dy2 + e11 * dy2 * dy2) / det

    # ---- reference sampling (effective background law mu_0) -----------------
    def sample_reference(self, n, rng):
        """iid draws from density prop. to a(u)^2 on the background diamond,
        coordinates (u, v, y1, y2).  Rejection sampling, exact."""
        out = []
        have = 0
        R = np.sqrt(0.5 * self.h * self.Mtot) * 1.0000001
        while have < n:
            k = max(4096, 6 * (n - have))
            u = rng.uniform(-0.5, 0.5, k)
            v = rng.uniform(0.0, self.h, k)
            y1 = rng.uniform(-R, R, k)
            y2 = rng.uniform(-R, R, k)
            r2 = y1 * y1 + y2 * y2
            mm = self.m(u)
            ok = (2.0 * v * mm >= r2) & (2.0 * (self.h - v) * (self.Mtot - mm) >= r2)
            ok &= rng.uniform(0.0, 1.0, k) < self.a(u) ** 2
            pts = np.stack([u[ok], v[ok], y1[ok], y2[ok]], axis=1)
            out.append(pts)
            have += len(pts)
        return np.concatenate(out)[:n]

    def to_minkowski(self, pts):
        """(u,v,y1,y2) -> (T,Z,y1,y2) of the background Minkowski chart."""
        mm = self.m(pts[:, 0])
        v = pts[:, 1]
        return np.stack(
            [(mm + v) / SQRT2, (mm - v) / SQRT2, pts[:, 2], pts[:, 3]], axis=1
        )

    def reference_volume(self):
        """V0 = integral a^2 over the background diamond = pi h^2 C2/K."""
        rq = self.rq
        c2 = (np.sin(rq) / rq - np.cos(rq)) / (2.0 * rq * rq)
        return np.pi * self.h**2 * c2 / (np.sin(rq) / rq)


# ----------------------------------------------------------------------------
# Exact polytope volume (Lasserre recursion, rational arithmetic)
# ----------------------------------------------------------------------------


def _normalise(cons, d):
    """Deduplicate same-direction constraints a.x <= b; None if infeasible."""
    best = {}
    for a, b in cons:
        nz = [i for i in range(d) if a[i] != 0]
        if not nz:
            if b < 0:
                return None
            continue
        sc = abs(a[nz[0]])
        key = tuple(ai / sc for ai in a)
        bb = b / sc
        if key not in best or bb < best[key]:
            best[key] = bb
    return list(best.items())


def polytope_volume(cons, d):
    """Volume of the bounded polytope {x in R^d : a.x <= b for (a,b) in cons}.

    Lasserre / divergence theorem with origin as base point:
        d * vol(P) = sum_i (b_i / |a_i|) vol_{d-1}(F_i),
        vol_{d-1}(F_i) = (|a_i| / |a_ij|) vol_{d-1}(projection of F_i along x_j).
    Exact with Fraction input.  Boundedness is the caller's responsibility
    (all uses include the cube constraints).
    """
    cons = _normalise(cons, d)
    if cons is None:
        return Fr(0)
    if d == 1:
        lo = None
        hi = None
        for (a,), b in cons:
            if a > 0:
                hi = b if hi is None else min(hi, b)
            else:
                lo = -b if lo is None else max(lo, -b)
        if lo is None or hi is None:
            raise ValueError("unbounded 1-d section")
        return hi - lo if hi > lo else Fr(0)
    total = Fr(0)
    for idx, (a, b) in enumerate(cons):
        if b == 0:
            continue
        total += b / abs(a[_pivot(a, d)]) * _facet_projected_volume(cons, idx, d)
    return total / d


def _pivot(a, d):
    for i in range(d):
        if a[i] != 0:
            return i
    raise ValueError("zero row")


def _facet_projected_volume(cons, idx, d):
    a, b = cons[idx]
    j = _pivot(a, d)
    sub = []
    for idx2, (a2, b2) in enumerate(cons):
        if idx2 == idx:
            continue
        coef = a2[j] / a[j]
        newa = tuple(a2[c] - coef * a[c] for c in range(d) if c != j)
        sub.append((newa, b2 - coef * b))
    return polytope_volume(sub, d - 1)


def cube_constraints(k):
    cons = []
    for i in range(k):
        e = [Fr(0)] * k
        e[i] = Fr(1)
        cons.append((tuple(e), Fr(1)))
        e2 = [Fr(0)] * k
        e2[i] = Fr(-1)
        cons.append((tuple(e2), Fr(0)))
    return cons


def pattern_constraints(k, rows, tau):
    """Unit cube [0,1]^k cut by sigma*(t_j - t_i) >= tau_e for rows (i,j,sigma)."""
    cons = cube_constraints(k)
    for (i, j, sg), th in zip(rows, tau, strict=False):
        a = [Fr(0)] * k
        a[j] = Fr(-sg)
        a[i] = Fr(sg)
        cons.append((tuple(a), -Fr(th)))
    return cons


def pattern_prob(k, rows, tau):
    """p(tau) for independent uniform times on [0,1] (exact Fraction)."""
    return polytope_volume(pattern_constraints(k, rows, tau), k)


def pattern_grad(k, rows, tau):
    """Exact gradient of p: d p / d tau_e = -(projected facet volume)."""
    cons = pattern_constraints(k, rows, tau)
    base = 2 * k
    grad = []
    for e in range(len(rows)):
        # the row constraint must not be merged with another one: rows are
        # distinct pairs and never parallel to a cube face, so index is stable
        grad.append(-_facet_projected_volume(cons, base + e, k))
    return grad
