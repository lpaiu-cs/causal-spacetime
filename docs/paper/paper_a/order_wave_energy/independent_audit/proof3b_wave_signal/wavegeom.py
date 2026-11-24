# ruff: noqa: E501, UP031
"""Referee-side exact geometry for the circular plane-wave family.

Written independently of experiments/positive_control/*.py:
  * the complex primitive is a panel Gauss-Legendre table (no ODE solver),
  * the background chain calibration uses its own 0F1 power series,
  * numpy only (runs under the system python and under the venv python).

Conventions (taken from full_order_resolution_probe.md and
null_coordinate_challenge.md, re-derived in report.md):

  a(u) = cos(k u), k = sqrt(q), u in [-1/2, 1/2], f = a^-2
  m(u) = int_{-1/2}^u f,  M = m(1/2)
  wave action matrix in the rescaled transverse chart
      Mt(s,t) = dm(s,t) I - B(s,t),
      B(s,t) <-> beta(s,t) = fac e^{i phi} int_s^t f e^{i u / lam} du,
      fac = 2 lam k / sqrt(1 + 4 q lam^2)
  a trace-free symmetric matrix [[b1, b2], [b2, -b1]] is stored as the
  complex number b1 + i b2; a symmetric matrix a I + B as the pair (a, b).
  With zeta = y1 + i y2:
      (a I + B) y      <->  a zeta + b conj(zeta)
      y^T (a I + B) y   =   a |zeta|^2 + Re(b conj(zeta)^2)
  causal relation P1 < P2 (u1 < u2):
      v2 - v1 >= (1/2) dy^T Mt(u1,u2)^-1 dy
  physical volume density proportional to a(u)^2 du dv d^2y.

ASCII only.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.polynomial.legendre import leggauss


def gauss01(n):
    """Gauss-Legendre nodes and weights on [0, 1]."""
    x, w = leggauss(int(n))
    return 0.5 * (x + 1.0), 0.5 * w


def hyp0f1_series(b, z, terms=24):
    """0F1(b; z) by its power series; used only for |z| <= 1/4."""
    z = np.asarray(z, dtype=float)
    total = np.ones_like(z)
    term = np.ones_like(z)
    for j in range(terms):
        term = term * z / ((b + j) * (j + 1.0))
        total = total + term
    return total


def sapply(a, b, zeta):
    """(a I + B) applied to the vector zeta."""
    return a * zeta + b * np.conj(zeta)


def sqf(a, b, zeta):
    """Quadratic form zeta^T (a I + B) zeta."""
    return a * np.abs(zeta) ** 2 + (b * np.conj(zeta) ** 2).real


class Background:
    """Effective conformally flat background and its order calibration."""

    def __init__(self, q):
        self.q = float(q)
        self.k = math.sqrt(q)
        self.M = 2.0 * math.tan(self.k / 2.0) / self.k
        self.gamma0 = float(self.alpha0(1.0))

    def a2(self, u):
        return np.cos(self.k * np.asarray(u, dtype=float)) ** 2

    def f(self, u):
        return 1.0 / self.a2(u)

    def fprime(self, u):
        u = np.asarray(u, dtype=float)
        return 2.0 * self.k * np.tan(self.k * u) / np.cos(self.k * u) ** 2

    def m(self, u):
        return (
            np.tan(self.k * np.asarray(u, dtype=float)) + math.tan(self.k / 2.0)
        ) / self.k

    def dm(self, s, t):
        s = np.asarray(s, dtype=float)
        t = np.asarray(t, dtype=float)
        return np.sin(self.k * (t - s)) / (
            self.k * np.cos(self.k * s) * np.cos(self.k * t)
        )

    def _z(self, length):
        return -self.q * np.asarray(length, dtype=float) ** 2 / 4.0

    def K0(self, length):
        length = np.asarray(length, dtype=float)
        return length * hyp0f1_series(1.5, self._z(length))

    def C2(self, length):
        length = np.asarray(length, dtype=float)
        return length**3 * hyp0f1_series(2.5, self._z(length)) / 6.0

    def C3(self, length):
        length = np.asarray(length, dtype=float)
        return length**5 * hyp0f1_series(3.5, self._z(length)) / 120.0

    def alpha0(self, length):
        length = np.asarray(length, dtype=float)
        z = self._z(length)
        return length**2 * hyp0f1_series(2.5, z) / (6.0 * hyp0f1_series(1.5, z))

    def p0(self, s):
        z = self._z(s)
        return (
            0.1
            * hyp0f1_series(1.5, z)
            * hyp0f1_series(3.5, z)
            / hyp0f1_series(2.5, z) ** 2
        )

    def invert_p0(self, c):
        """Solve p0(s) = c for s in (0, 1); p0 is strictly decreasing."""
        c = np.asarray(c, dtype=float)
        lo = np.zeros_like(c)
        hi = np.ones_like(c)
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            above = self.p0(mid) > c
            lo = np.where(above, mid, lo)
            hi = np.where(above, hi, mid)
        s = 0.5 * (lo + hi)
        for _ in range(3):  # Newton polish with a central difference slope
            h = 1e-5
            slope = (self.p0(s + h) - self.p0(s - h)) / (2 * h)
            s = s - (self.p0(s) - c) / slope
        return s

    def recon(self, p, f, c):
        """Observer-side state (u, m, v, R^2) from order probabilities and q."""
        s = self.invert_p0(c)
        u = s - 0.5
        m = self.m(u)
        t = m / self.M
        eta = np.sqrt(p * self.gamma0 / self.alpha0(s))
        zeta = np.sqrt(f * self.gamma0 / self.alpha0(1.0 - s))
        v = 1.0 - t + t * eta - (1.0 - t) * zeta
        rad2 = 2.0 * m * (self.M - m) / self.M * (1.0 - eta - zeta)
        return u, m, v, rad2

    def gram(self, state_a, state_x, b_ax):
        """Order-derived Gram entry g_AX of gram_fourier_challenge.md."""
        ua, ma, va, ra = state_a
        ux, mx, vx, rx = state_x
        slack0 = np.sqrt(b_ax * self.gamma0 / self.alpha0(ux - ua))
        return 0.5 * (ra + rx) - (mx - ma) * (vx - va - slack0)


class Wave(Background):
    """Exact circular vacuum wave, wavelength parameter lam = 1 / r."""

    def __init__(self, q, r, phi, ng=10):
        super().__init__(q)
        self.r = float(r)
        self.lam = 1.0 / self.r
        self.phi = float(phi)
        self.cfac = math.sqrt(1.0 + 4.0 * q * self.lam**2)
        self.fac = 2.0 * self.lam * self.k / self.cfac
        self.npan = int(math.ceil(2.0 * self.r))
        self.h = 1.0 / self.npan
        self.grid = -0.5 + self.h * np.arange(self.npan + 1)
        self.gx, self.gw = gauss01(ng)
        nodes = self.grid[:-1, None] + self.h * self.gx[None, :]
        vals = np.exp(1j * self.r * nodes) / np.cos(self.k * nodes) ** 2
        self.cum = np.concatenate(([0.0 + 0.0j], np.cumsum(vals @ self.gw * self.h)))

    def _direct(self, s, t):
        """int_s^t f e^{i r u} du by one Gauss panel (for |t - s| <= h)."""
        d = t - s
        nodes = s[..., None] + d[..., None] * self.gx
        vals = np.exp(1j * self.r * nodes) / np.cos(self.k * nodes) ** 2
        return (vals @ self.gw) * d

    def prim(self, u):
        """P(u) = int_{-1/2}^u f(x) e^{i r x} dx (metric phase not included)."""
        u = np.asarray(u, dtype=float)
        idx = np.clip(np.floor((u + 0.5) / self.h).astype(np.int64), 0, self.npan - 1)
        return self.cum[idx] + self._direct(self.grid[idx], u)

    def beta(self, s, t):
        """Complex representative of B(s, t)."""
        s, t = np.broadcast_arrays(
            np.asarray(s, dtype=float), np.asarray(t, dtype=float)
        )
        out = self.prim(t) - self.prim(s)
        close = np.abs(t - s) <= self.h
        if np.any(close):
            out = np.array(out)
            out[close] = self._direct(s[close], t[close])
        return self.fac * np.exp(1j * self.phi) * out

    def krel(self, s, t):
        """K_lam(s,t) / K_0(s,t) - 1 = sqrt(1 - |beta|^2 / dm^2) - 1."""
        z = np.abs(self.beta(s, t)) ** 2 / self.dm(s, t) ** 2
        return -z / (1.0 + np.sqrt(1.0 - z))

    def inv_action(self, s, t):
        """Mt(s,t)^-1 as the pair (a, b)."""
        d = self.dm(s, t)
        b = self.beta(s, t)
        den = d * d - np.abs(b) ** 2
        return d / den, b / den

    def margin(self, u1, v1, y1, u2, v2, y2):
        """Exact wave causal margin of physical points, u1 < u2."""
        a, b = self.inv_action(u1, u2)
        return v2 - v1 - 0.5 * sqf(a, b, y2 - y1)

    def slice(self, u):
        """Canonical slice data at u: y = T z, det T, left/right inverses."""
        u = np.asarray(u, dtype=float)
        m = self.m(u)
        mr = self.M - m
        al_, el_ = self.inv_action(-0.5 * np.ones_like(u), u)
        ar_, er_ = self.inv_action(u, 0.5 * np.ones_like(u))
        w0 = m * mr / self.M
        al = w0 * (al_ + ar_)
        e = w0 * (el_ + er_)
        ae = np.abs(e)
        lp, lm = al + ae, al - ae
        c0 = 0.5 * (lp**-0.5 + lm**-0.5)
        small = ae < 1e-150
        c1 = np.where(
            small,
            -0.5 * al**-1.5,
            (lp**-0.5 - lm**-0.5) / (2.0 * np.where(small, 1.0, ae)),
        )
        return {
            "u": u,
            "m": m,
            "mr": mr,
            "W0": w0,
            "aL": al_,
            "eL": el_,
            "aR": ar_,
            "eR": er_,
            "c0": c0,
            "tb": c1 * e,
            "detT": (lp * lm) ** -0.5,
            "relL": np.abs(el_) / al_,
            "relR": np.abs(er_) / ar_,
        }

    def to_physical(self, sl, vt, zeta):
        """Canonical (u, v_tilde, z) -> physical (v, y) on the slice sl."""
        y = sapply(sl["c0"], sl["tb"], zeta)
        eta = vt - np.abs(zeta) ** 2 / (2.0 * sl["m"])
        v = eta + 0.5 * sqf(sl["aL"], sl["eL"], y)
        return v, y

    def _krel_p(self, s, t, ps, pt):
        """krel with the primitive values at s and t supplied by the caller."""
        s, t = np.broadcast_arrays(
            np.asarray(s, dtype=float), np.asarray(t, dtype=float)
        )
        diff = np.array(np.broadcast_to(pt, s.shape) - np.broadcast_to(ps, s.shape))
        close = np.abs(t - s) <= self.h
        if np.any(close):
            diff[close] = self._direct(s[close], t[close])
        z = self.fac**2 * np.abs(diff) ** 2 / self.dm(s, t) ** 2
        return -z / (1.0 + np.sqrt(1.0 - z))

    def chain(self, s, t, pairs=False, n=None):
        """Exact wave/background ratios on [s, t].

        Returns (alpha_lam / alpha_0, c_lam / c_0) where alpha = C2 / K is
        the interval-volume coefficient and c the pair probability
        K C3 / (3 C2^2).  The second entry is None unless pairs is True.
        """
        s = float(s)
        t = float(t)
        length = t - s
        n = int(n or (math.ceil(self.r * length / 2.0) + 64))
        x, w = gauss01(n)
        xs = s + length * x
        ps, pt, pxs = self.prim(s), self.prim(t), self.prim(xs)
        k1 = self.K0(xs - s)
        k2 = self.K0(t - xs)
        d1 = self._krel_p(s, xs, ps, pxs)
        d2 = self._krel_p(xs, t, pxs, pt)
        e2 = float(np.sum(w * k1 * k2 * (d1 + d2 + d1 * d2)) * length / self.C2(length))
        dk = float(self._krel_p(s, t, ps, pt))
        alpha_ratio = (1.0 + e2) / (1.0 + dk)
        if not pairs:
            return alpha_ratio, None
        ys = np.broadcast_to(xs[None, :], (n, n))
        xx = s + (ys - s) * x[:, None]
        pxx = self.prim(xx)
        da = self._krel_p(s, xx, ps, pxx)
        db = self._krel_p(xx, ys, pxx, np.broadcast_to(pxs[None, :], (n, n)))
        dc = d2[None, :]
        change = da + db + dc + da * db + da * dc + db * dc + da * db * dc
        weight = (
            w[:, None]
            * w[None, :]
            * (ys - s)
            * self.K0(xx - s)
            * self.K0(ys - xx)
            * k2[None, :]
        )
        e3 = float(np.sum(weight * change) * length / self.C3(length))
        c_ratio = (1.0 + dk) * (1.0 + e3) / (1.0 + e2) ** 2
        return alpha_ratio, c_ratio


def bump(x):
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    inside = np.abs(x) < 1.0
    out[inside] = np.exp(-1.0 / (1.0 - x[inside] ** 2))
    return out


class Reporter:
    """PASS/FAIL bookkeeping; exit code nonzero on any FAIL."""

    def __init__(self):
        self.failed = 0

    def check(self, name, ok, detail=""):
        print(
            ("PASS " if ok else "FAIL ") + name + ((" | " + detail) if detail else ""),
            flush=True,
        )
        if not ok:
            self.failed += 1
        return ok

    def finish(self):
        print("SUMMARY failed=%d" % self.failed, flush=True)
        raise SystemExit(1 if self.failed else 0)
