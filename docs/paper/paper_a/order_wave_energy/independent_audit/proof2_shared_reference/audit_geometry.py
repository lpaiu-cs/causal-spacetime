# ruff: noqa: E501, UP031
"""Independent geometry helpers for the proof-2 (shared-reference) audit.

Everything here is derived from the Rosen metric

    ds^2 = -2 du dv + gamma_ij(u) dx^i dx^j,
    gamma = a(u)^2 R(theta) diag(e^{2 beta}, e^{-2 beta}) R(theta)^T,
    a = cos(sqrt(q) u), theta = (u/lambda + phi)/2, sinh(2 beta) = 2 lambda sqrt(q),

and from the referee's own derivation of the interval-volume hierarchy.  It does
not import or copy the author's sources.  ASCII only, numpy only.

Conventions (referee's own):
  * transverse coordinate y = x / sqrt(c), c = cosh(2 beta);
  * Mt(s,t) = m(s,t) I - eps S(I(s,t)), eps = sinh/cosh, S(w)=[[Re w, Im w],[Im w,-Re w]];
  * X < Y  iff  u_Y > u_X and v_Y - v_X >= 0.5 dy^T Mt(u_X,u_Y)^{-1} dy;
  * tips p=(u=-1/2, v=0, y=0), q=(u=+1/2, v=1, y=0);
  * Kbar(s,t) = a(s) a(t) sqrt(det Mt(s,t))   (Jacobi kernel / c);
  * Vol I(P,X) = pi eta^2 (Kbar*Kbar)/Kbar    (times a constant power of c);
  * P(two uniform events of I comparable) = Kbar C3 / (3 C2^2).
"""

from __future__ import annotations

import math

import numpy as np

A_TIP, B_TIP = -0.5, 0.5


def sym(w):
    """Trace-free symmetric matrix [[Re w, Im w], [Im w, -Re w]] (broadcasts)."""
    w = np.asarray(w, dtype=complex)
    out = np.empty(w.shape + (2, 2))
    out[..., 0, 0] = w.real
    out[..., 0, 1] = w.imag
    out[..., 1, 0] = w.imag
    out[..., 1, 1] = -w.real
    return out


def gauss(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


# ---------------------------------------------------------------------------
# Background closed forms (effective null-matter geometry), k = sqrt(q)
# ---------------------------------------------------------------------------


def c_series(r, z):
    """C_r(z) = 0F1(r+1/2; -z/4)/(2r-1)!: r-fold convolution of sin(k l)/k at l=1, z=k^2."""
    term = 1.0 / math.factorial(2 * r - 1)
    total = term
    n = 0
    while True:
        term *= -z / (2.0 * (n + 1) * (2 * r + 2 * n + 1))
        total += term
        n += 1
        if abs(term) < 1e-19 * abs(total) or n > 200:
            return total


def p0(z):
    """Background pair-comparability probability of a segment with q*length^2 = z."""
    return c_series(1, z) * c_series(3, z) / (3.0 * c_series(2, z) ** 2)


def alpha0(length, q):
    """C_2/K for a background segment of u-length `length`."""
    z = q * length * length
    return length * length * c_series(2, z) / c_series(1, z)


def t0(s, q):
    k = math.sqrt(q)
    return (math.tan(k * (s - 0.5)) + math.tan(k / 2)) / (2 * math.tan(k / 2))


def invert_p0(c, q):
    """Solve p0(q s^2) = c for s in (0,1) by bisection (p0 is decreasing)."""
    lo, hi = 0.0, q
    if not (p0(hi) < c < p0(lo)):
        raise ValueError("chain probability outside the calibrated range")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if p0(mid) > c:
            lo = mid
        else:
            hi = mid
    return math.sqrt(0.5 * (lo + hi) / q)


def calibration_v(p, f, c, q):
    """V(p,f,c;q) of null_coordinate_challenge.md eq. (3), referee's float version."""
    s = invert_p0(c, q)
    gamma = alpha0(1.0, q)
    t = t0(s, q)
    return (
        1
        - t
        + t * math.sqrt(p * gamma / alpha0(s, q))
        - (1 - t) * math.sqrt(f * gamma / alpha0(1 - s, q))
    )


def calibration_f(p, f, a, q):
    """F(p,f,a) = V(p,f,a/p^2;q): eq. (2) of local_null_challenge.md."""
    return calibration_v(p, f, a / (p * p), q)


# ---------------------------------------------------------------------------
# Exact circular-polarised vacuum wave
# ---------------------------------------------------------------------------


class Wave:
    def __init__(self, q, lam, phase=0.37, background=False, gl=12):
        self.q, self.k, self.lam, self.phase = q, math.sqrt(q), lam, phase
        self.c = math.sqrt(1 + 4 * q * lam * lam)
        self.eps = 0.0 if background else 2 * lam * self.k / self.c
        npan = max(64, int(math.ceil(1.0 / lam)))
        self.h = 1.0 / npan
        self.edges = A_TIP + self.h * np.arange(npan + 1)
        self.gx, self.gw = gauss(gl)
        nodes = self.edges[:-1, None] + self.h * self.gx[None, :]
        self.cum = np.concatenate(
            [[0.0], np.cumsum(self.h * (self._integrand(nodes) @ self.gw))]
        )
        self.Mtot = 2 * math.tan(self.k / 2) / self.k

    def _integrand(self, w):
        return np.exp(1j * w / self.lam) / np.cos(self.k * w) ** 2

    def prim(self, u):
        """int_{-1/2}^u sec^2(k w) exp(i w / lambda) dw, composite Gauss-Legendre."""
        u = np.asarray(u, dtype=float)
        flat = u.ravel()
        out = np.empty(flat.shape, dtype=complex)
        step = 200000
        for start in range(0, flat.size, step):
            part = flat[start : start + step]
            j = np.clip(
                np.floor((part - A_TIP) / self.h).astype(int), 0, len(self.edges) - 2
            )
            left = self.edges[j]
            d = part - left
            nodes = left[:, None] + d[:, None] * self.gx[None, :]
            out[start : start + step] = self.cum[j] + d * (
                self._integrand(nodes) @ self.gw
            )
        return out.reshape(u.shape)

    # scalar pieces -------------------------------------------------------
    def m(self, s, t):
        k = self.k
        return np.sin(k * (t - s)) / (k * np.cos(k * s) * np.cos(k * t))

    def osc(self, s, t):
        """Complex representative of J(s,t) (includes the global phase)."""
        return (self.prim(t) - self.prim(s)) * np.exp(1j * self.phase)

    def mt(self, s, t):
        """Effective action matrix Mt(s,t) in y coordinates (broadcasts)."""
        m = np.asarray(self.m(s, t), dtype=float)
        if self.eps == 0.0:
            return m[..., None, None] * np.eye(2)
        return m[..., None, None] * np.eye(2) - self.eps * sym(self.osc(s, t))

    def kbar(self, s, t):
        s, t = np.broadcast_arrays(np.asarray(s, float), np.asarray(t, float))
        m = self.m(s, t)
        det = m * m - self.eps**2 * np.abs(self.osc(s, t)) ** 2 if self.eps else m * m
        return np.cos(self.k * s) * np.cos(self.k * t) * np.sqrt(np.maximum(det, 0.0))

    def nodes_for(self, length, factor=1.0):
        return int(factor * (math.ceil(0.75 * length / self.lam) + 96))

    # chain hierarchy -------------------------------------------------------
    def c2(self, s0, s1, factor=1.0):
        x, w = gauss(self.nodes_for(s1 - s0, factor))
        u = s0 + (s1 - s0) * x
        return float((s1 - s0) * np.sum(w * self.kbar(s0, u) * self.kbar(u, s1)))

    def c3(self, s0, s1, factor=1.0):
        n = self.nodes_for(s1 - s0, factor)
        x, w = gauss(n)
        w1 = s0 + (s1 - s0) * x
        rest = s1 - w1
        w2 = w1[:, None] + rest[:, None] * x[None, :]
        inner = (self.kbar(w1[:, None], w2) * self.kbar(w2, s1)) @ w
        return float((s1 - s0) * np.sum(w * self.kbar(s0, w1) * rest * inner))

    def alpha(self, s0, s1, factor=1.0):
        return self.c2(s0, s1, factor) / float(self.kbar(s0, s1))

    def pair_probability(self, s0, s1, factor=1.0):
        return (
            float(self.kbar(s0, s1))
            * self.c3(s0, s1, factor)
            / (3.0 * self.c2(s0, s1, factor) ** 2)
        )

    # slice map and causal relation -------------------------------------------
    def w0(self, u):
        m = float(self.m(A_TIP, u))
        return m * (self.Mtot - m) / self.Mtot

    def slice_matrix(self, u):
        """(W0 Q_lambda)^(-1/2): canonical z -> actual y on the slice u."""
        q_mat = np.linalg.inv(self.mt(A_TIP, u)) + np.linalg.inv(self.mt(u, B_TIP))
        values, vectors = np.linalg.eigh(self.w0(u) * q_mat)
        return (vectors * values**-0.5) @ vectors.T

    def from_canonical(self, u, vt, z):
        """Canonical (u, v_tilde, z) -> wave coordinates (u, v, y1, y2)."""
        z = np.asarray(z, float)
        y = self.slice_matrix(u) @ z
        m = float(self.m(A_TIP, u))
        eta = vt - z @ z / (2 * m)
        v = eta + 0.5 * y @ np.linalg.solve(self.mt(A_TIP, u), y)
        return np.array([u, v, y[0], y[1]])

    def slacks(self, point):
        u, v, y = point[0], point[1], np.asarray(point[2:], float)
        eta = v - 0.5 * y @ np.linalg.solve(self.mt(A_TIP, u), y)
        zeta = 1 - v - 0.5 * y @ np.linalg.solve(self.mt(u, B_TIP), y)
        return float(eta), float(zeta)

    def v_tilde(self, point):
        """Background null coordinate of the slack-preserving image: 1-t+t eta-(1-t) zeta."""
        eta, zeta = self.slacks(point)
        t = t0(point[0] + 0.5, self.q)
        return 1 - t + t * eta - (1 - t) * zeta

    def margin(self, first, second):
        """Exact causal margin of first -> second (positive means first < second)."""
        dy = np.asarray(second[2:], float) - np.asarray(first[2:], float)
        return float(
            second[1]
            - first[1]
            - 0.5 * dy @ np.linalg.solve(self.mt(first[0], second[0]), dy)
        )

    def f_matrix_complex(self, u):
        """Complex representative of F_lambda(u) of null_coordinate_challenge.md eq. (4)."""
        m = float(self.m(A_TIP, u))
        e_u = -self.eps * complex(self.osc(A_TIP, u))
        e_b = -self.eps * complex(self.osc(A_TIP, B_TIP))
        return (e_u - (m / self.Mtot) * e_b) / (m * (self.Mtot - m))

    # exact population moments ---------------------------------------------
    def segment_data(self, u, factor=1.0):
        """(alpha(a,u), alpha(u,b), pair probability of [a,u]); depends on u only."""
        return (
            self.alpha(A_TIP, u, factor),
            self.alpha(u, B_TIP, factor),
            self.pair_probability(A_TIP, u, factor),
        )

    def gamma(self, factor=1.0):
        return self.alpha(A_TIP, B_TIP, factor)

    def moments(self, point, seg=None, gamma=None, factor=1.0):
        """Exact population (p, f, a) of an event: a = p^2 c."""
        eta, zeta = self.slacks(point)
        if min(eta, zeta) <= 0:
            raise ValueError("event outside the diamond")
        left, right, chain = (
            seg if seg is not None else self.segment_data(point[0], factor)
        )
        gamma = self.gamma(factor) if gamma is None else gamma
        p = left * eta**2 / gamma
        f = right * zeta**2 / gamma
        return np.array([p, f, p * p * chain])

    # iid physical-volume sampling (density a(u)^2 in (u,v,y); c cancels) ------
    def action(self, s, t, dy):
        """0.5 dy^T Mt(s,t)^{-1} dy, vectorised; Mt^{-1} = (m I + eps S(w))/(m^2-eps^2|w|^2)."""
        m = self.m(s, t)
        d1, d2 = dy[..., 0], dy[..., 1]
        if self.eps == 0.0:
            return 0.5 * (d1 * d1 + d2 * d2) / m
        w = self.eps * self.osc(s, t)
        num = (
            m * (d1 * d1 + d2 * d2)
            + w.real * (d1 * d1 - d2 * d2)
            + 2 * w.imag * d1 * d2
        )
        return 0.5 * num / (m * m - np.abs(w) ** 2)

    @staticmethod
    def _act(m, w, d1, d2):
        """0.5 d^T (m I - S(w))^{-1} d for arrays m (real), w (complex)."""
        num = (
            m * (d1 * d1 + d2 * d2)
            + w.real * (d1 * d1 - d2 * d2)
            + 2 * w.imag * d1 * d2
        )
        return 0.5 * num / (m * m - np.abs(w) ** 2)

    def sample(self, n, rng, batch=400000):
        """Rejection sampling from the normalised physical volume of the diamond."""
        half = math.sqrt(self.Mtot * (1 + self.eps) / 2.0)
        end = complex(self.prim(np.array(B_TIP))) * np.exp(1j * self.phase)
        out = []
        have = 0
        while have < n:
            u = rng.uniform(A_TIP, B_TIP, batch)
            v = rng.uniform(0.0, 1.0, batch)
            y = rng.uniform(-half, half, (batch, 2))
            keep = rng.uniform(0.0, 1.0, batch) < np.cos(self.k * u) ** 2
            u, v, y = u[keep], v[keep], y[keep]
            m1 = self.m(A_TIP, u)
            m2 = self.Mtot - m1
            r2 = np.sum(y * y, axis=1)
            # necessary condition (largest eigenvalue of Mt is at most m (1+eps)):
            pre = (r2 <= 2 * (1 + self.eps) * v * m1) & (
                r2 <= 2 * (1 + self.eps) * (1 - v) * m2
            )
            u, v, y, m1, m2 = u[pre], v[pre], y[pre], m1[pre], m2[pre]
            if self.eps == 0.0:
                w1 = np.zeros(len(u), dtype=complex)
                w2 = w1
            else:
                pu = self.prim(u) * np.exp(1j * self.phase)
                w1 = self.eps * pu
                w2 = self.eps * (end - pu)
            eta = v - self._act(m1, w1, y[:, 0], y[:, 1])
            zeta = 1 - v - self._act(m2, w2, y[:, 0], y[:, 1])
            inside = (eta > 0) & (zeta > 0)
            out.append(np.column_stack([u, v, y])[inside])
            have += int(inside.sum())
        return np.concatenate(out)[:n]

    def precedes(self, pts, point):
        """Boolean mask: pts[i] < point (exact relation)."""
        ok = pts[:, 0] < point[0]
        res = np.zeros(len(pts), dtype=bool)
        sub = pts[ok]
        dy = point[2:] - sub[:, 2:]
        act = self.action(sub[:, 0], np.full(len(sub), point[0]), dy)
        res[ok] = (point[1] - sub[:, 1]) >= act
        return res

    def follows(self, pts, point):
        """Boolean mask: point < pts[i]."""
        ok = pts[:, 0] > point[0]
        res = np.zeros(len(pts), dtype=bool)
        sub = pts[ok]
        dy = sub[:, 2:] - point[2:]
        act = self.action(np.full(len(sub), point[0]), sub[:, 0], dy)
        res[ok] = (sub[:, 1] - point[1]) >= act
        return res

    def comparable_pairs(self, pts, chunk=400):
        """Number of comparable unordered pairs among pts (exact relation)."""
        order = np.argsort(pts[:, 0])
        pts = pts[order]
        prim = self.prim(pts[:, 0]) * np.exp(1j * self.phase)
        n = len(pts)
        total = 0
        for start in range(0, n, chunk):
            a = pts[start : start + chunk]
            pa = prim[start : start + chunk]
            b = pts[start:]
            pb = prim[start:]
            m = self.m(a[:, 0][:, None], b[:, 0][None, :])
            w = self.eps * (pb[None, :] - pa[:, None])
            # Mt = m I - S(w): inverse = (m I + S(w)) / (m^2 - |w|^2)
            dy1 = b[None, :, 2] - a[:, None, 2]
            dy2 = b[None, :, 3] - a[:, None, 3]
            det = m * m - np.abs(w) ** 2
            quad = (
                m * (dy1 * dy1 + dy2 * dy2)
                + w.real * (dy1 * dy1 - dy2 * dy2)
                + 2 * w.imag * dy1 * dy2
            )
            dv = b[None, :, 1] - a[:, None, 1]
            later = b[None, :, 0] > a[:, None, 0]
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = later & (det > 0) & (dv * det >= 0.5 * quad)
            total += int(np.count_nonzero(rel))
        return total
