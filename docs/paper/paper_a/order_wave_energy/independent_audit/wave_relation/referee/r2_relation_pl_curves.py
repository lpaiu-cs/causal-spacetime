# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of wave_causal_relation.md.
"""Referee check R2: the causal relation tested directly on piecewise-linear curves.

(A) "only if": random piecewise-linear future causal curves, INCLUDING segments
    that pause in u (du = 0, dx = 0, dv >= 0) and segments with zero velocity,
    built so that every segment is certified causal (rigorous upper bound on
    max_u dx^T gamma(u) dx over the segment). The endpoint must satisfy
    Theorem 2. Also checks the sign of du along the curve.
(B) "if" / sharpness: minimise the v-cost over PL curves with K segments and
    fixed endpoints; the minimum must stay >= the bound and approach it.
(C) Lemma 3 (ii)-(iv) by a different quadrature: M = int gamma^-1 by adaptive
    quad_vec of the numerically inverted gamma; J and G by QUADPACK's QAWO
    (weight cos/sin). Compare M with c m I - s J and eig(E) with +-t|G|.
"""

import math

import numpy as np
from scipy import integrate, optimize

rng = np.random.default_rng(11)


def tm(psi):
    return np.array([[math.cos(psi), math.sin(psi)], [math.sin(psi), -math.cos(psi)]])


class Fam:
    def __init__(self, q, lam, phi, beta):
        self.q, self.lam, self.phi, self.beta = q, lam, phi, beta
        self.k = math.sqrt(q)
        self.c, self.s = math.cosh(2 * beta), math.sinh(2 * beta)

    def gamma(self, u):
        th = 0.5 * (u / self.lam + self.phi)
        R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
        D = np.diag([math.exp(2 * self.beta), math.exp(-2 * self.beta)])
        return math.cos(self.k * u) ** 2 * (R @ D @ R.T)

    def M(self, u1, u2):
        val, err = integrate.quad_vec(
            lambda w: np.linalg.inv(self.gamma(w)).ravel(),
            u1,
            u2,
            epsabs=1e-14,
            epsrel=1e-13,
            limit=20000,
        )
        return val.reshape(2, 2)

    def qf(self, us, dx):
        """dx^T gamma(u) dx, vectorised; gamma = a^2 (c I + s T(psi)) (Lemma 3 form),
        cross-checked against the R D R^T construction in check_qf below."""
        psi = us / self.lam + self.phi
        A = np.cos(self.k * us) ** 2
        return A * (self.c * (dx @ dx) + self.s * ((dx[0] ** 2 - dx[1] ** 2) * np.cos(psi)
                                                    + 2 * dx[0] * dx[1] * np.sin(psi)))

    def quadform_max(self, u_a, u_b, dx, ngrid=400):
        """Rigorous upper bound on max_{u in [u_a,u_b]} dx^T gamma(u) dx."""
        us = np.linspace(u_a, u_b, ngrid)
        vals = self.qf(us, dx)
        # |d/du dx^T gamma dx| <= |(a^2)'| (c+s)|dx|^2 + a^2 (s/lam) |dx|^2
        lip = (self.k * (self.c + self.s) + self.s / self.lam) * (dx @ dx)
        return float(vals.max() + lip * (u_b - u_a) / (ngrid - 1))


def theorem2(fam, X, Y):
    """Return (holds, margin) for X <= Y by Theorem 2 with M by quad_vec."""
    u1, v1, x1 = X[0], X[1], X[2:]
    u2, v2, x2 = Y[0], Y[1], Y[2:]
    dv, dx = v2 - v1, x2 - x1
    if u1 < u2:
        act = dx @ np.linalg.solve(fam.M(u1, u2), dx)
        return 2 * dv >= act, 2 * dv - act
    if u1 == u2:
        return bool(np.all(dx == 0)) and dv >= 0, None
    return False, None


# check the vectorised quadratic form against the R D R^T construction
_f = Fam(4.0, 0.1, 0.7, 0.9)
_dx = np.array([0.3, -1.2])
_us = np.linspace(-0.5, 0.5, 11)
assert np.max(np.abs(_f.qf(_us, _dx) - np.array([_dx @ _f.gamma(u) @ _dx for u in _us]))) < 1e-12

# ---------------- (A) random PL causal curves with pauses
fams = [
    Fam(1.0, 1 / 4, 0.3, 0.5 * math.asinh(2 * 0.25 * 1.0)),
    Fam(9.0, 1 / 8, 2.0, 0.5 * math.asinh(2 / 8 * 3.0)),
    Fam(4.0, 1 / 16, 5.0, 1.3),  # non-vacuum beta: Theorem 2 holds for any beta
    Fam(0.25, 1.0, 1.0, 0.0),  # background
]
n_curves = 0
n_pause = 0
n_zero = 0
min_ratio = np.inf
du_neg = 0
fails = 0
for fam in fams:
    for _ in range(150):
        nseg = rng.integers(2, 12)
        u = rng.uniform(-0.5, 0.0)
        X = np.array([u, rng.uniform(-1, 1), *rng.normal(size=2)])
        P = X.copy()
        for _k in range(nseg):
            kind = rng.uniform()
            if kind < 0.25:  # pause in u: velocity along d_v only
                P = P + np.array([0.0, rng.exponential(0.1), 0.0, 0.0])
                n_pause += 1
            elif kind < 0.35:  # zero-velocity segment
                n_zero += 1
            else:
                du = rng.uniform(0.0, 0.5 - P[0]) * rng.uniform(0.05, 0.4)
                if du <= 0:
                    continue
                dx = rng.normal(size=2) * rng.uniform(0.01, 1.0) * math.sqrt(du)
                qmax = fam.quadform_max(P[0], P[0] + du, dx)
                slack = rng.choice([0.0, rng.exponential(0.05)])
                dv = qmax / (2 * du) + slack  # 2 du dv >= dx^T gamma dx on the segment
                P = P + np.array([du, dv, *dx])
        n_curves += 1
        if P[0] < X[0]:
            du_neg += 1
        ok, margin = theorem2(fam, X, P)
        if not ok:
            fails += 1
        if margin is not None and P[0] > X[0]:
            dxt = P[2:] - X[2:]
            act = dxt @ np.linalg.solve(fam.M(X[0], P[0]), dxt)
            if act > 1e-12:
                min_ratio = min(min_ratio, 2 * (P[1] - X[1]) / act)
print(
    f"(A) {n_curves} PL causal curves ({n_pause} pause segments, {n_zero} zero-velocity "
    f"segments): {fails} violate Theorem 2; du<0 endpoints: {du_neg}; "
    f"min 2dv/act = {min_ratio:.6f}"
)


# ---------------- (B) optimise v-cost over PL curves with K segments
def pl_cost(fam, us, xs, ngrid=200):
    tot = 0.0
    for k in range(len(us) - 1):
        du = us[k + 1] - us[k]
        dx = xs[k + 1] - xs[k]
        uu = np.linspace(us[k], us[k + 1], ngrid)
        tot += fam.qf(uu, dx).max() / (2 * du)
    return tot  # grid max <= true max: this UNDER-estimates the true minimal dv


print("(B) min over PL curves of 2*dv_needed / (dx^T M^-1 dx):")
for fam, (u1, u2) in ((fams[0], (-0.45, 0.4)), (fams[1], (-0.3, 0.45)), (fams[2], (-0.2, 0.3))):
    x1 = np.array([0.2, -0.1])
    x2 = np.array([0.9, 0.5])
    act = (x2 - x1) @ np.linalg.solve(fam.M(u1, u2), x2 - x1)
    row = []
    for K in (1, 2, 4, 8, 16):
        us = np.linspace(u1, u2, K + 1)
        # start from the straight line in x
        x0 = np.linspace(x1, x2, K + 1)[1:-1].ravel()

        def obj(z):
            xs = np.vstack([x1, z.reshape(-1, 2), x2])
            return pl_cost(fam, us, xs)

        if K == 1:
            row.append(2 * obj(np.zeros(0)) / act)
            continue
        res = optimize.minimize(obj, x0, method="Powell",
                                options={"maxiter": 200000, "xtol": 1e-10, "ftol": 1e-13})
        row.append(2 * res.fun / act)
    print("   q, lam, beta =", fam.q, fam.lam, round(fam.beta, 4), "K=1,2,4,8,16:", ["%.5f" % r for r in row], flush=True)

# ---------------- (C) Lemma 3 by a different quadrature (QAWO)
worst_M = 0.0
worst_eig = 0.0
worst_tr = 0.0
for _ in range(60):
    q = rng.uniform(0.05, 9.5)
    lam = 10 ** rng.uniform(-1.8, 0.3)
    phi = rng.uniform(-6, 6)
    beta = rng.uniform(0, 1.5)
    fam = Fam(q, lam, phi, beta)
    u1, u2 = np.sort(rng.uniform(-0.5, 0.5, size=2))
    f = lambda w: 1 / math.cos(fam.k * w) ** 2  # noqa: E731
    m = integrate.quad(f, u1, u2, epsabs=1e-14, epsrel=1e-13)[0]
    Cc = integrate.quad(f, u1, u2, weight="cos", wvar=1 / lam, epsabs=1e-14, limit=500)[0]
    Ss = integrate.quad(f, u1, u2, weight="sin", wvar=1 / lam, epsabs=1e-14, limit=500)[0]
    G = np.exp(1j * phi) * (Cc + 1j * Ss)
    J = np.array([[G.real, G.imag], [G.imag, -G.real]])
    Mdir = fam.M(u1, u2)
    Mlem = fam.c * m * np.eye(2) - fam.s * J
    worst_M = max(worst_M, np.max(np.abs(Mdir - Mlem)) / np.max(np.abs(Mlem)))
    t = fam.s / fam.c
    E = Mdir / fam.c - m * np.eye(2)  # from the DIRECT M
    ev = np.sort(np.linalg.eigvalsh(E))
    e = t * abs(G)
    worst_eig = max(worst_eig, abs(ev[1] - e) / m, abs(ev[0] + e) / m)
    worst_tr = max(worst_tr, abs(np.trace(E)) / m, abs(E[0, 1] - E[1, 0]) / m)
print(f"(C) direct M vs c m I - s J: {worst_M:.1e}; eig(E) vs +-t|G| (rel to m): "
      f"{worst_eig:.1e}; trace/asym of E: {worst_tr:.1e}")
