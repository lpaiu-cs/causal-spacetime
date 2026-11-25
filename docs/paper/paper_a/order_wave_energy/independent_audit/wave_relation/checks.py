"""Deterministic checks for wave_causal_relation.md.

No sampling of orders. Seeded random draws are used only to test inequalities
and identities.

The curvature and geodesic checks do not use the note's Christoffel symbols.
They take finite differences of the 4x4 metric g(u, v, x1, x2) in all four
coordinates, and use the generic formulas
Gamma^a_bc = g^ad (d_b g_dc + d_c g_db - d_d g_bc)/2 and
R_bd = d_a Gamma^a_bd - d_d Gamma^a_ab + Gamma^a_ae Gamma^e_bd
       - Gamma^a_de Gamma^e_ab.

Pass criteria (written before the full run; one smoke case was run first to
debug the script):
  W1  Ricci. At 4 q x 3 lam x 2 phases x 3 values of u, with the vacuum s:
      every component is below 1e-5 * scale. The scale is
      2q + s^2/(2 lam^2) + sqrt(q)/lam. Positive control with s -> 1.1 s:
      R_uu agrees with 2q - s^2/(2 lam^2) to 1e-5 * scale, and the other
      components stay below it. Background: R_uu = 2q.
  W2  Family algebra. gamma^-1 and det gamma = a^4 against numerical
      inversion. M by quadrature of the numerical inverse against
      c m I - s J. y-chart identity and eigenvalues +-t|G|. Phase rotation
      identity. sqrt|det g| = a^2. All to relative 1e-9 (M to 1e-8).
  W3  Boundary by geodesic shooting. Null geodesics (finite-difference
      Christoffels, RK4) from random X reach points with
      2 Delta v = Delta x^T M^-1 Delta x to relative 1e-6, and u stays
      affine (du/dtau = 1 to 1e-8).
  W4  The inequality. For random fixed-endpoint paths x(u),
      int x'^T gamma x' du >= Delta x^T M^-1 Delta x (ratio >= 1 - 1e-9). The
      optimal path attains it (ratio 1 to 1e-8). The identity behind F' >= 0
      holds on random vectors to 1e-12. Lemma 1's sign conclusions hold on
      random future causal vectors.
  W5  Background. g_0 = Phi^*(a^2 eta) on random points to 1e-9; m from
      quadrature against the closed form. On random pairs, including
      Delta m = 0 and near-null pairs, Theorem 2 with gamma = a^2 I and
      Delta T >= |Delta z| classify every pair the same way. Density
      Jacobian a^4.
  W6  Chart. On random pairs, the x-chart relation with M and the y-chart
      relation with mI + E classify every pair the same way.

Run: python checks.py
"""

from __future__ import annotations

# The lambdas passed to Family.quad are evaluated inside the same loop iteration,
# so the late-binding warning B023 does not apply; zip pairs equal-length arrays.
# ruff: noqa: B023, B905
import math

import numpy as np

failed = 0
npass = 0


def check(name, ok, detail=""):
    global failed, npass
    line = ("PASS " if ok else "FAIL ") + name + (" | " + detail if detail else "")
    print(line, flush=True)
    if ok:
        npass += 1
    else:
        failed += 1


GX, GW = np.polynomial.legendre.leggauss(16)


def tmat(psi):
    return np.array([[math.cos(psi), math.sin(psi)], [math.sin(psi), -math.cos(psi)]])


def rot(th):
    return np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])


class Family:
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

    def quad(self, func, u1, u2, width):
        npan = max(8, math.ceil((u2 - u1) / width))
        ed = np.linspace(u1, u2, npan + 1)
        tot = 0.0
        for lo, hi in zip(ed[:-1], ed[1:]):
            half, mid = (hi - lo) / 2, (hi + lo) / 2
            for xg, wg in zip(GX, GW):
                tot = tot + half * wg * func(mid + half * xg)
        return tot


# ---------------------------------------------------------------- generic geometry
def d1(fun, p, k, h):
    e = np.zeros(4)
    e[k] = h
    return (-fun(p + 2 * e) + 8 * fun(p + e) - 8 * fun(p - e) + fun(p - 2 * e)) / (
        12 * h
    )


def christoffel(fam, p, hs):
    g = fam.metric(p)
    gi = np.linalg.inv(g)
    dg = np.array([d1(fam.metric, p, k, hs[k]) for k in range(4)])  # dg[k][a,b]
    # Gamma^a_bc = 1/2 g^ad (d_b g_dc + d_c g_db - d_d g_bc)
    t1 = np.einsum("bdc->dbc", dg)  # d_b g_dc  indexed [d,b,c]
    t2 = np.einsum("cdb->dbc", dg)  # d_c g_db
    t3 = dg  # d_d g_bc  indexed [d,b,c]
    return 0.5 * np.einsum("ad,dbc->abc", gi, t1 + t2 - t3)


def ricci(fam, p, hs):
    gam = christoffel(fam, p, hs)
    dgam = np.array(
        [d1(lambda pp: christoffel(fam, pp, hs), p, k, hs[k]) for k in range(4)]
    )
    # dgam[k][a,b,c] = d_k Gamma^a_bc
    term1 = np.einsum("aabd->bd", dgam)  # d_a Gamma^a_bd
    tr = np.einsum("aab->b", gam)  # Gamma^a_ab
    dtr = np.einsum("daab->bd", dgam)  # d_d Gamma^a_ab, indexed [b, d]
    term3 = np.einsum("e,ebd->bd", tr, gam)
    term4 = np.einsum("ade,eab->bd", gam, gam)
    return term1 - dtr + term3 - term4


# ================================================================ W1 Ricci
rng = np.random.default_rng(20261126)
w1_worst = 0.0
w1_ctrl_worst = 0.0
for q in (0.25, 1.0, 4.0, 9.0):
    for lam in (1.0, 1 / 8, 1 / 64):
        for phi in (0.0, 1.1):
            for u in (-0.4, 0.1, 0.45):
                hs = [lam / 40, 1e-2, 1e-2, 1e-2]
                p = np.array([u, 0.3, -0.2, 0.7])
                fam = Family(q, lam, phi)
                scale = 2 * q + fam.s**2 / (2 * lam**2) + math.sqrt(q) / lam
                ric = ricci(fam, p, hs)
                w1_worst = max(w1_worst, float(np.max(np.abs(ric))) / scale)
                ctrl = Family(q, lam, phi, s=1.1 * fam.s)
                scale_c = 2 * q + ctrl.s**2 / (2 * lam**2) + math.sqrt(q) / lam
                rc = ricci(ctrl, p, hs)
                want = 2 * q - ctrl.s**2 / (2 * lam**2)
                rest = rc.copy()
                rest[0, 0] = 0.0
                err = max(abs(rc[0, 0] - want), float(np.max(np.abs(rest))))
                w1_ctrl_worst = max(w1_ctrl_worst, err / scale_c)
check(
    "W1 Ricci = 0 for the vacuum family (all components)",
    w1_worst < 1e-5,
    f"max |R|/scale = {w1_worst:.2e}",
)
check(
    "W1 positive control s -> 1.1 s: R_uu = 2q - s^2/(2 lam^2), rest 0",
    w1_ctrl_worst < 1e-5,
    f"max err/scale = {w1_ctrl_worst:.2e}",
)
bg_worst = 0.0
for q in (0.25, 1.0, 4.0, 9.0):
    fam = Family(q, 1.0, 0.0, background=True)
    for u in (-0.4, 0.1, 0.45):
        rc = ricci(fam, np.array([u, 0.1, 0.2, -0.3]), [1e-3, 1e-2, 1e-2, 1e-2])
        rest = rc.copy()
        rest[0, 0] = 0.0
        bg_worst = max(
            bg_worst,
            abs(rc[0, 0] - 2 * q) / (2 * q),
            float(np.max(np.abs(rest))) / (2 * q),
        )
check(
    "W1 background R_uu = 2q, other components 0",
    bg_worst < 1e-6,
    f"max rel err = {bg_worst:.2e}",
)

# ================================================================ W2 family algebra
w2 = dict(inv=0.0, det=0.0, M=0.0, ychart=0.0, eig=0.0, rotid=0.0, vol=0.0)
for q in (0.25, 1.0, 4.0, 9.0):
    for lam in (1.0, 1 / 8, 1 / 64):
        for phi in (0.0, 1.1, 4.0):
            fam = Family(q, lam, phi)
            for u in (-0.45, 0.0, 0.3):
                gm = fam.gamma(u)
                closed = fam.a(u) ** -2 * (
                    fam.c * np.eye(2) - fam.s * tmat(u / lam + phi)
                )
                w2["inv"] = max(
                    w2["inv"],
                    float(np.max(np.abs(np.linalg.inv(gm) - closed)))
                    / float(np.max(np.abs(closed))),
                )
                w2["det"] = max(w2["det"], abs(np.linalg.det(gm) / fam.a(u) ** 4 - 1))
                gfull = fam.metric(np.array([u, 0, 0, 0]))
                w2["vol"] = max(
                    w2["vol"],
                    abs(math.sqrt(abs(np.linalg.det(gfull))) / fam.a(u) ** 2 - 1),
                )
            u1, u2 = -0.37, 0.41
            width = min(lam / 8, 0.02)
            m_num = fam.quad(lambda w: np.linalg.inv(fam.gamma(w)), u1, u2, width)
            f_int = fam.quad(lambda w: fam.a(w) ** -2, u1, u2, width)
            j_int = fam.quad(
                lambda w: fam.a(w) ** -2 * tmat(w / lam + phi), u1, u2, width
            )
            m_closed = fam.c * f_int * np.eye(2) - fam.s * j_int
            w2["M"] = max(
                w2["M"],
                float(np.max(np.abs(m_num - m_closed)))
                / float(np.max(np.abs(m_closed))),
            )
            g_c = fam.quad(
                lambda w: fam.a(w) ** -2 * np.exp(1j * w / lam), u1, u2, width
            )
            g_c = np.exp(1j * phi) * g_c
            e_mat = -fam.t * j_int
            ev = np.sort(np.linalg.eigvalsh(e_mat))
            e_val = fam.t * abs(g_c)
            w2["eig"] = max(
                w2["eig"],
                abs(ev[1] - e_val) / max(e_val, 1e-300),
                abs(ev[0] + e_val) / max(e_val, 1e-300),
            )
            for _ in range(20):
                dy = rng.normal(size=2)
                dx = math.sqrt(fam.c) * dy
                lhs = dx @ np.linalg.solve(m_num, dx)
                rhs = dy @ np.linalg.solve(f_int * np.eye(2) + e_mat, dy)
                w2["ychart"] = max(w2["ychart"], abs(lhs / rhs - 1))
for _ in range(200):
    psi, ph = rng.uniform(-10, 10, size=2)
    r = rot(ph / 2)
    w2["rotid"] = max(
        w2["rotid"], float(np.max(np.abs(r @ tmat(psi) @ r.T - tmat(psi + ph))))
    )
check("W2 gamma^-1 closed form", w2["inv"] < 1e-9, f"{w2['inv']:.1e}")
check("W2 det gamma = a^4", w2["det"] < 1e-9, f"{w2['det']:.1e}")
check("W2 sqrt|det g| = a^2", w2["vol"] < 1e-9, f"{w2['vol']:.1e}")
check(
    "W2 M = c m I - s J (quadrature of the numerical inverse)",
    w2["M"] < 1e-8,
    f"{w2['M']:.1e}",
)
check("W2 eigenvalues of E are +-t|G|", w2["eig"] < 1e-8, f"{w2['eig']:.1e}")
check("W2 x-chart action = y-chart action", w2["ychart"] < 1e-9, f"{w2['ychart']:.1e}")
check(
    "W2 R(phi/2) T(psi) R(phi/2)^T = T(psi + phi)",
    w2["rotid"] < 1e-12,
    f"{w2['rotid']:.1e}",
)


# ================================================================ W3 geodesic shooting
def geodesic(fam, p0, vel0, u_end, nstep, hs):
    def rhs(state):
        p, w = state[:4], state[4:]
        gam = christoffel(fam, p, hs)
        return np.concatenate([w, -np.einsum("abc,b,c->a", gam, w, w)])

    state = np.concatenate([p0, vel0])
    dtau = (u_end - p0[0]) / nstep
    for _ in range(nstep):
        k1 = rhs(state)
        k2 = rhs(state + 0.5 * dtau * k1)
        k3 = rhs(state + 0.5 * dtau * k2)
        k4 = rhs(state + dtau * k3)
        state = state + dtau / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return state


w3_worst = 0.0
w3_affine = 0.0
for q, lam in ((1.0, 1 / 4), (4.0, 1 / 8), (9.0, 1 / 16), (0.25, 1.0)):
    fam = Family(q, lam, 0.7)
    for _ in range(3):
        u1 = rng.uniform(-0.5, -0.1)
        u2 = rng.uniform(0.1, 0.5)
        p0 = np.array([u1, rng.uniform(-1, 1), *rng.normal(size=2)])
        xd = rng.normal(size=2)
        gm = fam.gamma(u1)
        vel0 = np.array([1.0, 0.5 * xd @ gm @ xd, *xd])
        nstep = int(max(800, 80 * (u2 - u1) / lam))
        st = geodesic(fam, p0, vel0, u2, nstep, [lam / 40, 1e-2, 1e-2, 1e-2])
        width = min(lam / 8, 0.02)
        m_num = fam.quad(lambda w: np.linalg.inv(fam.gamma(w)), u1, u2, width)
        dx = st[2:4] - p0[2:4]
        dv = st[1] - p0[1]
        act = dx @ np.linalg.solve(m_num, dx)
        w3_worst = max(w3_worst, abs(2 * dv - act) / max(abs(act), 1e-12))
        w3_affine = max(w3_affine, abs(st[4] - 1.0), abs(st[0] - u2))
check(
    "W3 null geodesics end on 2 Delta v = Delta x^T M^-1 Delta x",
    w3_worst < 1e-6,
    f"max rel diff {w3_worst:.1e}",
)
check("W3 u stays an affine parameter", w3_affine < 1e-8, f"max dev {w3_affine:.1e}")

# ================================================================ W4 the inequality
w4_min = np.inf
w4_opt = 0.0
for q, lam, phi in (
    (1.0, 1 / 4, 0.3),
    (4.0, 1 / 16, 2.0),
    (9.0, 1 / 8, 5.0),
    (0.25, 1.0, 1.0),
):
    fam = Family(q, lam, phi)
    u1, u2 = -0.45, 0.4
    npan = max(200, math.ceil((u2 - u1) / (lam / 10)))
    ed = np.linspace(u1, u2, npan + 1)
    us = ((ed[:-1] + ed[1:]) / 2)[:, None] + ((ed[1:] - ed[:-1]) / 2)[:, None] * GX[
        None, :
    ]
    ws = ((ed[1:] - ed[:-1]) / 2)[:, None] * GW[None, :]
    us, ws = us.ravel(), ws.ravel()
    gs = np.array([fam.gamma(u) for u in us])
    gis = np.linalg.inv(gs)
    m_cum = np.cumsum(ws[:, None, None] * gis, axis=0)
    m_tot = m_cum[-1]
    for _ in range(40):
        dx = rng.normal(size=2) * 2
        p = np.linalg.solve(m_tot, dx)
        act_min = dx @ p
        xp_opt = gis @ p  # x' = gamma^-1 p
        a_opt = float(np.sum(ws * np.einsum("ni,nij,nj->n", xp_opt, gs, xp_opt)))
        w4_opt = max(w4_opt, abs(a_opt / act_min - 1))
        for _ in range(10):
            ks = np.arange(1, 9)
            cs = rng.normal(size=(8, 2)) * rng.uniform(0.01, 2.0)
            arg = np.pi * (us - u1) / (u2 - u1)
            pert_d = (np.cos(np.outer(arg, ks)) * (np.pi * ks / (u2 - u1))) @ cs
            xp = xp_opt + pert_d
            a_path = float(np.sum(ws * np.einsum("ni,nij,nj->n", xp, gs, xp)))
            w4_min = min(w4_min, a_path / act_min)
check(
    "W4 random fixed-endpoint paths never beat Delta x^T M^-1 Delta x",
    w4_min >= 1 - 1e-9,
    f"min ratio {w4_min:.6f}",
)
check("W4 the optimal path attains it", w4_opt < 1e-8, f"max dev {w4_opt:.1e}")
w4_id = 0.0
w4_sign = True
for _ in range(2000):
    a_ = rng.normal(size=(2, 2))
    gm = a_ @ a_.T + 0.1 * np.eye(2)
    ev, evec = np.linalg.eigh(gm)
    gh = evec @ np.diag(np.sqrt(ev)) @ evec.T
    ghi = np.linalg.inv(gh)
    xd, p = rng.normal(size=2), rng.normal(size=2)
    ud = rng.uniform(1e-3, 3)
    lhs = xd @ gm @ xd / (2 * ud) - p @ xd + 0.5 * ud * p @ np.linalg.inv(gm) @ p
    rhs = np.sum((gh @ xd - ud * ghi @ p) ** 2) / (2 * ud)
    w4_id = max(w4_id, abs(lhs - rhs) / max(1.0, abs(lhs)))
    # Lemma 1: future causal vectors have w^u, w^v >= 0
    wx = rng.normal(size=2)
    wu = rng.uniform(-2, 2)
    wv = rng.uniform(-2, 2)
    causal = 2 * wu * wv >= wx @ gm @ wx
    future = wu + wv > 0
    if causal and future and (wu < 0 or wv < 0):
        w4_sign = False
check("W4 identity behind F' >= 0", w4_id < 1e-12, f"max rel err {w4_id:.1e}")
check("W4 Lemma 1: future causal vectors have w^u, w^v >= 0", w4_sign)

# ================================================================ W5 background
w5 = dict(pull=0.0, m=0.0, jac=0.0)
w5_mismatch = 0
w5_npairs = 0
for q in (0.25, 1.0, 4.0, 9.0):
    fam = Family(q, 1.0, 0.0, background=True)
    rq = math.sqrt(q)
    th = math.tan(rq / 2)

    def m_of(u, rq=rq, th=th):
        return (math.tan(rq * u) + th) / rq

    for u in (-0.45, -0.1, 0.2, 0.49):
        mq = fam.quad(lambda w: fam.a(w) ** -2, -0.5, u, 0.01)
        w5["m"] = max(w5["m"], abs(mq / m_of(u) - 1) if m_of(u) else 0.0)
        # Jacobian of (u, v, x) -> (T, Z, x)
        jac = np.zeros((4, 4))
        mp = fam.a(u) ** -2
        jac[0, :2] = [mp / math.sqrt(2), 1 / math.sqrt(2)]
        jac[1, :2] = [mp / math.sqrt(2), -1 / math.sqrt(2)]
        jac[2, 2] = jac[3, 3] = 1.0
        eta = np.diag([-1.0, 1.0, 1.0, 1.0])
        pull = fam.a(u) ** 2 * jac.T @ eta @ jac
        g0 = fam.metric(np.array([u, 0, 0, 0]))
        w5["pull"] = max(w5["pull"], float(np.max(np.abs(pull - g0))))
        w5["jac"] = max(
            w5["jac"], abs(abs(np.linalg.det(jac[:2, :2])) * fam.a(u) ** 2 - 1)
        )
    for _ in range(20000):
        u1, u2 = rng.uniform(-0.5, 0.5, size=2)
        if rng.uniform() < 0.1:
            u2 = u1
        x1, x2 = rng.normal(size=2), rng.normal(size=2)
        if rng.uniform() < 0.1:
            x2 = x1.copy()
        v1 = rng.uniform(-1, 1)
        dm = m_of(u2) - m_of(u1)
        dx = x2 - x1
        if dm > 0 and rng.uniform() < 0.3:  # near-null pair
            v2 = v1 + 0.5 * (dx @ dx) / dm * (1 + rng.normal() * 1e-6)
        else:
            v2 = v1 + rng.uniform(-1, 3)
        dv = v2 - v1
        if u1 < u2:
            thm = 2 * dv >= dx @ dx / dm
        elif u1 == u2:
            thm = bool(np.all(dx == 0)) and dv >= 0
        else:
            thm = False
        d_t = (dm + dv) / math.sqrt(2)
        d_z = (dm - dv) / math.sqrt(2)
        mink = d_t >= math.sqrt(d_z * d_z + dx @ dx)
        # skip only case-(a) pairs within 1e-12 of the null boundary, where
        # rounding could decide; case-(b) pairs (dm = 0) are always compared
        margin = 2 * dm * dv - dx @ dx
        if dm > 0 and abs(margin) < 1e-12 * (1 + abs(2 * dm * dv) + dx @ dx):
            continue
        w5_npairs += 1
        w5_mismatch += int(thm != mink)
check("W5 g_0 = Phi^*(a^2 eta)", w5["pull"] < 1e-9, f"{w5['pull']:.1e}")
check("W5 m(u) closed form against quadrature", w5["m"] < 1e-9, f"{w5['m']:.1e}")
check(
    "W5 |d(T,Z)/d(u,v)| = a^-2, so density a^4 in (T,Z,x)",
    w5["jac"] < 1e-12,
    f"{w5['jac']:.1e}",
)
check(
    "W5 Theorem 2 (gamma = a^2 I) and Delta T >= |Delta z| agree on random pairs",
    w5_mismatch == 0,
    f"{w5_mismatch} mismatches in {w5_npairs} pairs",
)

# ================================================================ W6 chart
w6_mismatch = 0
w6_n = 0
for q, lam, phi in ((1.0, 1 / 4, 0.3), (4.0, 1 / 16, 2.0), (9.0, 1 / 8, 5.0)):
    fam = Family(q, lam, phi)
    for _ in range(100):
        u1, u2 = np.sort(rng.uniform(-0.5, 0.5, size=2))
        width = min(lam / 8, 0.02)
        m_num = fam.quad(lambda w: np.linalg.inv(fam.gamma(w)), u1, u2, width)
        f_int = fam.quad(lambda w: fam.a(w) ** -2, u1, u2, width)
        j_int = fam.quad(lambda w: fam.a(w) ** -2 * tmat(w / lam + phi), u1, u2, width)
        mt = f_int * np.eye(2) - fam.t * j_int
        dy = rng.normal(size=2)
        dx = math.sqrt(fam.c) * dy
        qx = dx @ np.linalg.solve(m_num, dx)
        dv = 0.5 * qx * (1 + rng.normal() * 1e-3)
        qy = dy @ np.linalg.solve(mt, dy)
        if abs(2 * dv - qx) < 1e-10 * qx:
            continue
        w6_n += 1
        w6_mismatch += int((2 * dv >= qx) != (2 * dv >= qy))
check(
    "W6 x-chart and y-chart relations agree on near-null random pairs",
    w6_mismatch == 0,
    f"{w6_mismatch} mismatches in {w6_n} pairs",
)

print(f"SUMMARY pass={npass} failed={failed}")
raise SystemExit(1 if failed else 0)
