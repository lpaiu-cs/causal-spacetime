"""Referee check (a): the exact family, its causal relation, the background chart.

Verifies with SymPy (and a high-precision numeric fallback):
  A1  Ricci tensor of g = -2 du dv + gamma_AB(u) dx^A dx^B vanishes identically
      for gamma = a^2 R(theta) diag(e^{2b}, e^{-2b}) R(theta)^T,
      a = cos(sqrt(q) u), theta = (u/lam + phi)/2, sinh(2b) = 2 lam sqrt(q).
  A2  gamma^-1 = a^-2 [cosh(2b) I - sinh(2b) T(u/lam + phi)], hence
      M = c m I - s J with the J of the probe.
  A3  For general gamma(u), x' = gamma^-1 p, v' = p^T gamma^-1 p / 2 (u affine)
      is a null geodesic, so the null cone is 2 dv = dx^T M^-1 dx, M = int gamma^-1.
  A4  Numerical minimisation property behind the causal relation: for random
      paths x(u) with fixed endpoints, int x'^T gamma x' >= dx^T M^-1 dx.
  A5  Background: g0 = a^2 (-dT^2 + dZ^2 + dy^2) with T=(m+v)/sqrt2,
      Z=(m-v)/sqrt2, m' = a^-2; volume density a^4 in (T,Z,y); R_uu[g0] = 2q;
      Delta T^2 - Delta Z^2 = 2 Delta m Delta v.

Deterministic.  Run: python check_a_exact_family.py
"""

from __future__ import annotations

import mpmath as mp
import numpy as np
import sympy as sp
from audit_common import Reporter

rep = Reporter("check (a): exact family, causal relation, background chart")

u, v, x1, x2 = sp.symbols("u v x1 x2", real=True)
q, lam = sp.symbols("q lam", positive=True)
phi, beta = sp.symbols("phi beta", real=True)
coords = [u, v, x1, x2]


def christoffel(g):
    ginv = g.inv()
    n = 4
    gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                val = 0
                for d in range(n):
                    val += ginv[a, d] * (
                        sp.diff(g[d, b], coords[c])
                        + sp.diff(g[d, c], coords[b])
                        - sp.diff(g[b, c], coords[d])
                    )
                val = val / 2
                gam[a][b][c] = val
                gam[a][c][b] = val
    return gam


def ricci(g):
    n = 4
    gam = christoffel(g)
    ric = sp.zeros(n, n)
    for b in range(n):
        for c in range(b, n):
            expr = 0
            for a in range(n):
                expr += sp.diff(gam[a][b][c], coords[a]) - sp.diff(
                    gam[a][b][a], coords[c]
                )
                for d in range(n):
                    expr += gam[a][a][d] * gam[d][b][c] - gam[a][c][d] * gam[d][b][a]
            ric[b, c] = expr
            ric[c, b] = expr
    return ric, gam


def rosen_metric(gamma):
    g = sp.zeros(4, 4)
    g[0, 1] = g[1, 0] = -1
    g[2, 2] = gamma[0, 0]
    g[2, 3] = g[3, 2] = gamma[0, 1]
    g[3, 3] = gamma[1, 1]
    return g


# ---------------------------------------------------------------- A1 vacuum
a = sp.cos(sp.sqrt(q) * u)
theta = (u / lam + phi) / 2
R = sp.Matrix([[sp.cos(theta), -sp.sin(theta)], [sp.sin(theta), sp.cos(theta)]])
D = sp.diag(sp.exp(2 * beta), sp.exp(-2 * beta))
gamma = (a**2) * R * D * R.T
gamma = gamma.applyfunc(sp.simplify)
g_wave = rosen_metric(gamma)
ric, _ = ricci(g_wave)

# expected from the hand derivation: R_uu = 2q - sinh(2 beta)^2/(2 lam^2), rest zero
expected_uu = 2 * q - sp.sinh(2 * beta) ** 2 / (2 * lam**2)
rng = np.random.default_rng(20261124)
mp.mp.dps = 40
max_other = mp.mpf(0)
max_uu_dev = mp.mpf(0)
max_uu_vac = mp.mpf(0)
for _ in range(8):
    qv = mp.mpf(float(rng.uniform(0.05, 9.5)))
    lv = mp.mpf(float(rng.uniform(0.01, 0.7)))
    pv = mp.mpf(float(rng.uniform(-3, 3)))
    uv = mp.mpf(float(rng.uniform(-0.5, 0.5)))
    bfree = mp.mpf(float(rng.uniform(-0.8, 0.8)))
    bvac = mp.asinh(2 * lv * mp.sqrt(qv)) / 2
    sub_free = {q: qv, lam: lv, phi: pv, u: uv, beta: bfree}
    sub_vac = {q: qv, lam: lv, phi: pv, u: uv, beta: bvac}
    for i in range(4):
        for j in range(i, 4):
            val_free = mp.mpf(sp.N(ric[i, j].subs(sub_free), 40))
            val_vac = mp.mpf(sp.N(ric[i, j].subs(sub_vac), 40))
            if (i, j) == (0, 0):
                dev = abs(val_free - mp.mpf(sp.N(expected_uu.subs(sub_free), 40)))
                max_uu_dev = max(max_uu_dev, dev)
                max_uu_vac = max(max_uu_vac, abs(val_vac))
            else:
                max_other = max(max_other, abs(val_free), abs(val_vac))
rep.check(
    "A1a R_uu = 2q - sinh(2b)^2/(2 lam^2) for free constant beta (8 random points, 40 "
    "digits)",
    max_uu_dev < mp.mpf(10) ** -30,
    f"max dev {mp.nstr(max_uu_dev, 3)}",
)
rep.check(
    "A1b all other Ricci components vanish (free beta and vacuum beta)",
    max_other < mp.mpf(10) ** -30,
    f"max |R_ab| {mp.nstr(max_other, 3)}",
)
rep.check(
    "A1c R_uu = 0 at sinh(2 beta) = 2 lam sqrt(q): exactly vacuum",
    max_uu_vac < mp.mpf(10) ** -30,
    f"max |R_uu| {mp.nstr(max_uu_vac, 3)}",
)
# Symbolic route (fast): general Rosen identity, then the 2x2 trace formula.
h11, h12, h22 = [sp.Function(n)(u) for n in ("h11", "h12", "h22")]
gam_h = sp.Matrix([[h11, h12], [h12, h22]])
ric_h, _ = ricci(rosen_metric(gam_h))
gi_h = gam_h.inv()
trace_formula = (
    -sp.Rational(1, 2) * (gi_h * gam_h.diff(u, 2)).trace()
    + sp.Rational(1, 4) * ((gi_h * gam_h.diff(u)) * (gi_h * gam_h.diff(u))).trace()
)
gen_resid = sp.cancel(sp.together(ric_h[0, 0] - trace_formula))
others_zero = all(
    sp.cancel(sp.together(ric_h[i, j])) == 0
    for i in range(4)
    for j in range(i, 4)
    if (i, j) != (0, 0)
)
rep.check(
    "A1d symbolic, general gamma(u): R_uu = -tr(g^-1 g'')/2 + tr[(g^-1 g')^2]/4, "
    "other components 0",
    gen_resid == 0 and others_zero,
    f"residual {sp.sstr(gen_resid)[:60]}",
)
Ch, Sh = sp.symbols("Ch Sh", real=True)  # cosh(2 beta), sinh(2 beta), constant
al = u / lam + phi
Tm = sp.Matrix([[sp.cos(al), sp.sin(al)], [sp.sin(al), -sp.cos(al)]])
gam_cs = a**2 * (Ch * sp.eye(2) + Sh * Tm)
gi_cs = (Ch * sp.eye(2) - Sh * Tm) / a**2  # inverse when Ch^2 - Sh^2 = 1
inv_resid = (gam_cs * gi_cs).applyfunc(
    lambda z: sp.simplify(sp.expand(z).subs(Ch**2, 1 + Sh**2))
)
ruu_cs = (
    -sp.Rational(1, 2) * (gi_cs * gam_cs.diff(u, 2)).trace()
    + sp.Rational(1, 4) * ((gi_cs * gam_cs.diff(u)) * (gi_cs * gam_cs.diff(u))).trace()
)
ruu_cs = sp.expand(sp.trigsimp(sp.expand(ruu_cs)))
for _ in range(3):
    ruu_cs = sp.expand(
        ruu_cs.subs(Ch**4, (1 + Sh**2) ** 2)
        .subs(Ch**3, Ch * (1 + Sh**2))
        .subs(Ch**2, 1 + Sh**2)
    )
ruu_cs = sp.simplify(sp.trigsimp(ruu_cs))
resid_cs = sp.simplify(ruu_cs - (2 * q - Sh**2 / (2 * lam**2)))
rep.check(
    "A1e symbolic, circular family: R_uu = 2q - sinh(2b)^2/(2 lam^2)",
    inv_resid == sp.eye(2) and resid_cs == 0,
    f"residual {sp.sstr(resid_cs)[:60]}",
)
rep.check(
    "A1e2 R(theta) diag(e^{2b},e^{-2b}) R^T = cosh(2b) I + sinh(2b) T(2 theta)",
    (R * D * R.T - (sp.cosh(2 * beta) * sp.eye(2) + sp.sinh(2 * beta) * Tm)).applyfunc(
        lambda z: sp.simplify(z.rewrite(sp.exp))
    )
    == sp.zeros(2, 2),
)
vac_expr = sp.simplify((2 * q - Sh**2 / (2 * lam**2)).subs(Sh, 2 * lam * sp.sqrt(q)))
rep.check(
    "A1e3 that R_uu vanishes at sinh(2 beta) = 2 lam sqrt(q)",
    vac_expr == 0,
    sp.sstr(vac_expr),
)
# the scale factor: a'' = -q a exactly, same a as the background, for every lam
rep.check(
    "A1f a = cos(sqrt(q) u) solves a'' = -q a",
    sp.simplify(sp.diff(a, u, 2) + q * a) == 0,
)

# ---------------------------------------------------------------- A2 gamma^-1
T2 = sp.Matrix(
    [
        [sp.cos(u / lam + phi), sp.sin(u / lam + phi)],
        [sp.sin(u / lam + phi), -sp.cos(u / lam + phi)],
    ]
)
claim_inv = (sp.cosh(2 * beta) * sp.eye(2) - sp.sinh(2 * beta) * T2) / a**2
resid = (gamma * claim_inv - sp.eye(2)).applyfunc(
    lambda z: sp.simplify(z.rewrite(sp.exp))
)
rep.check(
    "A2 gamma^-1 = a^-2 [cosh(2b) I - sinh(2b) T(u/lam+phi)]  (so M = c m I - s J)",
    resid == sp.zeros(2, 2),
    "residual entries %s" % [sp.sstr(z) for z in resid],
)
det_resid = sp.simplify(sp.trigsimp(sp.expand((R * D * R.T).det())) - 1)
rep.check(
    "A2b det(R D R^T) = 1, so det(gamma) = a^4 and the volume density is a^2",
    det_resid == 0,
    f"residual {sp.sstr(det_resid)[:60]}",
)

# ---------------------------------------------------------------- A3 null geodesics
g11, g12, g22 = [sp.Function(n)(u) for n in ("g11", "g12", "g22")]
gam_gen = sp.Matrix([[g11, g12], [g12, g22]])
g_gen = rosen_metric(gam_gen)
chris = christoffel(g_gen)
p1, p2 = sp.symbols("p1 p2", real=True)
pvec = sp.Matrix([p1, p2])
xdot = gam_gen.inv() * pvec
vdot = (pvec.T * gam_gen.inv() * pvec)[0, 0] / 2
tangent = [sp.Integer(1), vdot, xdot[0], xdot[1]]
geo_ok = True
for mu in range(4):
    expr = sp.diff(tangent[mu], u)
    for al in range(4):
        for be in range(4):
            expr += chris[mu][al][be] * tangent[al] * tangent[be]
    if sp.simplify(expr) != 0:
        geo_ok = False
null_norm = sum(
    g_gen[i, j] * tangent[i] * tangent[j] for i in range(4) for j in range(4)
)
rep.check(
    "A3a x' = gamma^-1 p, v' = p^T gamma^-1 p/2 solves the geodesic equation (general "
    "gamma)",
    geo_ok,
)
rep.check("A3b that curve is null", sp.simplify(null_norm) == 0)

# ---------------------------------------------------------------- A4 minimisation
# Any causal curve from X to Y has u monotone and v' >= x'^T gamma x'/2, so
# Delta v >= (1/2) min int x'^T gamma x' = (1/2) dx^T M^-1 dx.  Spot-check the
# minimum on random perturbed paths for the actual family.
rngn = np.random.default_rng(7)
worst = np.inf
for _trial in range(40):
    qn = rngn.uniform(0.1, 9.0)
    ln = rngn.uniform(0.02, 0.5)
    ph = rngn.uniform(-3, 3)
    s_, t_ = np.sort(rngn.uniform(-0.5, 0.5, 2))
    nn = 4001
    uu = np.linspace(s_, t_, nn)
    an = np.cos(np.sqrt(qn) * uu)
    bn = 0.5 * np.arcsinh(2 * ln * np.sqrt(qn))
    th = 0.5 * (uu / ln + ph)
    c2, s2 = np.cos(2 * th), np.sin(2 * th)
    ch, sh = np.cosh(2 * bn), np.sinh(2 * bn)
    # gamma = a^2 [ch I + sh T(2 theta)],  gamma^-1 = a^-2 [ch I - sh T(2 theta)]
    gm = np.empty((nn, 2, 2))
    gm[:, 0, 0] = an**2 * (ch + sh * c2)
    gm[:, 1, 1] = an**2 * (ch - sh * c2)
    gm[:, 0, 1] = gm[:, 1, 0] = an**2 * sh * s2
    gi = np.linalg.inv(gm)
    w = np.full(nn, (t_ - s_) / (nn - 1))
    w[0] *= 0.5
    w[-1] *= 0.5
    M = (gi * w[:, None, None]).sum(axis=0)
    dx = rngn.normal(size=2)
    opt = dx @ np.linalg.solve(M, dx)
    # optimal path derivative gamma^-1 p, plus a perturbation with zero net shift
    p = np.linalg.solve(M, dx)
    xd = np.einsum("nij,j->ni", gi, p)
    pert = (
        rngn.normal(size=(3, 2))[:, None, :]
        * np.stack(
            [np.sin(2 * np.pi * (k + 1) * (uu - s_) / (t_ - s_)) for k in range(3)]
        )[:, :, None]
    )
    xd2 = xd + pert.sum(axis=0) * rngn.uniform(0.1, 2.0)
    xd2 = xd2 - ((xd2 - xd) * w[:, None]).sum(axis=0) / (t_ - s_)  # keep endpoints
    act_opt = (np.einsum("ni,nij,nj->n", xd, gm, xd) * w).sum()
    act_pert = (np.einsum("ni,nij,nj->n", xd2, gm, xd2) * w).sum()
    if abs(act_opt - opt) > 1e-6 * abs(opt):
        worst = -np.inf
    worst = min(worst, (act_pert - opt) / abs(opt))
rep.check(
    "A4 int x'^T gamma x' du >= dx^T M^-1 dx on 40 random perturbed paths (equality "
    "on gamma^-1 p)",
    worst > -1e-9,
    f"min relative excess {worst:.3e}",
)

# ---------------------------------------------------------------- A5 background
g0 = rosen_metric(a**2 * sp.eye(2))
ric0, _ = ricci(g0)
ric0 = ric0.applyfunc(sp.simplify)
target = sp.zeros(4, 4)
target[0, 0] = 2 * q
rep.check(
    "A5a Ricci[g0] = 2q du du only (null dust)", ric0 == target, sp.sstr(ric0[0, 0])
)
m_of_u = (sp.tan(sp.sqrt(q) * u) + sp.tan(sp.sqrt(q) / 2)) / sp.sqrt(q)
rep.check(
    "A5b m'(u) = a^-2 and m(-1/2) = 0",
    sp.simplify(sp.diff(m_of_u, u) - 1 / a**2) == 0
    and sp.simplify(m_of_u.subs(u, -sp.Rational(1, 2))) == 0,
)
Tc = (m_of_u + v) / sp.sqrt(2)
Zc = (m_of_u - v) / sp.sqrt(2)
new = [Tc, Zc, x1, x2]
Jac = sp.Matrix(4, 4, lambda i, j: sp.diff(new[i], coords[j]))
eta = sp.diag(-1, 1, 1, 1)
pull = (a**2 * Jac.T * eta * Jac).applyfunc(sp.simplify)
rep.check(
    "A5c g0 = a^2 (-dT^2 + dZ^2 + dy^2): background is conformally flat in (T,Z,y)",
    pull == g0,
)
dens = sp.simplify(sp.sqrt(sp.simplify(-g0.det())) / sp.Abs(Jac.det()))
dens_resid = sp.simplify(dens**2 - a**8)
rep.check(
    "A5d volume density in (T,Z,y) is a(u(m))^4",
    dens_resid == 0,
    f"density^2 - a^8 = {sp.sstr(dens_resid)}",
)
dm, dv = sp.symbols("dm dv", real=True)
dT = (dm + dv) / sp.sqrt(2)
dZ = (dm - dv) / sp.sqrt(2)
rep.check(
    "A5e Delta T^2 - Delta Z^2 = 2 Delta m Delta v",
    sp.simplify(dT**2 - dZ**2 - 2 * dm * dv) == 0,
)
# relation (2): Delta T >= |Delta z|  <=>  Delta m >= 0, Delta v >= 0,
# and 2 Delta m Delta v >= |Delta y|^2
rs = np.random.default_rng(11)
dmv = rs.normal(size=200000)
dvv = rs.normal(size=200000)
dyy = rs.normal(size=(200000, 2)) * rs.uniform(0, 1.5, size=(200000, 1))
lhs = (dmv + dvv) / np.sqrt(2) >= np.sqrt((dmv - dvv) ** 2 / 2 + (dyy**2).sum(axis=1))
rhs = (dmv >= 0) & (dvv >= 0) & (2 * dmv * dvv >= (dyy**2).sum(axis=1))
rep.check(
    "A5f eq.(2) equals {dm>=0, dv>=0, 2 dm dv >= |dy|^2} on 2e5 random displacements",
    bool(np.all(lhs == rhs)),
    f"mismatches {int(np.sum(lhs != rhs))}",
)

rep.finish()
