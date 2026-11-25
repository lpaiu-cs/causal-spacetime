# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of wave_causal_relation.md.
"""Referee check R1: Ricci tensor by exact derivatives (no finite differences).

Independent of checks.py: the metric's first and second u-derivatives are
written in closed form, Christoffels and their u-derivatives are assembled
from (g, g', g'') with explicit loops, the full Riemann tensor is built and
contracted R_bd = R^a_bad. The sign convention is calibrated on the unit
2-sphere (Ricci = +g).
"""

import math

import numpy as np

rng = np.random.default_rng(7)


def ricci_from_jets(g0, g1, g2):
    """Metric depending on coordinate 0 only; g1 = dg/dx0, g2 = d2g/dx0^2."""
    n = g0.shape[0]
    G0 = np.linalg.inv(g0)
    G1 = -G0 @ g1 @ G0
    dg = np.zeros((n, n, n))
    dg[0] = g1
    ddg = np.zeros((n, n, n))  # ddg[k] = d0 d_k g ; only k = 0
    ddg[0] = g2
    Gam = np.zeros((n, n, n))
    dGam = np.zeros((n, n, n, n))  # dGam[k, a, b, c] = d_k Gam^a_bc ; only k = 0
    for a in range(n):
        for b in range(n):
            for c in range(n):
                s0 = 0.0
                s1 = 0.0
                for d in range(n):
                    br = dg[b][d, c] + dg[c][d, b] - dg[d][b, c]
                    dbr = ddg[b][d, c] + ddg[c][d, b] - ddg[d][b, c]
                    s0 += 0.5 * G0[a, d] * br
                    s1 += 0.5 * G1[a, d] * br + 0.5 * G0[a, d] * dbr
                Gam[a, b, c] = s0
                dGam[0, a, b, c] = s1
    # Riemann R^a_bcd = d_c Gam^a_db - d_d Gam^a_cb + Gam^a_ce Gam^e_db - Gam^a_de Gam^e_cb
    Ric = np.zeros((n, n))
    for b in range(n):
        for d in range(n):
            tot = 0.0
            for a in range(n):
                c = a  # contraction R_bd = R^a_bad
                val = dGam[c, a, d, b] - dGam[d, a, c, b]
                for e in range(n):
                    val += Gam[a, c, e] * Gam[e, d, b] - Gam[a, d, e] * Gam[e, c, b]
                tot += val
            Ric[b, d] = tot
    return Ric, Gam


# ---- sign calibration: unit 2-sphere, coordinates (theta, phi)
th = 0.83
g0 = np.diag([1.0, math.sin(th) ** 2])
g1 = np.diag([0.0, 2 * math.sin(th) * math.cos(th)])
g2 = np.diag([0.0, 2 * math.cos(2 * th)])
Rs, _ = ricci_from_jets(g0, g1, g2)
print("sphere Ricci - g (should be 0):", np.max(np.abs(Rs - g0)))


def rosen_jets(gam, gamp, gampp):
    g0 = np.zeros((4, 4))
    g0[0, 1] = g0[1, 0] = -1.0
    g0[2:, 2:] = gam
    g1 = np.zeros((4, 4))
    g1[2:, 2:] = gamp
    g2 = np.zeros((4, 4))
    g2[2:, 2:] = gampp
    return g0, g1, g2


def tm(psi):
    return np.array([[math.cos(psi), math.sin(psi)], [math.sin(psi), -math.cos(psi)]])


def tmp(psi):  # d/dpsi T
    return np.array([[-math.sin(psi), math.cos(psi)], [math.cos(psi), math.sin(psi)]])


def circular(u, q, lam, phi, beta):
    """gamma = a^2 R(theta) diag(e^{2b}, e^{-2b}) R^T, built from R and D directly."""
    k = math.sqrt(q)
    th_ = 0.5 * (u / lam + phi)
    R = np.array([[math.cos(th_), -math.sin(th_)], [math.sin(th_), math.cos(th_)]])
    D = np.diag([math.exp(2 * beta), math.exp(-2 * beta)])
    Gm = R @ D @ R.T
    # closed-form derivatives via Gm = c I + s T(psi) (checked against Gm below)
    c, s = math.cosh(2 * beta), math.sinh(2 * beta)
    psi = u / lam + phi
    Gm_alt = c * np.eye(2) + s * tm(psi)
    assert np.max(np.abs(Gm - Gm_alt)) < 1e-12
    Gp = (s / lam) * tmp(psi)
    Gpp = -(s / lam**2) * tm(psi)
    A = math.cos(k * u) ** 2
    Ap = -k * math.sin(2 * k * u)
    App = -2 * q * math.cos(2 * k * u)
    gam = A * Gm
    gamp = Ap * Gm + A * Gp
    gampp = App * Gm + 2 * Ap * Gp + A * Gpp
    return gam, gamp, gampp


# ---- circular family: R_uu = 2q - s^2/(2 lam^2), all else 0, at random parameters
worst_formula = 0.0
worst_rest = 0.0
worst_chris = 0.0
for _ in range(400):
    q = rng.uniform(0.01, 9.8)
    lam = 10 ** rng.uniform(-2.5, 0.5)
    phi = rng.uniform(-7, 7)
    beta = rng.uniform(0.0, 2.0)
    umax = math.pi / (2 * math.sqrt(q))
    u = rng.uniform(-0.98 * umax, 0.98 * umax)
    gam, gamp, gampp = circular(u, q, lam, phi, beta)
    Ric, Gam = ricci_from_jets(*rosen_jets(gam, gamp, gampp))
    s = math.sinh(2 * beta)
    want = 2 * q - s**2 / (2 * lam**2)
    scale = 2 * q + s**2 / (2 * lam**2) + 1.0
    worst_formula = max(worst_formula, abs(Ric[0, 0] - want) / scale)
    rest = Ric.copy()
    rest[0, 0] = 0
    worst_rest = max(worst_rest, np.max(np.abs(rest)) / scale)
    # Lemma 4 Christoffels: Gam^v_ij = gam'_ij/2, Gam^i_uj = (gam^-1 gam')^i_j / 2
    L = np.zeros((4, 4, 4))
    L[1, 2:, 2:] = 0.5 * gamp
    P = np.linalg.solve(gam, gamp)
    L[2:, 0, 2:] = 0.5 * P
    L[2:, 2:, 0] = 0.5 * P
    worst_chris = max(worst_chris, np.max(np.abs(L - Gam)) / (1 + np.max(np.abs(Gam))))
print("circular R_uu vs 2q - s^2/(2lam^2), rel:", worst_formula)
print("circular other components, rel:", worst_rest)
print("Lemma 4 Christoffels vs generic, rel:", worst_chris)

# ---- vacuum at s = 2 lam sqrt q
worst_vac = 0.0
for _ in range(200):
    q = rng.uniform(0.01, 9.8)
    lam = 10 ** rng.uniform(-2.5, 0.5)
    phi = rng.uniform(-7, 7)
    beta = 0.5 * math.asinh(2 * lam * math.sqrt(q))
    umax = math.pi / (2 * math.sqrt(q))
    u = rng.uniform(-0.98 * umax, 0.98 * umax)
    Ric, _ = ricci_from_jets(*rosen_jets(*circular(u, q, lam, phi, beta)))
    worst_vac = max(worst_vac, np.max(np.abs(Ric)) / (4 * q + 1))
print("vacuum family max |Ric|/(4q+1):", worst_vac)

# ---- general Rosen gamma(u) = B B^T + eps I, B polynomial: Lemma 4's general formula
worst_gen = 0.0
worst_design = 0.0
worst_gen_rest = 0.0
for _ in range(300):
    B0, B1, B2 = rng.normal(size=(3, 2, 2))
    eps = 0.2
    u = rng.uniform(-1, 1)
    B = B0 + B1 * u + B2 * u**2
    Bp = B1 + 2 * B2 * u
    Bpp = 2 * B2
    gam = B @ B.T + eps * np.eye(2)
    gamp = Bp @ B.T + B @ Bp.T
    gampp = Bpp @ B.T + 2 * Bp @ Bp.T + B @ Bpp.T
    Ric, _ = ricci_from_jets(*rosen_jets(gam, gamp, gampp))
    gi = np.linalg.inv(gam)
    P = gi @ gamp
    Pp = -gi @ gamp @ gi @ gamp + gi @ gampp
    note = -0.5 * np.trace(Pp) - 0.25 * np.trace(P @ P)
    design = -0.5 * np.trace(gi @ gampp) + 0.25 * np.trace(P @ P)
    sc = 1 + abs(note) + np.trace(P @ P)
    worst_gen = max(worst_gen, abs(Ric[0, 0] - note) / sc)
    worst_design = max(worst_design, abs(Ric[0, 0] - design) / sc)
    rest = Ric.copy()
    rest[0, 0] = 0
    worst_gen_rest = max(worst_gen_rest, np.max(np.abs(rest)) / sc)
print("general Rosen R_uu vs -tr P'/2 - tr P^2/4, rel:", worst_gen)
print("general Rosen R_uu vs design.md form, rel:", worst_design)
print("general Rosen other components, rel:", worst_gen_rest)

# ---- tr Q^2 = 2 s^2/lam^2 for Gamma = R D R^T
worst_trq = 0.0
for _ in range(200):
    lam = 10 ** rng.uniform(-2, 0.5)
    beta = rng.uniform(0, 2)
    psi = rng.uniform(-7, 7)
    c, s = math.cosh(2 * beta), math.sinh(2 * beta)
    Gm = c * np.eye(2) + s * tm(psi)
    Gp = (s / lam) * tmp(psi)
    Q = np.linalg.solve(Gm, Gp)
    worst_trq = max(worst_trq, abs(np.trace(Q @ Q) - 2 * s**2 / lam**2) / (1 + 2 * s**2 / lam**2))
print("tr Q^2 vs 2 s^2/lam^2, rel:", worst_trq)
