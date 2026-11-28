# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of wave_causal_relation.md.
"""Referee check R4: can W4 (first part) and W4's Lemma-1 sign test fail?

W4 first part: x' = gamma^-1 p + delta', p = M^-1 dx with M the SAME quadrature
of gamma^-1, and delta' = sum c_k (k pi/L) cos(k pi (u-u1)/L). Then
  int x'^T gamma x' = act_min + 2 p^T int delta' + int delta'^T gamma delta'.
The middle term vanishes up to quadrature error, the last is >= 0. So the ratio
is >= 1 by construction, whatever gamma and M are. Here we reproduce the
construction (copied) and measure the cross term, and repeat the test with a
WRONG M (e.g. M scaled by 1.3, or the background M) to see whether W4's first
part would still pass. Also: how many draws of the Lemma-1 sign test are in
the tested class (causal and future)?
"""

import math

import numpy as np

GX, GW = np.polynomial.legendre.leggauss(16)
rng = np.random.default_rng(5)


def rot(th):
    return np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])


def gamma(u, q, lam, phi):
    s = 2 * lam * math.sqrt(q)
    beta = 0.5 * math.asinh(s)
    r = rot(0.5 * (u / lam + phi))
    d = np.diag([math.exp(2 * beta), math.exp(-2 * beta)])
    return math.cos(math.sqrt(q) * u) ** 2 * r @ d @ r.T


for q, lam, phi in ((1.0, 1 / 4, 0.3), (9.0, 1 / 8, 5.0)):
    u1, u2 = -0.45, 0.4
    npan = max(200, math.ceil((u2 - u1) / (lam / 10)))
    ed = np.linspace(u1, u2, npan + 1)
    us = ((ed[:-1] + ed[1:]) / 2)[:, None] + ((ed[1:] - ed[:-1]) / 2)[:, None] * GX[None, :]
    ws = ((ed[1:] - ed[:-1]) / 2)[:, None] * GW[None, :]
    us, ws = us.ravel(), ws.ravel()
    gs = np.array([gamma(u, q, lam, phi) for u in us])
    gis = np.linalg.inv(gs)
    m_tot = np.sum(ws[:, None, None] * gis, axis=0)
    worst_cross = 0.0
    min_ratio = {"true M": np.inf, "1.3*M (wrong)": np.inf, "M^-1 used as M (wrong)": np.inf}
    for _ in range(40):
        dx = rng.normal(size=2) * 2
        for label, mm in (
            ("true M", m_tot),
            ("1.3*M (wrong)", 1.3 * m_tot),
            ("M^-1 used as M (wrong)", np.linalg.inv(m_tot)),
        ):
            # W4 as written: the optimal path is built from the SAME gis, so for a
            # wrong M the 'optimal' path does not even reach dx. Rebuild it as W4 does.
            p = np.linalg.solve(mm, dx)
            act_min = dx @ p
            xp_opt = gis @ p
            for _ in range(10):
                ks = np.arange(1, 9)
                cs = rng.normal(size=(8, 2)) * rng.uniform(0.01, 2.0)
                arg = np.pi * (us - u1) / (u2 - u1)
                pert_d = (np.cos(np.outer(arg, ks)) * (np.pi * ks / (u2 - u1))) @ cs
                xp = xp_opt + pert_d
                a_path = float(np.sum(ws * np.einsum("ni,nij,nj->n", xp, gs, xp)))
                min_ratio[label] = min(min_ratio[label], a_path / act_min)
                if label == "true M":
                    cross = 2 * p @ np.sum(ws[:, None] * pert_d, axis=0)
                    worst_cross = max(worst_cross, abs(cross) / act_min)
    print(f"q={q} lam={lam}: |cross term|/act_min max = {worst_cross:.1e}; W4 min ratio:",
          {k: round(v, 4) for k, v in min_ratio.items()})

# Lemma-1 sign test: fraction of draws in the tested class
n_cf = 0
for _ in range(2000):
    a_ = rng.normal(size=(2, 2))
    gm = a_ @ a_.T + 0.1 * np.eye(2)
    wx = rng.normal(size=2)
    wu, wv = rng.uniform(-2, 2), rng.uniform(-2, 2)
    if 2 * wu * wv >= wx @ gm @ wx and wu + wv > 0:
        n_cf += 1
print(f"Lemma-1 sign test: {n_cf}/2000 draws are causal and future (the only ones tested)")
