# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of wave_causal_relation.md.
"""Referee check R5: is checks.py's einsum Ricci assembly right on a GENERIC metric?

The Rosen metric has so many zero Christoffels that an index transposition
in the assembly might not show. Here the author's christoffel/ricci (copied
verbatim from checks.py) are applied to a random Lorentzian metric depending
on all four coordinates, g(x) = eta + eps (B.x + x.C.x), and compared with an
exact computation (closed-form derivatives, explicit loops, full Riemann
contraction R_bd = R^a_bad, sign calibrated on S^2 in r1).
"""

import numpy as np

rng = np.random.default_rng(3)
n = 4
eta = np.diag([-1.0, 1.0, 1.0, 1.0])


def sym(a):
    return 0.5 * (a + a.T)


B = np.array([sym(rng.normal(size=(n, n))) for _ in range(n)])
C = np.zeros((n, n, n, n))
for k in range(n):
    for m in range(k, n):
        cc = sym(rng.normal(size=(n, n)))
        C[k, m] = cc
        C[m, k] = cc
eps = 0.15


def metric(x):
    return eta + eps * (np.einsum("k,kab->ab", x, B) + np.einsum("k,l,klab->ab", x, x, C))


def dmetric(x):  # dg[k] = eps (B_k + 2 sum_l C_kl x_l)
    return eps * (B + 2 * np.einsum("l,klab->kab", x, C))


def ddmetric():  # ddg[k, l] = 2 eps C_kl
    return 2 * eps * C


def exact_ricci(x):
    g = metric(x)
    G = np.linalg.inv(g)
    dg = dmetric(x)
    ddg = ddmetric()
    dG = np.array([-G @ dg[k] @ G for k in range(n)])
    Gam = np.zeros((n, n, n))
    dGam = np.zeros((n, n, n, n))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    br = dg[b][d, c] + dg[c][d, b] - dg[d][b, c]
                    Gam[a, b, c] += 0.5 * G[a, d] * br
                    for k in range(n):
                        dbr = ddg[k, b][d, c] + ddg[k, c][d, b] - ddg[k, d][b, c]
                        dGam[k, a, b, c] += 0.5 * dG[k][a, d] * br + 0.5 * G[a, d] * dbr
    Ric = np.zeros((n, n))
    for b in range(n):
        for d in range(n):
            for a in range(n):
                c = a
                val = dGam[c, a, d, b] - dGam[d, a, c, b]
                for e in range(n):
                    val += Gam[a, c, e] * Gam[e, d, b] - Gam[a, d, e] * Gam[e, c, b]
                Ric[b, d] += val
    return Ric


# ---- copied verbatim from checks.py (with fam.metric -> metric)
def d1(fun, p, k, h):
    e = np.zeros(4)
    e[k] = h
    return (-fun(p + 2 * e) + 8 * fun(p + e) - 8 * fun(p - e) + fun(p - 2 * e)) / (12 * h)


def christoffel(p, hs):
    g = metric(p)
    gi = np.linalg.inv(g)
    dg = np.array([d1(metric, p, k, hs[k]) for k in range(4)])
    t1 = np.einsum("bdc->dbc", dg)
    t2 = np.einsum("cdb->dbc", dg)
    t3 = dg
    return 0.5 * np.einsum("ad,dbc->abc", gi, t1 + t2 - t3)


def ricci(p, hs):
    gam = christoffel(p, hs)
    dgam = np.array([d1(lambda pp: christoffel(pp, hs), p, k, hs[k]) for k in range(4)])
    term1 = np.einsum("aabd->bd", dgam)
    tr = np.einsum("aab->b", gam)
    dtr = np.einsum("daab->bd", dgam)
    term3 = np.einsum("e,ebd->bd", tr, gam)
    term4 = np.einsum("ade,eab->bd", gam, gam)
    return term1 - dtr + term3 - term4


worst = 0.0
for _ in range(20):
    x = rng.uniform(-0.3, 0.3, size=4)
    Rex = exact_ricci(x)
    Rfd = ricci(x, [1e-3] * 4)
    worst = max(worst, np.max(np.abs(Rex - Rfd)) / np.max(np.abs(Rex)))
print(f"generic metric: max |Ric_checks.py - Ric_exact| / max|Ric| = {worst:.2e}")
print("typical max|Ric| =", np.max(np.abs(exact_ricci(np.zeros(4)))))
print("Ricci symmetric (exact)?", np.max(np.abs(Rex - Rex.T)))
