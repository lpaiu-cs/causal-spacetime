"""Referee check (d): the cube / uncertain-pair construction of Section 3.

Construction audited (referee's explicit version of the text):
  * cubes of side hc in the (T,Z,y1,y2) chart, random grid offset;
  * margins of cube centres  g+ = dT - |dz|,  g- = -dT - |dz|;
  * a pair of cubes is an EDGE iff min(|g+|, |g-|) <= L hc;
  * margins are Lipschitz: moving both points inside their cubes changes a
    margin by at most (1 + sqrt 3) hc;
  * hence the text's claim needs  (L - 1 - sqrt 3) hc >= c0 delta.
The test uses the TIGHTEST admissible cell, hc = c0 delta / (L - 1 - sqrt 3),
with L = 4, so any slack in the claim would show up as a violation.

Tests
  D0  the reference domain is the Minkowski diamond between the tips.
  D1  every non-edge has ONE relation, the same for the wave kernel, the
      background kernel and the cube-centre prediction (zero violations);
      the genuine wave/background disagreements exist and are all edges.
  D1n negative control: with L = 1 (below the Lipschitz constant) the same
      test must report violations.
  D2  orientation: for an edge whose two cube T-intervals do not overlap, the
      reverse relation is false for both kernels; an edge with ambiguous
      orientation has centres within (L+1) hc in every coordinate, and the
      number of such pairs obeys the explicit C N^2 hc^4 bound.
  D3  uniform uncertain-cube bound: for every sampled first event the fraction
      of other events in an uncertain cube is <= C_unc hc (+ binomial noise),
      C_unc = 4 (L + (1+sqrt3)/2) vol_3(proj D) / V_0;  E U <= C(N,2) C_unc hc.
  D4  tree counting without independence: the number of 3-vertex trees
      (paths i-j-k of uncertain pairs) is <= 3 C(N,3) (C_unc hc)^2, and the
      number of 4-vertex trees is <= 16 C(N,4) (C_unc hc)^3.
  D5  scaling: U / (N^2 hc) roughly constant when hc is halved at fixed L
      (uncertain pairs are a codimension-one event of width ~hc).
  D6  interior cubes: an event's cube touches the boundary of D with
      probability <= C_bd hc (measured; reported as the Nh term's constant).

Deterministic (fixed seeds).  Run: python check_d_bad_marks.py
"""

from __future__ import annotations

from math import comb

import numpy as np
from audit_common import Family, Reporter

rep = Reporter("check (d): cubes, uncertain pairs, bad marks")

L_EDGE = 4.0
LIP = 1.0 + np.sqrt(3.0)


def cube_data(X, hc, off):
    idx = np.floor((X - off) / hc)
    centre = off + (idx + 0.5) * hc
    return idx, centre


def in_diamond(fam, P):
    """Membership of Minkowski-chart points (T,Z,y1,y2) in the background diamond."""
    s2 = np.sqrt(2.0)
    Tt, Zt = (fam.Mtot + fam.h) / s2, (fam.Mtot - fam.h) / s2
    r_lo = np.sqrt(P[..., 1] ** 2 + P[..., 2] ** 2 + P[..., 3] ** 2)
    r_hi = np.sqrt((Zt - P[..., 1]) ** 2 + P[..., 2] ** 2 + P[..., 3] ** 2)
    return (P[..., 0] >= r_lo) & (Tt - P[..., 0] >= r_hi)


def analyse(fam, N, hc, rng, want_wave=True, L=L_EDGE):
    pts = fam.sample_reference(N, rng)
    X = fam.to_minkowski(pts)
    off = rng.uniform(0.0, hc, 4)
    idx, C = cube_data(X, hc, off)
    i, j = np.triu_indices(N, k=1)
    dTc = C[j, 0] - C[i, 0]
    dzc = np.sqrt(((C[j, 1:] - C[i, 1:]) ** 2).sum(axis=1))
    gp_c, gm_c = dTc - dzc, -dTc - dzc
    edge = np.minimum(np.abs(gp_c), np.abs(gm_c)) <= L * hc
    out = dict(N=N, hc=hc, i=i, j=j, edge=edge, C=C, X=X, idx=idx, off=off)
    dT = X[j, 0] - X[i, 0]
    dz = np.sqrt(((X[j, 1:] - X[i, 1:]) ** 2).sum(axis=1))
    rel0 = np.where(dT - dz >= 0, 1, np.where(-dT - dz >= 0, -1, 0))
    out["rel0"] = rel0
    out["pred"] = np.where(gp_c > 0, 1, np.where(gm_c > 0, -1, 0))
    out["dTc"], out["dzc"] = dTc, dzc
    if want_wave:
        P = fam.prim(pts[:, 0])
        mm = fam.m(pts[:, 0])
        fwd = pts[j, 0] > pts[i, 0]
        lo = np.where(fwd, i, j)
        hi = np.where(fwd, j, i)
        b = mm[hi] - mm[lo]
        G = fam.pair_G(P[lo], P[hi])
        dv = pts[hi, 1] - pts[lo, 1]
        qf = fam.wave_form(b, G, pts[hi, 2] - pts[lo, 2], pts[hi, 3] - pts[lo, 3])
        related = 2.0 * dv >= qf
        out["relw"] = np.where(related, np.where(fwd, 1, -1), 0)
    return out


def reference_constants(fam):
    """Explicit constants of the referee's version of the Section 3 bounds."""
    V0 = fam.reference_volume()
    A = (fam.Mtot + fam.h) / (2.0 * np.sqrt(2.0))  # semi-major axis of proj_z D
    B2 = fam.Mtot * fam.h / 2.0  # squared semi-minor axis
    vol3 = 4.0 / 3.0 * np.pi * A * B2
    c_unc = 4.0 * (L_EDGE + 0.5 * LIP) * vol3 / V0
    # ambiguous orientation: |dT_c| <= hc and |dz_c| <= (L+1) hc
    c_amb = (
        (2.0 * (1.0 + 0.5) + 0.0)
        * (4.0 / 3.0 * np.pi * (L_EDGE + 1.0 + 0.5 * np.sqrt(3.0)) ** 3)
        / V0
    )
    return V0, vol3, c_unc, c_amb


rng = np.random.default_rng(271828)

# ------------------------------------------------------------------ D0
# The reference domain is a Minkowski diamond in the (T,Z,y) chart: samples
# drawn in (u,v,y) from the a^2-weighted background diamond land exactly in
# J+(tip-) cap J-(tip+) of the flat chart, and the chart map is onto it.
fam0 = Family(4.0, 1 / 16, phi=0.1)
pts0 = fam0.sample_reference(200000, rng)
X0 = fam0.to_minkowski(pts0)
inside = in_diamond(fam0, X0)
s2 = np.sqrt(2.0)
box_lo = np.array([0.0, -fam0.h / s2, -1.0, -1.0])
box_hi = np.array([(fam0.Mtot + fam0.h) / s2, fam0.Mtot / s2, 1.0, 1.0])
Y0 = box_lo + (box_hi - box_lo) * rng.uniform(size=(400000, 4))
m_ = (Y0[:, 0] + Y0[:, 1]) / s2
v_ = (Y0[:, 0] - Y0[:, 1]) / s2
r2_ = Y0[:, 2] ** 2 + Y0[:, 3] ** 2
uvy = (
    (m_ >= 0)
    & (m_ <= fam0.Mtot)
    & (2 * v_ * m_ >= r2_)
    & (v_ >= 0)
    & (2 * (fam0.h - v_) * (fam0.Mtot - m_) >= r2_)
    & (v_ <= fam0.h)
)
rep.check(
    "D0 reference domain = Minkowski diamond between the two tips in the (T,Z,y) chart",
    bool(inside.all()) and bool(np.all(uvy == in_diamond(fam0, Y0))),
    "200000 samples inside; 400000 box points: memberships agree, "
    f"{int(uvy.sum())} inside",
)

# ------------------------------------------------------------------ D1, D2
for q, lam in (
    (1.0, 1 / 16),
    (1.0, 1 / 32),
    (1.0, 1 / 64),
    (0.25, 1 / 32),
    (4.0, 1 / 64),
):
    fam = Family(q, lam, phi=float(rng.uniform(0, 2 * np.pi)))
    hc = fam.c0 * fam.delta / (L_EDGE - LIP)
    N = 1400
    viol = 0
    ndis = 0
    dis_nonedge = 0
    rev_bad = 0
    amb = 0
    amb_far = 0
    nedge = 0
    npair = 0
    for _rpt in range(3):
        d = analyse(fam, N, hc, rng)
        ne = ~d["edge"]
        viol += int(np.sum(ne & ((d["rel0"] != d["relw"]) | (d["rel0"] != d["pred"]))))
        dis = d["rel0"] != d["relw"]
        ndis += int(dis.sum())
        dis_nonedge += int(np.sum(dis & ne))
        ed = d["edge"]
        nedge += int(ed.sum())
        npair += len(ed)
        fixed = ed & (np.abs(d["dTc"]) >= hc)  # cube T-intervals do not overlap
        sgn = np.sign(d["dTc"]).astype(int)
        # reverse relation must be absent in both kernels
        rev_bad += int(np.sum(fixed & ((d["rel0"] == -sgn) | (d["relw"] == -sgn))))
        ambiguous = ed & ~fixed
        amb += int(ambiguous.sum())
        amb_far += int(np.sum(ambiguous & (d["dzc"] > (L_EDGE + 1.0) * hc)))
    V0, vol3, c_unc, c_amb = reference_constants(fam)
    tag = f"q={q:.2f} lam=1/{round(1 / lam)} hc={hc:.5f}"
    rep.check(
        "D1 non-edges deterministic and shared by both kernels " + tag,
        viol == 0 and ndis > 0 and dis_nonedge == 0,
        f"pairs={npair} edges={nedge} violations={viol} "
        f"wave/background disagreements={ndis} (all on edges)",
    )
    amb_bound = 3 * comb(N, 2) * c_amb * hc**4
    rep.check(
        "D2 orientation fixed off the O(hc)-coincidence set " + tag,
        rev_bad == 0 and amb_far == 0 and amb <= max(3.0 * amb_bound, amb_bound + 10),
        f"reverse-true={rev_bad} ambiguous={amb} "
        f"(bound on mean {amb_bound:.2f}) ambiguous-but-far={amb_far}",
    )

# ------------------------------------------------------------------ D1n
# Negative control: with L below the Lipschitz constant 1 + sqrt 3 the claim
# must FAIL (some non-edge pairs have a relation that varies inside the cubes).
# This shows that D1 has power and that "fixed large L" is a real requirement.
fam = Family(1.0, 1 / 32, phi=2.2)
hc = fam.c0 * fam.delta / (L_EDGE - LIP)
d = analyse(fam, 1400, hc, rng, L=1.0)
ne = ~d["edge"]
bad = int(np.sum(ne & ((d["rel0"] != d["relw"]) | (d["rel0"] != d["pred"]))))
rep.check(
    "D1n negative control: L = 1 < 1 + sqrt 3 breaks non-edge determinism (test has "
    "power)",
    bad > 0,
    f"violations with L=1: {bad}",
)

# ------------------------------------------------------------------ D3, D4, D5
fam = Family(1.0, 1 / 64, phi=1.234)
V0, vol3, c_unc, c_amb = reference_constants(fam)
N = 1600
prev_ratio = None
for hc in (0.004, 0.002, 0.001):
    U = 0
    W3 = 0
    W4 = 0
    maxfrac = 0.0
    reps = 3
    for _rpt in range(reps):
        d = analyse(fam, N, hc, rng, want_wave=False)
        ed = d["edge"]
        adj = np.zeros((N, N), dtype=np.float64)
        adj[d["i"][ed], d["j"][ed]] = 1.0
        adj = adj + adj.T
        deg = np.rint(adj.sum(axis=1)).astype(np.int64)
        U += int(ed.sum())
        maxfrac = max(maxfrac, float(deg.max()) / (N - 1))
        W3 += int((deg * (deg - 1) // 2).sum())
        stars = int((deg * (deg - 1) * (deg - 2) // 6).sum())
        tri = int(round(float(np.einsum("ij,ji->", adj @ adj, adj)))) // 6
        ii, jj = d["i"][ed], d["j"][ed]
        paths4 = int(((deg[ii] - 1) * (deg[jj] - 1)).sum()) - 3 * tri
        W4 += stars + paths4
    p = c_unc * hc
    noise = 5.0 * np.sqrt(p * (1 - min(p, 1.0)) / (N - 1) + 1e-12)
    ratio = U / reps / (N * N * hc)
    tag = f"q=1 hc={hc:.4f} N={N}"
    rep.check(
        "D3 uncertain-cube probability <= C_unc hc uniformly in the first event " + tag,
        maxfrac <= p + noise and U / reps <= comb(N, 2) * p,
        f"max_x fraction={maxfrac:.4f} mean "
        f"fraction={2.0 * U / reps / (N * (N - 1)):.4f} "
        f"bound C_unc*hc={p:.4f} (C_unc={c_unc:.1f})",
    )
    rep.check(
        "D4 tree counts <= (number of labelled trees) C(N,k) (C_unc hc)^(k-1) " + tag,
        W3 / reps <= 3 * comb(N, 3) * p**2 and W4 / reps <= 16 * comb(N, 4) * p**3,
        f"3-trees {W3 / reps:.3e} <= {3 * comb(N, 3) * p**2:.3e} ; 4-trees "
        f"{W4 / reps:.3e} <= {16 * comb(N, 4) * p**3:.3e}",
    )
    if prev_ratio is not None:
        rep.check(
            "D5 U/(N^2 hc) stable under halving hc " + tag,
            0.8 < ratio / prev_ratio < 1.25,
            f"ratio {ratio:.4f} vs {prev_ratio:.4f}",
        )
    else:
        rep.info(f"D5 baseline U/(N^2 hc) = {ratio:.4f} at hc={hc:.4f}")
    prev_ratio = ratio

# ------------------------------------------------------------------ D6
fam = Family(1.0, 1 / 64, phi=0.5)
corners = np.array([[(k >> b) & 1 for b in range(4)] for k in range(16)], dtype=float)
prev = None
for hc in (0.02, 0.01, 0.005):
    pts = fam.sample_reference(200000, rng)
    X = fam.to_minkowski(pts)
    off = rng.uniform(0.0, hc, 4)
    idx = np.floor((X - off) / hc)
    low = off + idx * hc
    allin = np.ones(len(X), dtype=bool)
    for cvec in corners:
        allin &= in_diamond(fam, low + cvec * hc)
    frac = 1.0 - allin.mean()
    ok = frac / hc < 60.0 and (prev is None or 0.7 < (frac / hc) / prev < 1.3)
    rep.check(
        f"D6 boundary-cube probability is O(hc) hc={hc:.3f}",
        ok,
        f"P(cube not inside D)={frac:.4f}, /hc={frac / hc:.2f}",
    )
    prev = frac / hc

rep.finish()
