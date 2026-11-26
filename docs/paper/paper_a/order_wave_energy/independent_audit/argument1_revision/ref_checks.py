# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of quantized_order_bound_r2.md.
"""Referee checks for quantized_order_bound_r2.md (argument 1, revision 2).

Independent of the author's checks.py: nothing is imported from the repository.
Deterministic, seeded draws only to test inequalities. No causal order is
sampled. CPU only. Run: python ref_checks.py [path-to-worktree]

Sections
  R0  provenance: SHA-256 of the original, section identity of 6 and 7
  R1  Lemma 2.1 (operator bound) on random endpoints, incl. tiny l and lam=1
  R2  Lemma 2.2 (sandwich) with the exact E matrix and phase, all u-orderings
  R3  Section 3 non-edge claim end-to-end with the exact wave kernel at h=h_min
  R4  Section 3 geometric constants (1+sqrt3)h and (L+(1+sqrt3)/2)h, tightness
  R5  Lemma 4.1 exact (piecewise-polynomial) G = -grad p for k=3, non-uniform
  R6  Lemmas 4.2, 4.3 on the same exact p
  R7  Section 5 algebra: (8) for every N>=1, 0<lam<=1; coverage; margins
"""

from __future__ import annotations

import hashlib
import itertools
import math
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

WT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    r"E:\lab\causal-spacetime\.claude\worktrees\upbeat-newton-3e92fe")
REL = "docs/paper/paper_a/order_wave_energy/"
failed = 0
npass = 0


def check(name, ok, detail=""):
    global failed, npass
    print(("PASS " if ok else "FAIL ") + name + (" | " + detail if detail else ""))
    if ok:
        npass += 1
    else:
        failed += 1


def git_show(rev_path):
    return subprocess.run(["git", "-C", str(WT), "show", rev_path],
                          capture_output=True, check=True).stdout


# ---------------------------------------------------------------- R0
orig_bytes = git_show("a032f0e:" + REL + "quantized_order_bound.md")
r2_bytes = git_show("a032f0e:" + REL + "quantized_order_bound_r2.md")
sha = hashlib.sha256(orig_bytes).hexdigest()
check("R0 SHA-256 of original at a032f0e matches r2's record",
      sha == "39217a16c9770ff9b56b326984262b1424f4aaa4e1a2bd934238d658bed0564b", sha)
frozen = git_show("5dda03d:" + REL + "quantized_order_bound.md").decode("utf-8")


def split_sections(text):
    out, key, buf = {}, "pre", []
    for ln in text.splitlines():
        if ln.startswith("## ") and ln[3:4].isdigit():
            out[key] = buf
            key, buf = ln.split(".")[0][3:], []
        buf.append(ln)
    out[key] = buf
    return out


so = split_sections(orig_bytes.decode("utf-8"))
sf = split_sections(frozen)
sr = split_sections(r2_bytes.decode("utf-8"))
check("R0 original Sections 1-5 identical to the 5dda03d freeze",
      all(so[k] == sf[k] for k in "12345"))
s6 = [ln for ln in sr["6"]]
note_end = next(i for i, ln in enumerate(s6) if ln.endswith("checks.py).*"))
s6_wo_note = s6[:1] + s6[note_end + 1:]
check("R0 r2 Section 6 = original Section 6 apart from the italic note",
      s6_wo_note == so["6"])
s7 = sr["7"]
while s7 and s7[-1] == "":
    s7 = s7[:-1]
o7 = so["7"]
while o7 and o7[-1] == "":
    o7 = o7[:-1]
check("R0 r2 Section 7 = original Section 7 (trailing blank lines ignored)", s7 == o7)
txt = r2_bytes.decode("utf-8")
check("R0 r2 contains no sentence deriving (8)->0 from (1) for general sequences "
      "(Section 8 says it is 'written above')",
      "This proves (1)" not in txt and "Under (1)" not in txt,
      "neither 'This proves (1)' nor 'Under (1)' occurs; reported as finding")

# ---------------------------------------------------------------- R1
gx, gw = np.polynomial.legendre.leggauss(16)


def osc_int(u1, u2, lam, q, phase=0.0):
    """int_u1^u2 sec^2(sqrt q u) exp(i(u/lam+phase)) du, panels <= lam/6."""
    n = max(2, math.ceil((u2 - u1) / (lam / 6)))
    edges = np.linspace(u1, u2, n + 1)
    mid, half = (edges[1:] + edges[:-1]) / 2, (edges[1:] - edges[:-1]) / 2
    u = (mid[:, None] + half[:, None] * gx[None, :]).ravel()
    w = (half[:, None] * gw[None, :]).ravel()
    return np.sum(w * np.exp(1j * (u / lam + phase)) / np.cos(math.sqrt(q) * u) ** 2)


def fmax(q):
    return 1 / math.cos(math.sqrt(q) / 2) ** 2


def CE(q):
    return 2 * math.sqrt(q) * (4 * fmax(q) - 2)


def bexact(u1, u2, q):
    r = math.sqrt(q)
    return (math.tan(r * u2) - math.tan(r * u1)) / r


rng = np.random.default_rng(20261126)
worst_r, worst_b, ntest = 0.0, 0.0, 0
for q in (0.25, 1.0, 4.0, 9.0, 9.8):
    for lam in (1.0, 0.5, 1 / 7, 1 / 33, 1 / 200):
        tpar = 2 * lam * math.sqrt(q) / math.sqrt(1 + 4 * q * lam * lam)
        for _ in range(40):
            u1 = rng.uniform(-0.5, 0.5)
            ell = 10 ** rng.uniform(-7, 0)
            u2 = min(0.5, u1 + ell)
            if u2 <= u1:
                continue
            e = tpar * abs(osc_int(u1, u2, lam, q))
            worst_r = max(worst_r, e / (CE(q) * lam * min(u2 - u1, lam)))
            worst_b = max(worst_b, e / bexact(u1, u2, q))
            ntest += 1
        # the end-to-end segment, where TV f is largest
        e = tpar * abs(osc_int(-0.5, 0.5, lam, q))
        worst_r = max(worst_r, e / (CE(q) * lam * min(1.0, lam)))
check("R1 Lemma 2.1 e <= C_E lam min(l,lam) (q up to 9.8, lam up to 1, l down to 1e-7)",
      worst_r <= 1.0, f"{ntest} segments, max ratio {worst_r:.4f}")
check("R1 Lemma 2.1 e < b", worst_b < 1.0, f"max e/b {worst_b:.4f}")

# ---------------------------------------------------------------- R2


def wave_E(u1, u2, lam, q, phase):
    tpar = 2 * lam * math.sqrt(q) / math.sqrt(1 + 4 * q * lam * lam)
    z = osc_int(u1, u2, lam, q, phase)
    return -tpar * np.array([[z.real, z.imag], [z.imag, -z.real]])


def wave_rel(u1, u2, dv, dy, E, b):
    if u1 < u2:
        return 2 * dv >= dy @ np.linalg.solve(b * np.eye(2) + E, dy)
    if u1 == u2:
        return bool(np.all(dy == 0)) and dv >= 0
    return False


def gmarg(b, dv, dy):
    dT = (b + dv) / math.sqrt(2)
    dZ = (b - dv) / math.sqrt(2)
    return dT - math.sqrt(dZ * dZ + dy @ dy)


viol = {"first": 0, "second": 0, "first_tight": 0, "second_tight": 0}
worst_ratio_pos, worst_ratio_neg = 0.0, 0.0
ntest = 0
for q in (0.25, 1.0, 4.0, 9.0):
    for lam in (0.7, 1 / 9, 1 / 40):
        delta = lam * lam
        c0 = 2 * math.sqrt(2) * CE(q)
        for _ in range(60):
            u1 = rng.uniform(-0.5, 0.45)
            ell = 10 ** rng.uniform(-6, 0)
            u2 = min(0.5, u1 + ell)
            phase = rng.uniform(0, 2 * math.pi)
            b = bexact(u1, u2, q)
            E = wave_E(u1, u2, lam, q, phase)
            e = abs(np.linalg.eigvalsh(E)).max()
            ctight = 2 * math.sqrt(2) * e / delta
            for _ in range(30):
                r = math.sqrt(b) * 10 ** rng.uniform(-3, 0.5)
                ang = rng.uniform(0, 2 * math.pi)
                dy = r * np.array([math.cos(ang), math.sin(ang)])
                vstar = 0.5 * dy @ np.linalg.solve(b * np.eye(2) + E, dy)
                for s in (-1.0, -1e-3, -1e-9, 0.0, 1e-9, 1e-3, 1.0):
                    dv = vstar * (1 + s * e / b) if s != 0 else vstar
                    rel = wave_rel(u1, u2, dv, dy, E, b)
                    g = gmarg(b, dv, dy)
                    ntest += 1
                    tol = 1e-9 * (abs(b) + abs(dv) + r * r) / delta
                    if g >= c0 * delta and not rel:
                        viol["first"] += 1
                    if rel and g < -c0 * delta:
                        viol["second"] += 1
                    if g >= (ctight + tol) * delta and not rel:
                        viol["first_tight"] += 1
                    if rel and g < -(ctight + tol) * delta:
                        viol["second_tight"] += 1
                    if abs(s) <= 1e-9:
                        if g > 0:
                            worst_ratio_pos = max(worst_ratio_pos, g / (c0 * delta))
                        else:
                            worst_ratio_neg = max(worst_ratio_neg, -g / (c0 * delta))
check("R2 Lemma 2.2, u1<u2, both inclusions with c0 (exact E incl. phase)",
      viol["first"] == 0 and viol["second"] == 0, f"{ntest} pairs, {viol}")
check("R2 Lemma 2.2 with per-pair c = 2 sqrt2 e/delta (rounding tol 1e-9 relative)",
      viol["first_tight"] == 0 and viol["second_tight"] == 0,
      f"max |g|/(c0 delta) on wave null cone: +{worst_ratio_pos:.3f} / -{worst_ratio_neg:.3f}")
# u1 = u2 and u1 > u2
bad = 0
for _ in range(20000):
    b = -abs(rng.normal()) if rng.random() < 0.7 else 0.0
    dv = rng.normal()
    dy = rng.normal(size=2) * (0 if (b == 0 and rng.random() < 0.5) else 1)
    g = gmarg(b, dv, dy)
    rel = (b == 0 and bool(np.all(dy == 0)) and dv >= 0)
    if g >= 1e-12:  # first set must be empty
        bad += 1
    if rel and g < -1e-12:
        bad += 1
check("R2 Lemma 2.2, u1=u2 and u1>u2: first set empty, relation => g=0", bad == 0,
      "20000 draws")

# ---------------------------------------------------------------- R3
# End-to-end: points in the reference diamond, (T,z) -> (u,v,y), exact wave kernel,
# h = h_min with L = 4, cubes on the grid. Every NON-edge pair must have the same
# relation under both kernels; report how many genuine disagreements sit on edges.
q, lam, phase, Lc = 1.0, 1 / 12, 0.7, 4.0
delta = lam * lam
c0 = 2 * math.sqrt(2) * CE(q)
hmin = c0 * delta / (Lc - 1 - math.sqrt(3))
h = hmin
rq = math.sqrt(q)
M = (math.tan(rq / 2) * 2) / rq
# primitive F(u) = int_{-1/2}^u f e^{iu/lam} on a fine grid, then local Gauss
ugrid = np.linspace(-0.5, 0.5, 4001)
Fgrid = np.zeros(len(ugrid), complex)
for k in range(1, len(ugrid)):
    Fgrid[k] = Fgrid[k - 1] + osc_int(ugrid[k - 1], ugrid[k], lam, q)
g8x, g8w = np.polynomial.legendre.leggauss(10)


def Fprim(u):
    k = min(int((u + 0.5) / (1 / 4000)), 3999)
    a = ugrid[k]
    if u <= a:
        return Fgrid[k]
    mid, half = (a + u) / 2, (u - a) / 2
    x = mid + half * g8x
    return Fgrid[k] + np.sum(half * g8w * np.exp(1j * x / lam) / np.cos(rq * x) ** 2)


def to_uvy(P):
    T, Z, y1, y2 = P
    m = (T + Z) / math.sqrt(2)
    v = (T - Z) / math.sqrt(2)
    u = math.atan(rq * m - math.tan(rq / 2)) / rq
    return u, v, np.array([y1, y2])


tpar = 2 * lam * rq / math.sqrt(1 + 4 * q * lam * lam)
top = np.array([(M + 1) / math.sqrt(2), (M - 1) / math.sqrt(2), 0.0, 0.0])


def in_D(P):
    return P[0] >= np.linalg.norm(P[1:]) and top[0] - P[0] >= np.linalg.norm(top[1:] - P[1:])


pts = []
while len(pts) < 1800:
    P = np.array([rng.uniform(0, top[0]), rng.uniform(-0.9, 0.9 + top[1]),
                  rng.uniform(-1.1, 1.1), rng.uniform(-1.1, 1.1)])
    if in_D(P):
        pts.append(P)
pts = np.array(pts)
uvy = [to_uvy(P) for P in pts]
U = np.array([x[0] for x in uvy])
V = np.array([x[1] for x in uvy])
Yc = np.array([x[2] for x in uvy])
Fu = np.array([Fprim(u) for u in U])
cube = np.floor(pts / h)
centre = (cube + 0.5) * h
I, J = np.triu_indices(len(pts), 1)
d = pts[J] - pts[I]
dT, dz = d[:, 0], np.linalg.norm(d[:, 1:], axis=1)
dc = centre[J] - centre[I]
dTc, dzc = dc[:, 0], np.linalg.norm(dc[:, 1:], axis=1)
# orient each pair so that u_a < u_b (ties have probability zero)
swap = U[I] > U[J]
A_ = np.where(swap, J, I)
B_ = np.where(swap, I, J)
bvec = (np.tan(rq * U[B_]) - np.tan(rq * U[A_])) / rq
zc = np.exp(1j * phase) * (Fu[B_] - Fu[A_])
e11, e12 = -tpar * zc.real, -tpar * zc.imag      # E = [[e11, e12], [e12, -e11]]
dy = Yc[B_] - Yc[A_]
ee2 = e11 ** 2 + e12 ** 2
quad_E = e11 * (dy[:, 0] ** 2 - dy[:, 1] ** 2) + 2 * e12 * dy[:, 0] * dy[:, 1]
Qv = (bvec * np.sum(dy * dy, axis=1) - quad_E) / (bvec ** 2 - ee2)   # (bI+E)^-1 = (bI-E)/(b^2-e^2)
wav = 2 * (V[B_] - V[A_]) >= Qv
dd = pts[B_] - pts[A_]
bgd = dd[:, 0] >= np.linalg.norm(dd[:, 1:], axis=1)
dis = wav != bgd
gpair = np.abs(np.abs(dT) - dz)


def classify(Lval):
    edge = (np.abs(dTc - dzc) <= Lval * h) | (np.abs(-dTc - dzc) <= Lval * h)
    return int(np.sum(dis & edge)), int(np.sum(dis & ~edge))


on_e, off_e = classify(Lc)
check("R3 non-edge claim end-to-end at h = h_min (q=1, lam=1/12, L=4)",
      off_e == 0,
      f"{len(I)} pairs, {int(np.sum(gpair < c0 * delta))} with |margin|<c0 delta, "
      f"disagreements on edges {on_e}, on non-edges {off_e}, h_min={hmin:.4f}")
check("R3 sandwich end-to-end: every disagreement has |margin| <= c0 delta",
      bool(np.all(gpair[dis] <= c0 * delta)),
      f"max |margin|/(c0 delta) among {int(dis.sum())} disagreements "
      f"{(gpair[dis].max() / (c0 * delta)) if dis.any() else 0:.3f}")
on_c, off_c = classify(1.0)
check("R3 control: classifying with L = 1 at the same h puts disagreements on non-edges",
      off_c > 0, f"non-edge disagreements with L=1: {off_c} (power of the test above)")

# ---------------------------------------------------------------- R4
worst_shift, worst_unif, worst_edge_imp = 0.0, 0.0, 0.0
for _ in range(200000):
    hh = 1.0
    nc = rng.integers(-12, 13, size=4)
    cX = np.zeros(4) + 0.5
    cY = nc + 0.5
    dTc = cY[0] - cX[0]
    dzc = np.linalg.norm(cY[1:] - cX[1:])
    # corners and interior
    oX = rng.choice([-0.5, 0.5], size=4) if rng.random() < 0.5 else rng.uniform(-0.5, 0.5, 4)
    oY = rng.choice([-0.5, 0.5], size=4) if rng.random() < 0.5 else rng.uniform(-0.5, 0.5, 4)
    X, Y = cX + oX, cY + oY
    dT, dz = Y[0] - X[0], np.linalg.norm(Y[1:] - X[1:])
    for sgn in (1, -1):
        worst_shift = max(worst_shift, abs((sgn * dT - dz) - (sgn * dTc - dzc)) / (1 + math.sqrt(3)))
    Lt = 4.0
    edge = abs(dTc - dzc) <= Lt or abs(-dTc - dzc) <= Lt
    if edge:
        worst_edge_imp = max(worst_edge_imp, abs(abs(dTc) - dzc) / Lt)
        val = abs(abs(Y[0] - cX[0]) - np.linalg.norm(Y[1:] - cX[1:]))
        worst_unif = max(worst_unif, val / (Lt + (1 + math.sqrt(3)) / 2))
check("R4 margin moves by at most (1+sqrt3)h inside two cubes", worst_shift <= 1 + 1e-12,
      f"max ratio {worst_shift:.4f} (1 = tight)")
check("R4 edge => ||dT_c|-|dz_c|| <= L h", worst_edge_imp <= 1 + 1e-12, f"{worst_edge_imp:.4f}")
check("R4 edge => ||T_Y-T_cX|-|z_Y-z_cX|| <= (L+(1+sqrt3)/2)h", worst_unif <= 1 + 1e-12,
      f"max ratio {worst_unif:.4f}")

# ---------------------------------------------------------------- R5, R6
# k = 3, rows (head, tail): a.s = s_head - s_tail >= tau. Densities rho_i(s) = al+be s.
G5x, G5w = np.polynomial.legendre.leggauss(6)


def G_row(e, rows, tau, dens):
    """Exact slice integral G_e(tau) = -d p / d tau_e, piecewise-polynomial quadrature."""
    j, i = rows[e]
    c = ({0, 1, 2} - {i, j}).pop()
    te = tau[e]
    lows, ups = [(0.0, 0.0)], [(0.0, 1.0)]  # (slope in s_i, intercept)
    for f, (x, y) in enumerate(rows):
        if f == e:
            continue
        tf = tau[f]
        if (x, y) == (c, i):
            lows.append((1.0, tf))
        elif (x, y) == (i, c):
            ups.append((1.0, -tf))
        elif (x, y) == (c, j):
            lows.append((1.0, te + tf))
        elif (x, y) == (j, c):
            ups.append((1.0, te - tf))
        else:
            raise ValueError
    a0, a1 = max(0.0, -te), min(1.0, 1.0 - te)
    if a1 <= a0:
        return 0.0
    fns = lows + ups
    bps = {a0, a1}
    for (s1, i1), (s2, i2) in itertools.combinations(fns, 2):
        if s1 != s2:
            x = (i2 - i1) / (s1 - s2)
            if a0 < x < a1:
                bps.add(x)
    bps = sorted(bps)
    (ai, bi), (aj, bj), (ac, bc) = dens[i], dens[j], dens[c]
    tot = 0.0
    for lo_, hi_ in zip(bps[:-1], bps[1:]):
        half, mid = (hi_ - lo_) / 2, (hi_ + lo_) / 2
        for xg, wg in zip(G5x, G5w):
            s = mid + half * xg
            lo = max(sl * s + ic for sl, ic in lows)
            hi = min(sl * s + ic for sl, ic in ups)
            if hi <= lo:
                continue
            inner = ac * (hi - lo) + bc * (hi * hi - lo * lo) / 2
            tot += half * wg * (ai + bi * s) * (aj + bj * (s + te)) * inner
    return tot


def p_exact(rows, tau, dens):
    """p(tau) = int_{max(tau_0,-1)}^{1} G_0(sigma) d sigma, exact on sigma-pieces."""
    j, i = rows[0]
    c = ({0, 1, 2} - {i, j}).pop()
    # candidate s_i breakpoints are A + B sigma with B in {0,-1}; collect A's
    A0, A1 = [0.0, 1.0], [0.0, 1.0]  # B=0 / B=-1 constants
    for f, (x, y) in enumerate(rows):
        if f == 0:
            continue
        tf = tau[f]
        for k in (0.0, 1.0):
            if c in (x, y) and i in (x, y):
                A0 += [k - tf, k + tf]
            if c in (x, y) and j in (x, y):
                A1 += [k - tf, k + tf]
    sig_bp = {max(tau[0], -1.0), 1.0}
    for a in A0:
        for a2 in A1:
            s = a2 - a  # a = a2 - sigma
            if max(tau[0], -1.0) < s < 1.0:
                sig_bp.add(s)
    sig_bp = sorted(sig_bp)
    tot = 0.0
    for lo_, hi_ in zip(sig_bp[:-1], sig_bp[1:]):
        half, mid = (hi_ - lo_) / 2, (hi_ + lo_) / 2
        for xg, wg in zip(G5x, G5w):
            t2 = list(tau)
            t2[0] = mid + half * xg
            tot += half * wg * G_row(0, rows, t2, dens)
    return tot


def grad_p(rows, tau, dens):
    return -np.array([G_row(e, rows, tau, dens) for e in range(len(rows))])


def make_dens(r):
    out = []
    for _ in range(3):
        be = r.uniform(-1.8, 1.8)
        al = 1.0 - be / 2 + r.uniform(0.0, 0.3)  # min over [0,1] stays >= ~0.1
        al = max(al, -be + 0.1, 0.1)
        out.append((al, be))
    return out


def Lam(rows, dens):
    R = [max(al, al + be) for al, be in dens]
    Rp = [abs(be) for al, be in dens]
    m = len(rows)
    return (m + 1) * np.prod(R) + max(Rp[j] * np.prod([R[l] for l in range(3) if l != j])
                                      for j in range(3))


# validate the exact G/p code against a brute-force 3D midpoint grid once
rows_t = [(1, 0), (2, 1), (0, 2)]
dens_t = [(1.0, 0.0), (0.6, 0.8), (1.4, -0.8)]
tau_t = [-0.15, -0.1, -0.2]
def p_brute(nb):
    s = (np.arange(nb) + 0.5) / nb
    tot = 0.0
    S1, S2 = np.meshgrid(s, s, indexing="ij")
    for s0 in s:
        Sv = [np.full_like(S1, s0), S1, S2]
        ind = np.ones_like(S1)
        for (hd, tl), tv in zip(rows_t, tau_t):
            ind *= (Sv[hd] - Sv[tl] >= tv)
        rho = np.ones_like(S1)
        for k, (al, be) in enumerate(dens_t):
            rho *= al + be * Sv[k]
        tot += float(np.sum(rho * ind))
    return tot / nb ** 3


p_ex = p_exact(rows_t, tau_t, dens_t)
pb1, pb2 = p_brute(150), p_brute(300)
check("R5 exact p code agrees with brute-force grids (error shrinks ~1/n)",
      abs(p_ex - pb2) < 2e-3 and abs(p_ex - pb2) < abs(p_ex - pb1) + 1e-4,
      f"exact {p_ex:.6f}, grid150 {pb1:.6f}, grid300 {pb2:.6f}")
# also check grad p = -G against a central difference of p_exact
eps = 1e-5
fd = []
for e in range(3):
    tp = list(tau_t); tm = list(tau_t)
    tp[e] += eps; tm[e] -= eps
    fd.append((p_exact(rows_t, tp, dens_t) - p_exact(rows_t, tm, dens_t)) / (2 * eps))
gp = grad_p(rows_t, tau_t, dens_t)
check("R5 grad p = -G (central difference of exact p)", np.max(np.abs(np.array(fd) - gp)) < 1e-6,
      f"max diff {np.max(np.abs(np.array(fd) - gp)):.2e}")

# all orientation patterns on the triangle, and two-row paths sharing a vertex
configs = []
for signs in itertools.product((0, 1), repeat=3):
    base = [(1, 0), (2, 1), (2, 0)]
    configs.append([(a, b) if s_ == 0 else (b, a) for (a, b), s_ in zip(base, signs)])
for signs in itertools.product((0, 1), repeat=2):
    base = [(1, 0), (2, 1)]
    configs.append([(a, b) if s_ == 0 else (b, a) for (a, b), s_ in zip(base, signs)])
    base = [(1, 0), (2, 0)]
    configs.append([(a, b) if s_ == 0 else (b, a) for (a, b), s_ in zip(base, signs)])
worst_lip, worst_row, worst_gp, worst_sq, worst_squeeze = 0.0, 0.0, 0.0, 0.0, 0.0
n_lip = 0
for ci, rows in enumerate(configs):
    for _ in range(3):
        dens = make_dens(rng)
        Lm = Lam(rows, dens)
        Rs = [max(al, al + be) for al, be in dens]
        Rps = [abs(be) for al, be in dens]
        m = len(rows)
        for _ in range(8):
            tau = list(rng.uniform(-1.1, 0.9, size=m))
            G0 = grad_p(rows, tau, dens)
            # Lipschitz of the gradient along random directions and steps
            for step in (1e-3, 0.05, 0.4):
                d = rng.normal(size=m)
                d *= step / np.linalg.norm(d)
                t2 = list(np.array(tau) + d)
                G1 = grad_p(rows, t2, dens)
                worst_lip = max(worst_lip, np.linalg.norm(G1 - G0) / step / Lm)
                n_lip += 1
            # row-wise constants c_ef
            for f in range(m):
                eta = 0.03
                t2 = list(tau); t2[f] += eta
                G1 = grad_p(rows, t2, dens)
                for e in range(m):
                    jh = rows[e][0]
                    cef = (np.prod(Rs) if f != e else
                           2 * np.prod(Rs) + Rps[jh] * np.prod([Rs[l] for l in range(3) if l != jh]))
                    worst_row = max(worst_row, abs(G1[e] - G0[e]) / eta / cef)
            # Lemma 4.2: |grad p|^2 <= 2 Lam p, and sqrt p Lipschitz
            p0 = p_exact(rows, tau, dens)
            if p0 > 1e-14:
                worst_gp = max(worst_gp, (G0 @ G0) / (2 * Lm * p0))
            d = rng.normal(size=m)
            d *= rng.choice([1e-3, 0.1, 0.6]) / np.linalg.norm(d)
            p1 = p_exact(rows, list(np.array(tau) + d), dens)
            worst_sq = max(worst_sq, abs(math.sqrt(max(p1, 0)) - math.sqrt(max(p0, 0)))
                           / np.linalg.norm(d) / math.sqrt(Lm / 2))
            # Lemma 4.3 squeeze width
            for epsq in (1e-3, 0.02):
                pp = p_exact(rows, list(np.array(tau) - epsq), dens)
                pm = p_exact(rows, list(np.array(tau) + epsq), dens)
                worst_squeeze = max(worst_squeeze, (math.sqrt(max(pp, 0)) - math.sqrt(max(pm, 0)))
                                    / (2 * math.sqrt(Lm / 2) * epsq * math.sqrt(m)))
check("R5 Lemma 4.1: grad p is Lambda-Lipschitz (k=3, all orientations, non-uniform)",
      worst_lip <= 1 + 1e-9, f"{len(configs)} row sets, {n_lip} pairs, max ratio {worst_lip:.4f}")
check("R5 Lemma 4.1 (iii): row-wise c_ef bounds", worst_row <= 1 + 1e-9, f"max ratio {worst_row:.4f}")
check("R6 Lemma 4.2: |grad p|^2 <= 2 Lambda p", worst_gp <= 1 + 1e-9, f"max ratio {worst_gp:.4f}")
check("R6 Lemma 4.2: sqrt p is sqrt(Lambda/2)-Lipschitz", worst_sq <= 1 + 1e-9,
      f"max ratio {worst_sq:.4f}")
check("R6 Lemma 4.3 squeeze width <= 2 sqrt(Lambda/2) eps sqrt(m)", worst_squeeze <= 1 + 1e-9,
      f"max ratio {worst_squeeze:.4f}")

# targeted near the zero set of the cyclic triangle with non-uniform densities
rows = [(1, 0), (2, 1), (0, 2)]
dens = [(0.2, 1.6), (1.8, -1.6), (1.0, 0.0)]
Lm = Lam(rows, dens)
worst_gp_c = 0.0
for epsv in (1e-4, 1e-3, 1e-2, 0.05):
    for tilt in itertools.product((-1, 0, 1), repeat=3):
        tau = [-epsv * (1 + 0.3 * t) for t in tilt]
        p0 = p_exact(rows, tau, dens)
        G0 = grad_p(rows, tau, dens)
        if p0 > 0:
            worst_gp_c = max(worst_gp_c, (G0 @ G0) / (2 * Lm * p0))
check("R6 Lemma 4.2 near the cyclic zero set, non-uniform densities", worst_gp_c <= 1,
      f"max |grad p|^2/(2 Lambda p) = {worst_gp_c:.4f}, Lambda = {Lm:.3f}")

# ---------------------------------------------------------------- R7


def rhs7_minus_rhs8(lnN, lnd, K, cst, CK):
    """Return (case, ok). Logs throughout."""
    lnh = 2 * (lnd - K * lnN) / (2 * K + 1)
    c1 = cst["c1"]
    t2 = (lnN + 2 * lnd) / (2 * K + 1)          # ln (N delta^2)^(1/(2K+1)) = ln N h
    t3 = ((3 * K + 1) * lnN + 2 * K * lnd) / (2 * K + 1)
    rhs8 = CK * (math.exp(min(lnN + lnd, 50)) + math.exp(min(t2, 50)) + math.exp(min(t3, 50)))
    if lnh > 0:
        return "h>1", rhs8 >= 1
    if lnh < math.log(c1) + lnd:
        return "h<hmin", rhs8 >= 1
    Nh = math.exp(lnN + lnh)
    if Nh > 1:
        return "Nh>1", rhs8 >= 1
    if rhs8 >= 1:
        return "in(3a),trivial", True
    tv7 = (cst["C1"] * math.exp(lnN + lnd) + cst["Cbd"] * Nh
           + cst["Camb"] * math.exp(2 * lnN + 4 * lnh)
           + cst["Ctree"] * math.exp((K + 1) * lnN + K * lnh)
           + cst["CHell"] * math.exp(lnN + lnd - lnh / 2))
    return "in(3a)", tv7 <= rhs8 * (1 + 1e-12)


okall, cases = True, {}
for K in (1, 2, 3, 5, 14):
    for c1 in (0.5, 5.04, 14.25, 104.4, 1.07e4):
        cst = dict(c1=c1, C1=rng.uniform(0.1, 10), Cbd=rng.uniform(0.1, 10),
                   Camb=rng.uniform(0.1, 10), Ctree=rng.uniform(0.1, 10),
                   CHell=rng.uniform(0.1, 10))
        CK = max(1.0, c1 ** ((3 * K + 1) / (2 * K)), cst["C1"], cst["Cbd"] + cst["Camb"],
                 cst["Ctree"] + cst["CHell"])
        for lnd in np.linspace(math.log(1e-14), 0.0, 41):
            for lnN in np.linspace(0.0, math.log(1e40), 121):
                case, ok = rhs7_minus_rhs8(lnN, lnd, K, cst, CK)
                cases[case] = cases.get(case, 0) + 1
                okall &= ok
check("R7 (8) holds on the whole grid N in [1,1e40], delta in [1e-14,1]", okall, str(cases))
# N^2 h^4 <= (N h)^3 needs h <= N: in (3a) h <= 1 <= N.
# exponent coverage, exact rationals, and the margin as a function of K at alpha = 1.3
ok = True
for jn in range(1, 20000):
    a = Fraction(jn, 15000)
    if a >= Fraction(4, 3):
        break
    K = math.floor(a / (4 - 3 * a)) + 1
    ok &= (4 * K - (3 * K + 1) * a) > 0
    if K > 1:
        Km = K - 1
        ok &= (4 * Km - (3 * Km + 1) * a) <= 0  # minimal K
check("R7 K = floor(alpha/(4-3alpha))+1 covers and is minimal, alpha = j/15000 < 4/3", ok)
a = Fraction(13, 10)
marg = {K: float((4 * K - (3 * K + 1) * a) / (2 * K + 1)) for K in (14, 20, 50, 1000)}
check("R7 alpha=1.3 margin: 0.00345 at K=14 but grows with K (to 0.05)",
      marg[14] < 0.0035 and marg[1000] > 0.0035, str({k: round(v, 5) for k, v in marg.items()}))
# C_K's dependence on 2^(M+1): it enters through sqrt(C_H)
for K in (2, 5, 14):
    Mk = K * (K - 1) // 2
    print(f"INFO K={K}: M={Mk}, 2^(M+1)={2.0 ** (Mk + 1):.3g}, sqrt=2^((M+1)/2)="
          f"{2.0 ** ((Mk + 1) / 2):.3g}, (K+1)^(K-1)/(K+1)!={(K + 1) ** (K - 1) / math.factorial(K + 1):.3g}")
for q in (0.25, 1.0, 4.0, 9.0):
    c0v = 2 * math.sqrt(2) * CE(q)
    print(f"INFO q={q}: c0={c0v:.4g}, c1 (L=4)={c0v / (3 - math.sqrt(3)):.4g}")

# ---------------------------------------------------------------- R8
# The imported one-point bound TV(mu_lam, mu_0) <= C_1 lam^2, deterministically.
# Both laws: density a^2 in (u,v,y) on their diamonds between tips (u,v,y) =
# (-1/2,0,0) and (1/2,1,0). For fixed (u, direction n), the y-plane integral of the
# v-length is exact: int_0^inf max(0, 1 - rho (A_+ + A_-)/2) d rho / 2 = 1/(A_+ + A_-),
# and for the intersection replace A_+- by max(A_+-, B_+-). TV = 1 - mu(int)/max(mu).


def one_point_tv(q, lam, phase=0.3, nth=96):
    rq = math.sqrt(q)
    tpar = 2 * lam * rq / math.sqrt(1 + 4 * q * lam * lam)
    npan = max(64, math.ceil(1 / (lam / 6)))
    edges = np.linspace(-0.5, 0.5, npan + 1)
    xg, wg = np.polynomial.legendre.leggauss(10)
    us, ws, Fs = [], [], []
    Fcum = 0j
    for a, b in zip(edges[:-1], edges[1:]):
        mid, half = (a + b) / 2, (b - a) / 2
        x = mid + half * xg
        for xx, ww in zip(x, half * wg):
            # primitive at node: F(a) + int_a^xx, by a 10-point rule on [a, xx]
            m2, h2 = (a + xx) / 2, (xx - a) / 2
            yy = m2 + h2 * xg
            Fn = Fcum + np.sum(h2 * wg * np.exp(1j * yy / lam) / np.cos(rq * yy) ** 2)
            us.append(xx); ws.append(ww); Fs.append(Fn)
        Fcum = Fcum + np.sum(half * wg * np.exp(1j * x / lam) / np.cos(rq * x) ** 2)
    us, ws, Fs = np.array(us), np.array(ws), np.array(Fs)
    Ftot = Fcum
    th = (np.arange(nth) + 0.5) * 2 * math.pi / nth
    c2, s2 = np.cos(2 * th), np.sin(2 * th)
    mminus = (np.tan(rq * us) + math.tan(rq / 2)) / rq
    mplus = (math.tan(rq / 2) - np.tan(rq * us)) / rq
    zm = np.exp(1j * phase) * Fs
    zp = np.exp(1j * phase) * (Ftot - Fs)

    def A(mv, zv):
        a_, b_ = -tpar * zv.real, -tpar * zv.imag
        nEn = a_[:, None] * c2[None, :] + b_[:, None] * s2[None, :]
        e2 = a_ ** 2 + b_ ** 2
        return (mv[:, None] - nEn) / (mv[:, None] ** 2 - e2[:, None])

    Am, Ap = A(mminus, zm), A(mplus, zp)
    Bm, Bp = (1 / mminus)[:, None], (1 / mplus)[:, None]
    a2 = np.cos(rq * us) ** 2
    wgt = (ws * a2)[:, None] * (2 * math.pi / nth)
    mu_l = np.sum(wgt / (Ap + Am))
    mu_0 = np.sum(wgt / (Bp + Bm) * np.ones_like(Am))
    mu_i = np.sum(wgt / (np.maximum(Ap, Bp) + np.maximum(Am, Bm)))
    return 1 - mu_i / max(mu_l, mu_0)


rows8 = []
for q in (0.25, 1.0, 4.0):
    for lam in (1 / 8, 1 / 16, 1 / 32, 1 / 64):
        rows8.append((q, lam, one_point_tv(q, lam) / lam ** 2))
for q in (0.25, 1.0, 4.0):
    vals = [r[2] for r in rows8 if r[0] == q]
    print(f"INFO one-point TV/lam^2 at q={q}, lam=1/8..1/64: " + ", ".join(f"{v:.3f}" for v in vals))
ok8 = all(max(r[2] for r in rows8 if r[0] == q) < 2 * min(r[2] for r in rows8 if r[0] == q)
          for q in (0.25, 1.0, 4.0))
check("R8 imported one-point bound: TV(mu_lam,mu_0)/lam^2 stays bounded (lam 1/8..1/64)", ok8)

print(f"SUMMARY pass={npass} failed={failed}")
raise SystemExit(1 if failed else 0)
