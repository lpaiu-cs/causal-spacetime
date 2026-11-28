"""Deterministic checks for local_null_challenge_r2.md (argument 2 revised).

No order is sampled. Seeded random draws are used only to place test points.
Wave integrals use a composite 16-point Gauss-Legendre rule with panels at
most lambda/4 wide (the helpers are copied from argument3_lemmas/checks.py).

Pass criteria (written before the full run):
  N1  p0(z) = C1 C3/(3 C2^2), C_r(z) = 0F1(r + 1/2; -z/4)/(2r - 1)!, by
      mpmath. |p0(1) - 0.095919| < 2e-6, |p0(0.25) - 0.099031| < 2e-6, and
      the window widths (1/10 - p0(q))/(1/10) are within 5e-4 of 0.041
      (q = 1) and 0.010 (q = 0.25). The fast series used below agrees with
      mpmath to 1e-14 at 9 points of [0, 1].
  N2  The 12 archived rows (local_null_evidence/rows.jsonl, moments are
      (p, f, c), a = p^2 c). (a) Every row has p_Y >= p_X, f_Y <= f_X and
      a_Y >= a_X. (b) The limit range [a_X/p_Y^2, a_Y/p_X^2] reproduces the
      audit referee's s6 X4 intervals to 1e-5 at its three rows (q = 0.25,
      r = 1024: [0.097688, 0.101401], outside; q = 0.25, r = 2048:
      [0.098604, 0.100462], outside; q = 1, r = 2048: [0.097139, 0.099002],
      inside (p0(q), 1/10)). The verdict for all 12 rows is printed.
  N3  For N = 4..3000, N = 2^j and 2^j +/- 1 (j = 2..40), and N = 10^j
      (j = 1..12): with B_least = (32N^4 - 1).bit_length() and
      B_impl = (32N^4).bit_length(), both satisfy 32N^4 <= 2^B <= 64N^4,
      B <= 6 + 4 log2 N, and (6N + 3N^2) 2 e^-B <= (3/4) N^-2. B_impl =
      B_least + 1 exactly when N is a power of two.
  N4  Lemma 2.1 by exact binomial sums: for kk in {10, 40, 160}, d in
      {0.001, 0.02, 0.1, 0.3, 0.5, 0.8}, B in {1, 3, 6}, the probability
      that |xi_bar - d| > r_D(xi_bar) is at most 4 e^-B in all 54 cases.
  N5  Lemma 2.2. (a) For n = 2..7 and seeded integer kernels, the average of
      V_pi over all n! permutations equals U exactly (Fractions). (b) For a
      three-valued R with probabilities (1/2, 1/3, 1/6), a fixed {0,1}
      symmetric kernel, n in {4, 5} and s in {-8, -2, 2, 8}:
      E e^(sU) <= E e^(sV) (relative 1e-12), by full enumeration.
  N6  Lemma 2.4 on 20000 seeded configurations, with an explicit function
      whose derivative bounds on R are exact upper bounds: (4) holds in
      every case (relative 1e-12), and the largest ratio of the two sides
      is printed.
  N7  Lemma 3.1. (a) q in {0.25, 1}, phase in {0, 0.37, 2.1}, lambda in
      {1/256, 1/1024, 1/4096}, u on 9 points of [-1/4, 1/4], tau in
      {0.3, pi/2, pi, 3pi/2, 6.0}: the residual |Phi(u + tau lam) - Phi(u)
      - lam^2 l| / lam^3 is at most 3000 at every lambda, and its maximum
      at 1/4096 is at most twice its maximum at 1/256. On tau in [pi/2,
      3pi/2], |l| >= c_L. (b) The sign of (5'): q in {0.25, 1}, phase
      0.37, lambda in {1/256, 1/1024}, u in {-1/4, 0, 1/4}, |z| = 0.2 at 6
      angles: |dv + z^T Phi z / 2| / lam^4 <= 500 at both lambda.
  N8  Lemma 3.3. q in {0.25, 1}, phase in {0.37, 2.1}, lambda in {1/512,
      1/2048, 1/8192}; 300 seeded pairs from W_lam plus 16 corner pairs.
      At lambda = 1/2048 and 1/8192: every exact wave slack is positive,
      and both canonical points are in K. At 1/8192: min kappa/lam^2 >=
      (3/64) c_L rho_1^2. At 1/512: information only.
  N9  Lemma 4.2, with background moments at the canonical points (the wave
      correction O(lam^4 log) is not included). q in {0.25, 1}; lambda in
      {2^-11, 2^-14, 2^-17, 2^-20}; 2000 seeded pairs from W_lam (phase
      0.37). (i) max |theta_0(Y) - theta_0(X)|_inf / lam <= 20 at every
      lambda, and its value at 2^-20 is within a factor 1.5 of its value at
      2^-14. (ii) [p0(0.64q), p0(0.04q)] lies strictly inside
      (p0(q), 1/10). (iii) At 2^-20 the corner-box c range of every pair
      lies in the window. The fractions at the other lambda are printed.
  N10 Information only: for the 12 archived rows,
      |population decrease - wave margin| / (lam^4 log(1/lam)).

Before the full run (debugging, recorded here): the helpers were run on
single cases and not compared with any criterion. One archived row (q =
0.25, r = 512) was read by hand before N2 was written; its limit range
leaves the window.

Run: python checks.py   (system Python 3.11 with numpy and mpmath)
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
ROWS = HERE.parent.parent / "local_null_evidence" / "rows.jsonl"
NODES, WEIGHTS = np.polynomial.legendre.leggauss(16)
RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, info: str) -> None:
    RESULTS.append((name, bool(ok), info))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {info}", flush=True)


def info(name: str, text: str) -> None:
    print(f"[INFO] {name}: {text}", flush=True)


# ---------------------------------------------------------------- geometry


class Geo:
    """Background quantities for fixed q."""

    def __init__(self, q: float) -> None:
        self.q = q
        self.k = math.sqrt(q)
        self.fmax = 1.0 / math.cos(self.k / 2) ** 2
        self.CE = 2 * self.k * (4 * self.fmax - 2)
        self.M = float(self.m(-0.5, 0.5))
        self.gamma = float(self.alpha(1.0))
        self.lam1 = 1.0 / (4 * self.CE)

    def a(self, u):
        return np.cos(self.k * u)

    def f(self, u):
        return 1.0 / np.cos(self.k * u) ** 2

    def m(self, s, t):
        return (np.tan(self.k * t) - np.tan(self.k * s)) / self.k

    def alpha(self, ell):
        x = self.k * np.asarray(ell, dtype=float)
        return (1 - x / np.tan(x)) / (2 * self.q)


def panels(s: float, t: float, width: float):
    n = max(1, math.ceil((t - s) / width - 1e-12))
    edges = np.linspace(s, t, n + 1)
    lo, hi = edges[:-1, None], edges[1:, None]
    x = 0.5 * (hi - lo) * NODES[None, :] + 0.5 * (hi + lo)
    w = 0.5 * (hi - lo) * WEIGHTS[None, :]
    return edges, x.ravel(), w.ravel()


class Wave:
    """beta(s, t) through an accurate primitive of f e^(iu/lam)."""

    def __init__(self, geo: Geo, lam: float, phi: float) -> None:
        self.g = geo
        self.lam = lam
        self.phi = phi
        self.c = math.sqrt(1 + 4 * geo.q * lam**2)
        self.pref = (2 * lam * geo.k / self.c) * np.exp(1j * phi)
        self.width = min(lam / 4, 1 / 64)
        edges, x, w = panels(-0.5, 0.5, self.width)
        self.edges = edges
        vals = (w * self.integrand(x)).reshape(-1, 16).sum(axis=1)
        self.Gedge = np.concatenate([[0.0], np.cumsum(vals)])

    def integrand(self, x):
        return self.g.f(x) * np.exp(1j * x / self.lam)

    def G(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        j = np.clip(np.searchsorted(self.edges, x, side="right") - 1, 0, None)
        j = np.minimum(j, len(self.edges) - 2)
        a0 = self.edges[j]
        span = (x - a0)[:, None]
        nodes = a0[:, None] + span * (NODES[None, :] + 1) / 2
        part = (span * WEIGHTS[None, :] / 2 * self.integrand(nodes)).sum(axis=1)
        return self.Gedge[j] + part

    def beta(self, s, t):
        return self.pref * (self.G(t) - self.G(s))


def bmat(b: complex) -> np.ndarray:
    return np.array([[b.real, b.imag], [b.imag, -b.real]])


def sym_pow(A: np.ndarray, p: float) -> np.ndarray:
    w, V = np.linalg.eigh(A)
    return V @ np.diag(w**p) @ V.T


def frame(wv: Wave, u: float):
    g = wv.g
    m = float(g.m(-0.5, u))
    bL = complex(wv.beta(-0.5, u)[0])
    bR = complex(wv.beta(u, 0.5)[0])
    ML = m * np.eye(2) - bmat(bL)
    MR = (g.M - m) * np.eye(2) - bmat(bR)
    W0 = m * (g.M - m) / g.M
    Q = np.linalg.inv(ML) + np.linalg.inv(MR)
    T = sym_pow(W0 * Q, -0.5)
    return dict(m=m, bL=bL, bR=bR, ML=ML, MR=MR, W0=W0, Q=Q, T=T)


def to_phys(wv: Wave, u: float, vt: float, z: np.ndarray):
    fr = frame(wv, u)
    y = fr["T"] @ z
    v = vt - z @ z / (2 * fr["m"]) + y @ np.linalg.solve(fr["ML"], y) / 2
    return v, y, fr


def phi_rep(wv: Wave, u: float) -> complex:
    """Complex representative of Phi_lam(u) = -(B_L - t B_tot)/(m (M - m))."""
    g = wv.g
    m = float(g.m(-0.5, u))
    bL = complex(wv.beta(-0.5, u)[0])
    btot = complex(wv.beta(-0.5, 0.5)[0])
    return -(bL - (m / g.M) * btot) / (m * (g.M - m))


def ell(geo: Geo, lam: float, phi: float, u: float, tau: float) -> complex:
    c = math.sqrt(1 + 4 * geo.q * lam**2)
    m = float(geo.m(-0.5, u))
    pref = 2j * geo.k * float(geo.f(u)) / (c * m * (geo.M - m))
    return pref * np.exp(1j * (u / lam + phi)) * (np.exp(1j * tau) - 1)


# ------------------------------------------------------------ calibration


def C_r_mp(r: int, z):
    return mp.hyp0f1(r + mp.mpf(1) / 2, -mp.mpf(z) / 4) / mp.factorial(2 * r - 1)


def p0_mp(z) -> float:
    return float(C_r_mp(1, z) * C_r_mp(3, z) / (3 * C_r_mp(2, z) ** 2))


def C_r_fast(r: int, z):
    """0F1(r + 1/2; -z/4)/(2r - 1)! by its series; z in [0, 1]."""
    z = np.asarray(z, dtype=float)
    b = r + 0.5
    term = np.ones_like(z)
    total = np.ones_like(z)
    for n in range(1, 25):
        term = term * (-z / 4) / ((b + n - 1) * n)
        total = total + term
    return total / math.factorial(2 * r - 1)


def p0_fast(z):
    return C_r_fast(1, z) * C_r_fast(3, z) / (3 * C_r_fast(2, z) ** 2)


# ------------------------------------------------------------------ N1


def n1() -> None:
    p1, p25 = p0_mp(1), p0_mp(0.25)
    w1, w25 = (0.1 - p1) / 0.1, (0.1 - p25) / 0.1
    ok = abs(p1 - 0.095919) < 2e-6 and abs(p25 - 0.099031) < 2e-6
    ok = ok and abs(w1 - 0.041) < 5e-4 and abs(w25 - 0.010) < 5e-4
    record(
        "N1 window",
        ok,
        f"p0(1) = {p1:.7f} (width {100 * w1:.3f}%), "
        f"p0(0.25) = {p25:.7f} (width {100 * w25:.3f}%)",
    )
    zs = np.linspace(0, 1, 9)
    err = max(abs(float(p0_fast(z)) - p0_mp(z)) for z in zs)
    record("N1 fast series", err < 1e-14, f"max |fast - mpmath| = {err:.2e}")


# ------------------------------------------------------------------ N2


def n2() -> list[dict]:
    rows = [json.loads(line) for line in ROWS.read_text().splitlines() if line]
    nested = True
    table = []
    for row in rows:
        (pX, fX, cX), (pY, fY, cY) = row["moments"]
        aX, aY = pX * pX * cX, pY * pY * cY
        nested = nested and pY >= pX and fY <= fX and aY >= aX
        lo, hi = aX / pY**2, aY / pX**2
        p0q = p0_mp(row["q"])
        inside = p0q < lo and hi < 0.1
        table.append(dict(q=row["q"], r=row["r"], lo=lo, hi=hi, inside=inside))
        info(
            "N2 row",
            f"q={row['q']}, r={row['r']}: [{lo:.6f}, {hi:.6f}] vs "
            f"({p0q:.6f}, 0.1) -> {'evaluated' if inside else 'SKIPPED'}",
        )
    record("N2 nestedness of the 12 rows", nested, "p up, f down, a up in every row")
    want = {
        (0.25, 1024): (0.097688, 0.101401, False),
        (0.25, 2048): (0.098604, 0.100462, False),
        (1.0, 2048): (0.097139, 0.099002, True),
    }
    ok, worst = True, 0.0
    for t in table:
        key = (t["q"], t["r"])
        if key in want:
            lo, hi, inside = want[key]
            worst = max(worst, abs(t["lo"] - lo), abs(t["hi"] - hi))
            ok = ok and t["inside"] == inside
    record(
        "N2 reproduces the referee's X4",
        ok and worst < 1e-5,
        f"max endpoint difference {worst:.2e}",
    )
    skipped = [(t["q"], t["r"]) for t in table if not t["inside"]]
    info("N2 skipped at every reference count", str(skipped))
    return rows


# ------------------------------------------------------------------ N3


def n3() -> None:
    Ns = set(range(4, 3001))
    for j in range(2, 41):
        Ns.update({2**j - 1, 2**j, 2**j + 1})
    Ns.update(10**j for j in range(1, 13))
    ok_box = ok_log = ok_union = ok_pow = True
    for N in sorted(Ns):
        x = 32 * N**4
        bl, bi = (x - 1).bit_length(), x.bit_length()
        for B in (bl, bi):
            ok_box = ok_box and x <= 2**B <= 2 * x
            ok_log = ok_log and B <= 6 + 4 * math.log2(N) + 1e-9
            lhs = math.log(6 * N + 3 * N * N) + math.log(2) - B
            ok_union = ok_union and lhs <= math.log(0.75) - 2 * math.log(N)
        is_pow = N & (N - 1) == 0
        ok_pow = ok_pow and ((bi == bl + 1) == is_pow) and bi in (bl, bl + 1)
    record("N3 32N^4 <= 2^B <= 64N^4", ok_box, f"{len(Ns)} values of N, both B")
    record("N3 B <= 6 + 4 log2 N", ok_log, "both B")
    record("N3 union bound (3/4) N^-2", ok_union, "both B")
    record("N3 implementation's B", ok_pow, "least B, plus one iff N = 2^j")


# ------------------------------------------------------------------ N4


def n4() -> None:
    worst, ok = 0.0, True
    for kk in (10, 40, 160):
        for d in (0.001, 0.02, 0.1, 0.3, 0.5, 0.8):
            for B in (1, 3, 6):
                rH = math.sqrt(B / (2 * kk))
                fail = 0.0
                for j in range(kk + 1):
                    xbar = j / kk
                    dU = min(1.0, xbar + rH)
                    rD = min(rH, math.sqrt(2 * dU * B / kk) + 2 * B / (3 * kk))
                    if abs(xbar - d) > rD:
                        logp = (
                            math.lgamma(kk + 1)
                            - math.lgamma(j + 1)
                            - math.lgamma(kk - j + 1)
                            + j * math.log(d)
                            + (kk - j) * math.log1p(-d)
                        )
                        fail += math.exp(logp)
                ratio = fail / (4 * math.exp(-B))
                worst = max(worst, ratio)
                ok = ok and ratio <= 1.0
    record("N4 Lemma 2.1 exact binomial", ok, f"worst failure / 4e^-B = {worst:.3f}")


# ------------------------------------------------------------------ N5


def n5() -> None:
    rng = np.random.default_rng(5)
    ok_id = True
    for n in range(2, 8):
        h = {}
        for j in range(n):
            for jj in range(j + 1, n):
                h[(j, jj)] = h[(jj, j)] = int(rng.integers(0, 10))
        U = Fr(
            sum(h[(j, jj)] for j in range(n) for jj in range(j + 1, n)), math.comb(n, 2)
        )
        kk = n // 2
        total, count = Fr(0), 0
        for perm in itertools.permutations(range(n)):
            total += Fr(sum(h[(perm[2 * i], perm[2 * i + 1])] for i in range(kk)), kk)
            count += 1
        ok_id = ok_id and total / count == U
    record("N5 averaging identity", ok_id, "n = 2..7, exact")

    probs = [0.5, 1 / 3, 1 / 6]
    H = [[1, 0, 1], [0, 0, 1], [1, 1, 1]]
    d = sum(probs[i] * probs[j] * H[i][j] for i in range(3) for j in range(3))
    ok_mgf, worst = True, 0.0
    for n in (4, 5):
        kk = n // 2
        outcomes = list(itertools.product(range(3), repeat=n))
        weights = [math.prod(probs[i] for i in o) for o in outcomes]
        Us = [
            sum(H[o[j]][o[jj]] for j in range(n) for jj in range(j + 1, n))
            / math.comb(n, 2)
            for o in outcomes
        ]
        for s in (-8, -2, 2, 8):
            eu = sum(w * math.exp(s * u) for w, u in zip(weights, Us, strict=True))
            ev = (1 - d + d * math.exp(s / kk)) ** kk
            worst = max(worst, eu / ev)
            ok_mgf = ok_mgf and eu <= ev * (1 + 1e-12)
    record("N5 E e^(sU) <= E e^(sV)", ok_mgf, f"largest ratio {worst:.4f}")


# ------------------------------------------------------------------ N6


def n6() -> None:
    rng = np.random.default_rng(6)
    al = np.array([1.3, -0.7, 2.1])
    be = np.array([[0, 0.8, -1.1], [0.8, 0, 0.5], [-1.1, 0.5, 0]])
    ga = 1.7

    def F(x):
        return (
            np.sum(np.exp(al * x))
            + sum(be[i, j] * x[i] * x[j] for i in range(3) for j in range(i + 1, 3))
            + ga * x[0] * x[1] * x[2]
        )

    def absmax(lo, hi):
        return np.maximum(np.abs(lo), np.abs(hi))

    ok, worst = True, 0.0
    for _ in range(20000):
        b = rng.uniform(0.1, 0.9, 3)
        d = rng.uniform(0.0, 0.2, 3) * np.array([1, -1, 1])
        rb = rng.uniform(0, 0.05, 3)
        rd = rng.uniform(0, 0.05, 3)
        bh = b + rng.uniform(-1, 1, 3) * rb
        dh = d + rng.uniform(-1, 1, 3) * rd
        blo, bhi = bh - rb, bh + rb
        dlo, dhi = dh - rd, dh + rd
        lo = blo + np.minimum(0, dlo)
        hi = bhi + np.maximum(0, dhi)
        X = absmax(lo, hi)
        L = np.empty(3)
        for i in range(3):
            o = [j for j in range(3) if j != i]
            L[i] = (
                np.max(np.abs(al[i] * np.exp(al[i] * np.array([lo[i], hi[i]]))))
                + sum(abs(be[i, j]) * X[j] for j in o)
                + abs(ga) * X[o[0]] * X[o[1]]
            )
        Hm = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                if i == j:
                    Hm[i, j] = np.max(
                        al[i] ** 2 * np.exp(al[i] * np.array([lo[i], hi[i]]))
                    )
                else:
                    kx = 3 - i - j
                    Hm[i, j] = abs(be[i, j]) + abs(ga) * X[kx]
        Dj = absmax(dlo, dhi)
        lhs = abs((F(bh + dh) - F(bh)) - (F(b + d) - F(b)))
        rhs = float(L @ rd + rb @ Hm @ Dj)
        worst = max(worst, lhs / rhs)
        ok = ok and lhs <= rhs * (1 + 1e-12)
    record("N6 Lemma 2.4 on 20000 boxes", ok, f"largest lhs/rhs = {worst:.3f}")


# ------------------------------------------------------------------ N7


def c_L(geo: Geo) -> float:
    return 8 * math.sqrt(2) * geo.k / (1.008 * geo.M**2)


def n7() -> None:
    us = np.linspace(-0.25, 0.25, 9)
    taus = (0.3, math.pi / 2, math.pi, 3 * math.pi / 2, 6.0)
    ok_all = True
    for q in (0.25, 1.0):
        geo = Geo(q)
        cl = c_L(geo)
        mx = {}
        lmin = math.inf
        for lam in (1 / 256, 1 / 1024, 1 / 4096):
            worst = 0.0
            for phi in (0.0, 0.37, 2.1):
                wv = Wave(geo, lam, phi)
                for u in us:
                    base = phi_rep(wv, float(u))
                    for tau in taus:
                        diff = phi_rep(wv, float(u) + tau * lam) - base
                        el = ell(geo, lam, phi, float(u), tau)
                        worst = max(worst, abs(diff - lam**2 * el) / lam**3)
                        if math.pi / 2 - 1e-12 <= tau <= 3 * math.pi / 2 + 1e-12:
                            lmin = min(lmin, abs(el) / cl)
            mx[lam] = worst
        ok = all(v <= 3000 for v in mx.values()) and mx[1 / 4096] <= 2 * mx[1 / 256]
        ok = ok and lmin >= 1.0
        ok_all = ok_all and ok
        record(
            f"N7a Lemma 3.1 q={q}",
            ok,
            "max residual/lam^3: "
            + ", ".join(f"1/{round(1 / k)}: {v:.1f}" for k, v in mx.items())
            + f"; min |l|/c_L = {lmin:.3f} (c_L = {cl:.4f})",
        )
    for q in (0.25, 1.0):
        geo = Geo(q)
        res = {}
        for lam in (1 / 256, 1 / 1024):
            wv = Wave(geo, lam, 0.37)
            worst = 0.0
            for u in (-0.25, 0.0, 0.25):
                ph = bmat(phi_rep(wv, u))
                for th in np.linspace(0, math.pi, 6, endpoint=False):
                    z = 0.2 * np.array([math.cos(th), math.sin(th)])
                    v, _, _ = to_phys(wv, u, 0.5, z)
                    worst = max(worst, abs((v - 0.5) + z @ ph @ z / 2) / lam**4)
            res[lam] = worst
        record(
            f"N7b sign of (5') q={q}",
            all(v <= 500 for v in res.values()),
            "max |dv + z Phi z/2|/lam^4: "
            + ", ".join(f"1/{round(1 / k)}: {v:.2f}" for k, v in res.items()),
        )


# ------------------------------------------------------------------ N8


RHO1, RHO2 = 0.1, 0.2


def constants(geo: Geo):
    cl = c_L(geo)
    g1, g2 = cl * RHO1**2 / 16, cl * RHO1**2 / 8
    cp = math.sqrt(3 * math.pi * cl * RHO1**2 / 128)
    return cl, g1, g2, cp


def sample_pairs(geo: Geo, lam: float, phi: float, n: int, rng, corners: bool):
    """Canonical witness pairs (u1, vt1, z1, u2, vt2, z2)."""
    cl, g1, g2, cp = constants(geo)
    out = []
    specs = []
    for _ in range(n):
        specs.append(
            dict(
                u1=rng.uniform(-0.25, 0.25),
                tau=rng.uniform(math.pi / 2, 3 * math.pi / 2),
                eta=rng.uniform(0.3, 0.4),
                r=math.sqrt(rng.uniform(RHO1**2, RHO2**2)),
                off=rng.uniform(-math.pi / 6, math.pi / 6)
                + math.pi * rng.integers(0, 2),
                wr=cp * lam**1.5 * math.sqrt(rng.uniform()),
                wth=rng.uniform(0, 2 * math.pi),
                dv=-rng.uniform(g1, g2) * lam**2,
            )
        )
    if corners:
        for u1, tau, off in itertools.product(
            (-0.25, 0.25), (math.pi / 2, 3 * math.pi / 2), (-math.pi / 6, math.pi / 6)
        ):
            for wth_rel in (0.0, math.pi / 2):
                specs.append(
                    dict(
                        u1=u1,
                        tau=tau,
                        eta=0.4,
                        r=RHO1,
                        off=off,
                        wr=cp * lam**1.5,
                        wth=None,
                        wth_rel=wth_rel,
                        dv=-g2 * lam**2,
                    )
                )
    for sp in specs:
        u1, tau = sp["u1"], sp["tau"]
        el = ell(geo, lam, phi, u1, tau)
        psi_minus = np.angle(el) / 2 + math.pi / 2
        psi = psi_minus + sp["off"]
        z1 = sp["r"] * np.array([math.cos(psi), math.sin(psi)])
        wth = sp["wth"] if sp["wth"] is not None else psi + sp["wth_rel"]
        z2 = z1 + sp["wr"] * np.array([math.cos(wth), math.sin(wth)])
        m1 = float(geo.m(-0.5, u1))
        vt1 = sp["eta"] + z1 @ z1 / (2 * m1)
        out.append((u1, vt1, z1, u1 + tau * lam, vt1 + sp["dv"], z2))
    return out


def in_K(geo: Geo, u, vt, z) -> bool:
    m = float(geo.m(-0.5, u))
    eta = vt - z @ z / (2 * m)
    zeta = 1 - vt - z @ z / (2 * (geo.M - m))
    return abs(u) <= 0.3 and eta >= 0.25 and zeta >= 0.4


def n8() -> None:
    for q in (0.25, 1.0):
        geo = Geo(q)
        cl, g1, g2, cp = constants(geo)
        info(
            f"N8 constants q={q}",
            f"c_L = {cl:.5f}, g1 = {g1:.6f}, g2 = {g2:.6f}, c' = {cp:.5f}, "
            f"floor (3/64) c_L rho1^2 = {3 / 64 * cl * RHO1**2:.6f}",
        )
        rng = np.random.default_rng(8 + int(4 * q))
        for lam in (1 / 512, 1 / 2048, 1 / 8192):
            related, kmin, inK, total = 0, math.inf, True, 0
            for phi in (0.37, 2.1):
                wv = Wave(geo, lam, phi)
                for u1, vt1, z1, u2, vt2, z2 in sample_pairs(
                    geo, lam, phi, 300, rng, True
                ):
                    v1, y1, _ = to_phys(wv, u1, vt1, z1)
                    v2, y2, _ = to_phys(wv, u2, vt2, z2)
                    b12 = complex(wv.beta(u1, u2)[0])
                    M12 = float(geo.m(u1, u2)) * np.eye(2) - bmat(b12)
                    dy = y2 - y1
                    kap = (v2 - v1) - dy @ np.linalg.solve(M12, dy) / 2
                    related += kap > 0
                    kmin = min(kmin, kap / lam**2)
                    inK = inK and in_K(geo, u1, vt1, z1) and in_K(geo, u2, vt2, z2)
                    total += 1
            text = (
                f"lam=1/{round(1 / lam)}: related {related}/{total}, "
                f"min kappa/lam^2 = {kmin:.6f}, both in K: {inK}"
            )
            if lam == 1 / 512:
                info(f"N8 q={q}", text)
                continue
            ok = related == total and inK
            if lam == 1 / 8192:
                ok = ok and kmin >= 3 / 64 * cl * RHO1**2
            record(f"N8 Lemma 3.3 q={q}", ok, text)


# ------------------------------------------------------------------ N9


def theta0(geo: Geo, u, vt, z):
    m = geo.m(-0.5, u)
    s = u + 0.5
    zz = np.sum(z * z, axis=-1)
    eta = vt - zz / (2 * m)
    zeta = 1 - vt - zz / (2 * (geo.M - m))
    p = eta**2 * geo.alpha(s) / geo.gamma
    f = zeta**2 * geo.alpha(1 - s) / geo.gamma
    a = p**2 * p0_fast(geo.q * s**2)
    return np.stack([p, f, a], axis=-1)


def n9() -> None:
    for q in (0.25, 1.0):
        geo = Geo(q)
        p0q = float(p0_fast(q))
        lo_c, hi_c = float(p0_fast(0.64 * q)), float(p0_fast(0.04 * q))
        record(
            f"N9 c_0(K) inside the window q={q}",
            p0q < lo_c and hi_c < 0.1,
            f"[{lo_c:.6f}, {hi_c:.6f}] in ({p0q:.6f}, 0.1); margins "
            f"{lo_c - p0q:.2e}, {0.1 - hi_c:.2e}",
        )
        rng = np.random.default_rng(90 + int(4 * q))
        ratios, fracs = {}, {}
        for lam in (2.0**-11, 2.0**-14, 2.0**-17, 2.0**-20):
            pairs = sample_pairs(geo, lam, 0.37, 2000, rng, False)
            X = np.array([theta0(geo, u1, vt1, z1) for u1, vt1, z1, *_ in pairs])
            Y = np.array([theta0(geo, u2, vt2, z2) for *_, u2, vt2, z2 in pairs])
            ratios[lam] = float(np.max(np.abs(Y - X)) / lam)
            amin, amax = np.minimum(X[:, 2], Y[:, 2]), np.maximum(X[:, 2], Y[:, 2])
            pmin, pmax = np.minimum(X[:, 0], Y[:, 0]), np.maximum(X[:, 0], Y[:, 0])
            inside = (amin / pmax**2 > p0q) & (amax / pmin**2 < 0.1)
            fracs[lam] = float(np.mean(inside))
        ok = all(v <= 20 for v in ratios.values())
        ok = ok and ratios[2.0**-20] <= 1.5 * ratios[2.0**-14]
        ok = ok and ratios[2.0**-14] <= 1.5 * ratios[2.0**-20]
        record(
            f"N9 increments q={q}",
            ok,
            "max |Delta theta_0|/lam: "
            + ", ".join(f"2^{round(math.log2(k))}: {v:.3f}" for k, v in ratios.items()),
        )
        record(
            f"N9 limit boxes in the window q={q}",
            fracs[2.0**-20] == 1.0,
            "fraction inside: "
            + ", ".join(f"2^{round(math.log2(k))}: {v:.3f}" for k, v in fracs.items()),
        )


# ------------------------------------------------------------------ N10


def n10(rows: list[dict]) -> None:
    for row in rows:
        lam = 1 / row["r"]
        val = abs(row["population_v_decrease"] - row["wave_margin"]) / (
            lam**4 * math.log(1 / lam)
        )
        info("N10", f"q={row['q']}, r={row['r']}: {val:.1f}")


def main() -> int:
    n1()
    rows = n2()
    n3()
    n4()
    n5()
    n6()
    n7()
    n8()
    n9()
    n10(rows)
    fails = [r for r in RESULTS if not r[1]]
    print(f"\n{len(RESULTS) - len(fails)} PASS, {len(fails)} FAIL")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
