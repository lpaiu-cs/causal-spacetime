"""Referee check (b): ||E|| <= C lam min(ell, lam) and the uniform sandwich (4).

Referee's explicit constants (derived by hand in the referee report):
    |int_s^t f e^{iu/lam} du| <= min( m(s,t), lam (f(s)+f(t)+TV_[s,t] f) )
    e = ||E|| <= C_E lam min(ell, lam),  C_E = 2 sqrt(q) (4 fmax - 2),
    e/b <= tanh(2 beta) < 1,
    c0 = 2 sqrt(2) C_E,  and pointwise on the wave null cone |margin| <= 2 sqrt(2) e.

Tests
  B0  helper vs direct quadrature of int gamma^-1 du / c from the definition.
  B1  oscillatory-integral bound, including ell down to 1e-7.
  B2  e <= C_E lam min(ell,lam) and e/b <= tanh(2 beta).
  B3  boundary-targeted: put Y exactly on the WAVE null cone of X and measure the
      background null margin g = dT - |dz|.  (4) holds with constant c0 iff
      sup |g| <= c0 delta, because g is monotone in dv.  Reported per decade
      of u-separation: the bound must not degrade for short segments.
  B4  iid reference pairs in the diamond, both time orientations: count pairs
      violating (4) (must be exactly zero); count genuine disagreements (must be
      positive, otherwise the test is vacuous); T is a common time function.

Deterministic (fixed seeds).  Run: python check_b_sandwich.py
"""

from __future__ import annotations

import numpy as np
from audit_common import SQRT2, Family, Reporter

rep = Reporter("check (b): operator bound and uniform sandwich (4)")

QS = (0.25, 1.0, 4.0, 9.0)
LAMS = (1 / 4, 1 / 16, 1 / 64, 1 / 256, 1 / 1024, 1 / 4096)


def tv_f(fam, s, t):
    """Total variation of f = sec^2(sqrt(q) u) on [s,t] (f decreases then increases)."""
    fs, ft = fam.f(s), fam.f(t)
    straddle = (s < 0) & (t > 0)
    return np.where(straddle, fs + ft - 2.0, np.abs(ft - fs))


# ---------------------------------------------------------------- B0
rng = np.random.default_rng(1)
worst = 0.0
for q in (0.25, 4.0, 9.0):
    for lam in (1 / 4, 1 / 64):
        fam = Family(q, lam, phi=0.37)
        for _ in range(4):
            s, t = np.sort(rng.uniform(-0.5, 0.5, 2))
            nn = 200001
            uu = np.linspace(s, t, nn)
            w = np.full(nn, (t - s) / (nn - 1))
            w[0] *= 0.5
            w[-1] *= 0.5
            beta = 0.5 * np.arcsinh(2 * lam * np.sqrt(q))
            th = 0.5 * (uu / lam + fam.phi)
            a2 = np.cos(np.sqrt(q) * uu) ** 2
            ct, st = np.cos(th), np.sin(th)
            e2b, em2b = np.exp(2 * beta), np.exp(-2 * beta)
            # gamma = a^2 R diag(e2b, em2b) R^T straight from the definition
            g11 = a2 * (ct * ct * e2b + st * st * em2b)
            g22 = a2 * (st * st * e2b + ct * ct * em2b)
            g12 = a2 * (ct * st * (e2b - em2b))
            det = g11 * g22 - g12 * g12
            M11 = (g22 / det * w).sum() / fam.c
            M22 = (g11 / det * w).sum() / fam.c
            M12 = (-g12 / det * w).sum() / fam.c
            b = fam.m(t) - fam.m(s)
            G = fam.pair_G(fam.prim(np.array([s]))[0], fam.prim(np.array([t]))[0])
            k = fam.s / fam.c
            H11, H22, H12 = b - k * G.real, b + k * G.real, -k * G.imag
            worst = max(worst, abs(M11 - H11), abs(M22 - H22), abs(M12 - H12))
rep.check(
    "B0 helper Mtilde = b I + E matches direct trapezoid of gamma^-1/c from the "
    "definition",
    worst < 2e-8,
    f"max abs diff {worst:.2e}",
)

# ---------------------------------------------------------------- B1, B2, B3
rng = np.random.default_rng(20261124)
NPAIR = 120000
print(
    "INFO B3 table: sup over wave-null pairs of |background margin| / delta, by "
    "decade of u-separation"
)
for q in QS:
    for lam in LAMS:
        fam = Family(q, lam, phi=float(rng.uniform(0, 2 * np.pi)))
        ell = 10.0 ** rng.uniform(-7.0, 0.0, NPAIR)
        ell = np.minimum(ell, 1.0)
        u1 = -0.5 + rng.uniform(0.0, 1.0, NPAIR) * (1.0 - ell)
        u2 = np.minimum(u1 + ell, 0.5)
        ell = u2 - u1
        seg = fam.seg(u1, u2)
        G = np.exp(1j * fam.phi) * seg
        # b = m(u1,u2) from the tangent subtraction formula (no cancellation)
        b = np.sin(fam.rq * ell) / (fam.rq * np.cos(fam.rq * u1) * np.cos(fam.rq * u2))
        absI = np.abs(seg)
        bound_I = np.minimum(b, lam * (fam.f(u1) + fam.f(u2) + tv_f(fam, u1, u2)))
        tol = 1e-9 * bound_I + 1e-13 * fam.fmax
        ok1 = bool(np.all(absI <= bound_I + tol))
        e = fam.tanh2b * absI
        ok2 = bool(
            np.all(e <= fam.C_E * lam * np.minimum(ell, lam) * (1 + 1e-9))
        ) and bool(np.all(e <= fam.tanh2b * (b + tol)))
        # ---- B3: Y on the wave null cone
        ang = rng.uniform(0, 2 * np.pi, NPAIR)
        n1, n2 = np.cos(ang), np.sin(ang)
        dv = 10.0 ** rng.uniform(-6.0, 0.0, NPAIR)  # dv <= slack h = 1
        qn = fam.wave_form(b, G, n1, n2)  # n^T (bI+E)^-1 n > 0
        r2 = 2.0 * dv / qn
        dT = (b + dv) / SQRT2
        dz = np.sqrt((b - dv) ** 2 / 2.0 + r2)
        g = (2.0 * b * dv - r2) / (dT + dz)
        ratio = np.abs(g) / fam.delta
        ok3a = bool(np.all(np.abs(g) <= 2.0 * SQRT2 * e * (1 + 1e-7) + 1e-18))
        ok3b = bool(ratio.max() <= fam.c0)
        dec = np.clip(np.floor(np.log10(ell)).astype(int), -7, -1)
        cells = []
        for d in range(-7, 0):
            sel = dec == d
            top = ratio[sel].max() if sel.any() else float("nan")
            cells.append(f"1e{d}:{top:.3g}")
        tag = f"q={q:.2f} lam=1/{round(1 / lam)}"
        rep.check(
            "B1 oscillatory bound " + tag,
            ok1,
            f"max |I|/bound {float((absI / bound_I).max()):.4f}",
        )
        rep.check(
            "B2 e <= C_E lam min(ell,lam), e/b <= tanh(2b) " + tag,
            ok2,
            f"C_E={fam.C_E:.3f} max e/(lam min) "
            f"{float((e / (lam * np.minimum(ell, lam))).max()):.4f}",
        )
        rep.check(
            "B3 wave-null margin <= 2 sqrt2 e and <= c0 delta " + tag,
            ok3a and ok3b,
            "c0={:.2f} sup|g|/delta={:.4f} [{}]".format(
                fam.c0, float(ratio.max()), " ".join(cells)
            ),
        )

# ---------------------------------------------------------------- B4
rng = np.random.default_rng(424242)
NPTS = 1500
for q in (0.25, 1.0, 4.0):
    for lam in (1 / 8, 1 / 32, 1 / 128):
        fam = Family(q, lam, phi=float(rng.uniform(0, 2 * np.pi)))
        pts = fam.sample_reference(NPTS, rng)
        X = fam.to_minkowski(pts)
        P = fam.prim(pts[:, 0])
        mm = fam.m(pts[:, 0])
        i, j = np.triu_indices(NPTS, k=1)
        # orient each pair so that u_i < u_j ("forward" candidate i -> j)
        swap = pts[i, 0] > pts[j, 0]
        i, j = np.where(swap, j, i), np.where(swap, i, j)
        b = mm[j] - mm[i]
        G = fam.pair_G(P[i], P[j])
        dv = pts[j, 1] - pts[i, 1]
        dy1 = pts[j, 2] - pts[i, 2]
        dy2 = pts[j, 3] - pts[i, 3]
        wave_fwd = 2.0 * dv >= fam.wave_form(b, G, dy1, dy2)
        dT = X[j, 0] - X[i, 0]
        dzn = np.sqrt(((X[j, 1:] - X[i, 1:]) ** 2).sum(axis=1))
        g_fwd = dT - dzn
        g_bwd = -dT - dzn
        back_fwd = g_fwd >= 0.0
        c0d = fam.c0 * fam.delta
        # (4) forward: {g >= c0 d} subset wave subset {g >= -c0 d}
        viol_fwd = int(
            np.sum((g_fwd >= c0d) & ~wave_fwd) + np.sum(wave_fwd & (g_fwd < -c0d))
        )
        # (4) backward orientation j -> i: wave relation is impossible since u_j > u_i,
        # so the sandwich requires g_bwd < c0 d for every such pair.
        viol_bwd = int(np.sum(g_bwd >= c0d))
        disagree = wave_fwd != back_fwd
        ndis = int(disagree.sum())
        emp_c0 = float((np.abs(g_fwd[disagree]) / fam.delta).max()) if ndis else 0.0
        time_fn = bool(np.all(dT[wave_fwd] > 0.0))
        tag = f"q={q:.2f} lam=1/{round(1 / lam)}"
        rep.check(
            "B4 no pair violates (4), either orientation, T is a time function " + tag,
            viol_fwd == 0 and viol_bwd == 0 and time_fn and ndis > 0,
            f"pairs={len(i)} disagreements={ndis} "
            f"(frac/lam^2={ndis / len(i) / fam.delta:.3f}) max|g|/delta on "
            f"disagreements={emp_c0:.3f} vs c0={fam.c0:.2f}",
        )

rep.finish()
