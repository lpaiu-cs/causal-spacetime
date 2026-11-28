"""Deterministic checks for one_point_strip_bound.md (the one-point bound written out).

No sampling. Every number below is a deterministic quadrature.

Setting (one_point_strip_bound.md Section 1). a = cos(sqrt(q) u), f = a^-2,
m_-(u) = int_{-1/2}^u f, m_+(u) = int_u^{1/2} f, M = m_- + m_+. Tips
(u, v, y) = (-1/2, 0, 0) and (1/2, h, 0), h = 1. E_- = E_lam(-1/2, u) and
E_+ = E_lam(u, 1/2) have eigenvalues +-e_-, +-e_+, with
e = t |G|, G = e^{i phi} int f e^{i w / lam} dw, t = 2 lam sqrt(q)/sqrt(1 + 4 q lam^2).
For a unit vector n(theta): A = n^T (m I + E)^-1 n, B = 1/m, dA = A - B,
S_lam = A_+ + A_-, S_0 = B_+ + B_-, S_cap = max(A_+, B_+) + max(A_-, B_-).

Pass criteria (written before the full run; one smoke case, q = 1,
lam = 1/8, was run first to debug the script):
  K1  At every quadrature node (u, theta) of every case:
        e <= t m and e <= t lam (4 f_max - 2) <= C_E lam^2;
        |dA| <= e / (m (m - e));
        A >= 1 / (m + e).
      All up to a relative slack of 1e-9.
  K2  Polar reduction. At 5 values of u, for 2 cases, the y-plane integrals of
      the v-section lengths of D_lam, D_0 and D_lam cap D_0 on a Cartesian grid
      agree with (h^2/2) int dtheta 1/S to relative 1e-3.
  K3  theta quadrature. V_lam - V_0 from the theta grid agrees with the closed
      form pi h^2 int a^2 [1/sqrt(P^2 - |Z|^2) - 1/S_0] du, to within
      1e-6 |V_lam - V_0| + 1e-12 V_0.
  K4  V_0 from quadrature agrees with the closed form pi h^2 Vhat to 1e-9
      relative, and W, Vhat agree with their closed forms to 1e-9.
  K5  Main. For every case:
        TV <= B(lam, phi) <= C_1(q, t) lam^2 <= C_1(q) lam^2,
      and Delta_lam, Delta_0 each lie below their separate bounds.
  K6  Resolution. Doubling the panels and the theta grid changes TV by less
      than 1e-3 relative (q in the case list, lam in {2^-3, 2^-8}, phi = 0.3).
  K7  Agreement with the previous referee's R8 values at phi = 0.3
      (argument1_revision/review.md): TV/lam^2 within 1 percent at
      lam = 1/8 and 1/64 for q = 0.25, 1 and 4.
Cases for K1, K3 and K5: q in {0.25, 1, 4, 9, 9.8}, phi in {k pi/4, k = 0..7},
lam = 2^(-k/2) for k = 0..20 (1 down to 1/1024).

Run: python checks.py
"""

from __future__ import annotations

import math

import numpy as np

failed = 0
npass = 0


def check(name, ok, detail=""):
    global failed, npass
    print(
        ("PASS " if ok else "FAIL ") + name + (" | " + detail if detail else ""),
        flush=True,
    )
    if ok:
        npass += 1
    else:
        failed += 1


GX, GW = np.polynomial.legendre.leggauss(12)
H = 1.0


# ---------------------------------------------------------------- constants
def consts(q):
    s = math.sqrt(q)
    th = math.tan(s / 2)
    fmax = 1.0 / math.cos(s / 2) ** 2
    ce = 2 * s * (4 * fmax - 2)
    i_c = 0.5 + math.sin(s) / (2 * s)  # int a^2
    i_s = 0.5 - math.sin(s) / (2 * s)  # int sin^2(s u)
    big_m = 2 * th / s
    vhat = (th * th * i_c - i_s) / (2 * s * th)  # int a^2 m_+ m_- / M
    w = (i_s + th * th * i_c) / (2 * th * th)  # int a^2 (m_+^2 + m_-^2) / M^2
    t1 = 2 * s / math.sqrt(1 + 4 * q)
    return dict(s=s, th=th, fmax=fmax, ce=ce, M=big_m, vhat=vhat, w=w, t1=t1)


def t_of(q, lam):
    return 2 * lam * math.sqrt(q) / math.sqrt(1 + 4 * q * lam * lam)


def c1_of(cst, t):
    return (1 + t) / (1 - t) * cst["ce"] * cst["w"] / cst["vhat"]


# ---------------------------------------------------------------- quadrature in u
def u_rule(q, lam, factor=1):
    """Composite 12-point Gauss nodes on [-1/2, 1/2] and the primitive
    F(u) = int_{-1/2}^u f(w) e^{i w/lam} dw at every node."""
    s = math.sqrt(q)
    d = (math.pi / 2 - s / 2) / s  # distance from u = 1/2 to the pole of f
    width = min(lam / 6, d / 2, 1 / 32) / factor
    npan = math.ceil(1 / width)
    edges = np.linspace(-0.5, 0.5, npan + 1)
    a, b = edges[:-1], edges[1:]
    half, mid = (b - a) / 2, (a + b) / 2
    x = mid[:, None] + half[:, None] * GX[None, :]
    wx = half[:, None] * GW[None, :]

    def g(w):
        return np.exp(1j * w / lam) / np.cos(s * w) ** 2

    pan = np.sum(wx * g(x), axis=1)
    fstart = np.concatenate([[0.0 + 0.0j], np.cumsum(pan)[:-1]])
    hh = (x - a[:, None]) / 2
    inner_nodes = a[:, None, None] + hh[:, :, None] * (1 + GX[None, None, :])
    inner = np.sum(hh[:, :, None] * GW[None, None, :] * g(inner_nodes), axis=2)
    f_at = fstart[:, None] + inner
    return x.ravel(), wx.ravel(), f_at.ravel(), np.sum(pan)


def pos(z):
    return np.maximum(z, 0.0)


# ---------------------------------------------------------------- one case
def one_case(q, lam, phi, factor=1, nth=256, pointwise=False, chunk=8192):
    cst = consts(q)
    s, th, big_m = cst["s"], cst["th"], cst["M"]
    t = t_of(q, lam)
    x, wx, f_at, ftot = u_rule(q, lam, factor)
    a2 = np.cos(s * x) ** 2
    m_lo = (np.tan(s * x) + th) / s
    m_hi = (th - np.tan(s * x)) / s
    g_lo = np.exp(1j * phi) * f_at
    g_hi = np.exp(1j * phi) * (ftot - f_at)
    e_lo, e_hi = t * np.abs(g_lo), t * np.abs(g_hi)
    s0 = big_m / (m_lo * m_hi)
    ang = (np.arange(nth) + 0.5) * math.pi / nth  # theta in [0, pi); period pi
    c2, s2 = np.cos(2 * ang), np.sin(2 * ang)
    wth = 2 * math.pi / nth

    d_lam = d_0 = 0.0
    worst = dict(e_tm=0.0, e_ce=0.0, da=0.0, alow=0.0)
    ibp = t * lam * (4 * cst["fmax"] - 2)
    for i0 in range(0, x.size, chunk):
        sl = slice(i0, i0 + chunk)
        das = []
        for m, g, e in ((m_lo[sl], g_lo[sl], e_lo[sl]), (m_hi[sl], g_hi[sl], e_hi[sl])):
            nen = -t * (g.real[:, None] * c2[None, :] + g.imag[:, None] * s2[None, :])
            da = (e[:, None] ** 2 - m[:, None] * nen) / (
                m[:, None] * (m[:, None] ** 2 - e[:, None] ** 2)
            )
            das.append(da)
            if pointwise:
                worst["e_tm"] = max(worst["e_tm"], float(np.max(e / (t * m))))
                worst["e_ce"] = max(worst["e_ce"], float(np.max(e)) / ibp)
                bound = e / (m * (m - e))
                with np.errstate(invalid="ignore", divide="ignore"):
                    r = np.where(bound[:, None] > 0, np.abs(da) / bound[:, None], 0.0)
                worst["da"] = max(worst["da"], float(np.max(r)))
                alow = 1.0 / (m + e)
                worst["alow"] = max(
                    worst["alow"],
                    float(np.max(alow[:, None] / (1.0 / m[:, None] + da))),
                )
        da_lo, da_hi = das
        sz = s0[sl][:, None]
        s_lam = sz + da_lo + da_hi
        s_cap = sz + pos(da_lo) + pos(da_hi)
        neg = pos(-da_lo) + pos(-da_hi)
        pp = pos(da_lo) + pos(da_hi)
        wgt = (wx[sl] * a2[sl])[:, None] * wth * (H * H / 2)
        d_lam += float(np.sum(wgt * neg / (s_lam * s_cap)))
        d_0 += float(np.sum(wgt * pp / (sz * s_cap)))

    v0_num = float(np.sum(wx * a2 / s0)) * 2 * math.pi * (H * H / 2)
    v0 = math.pi * H * H * cst["vhat"]
    v_lam = v0 + d_lam - d_0
    tv = max(d_lam, d_0) / max(v_lam, v0)

    # closed-form theta integral of 1/S_lam (K3)
    den_lo = m_lo**2 - e_lo**2
    den_hi = m_hi**2 - e_hi**2
    p_ = m_lo / den_lo + m_hi / den_hi
    z_ = t * (g_lo / den_lo + g_hi / den_hi)
    dv_closed = (
        math.pi
        * H
        * H
        * float(np.sum(wx * a2 * (1 / np.sqrt(p_ * p_ - np.abs(z_) ** 2) - 1 / s0)))
    )

    jint = float(np.sum(wx * a2 * (e_lo * m_hi**2 + e_hi * m_lo**2))) / big_m**2
    b_lam = math.pi * H * H * (1 + t) / (1 - t) * jint
    b_0 = math.pi * H * H / (1 - t) * jint
    big_b = (1 + t) / (1 - t) * jint / cst["vhat"]
    w_num = float(np.sum(wx * a2 * (m_lo**2 + m_hi**2))) / big_m**2
    vhat_num = float(np.sum(wx * a2 * m_lo * m_hi)) / big_m
    return dict(
        tv=tv,
        d_lam=d_lam,
        d_0=d_0,
        v0=v0,
        v0_num=v0_num,
        v_lam=v_lam,
        dv_closed=dv_closed,
        b_lam=b_lam,
        b_0=b_0,
        B=big_b,
        t=t,
        w_num=w_num,
        vhat_num=vhat_num,
        worst=worst,
        nu=x.size,
    )


# ---------------------------------------------------------------- K2 polar reduction
def polar_check(q, lam, phi, us, ngrid=1601, nth=4096):
    cst = consts(q)
    s, th = cst["s"], cst["th"]
    t = t_of(q, lam)
    worst = 0.0
    ang = (np.arange(nth) + 0.5) * math.pi / nth
    c2, s2 = np.cos(2 * ang), np.sin(2 * ang)
    for u in us:

        def prim(lo, hi):
            npan = max(64, math.ceil((hi - lo) / (lam / 8)))
            ed = np.linspace(lo, hi, npan + 1)
            hlf, md = (ed[1:] - ed[:-1]) / 2, (ed[1:] + ed[:-1]) / 2
            xx = md[:, None] + hlf[:, None] * GX[None, :]
            return np.sum(
                hlf[:, None] * GW[None, :] * np.exp(1j * xx / lam) / np.cos(s * xx) ** 2
            )

        mats, ms, gs = [], [], []
        for lo, hi, m in (
            (-0.5, u, (math.tan(s * u) + th) / s),
            (u, 0.5, (th - math.tan(s * u)) / s),
        ):
            gg = np.exp(1j * phi) * prim(lo, hi)
            e_mat = -t * np.array([[gg.real, gg.imag], [gg.imag, -gg.real]])
            mats.append(np.linalg.inv(m * np.eye(2) + e_mat))
            ms.append(m)
            gs.append(gg)
        r = 1.01 * math.sqrt(2 * H * (max(ms) + t * max(abs(gs[0]), abs(gs[1]))))
        yy = -r + (np.arange(ngrid) + 0.5) * 2 * r / ngrid
        y1, y2 = np.meshgrid(yy, yy, indexing="ij")
        cell = (2 * r / ngrid) ** 2

        def qf(mi, y1=y1, y2=y2):
            return mi[0, 0] * y1 * y1 + 2 * mi[0, 1] * y1 * y2 + mi[1, 1] * y2 * y2

        ql, qh = qf(mats[0]), qf(mats[1])
        r2 = y1 * y1 + y2 * y2
        q0l, q0h = r2 / ms[0], r2 / ms[1]
        grid = {
            "lam": cell * np.sum(pos(H - ql / 2 - qh / 2)),
            "0": cell * np.sum(pos(H - q0l / 2 - q0h / 2)),
            "cap": cell
            * np.sum(pos(H - np.maximum(ql, q0l) / 2 - np.maximum(qh, q0h) / 2)),
        }
        a_lo = (ms[0] + t * (gs[0].real * c2 + gs[0].imag * s2)) / (
            ms[0] ** 2 - (t * abs(gs[0])) ** 2
        )
        a_hi = (ms[1] + t * (gs[1].real * c2 + gs[1].imag * s2)) / (
            ms[1] ** 2 - (t * abs(gs[1])) ** 2
        )
        b_lo, b_hi = 1 / ms[0], 1 / ms[1]
        formula = {
            "lam": H * H / 2 * np.sum(1 / (a_lo + a_hi)) * 2 * math.pi / nth,
            "0": H * H / 2 * 2 * math.pi / (b_lo + b_hi),
            "cap": H
            * H
            / 2
            * np.sum(1 / (np.maximum(a_lo, b_lo) + np.maximum(a_hi, b_hi)))
            * 2
            * math.pi
            / nth,
        }
        for k in grid:
            worst = max(worst, abs(grid[k] / formula[k] - 1))
    return worst


# ================================================================ run
QS = (0.25, 1.0, 4.0, 9.0, 9.8)
PHIS = [k * math.pi / 4 for k in range(8)]
LAMS = [2.0 ** (-k / 2) for k in range(21)]

print(
    "INFO constants: q, f_max, C_E, W/Vhat, C_1^0 = C_E W/Vhat (lam -> 0), "
    "t_1, C_1(q) (lam <= 1)"
)
for q in QS:
    c = consts(q)
    print(
        f"INFO q={q}: f_max={c['fmax']:.6g} C_E={c['ce']:.6g} "
        f"W/Vhat={c['w'] / c['vhat']:.6g} "
        f"C_1^0={c['ce'] * c['w'] / c['vhat']:.6g} t_1={c['t1']:.6f} "
        f"C_1={c1_of(c, c['t1']):.6g}",
        flush=True,
    )

# K2
for q, lam, phi in ((1.0, 0.25, 0.3), (4.0, 0.5, 1.1)):
    wst = polar_check(q, lam, phi, (-0.4, -0.2, 0.0, 0.25, 0.45))
    check(
        f"K2 polar reduction q={q} lam={lam} phi={phi}",
        wst < 1e-3,
        f"max rel diff {wst:.2e}",
    )

rows = []
k1_worst = dict(e_tm=0.0, e_ce=0.0, da=0.0, alow=0.0)
k3_ok = k4_ok = k5_ok = True
k3_worst = k4_worst = 0.0
k5_fail = []
for q in QS:
    cst = consts(q)
    for lam in LAMS:
        for phi in PHIS:
            r = one_case(q, lam, phi, pointwise=True)
            for k in k1_worst:
                k1_worst[k] = max(k1_worst[k], r["worst"][k])
            dvq = r["v_lam"] - r["v0"]
            tol = 1e-6 * abs(dvq) + 1e-12 * r["v0"]
            k3_worst = max(
                k3_worst, abs(dvq - r["dv_closed"]) / (abs(dvq) + 1e-12 * r["v0"])
            )
            if abs(dvq - r["dv_closed"]) > tol:
                k3_ok = False
            k4e = max(
                abs(r["v0_num"] / r["v0"] - 1),
                abs(r["w_num"] / cst["w"] - 1),
                abs(r["vhat_num"] / cst["vhat"] - 1),
            )
            k4_worst = max(k4_worst, k4e)
            if k4e > 1e-9:
                k4_ok = False
            c1t = c1_of(cst, r["t"]) * lam * lam
            c1u = c1_of(cst, cst["t1"]) * lam * lam
            ok = (
                r["tv"] <= r["B"] <= c1t * (1 + 1e-12) <= c1u * (1 + 1e-12)
                and r["d_lam"] <= r["b_lam"]
                and r["d_0"] <= r["b_0"]
            )
            if not ok:
                k5_ok = False
                k5_fail.append((q, lam, phi, r["tv"], r["B"], c1t))
            rows.append((q, lam, phi, r["tv"], r["B"], c1t, r["d_lam"], r["d_0"]))
        sub = [x for x in rows if x[0] == q and x[1] == lam]
        tvs = [x[3] / lam**2 for x in sub]
        slack = max(x[3] / x[5] for x in sub)
        print(
            f"INFO q={q} lam=2^-{round(-2 * math.log2(lam))}/2: "
            f"TV/lam^2 min {min(tvs):.5g} "
            f"max {max(tvs):.5g}; max TV/(C_1(q,t) lam^2) {slack:.3g}; nu={r['nu']}",
            flush=True,
        )

check(
    "K1 e <= t m at every node",
    k1_worst["e_tm"] <= 1 + 1e-9,
    f"max e/(t m) = {k1_worst['e_tm']:.6f}",
)
check(
    "K1 e <= t lam (4 f_max - 2) at every node",
    k1_worst["e_ce"] <= 1 + 1e-9,
    f"max ratio = {k1_worst['e_ce']:.4f}",
)
check(
    "K1 |dA| <= e/(m(m-e)) at every node",
    k1_worst["da"] <= 1 + 1e-9,
    f"max ratio = {k1_worst['da']:.6f}",
)
check(
    "K1 A >= 1/(m+e) at every node",
    k1_worst["alow"] <= 1 + 1e-9,
    f"max ratio = {k1_worst['alow']:.12f}",
)
check(
    "K3 theta grid vs closed form for V_lam - V_0",
    k3_ok,
    f"max rel diff {k3_worst:.2e}",
)
check("K4 V_0, W, Vhat vs closed forms", k4_ok, f"max rel diff {k4_worst:.2e}")
check(
    f"K5 TV <= B <= C_1(q,t) lam^2 <= C_1(q) lam^2 in all {len(rows)} cases",
    k5_ok,
    f"failures {k5_fail[:3]}",
)

# K6
k6_worst = 0.0
for q in QS:
    for lam in (2.0**-3, 2.0**-8):
        a = one_case(q, lam, 0.3)["tv"]
        b = one_case(q, lam, 0.3, factor=2, nth=512)["tv"]
        k6_worst = max(k6_worst, abs(b / a - 1))
check(
    "K6 doubling panels and theta grid changes TV < 1e-3",
    k6_worst < 1e-3,
    f"max rel change {k6_worst:.2e}",
)

# K7
ref = {
    (0.25, 1 / 8): 1.644,
    (0.25, 1 / 64): 1.623,
    (1.0, 1 / 8): 3.440,
    (1.0, 1 / 64): 3.283,
    (4.0, 1 / 8): 7.08,
    (4.0, 1 / 64): 8.05,
}
k7_worst = 0.0
for (q, lam), want in ref.items():
    got = one_case(q, lam, 0.3)["tv"] / lam**2
    print(
        f"INFO K7 q={q} lam={lam}: TV/lam^2 = {got:.5f} (referee R8: {want})",
        flush=True,
    )
    k7_worst = max(k7_worst, abs(got / want - 1))
check(
    "K7 agreement with the referee's R8 values (1 percent)",
    k7_worst < 1e-2,
    f"max rel diff {k7_worst:.2e}",
)

print(f"SUMMARY pass={npass} failed={failed}")
raise SystemExit(1 if failed else 0)
