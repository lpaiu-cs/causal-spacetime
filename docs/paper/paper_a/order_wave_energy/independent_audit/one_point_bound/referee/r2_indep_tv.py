# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of one_point_strip_bound.md.
"""Referee r2: independent deterministic computation of TV(mu_lam, mu_0).

Different code path from checks.py / ref_checks.py R8:
  * m and G from cumulative GL-16 on the intervals BETWEEN consecutive u-nodes
    (no closed form for m, no panel primitive), spot-checked against scipy QAWO;
  * the 2x2 inverse by adjugate/determinant (no Step-1 formula);
  * V_lam(u) by the ellipse formula pi h^2 / sqrt(det(K_- + K_+)) (no theta grid);
  * V_cap(u) by piecewise GL-20 in theta with the kinks located from the
    matrix entries (no midpoint grid);
  * TV = (1/2) int |p_lam - p_0| assembled from V_cap, V_lam, V_0 (not formula (4)).
Also evaluates the note's bound (1+t)/(1-t) J/Vhat with the referee's own e, m.
"""
import math
import sys
import numpy as np
from scipy.integrate import quad

GX10, GW10 = np.polynomial.legendre.leggauss(10)
GX16, GW16 = np.polynomial.legendre.leggauss(16)
GX20, GW20 = np.polynomial.legendre.leggauss(20)
H = 1.0


def gl_on(a, b, func, gx, gw):
    half = (b - a) / 2
    mid = (a + b) / 2
    x = mid[:, None] + half[:, None] * gx[None, :]
    return np.sum(half[:, None] * gw[None, :] * func(x), axis=1)


def case(q, lam, phi, spot=False):
    s = math.sqrt(q)
    d = (math.pi / 2 - s / 2) / s
    width = min(lam / 8, d / 3, 1 / 40)
    npan = math.ceil(1 / width)
    ed = np.linspace(-0.5, 0.5, npan + 1)
    half = (ed[1:] - ed[:-1]) / 2
    mid = (ed[1:] + ed[:-1]) / 2
    u = (mid[:, None] + half[:, None] * GX10[None, :]).ravel()
    wu = (half[:, None] * GW10[None, :]).ravel()
    pts = np.concatenate([[-0.5], u, [0.5]])
    lo, hi = pts[:-1], pts[1:]
    f = lambda x: 1.0 / np.cos(s * x) ** 2
    g = lambda x: np.exp(1j * x / lam) / np.cos(s * x) ** 2
    inc_f = gl_on(lo, hi, f, GX16, GW16)
    inc_g = gl_on(lo, hi, g, GX16, GW16)
    m_lo = np.cumsum(inc_f)[:-1]
    m_hi = np.cumsum(inc_f[::-1])[::-1][1:]
    G_lo = np.exp(1j * phi) * np.cumsum(inc_g)[:-1]
    G_hi = np.exp(1j * phi) * np.cumsum(inc_g[::-1])[::-1][1:]
    M = m_lo[-1] + m_hi[-1]
    if spot:
        idx = np.linspace(0, u.size - 1, 7).astype(int)
        worst = 0.0
        for k in idx:
            re = quad(lambda x: 1 / math.cos(s * x) ** 2, -0.5, u[k], weight="cos", wvar=1 / lam, limit=2000, epsabs=1e-15)[0]
            im = quad(lambda x: 1 / math.cos(s * x) ** 2, -0.5, u[k], weight="sin", wvar=1 / lam, limit=2000, epsabs=1e-15)[0]
            ref = np.exp(1j * phi) * (re + 1j * im)
            worst = max(worst, abs(ref - G_lo[k]) / max(abs(ref), 1e-300))
        print(f"  spot QAWO vs cumulative GL16 for G_-: max rel diff {worst:.2e}")
    t = 2 * lam * s / math.sqrt(1 + 4 * q * lam * lam)
    a2 = np.cos(s * u) ** 2

    def kmat(m, G):
        e11 = m - t * G.real
        e22 = m + t * G.real
        e12 = -t * G.imag
        det = e11 * e22 - e12 * e12
        return e22 / det, -e12 / det, e11 / det  # K11, K12, K22

    K_lo = kmat(m_lo, G_lo)
    K_hi = kmat(m_hi, G_hi)
    # V_lam(u): ellipse formula
    S11 = K_lo[0] + K_hi[0]
    S12 = K_lo[1] + K_hi[1]
    S22 = K_lo[2] + K_hi[2]
    vlam_u = math.pi * H * H / np.sqrt(S11 * S22 - S12 * S12)
    B_lo, B_hi = 1 / m_lo, 1 / m_hi
    v0_u = math.pi * H * H / (B_lo + B_hi)

    # V_cap(u): theta in [0, pi), kinks from matrix entries
    def coeffs(K):
        al = (K[0] + K[2]) / 2
        be = (K[0] - K[2]) / 2
        ga = K[1]
        return al, be, ga

    cl = coeffs(K_lo)
    ch = coeffs(K_hi)
    brk = [np.zeros_like(u), np.full_like(u, math.pi)]
    for (al, be, ga), B in ((cl, B_lo), (ch, B_hi)):
        R = np.hypot(be, ga)
        psi = np.arctan2(ga, be)
        with np.errstate(invalid="ignore", divide="ignore"):
            c = (B - al) / R
        ok = np.abs(c) < 1
        ac = np.arccos(np.clip(c, -1, 1))
        for sg in (1, -1):
            th = np.mod((psi + sg * ac) / 2, math.pi)
            brk.append(np.where(ok, th, 0.0))
    brk = np.sort(np.stack(brk, axis=1), axis=1)  # nu x 6
    a_ = brk[:, :-1]
    b_ = brk[:, 1:]
    hh = (b_ - a_) / 2
    mm = (b_ + a_) / 2
    th = mm[:, :, None] + hh[:, :, None] * GX20[None, None, :]
    wth = hh[:, :, None] * GW20[None, None, :]
    c2, s2 = np.cos(2 * th), np.sin(2 * th)

    def A(cf):
        al, be, ga = cf
        return al[:, None, None] + be[:, None, None] * c2 + ga[:, None, None] * s2

    A_lo, A_hi = A(cl), A(ch)
    Scap = np.maximum(A_lo, B_lo[:, None, None]) + np.maximum(A_hi, B_hi[:, None, None])
    vcap_u = H * H * np.sum(wth / Scap, axis=(1, 2))  # (h^2/2) * 2 * int_0^pi
    Vl = float(np.sum(wu * a2 * vlam_u))
    V0 = float(np.sum(wu * a2 * v0_u))
    Vc = float(np.sum(wu * a2 * vcap_u))
    tv = 0.5 * (Vc * abs(1 / Vl - 1 / V0) + (Vl - Vc) / Vl + (V0 - Vc) / V0)
    # note's bound with referee's e and m
    e_lo = t * np.abs(G_lo)
    e_hi = t * np.abs(G_hi)
    J = float(np.sum(wu * a2 * (e_lo * m_hi**2 + e_hi * m_lo**2))) / M**2
    vhat = V0 / (math.pi * H * H)
    bound = (1 + t) / (1 - t) * J / vhat
    fmax = 1 / math.cos(s / 2) ** 2
    CE = 2 * s * (4 * fmax - 2)
    W = float(np.sum(wu * a2 * (m_lo**2 + m_hi**2))) / M**2
    c1t = (1 + t) / (1 - t) * CE * W / vhat
    return dict(tv=tv, bound=bound, c1t=c1t, Vl=Vl, V0=V0, Vc=Vc, nu=u.size,
                emax_over_tm=float(max(np.max(e_lo / (t * m_lo)), np.max(e_hi / (t * m_hi)))))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "grid"
    if mode == "grid":
        note = {  # TV/lam^2 printed in checks_output.txt (5 s.f.)
            (0.25, 0): 0.39828, (0.25, 4): 1.3391, (0.25, 8): 1.6372, (0.25, 12): 1.623, (0.25, 16): 1.6201,
            (1.0, 0): 0.63365, (1.0, 4): 2.9429, (1.0, 8): 3.358, (1.0, 12): 3.2828, (1.0, 16): 3.2705,
            (4.0, 0): 0.8113, (4.0, 4): 5.9759, (4.0, 8): 7.7029, (4.0, 12): 8.047, (4.0, 16): 8.0638,
            (9.0, 0): 0.91691, (9.0, 4): 9.5699, (9.0, 8): 27.822, (9.0, 12): 53.007, (9.0, 16): 75.906,
            (9.8, 0): 0.94663, (9.8, 4): 10.983, (9.8, 8): 43.248, (9.8, 12): 120.06, (9.8, 16): 354.04,
        }
        worst = 0.0
        first = True
        for (q, k), want in note.items():
            lam = 2.0 ** (-k / 2)
            for phi in (0.0, 1.234):
                r = case(q, lam, phi, spot=first)
                first = False
                got = r["tv"] / lam**2
                rel = abs(got / want - 1)
                worst = max(worst, rel)
                print(f"q={q} lam=2^-{k}/2 phi={phi}: TV/lam^2 referee {got:.6f} note {want} rel {rel:.1e} | "
                      f"TV/bound {r['tv']/r['bound']:.4f} TV/(C1t lam^2) {r['tv']/(r['c1t']*lam*lam):.4g} "
                      f"(V_lam-V_0)/V_0 {(r['Vl']-r['V0'])/r['V0']:.3e} nu={r['nu']}", flush=True)
        print(f"max rel diff vs note (5 s.f. rounding ~ 1e-5): {worst:.2e}")
    elif mode == "random":
        rng = np.random.default_rng(20261126)
        nbad = 0
        worst = 0.0
        for i in range(24):
            q = float(rng.uniform(0.05, 9.8))
            lam = float(2.0 ** rng.uniform(-8, 0))
            phi = float(rng.uniform(0, 2 * math.pi))
            r = case(q, lam, phi)
            ok = r["tv"] <= r["bound"] <= r["c1t"] * lam * lam * (1 + 1e-12)
            nbad += not ok
            worst = max(worst, r["tv"] / (r["c1t"] * lam * lam))
            print(f"q={q:.4f} lam={lam:.5f} phi={phi:.3f}: TV/lam^2={r['tv']/lam**2:.5g} TV/bound={r['tv']/r['bound']:.4f} "
                  f"TV/(C1t lam^2)={r['tv']/(r['c1t']*lam*lam):.4g} max e/(t m)={r['emax_over_tm']:.9f} ok={ok}", flush=True)
        print(f"random cases violating Prop 1 chain: {nbad}; max TV/(C1(q,t) lam^2) = {worst:.4f}")
    elif mode == "small":
        for q in (0.25, 1.0, 4.0):
            for k in (18, 20):
                lam = 2.0 ** (-k / 2)
                r = case(q, lam, 0.7)
                print(f"q={q} lam=2^-{k}/2: TV/lam^2 {r['tv']/lam**2:.6f} (V_lam-V_0)/V_0 {(r['Vl']-r['V0'])/r['V0']:.3e} nu={r['nu']}", flush=True)
