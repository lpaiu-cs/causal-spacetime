# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of one_point_strip_bound.md.
"""Referee r3: convergence study at q = 9.8 (and q = 1 control), where the
referee's r2 value and the note's printed value differ by ~5e-5 relative.

Separates the u-panel refinement from the theta rule:
  * theta: kink-aware piecewise GL-n (n = 20, 40) versus the note's 256 / 512 / 2048
    midpoints on a half period (the note's rule, re-implemented here);
  * u: panel width divisor 1, 2, 4 (relative to r2's choice).
"""
import math
import numpy as np

GX10, GW10 = np.polynomial.legendre.leggauss(10)
GX16, GW16 = np.polynomial.legendre.leggauss(16)
H = 1.0


def setup(q, lam, phi, div):
    s = math.sqrt(q)
    d = (math.pi / 2 - s / 2) / s
    width = min(lam / 8, d / 3, 1 / 40) / div
    npan = math.ceil(1 / width)
    ed = np.linspace(-0.5, 0.5, npan + 1)
    half = (ed[1:] - ed[:-1]) / 2
    mid = (ed[1:] + ed[:-1]) / 2
    u = (mid[:, None] + half[:, None] * GX10[None, :]).ravel()
    wu = (half[:, None] * GW10[None, :]).ravel()
    pts = np.concatenate([[-0.5], u, [0.5]])
    lo, hi = pts[:-1], pts[1:]
    hl, md = (hi - lo) / 2, (hi + lo) / 2
    x = md[:, None] + hl[:, None] * GX16[None, :]
    inc_f = np.sum(hl[:, None] * GW16 / np.cos(s * x) ** 2, axis=1)
    inc_g = np.sum(hl[:, None] * GW16 * np.exp(1j * x / lam) / np.cos(s * x) ** 2, axis=1)
    m_lo = np.cumsum(inc_f)[:-1]
    m_hi = np.cumsum(inc_f[::-1])[::-1][1:]
    G_lo = np.exp(1j * phi) * np.cumsum(inc_g)[:-1]
    G_hi = np.exp(1j * phi) * np.cumsum(inc_g[::-1])[::-1][1:]
    t = 2 * lam * s / math.sqrt(1 + 4 * q * lam * lam)
    a2 = np.cos(s * u) ** 2

    def kmat(m, G):
        e11 = m - t * G.real
        e22 = m + t * G.real
        e12 = -t * G.imag
        det = e11 * e22 - e12 * e12
        K11, K12, K22 = e22 / det, -e12 / det, e11 / det
        return (K11 + K22) / 2, (K11 - K22) / 2, K12

    return u, wu, a2, kmat(m_lo, G_lo), kmat(m_hi, G_hi), 1 / m_lo, 1 / m_hi


def deltas_gl(u, wu, a2, cl, ch, Blo, Bhi, n, chunk=4000):
    gx, gw = np.polynomial.legendre.leggauss(n)
    dl = d0 = 0.0
    for i0 in range(0, u.size, chunk):
        sl = slice(i0, i0 + chunk)
        brk = [np.zeros(u[sl].size), np.full(u[sl].size, math.pi)]
        for (al, be, ga), B in ((cl, Blo), (ch, Bhi)):
            al, be, ga, B = al[sl], be[sl], ga[sl], B[sl]
            R = np.hypot(be, ga)
            psi = np.arctan2(ga, be)
            c = (B - al) / R
            ok = np.abs(c) < 1
            ac = np.arccos(np.clip(c, -1, 1))
            for sg in (1, -1):
                brk.append(np.where(ok, np.mod((psi + sg * ac) / 2, math.pi), 0.0))
        brk = np.sort(np.stack(brk, axis=1), axis=1)
        a_, b_ = brk[:, :-1], brk[:, 1:]
        hh, mm = (b_ - a_) / 2, (b_ + a_) / 2
        th = mm[:, :, None] + hh[:, :, None] * gx[None, None, :]
        wth = hh[:, :, None] * gw[None, None, :]
        dl_, d0_ = integrand(th, cl, ch, Blo, Bhi, sl)
        w = (wu[sl] * a2[sl])[:, None, None] * wth * H * H  # (h^2/2)*2 over [0,pi)
        dl += float(np.sum(w * dl_))
        d0 += float(np.sum(w * d0_))
    return dl, d0


def integrand(th, cl, ch, Blo, Bhi, sl):
    c2, s2 = np.cos(2 * th), np.sin(2 * th)
    sh = (slice(None),) + (None,) * (th.ndim - 1)
    Al = cl[0][sl][sh] + cl[1][sl][sh] * c2 + cl[2][sl][sh] * s2
    Ah = ch[0][sl][sh] + ch[1][sl][sh] * c2 + ch[2][sl][sh] * s2
    bl, bh = Blo[sl][sh], Bhi[sl][sh]
    Sl = Al + Ah
    S0 = bl + bh
    Sc = np.maximum(Al, bl) + np.maximum(Ah, bh)
    return (Sc - Sl) / (Sl * Sc), (Sc - S0) / (S0 * Sc)


def deltas_mid(u, wu, a2, cl, ch, Blo, Bhi, nth, chunk=4000):
    th1 = (np.arange(nth) + 0.5) * math.pi / nth
    dl = d0 = 0.0
    for i0 in range(0, u.size, chunk):
        sl = slice(i0, i0 + chunk)
        th = np.broadcast_to(th1[None, :], (u[sl].size, nth))
        dl_, d0_ = integrand(th, cl, ch, Blo, Bhi, sl)
        w = (wu[sl] * a2[sl])[:, None] * (math.pi / nth) * H * H
        dl += float(np.sum(w * dl_))
        d0 += float(np.sum(w * d0_))
    return dl, d0


def tv_from(dl, d0, u, wu, a2, Blo, Bhi):
    V0 = float(np.sum(wu * a2 * math.pi * H * H / (Blo + Bhi)))
    Vl = V0 + dl - d0
    return max(dl, d0) / max(Vl, V0)


for q, lam in ((9.8, 1 / 64), (9.8, 1 / 256), (9.0, 1 / 64), (1.0, 1 / 64), (4.0, 1 / 8)):
    print(f"--- q={q} lam={lam}")
    for div in (1, 2, 4):
        u, wu, a2, cl, ch, Blo, Bhi = setup(q, lam, 0.3, div)
        out = []
        for n in (20, 40):
            dl, d0 = deltas_gl(u, wu, a2, cl, ch, Blo, Bhi, n)
            out.append(f"GL{n}: {tv_from(dl, d0, u, wu, a2, Blo, Bhi)/lam**2:.7f}")
        for nth in (256, 512, 2048):
            dl, d0 = deltas_mid(u, wu, a2, cl, ch, Blo, Bhi, nth)
            out.append(f"mid{nth}: {tv_from(dl, d0, u, wu, a2, Blo, Bhi)/lam**2:.7f}")
        print(f"  u-div {div} (nu={u.size}): " + "  ".join(out), flush=True)
