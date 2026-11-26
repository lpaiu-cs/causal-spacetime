# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of gram_fourier_lemmas.md.
"""Referee check r2: Lemma F (7), Lemma G density bound, Corollary K2(ii), Remark K1.

Independent code. Deterministic quadrature only.
 (a) Exact TV between normalised densities a^2 det T / Z_lam and a^2 / Z_0 on D_0
     (both u-only times slice volume pi W0) against (7); ratio TV/(lam^4 log(1/lam)).
 (b) sup of the nu_lam density a^2 det T / Z_lam against 2/Z_0 (Lemma G).
 (c) K2(ii): |Lambda_lam/c - Lambda_0| against 3 fmax C_E^2 lam^2 Phi(l), including
     SHORT intervals l < lam (the note claims it 'for any interval of length l'),
     and the growth of dev/lam^4 with log(1/lam) (is the log genuine?).
 (d) Remark K1: Lambda_0 = alpha; Lambda_20/(3 Lambda_0^2) = C1 C3/(3 C2^2).
"""
import math
import numpy as np
import mpmath as mp

GN, GW = np.polynomial.legendre.leggauss(24)


class Geo:
    def __init__(self, q):
        self.q = q
        self.k = math.sqrt(q)
        self.fmax = 1 / math.cos(self.k / 2) ** 2
        self.CE = 2 * self.k * (4 * self.fmax - 2)
        self.M = self.m(-0.5, 0.5)
        self.lam1 = 1 / (4 * self.CE)
        self.gamma = self.alpha(1.0)

    def a2(self, w):
        return np.cos(self.k * w) ** 2

    def f(self, w):
        return 1 / np.cos(self.k * w) ** 2

    def m(self, s, t):
        k = self.k
        return np.sin(k * (t - s)) / (k * np.cos(k * s) * np.cos(k * t))

    def alpha(self, ell):
        x = self.k * ell
        return (1 - x / math.tan(x)) / (2 * self.q)


def grid(s, t, h):
    n = max(1, math.ceil((t - s) / h))
    e = np.linspace(s, t, n + 1)
    lo, hi = e[:-1, None], e[1:, None]
    x = (0.5 * (hi - lo) * GN + 0.5 * (hi + lo)).ravel()
    w = (0.5 * (hi - lo) * GW).ravel()
    return x, w


class Wave:
    def __init__(self, g, lam, phi):
        self.g, self.lam = g, lam
        self.c = math.sqrt(1 + 4 * g.q * lam**2)
        self.pref = 2 * lam * g.k / self.c * np.exp(1j * phi)
        self.h = min(lam / 8, 1 / 256)

    def I(self, s, t):
        """int_s^t f e^{iw/lam} dw, vectorised over t (s scalar)."""
        t = np.atleast_1d(t)
        out = np.zeros(t.shape, complex)
        for i, ti in enumerate(t):
            x, w = grid(s, ti, self.h)
            out[i] = np.sum(w * self.g.f(x) * np.exp(1j * x / self.lam))
        return out


def lam_ratio(wv, s, t):
    """Lambda_lam(s,t)/c with the exact integrand, and Lambda_0."""
    g = wv.g
    x, w = grid(s, t, wv.h)
    # primitive at the nodes via cumulative panels
    xs = np.sort(x)
    Ix = np.array([0j] * len(x))
    # cumulative: integrate between consecutive sorted nodes
    pts = np.concatenate([[s], xs])
    seg = np.array([wv.I(pts[i], pts[i + 1])[0] for i in range(len(xs))])
    Ix = np.cumsum(seg)
    Itot = Ix[-1] + wv.I(xs[-1], t)[0]
    b1 = np.abs(wv.pref * Ix) ** 2
    b2 = np.abs(wv.pref * (Itot - Ix)) ** 2
    b = abs(wv.pref * Itot) ** 2
    m1, m2, mm = g.m(s, xs), g.m(xs, t), g.m(s, t)
    integ_l = g.a2(xs) * np.sqrt((m1**2 - b1) * (m2**2 - b2) / (mm**2 - b))
    integ_0 = g.a2(xs) * m1 * m2 / mm
    order = np.argsort(x)
    ww = w[order]
    return float(np.sum(ww * integ_l)), float(np.sum(ww * integ_0))


def Phi(l, lam):
    return l * l / 2 if l <= lam else lam**2 / 2 + lam**2 * math.log(l / lam)


def detT_u(wv, u_nodes):
    g = wv.g
    out = []
    # cumulative primitive along sorted nodes
    xs = np.sort(u_nodes)
    pts = np.concatenate([[-0.5], xs, [0.5]])
    seg = np.array([wv.I(pts[i], pts[i + 1])[0] for i in range(len(pts) - 1)])
    cum = np.concatenate([[0], np.cumsum(seg)])
    Itot = cum[-1]
    IL = cum[1:-1]
    bL = wv.pref * IL
    bR = wv.pref * (Itot - IL)
    m = g.m(-0.5, xs)
    mr = g.m(xs, 0.5)
    eL, eR = np.abs(bL) / m, np.abs(bR) / mr
    sig = m / g.M
    a0 = (1 - sig) / (1 - eL**2) + sig / (1 - eR**2)
    Ec = (1 - sig) * (bL / m) / (1 - eL**2) + sig * (bR / mr) / (1 - eR**2)
    dT = (a0**2 - np.abs(Ec) ** 2) ** -0.5
    return xs, dT, m, mr


def main():
    out = []
    # (a), (b)
    for q in (0.25, 1.0):
        g = Geo(q)
        Z0 = math.pi * g.gamma
        for lam in (1 / 32, 1 / 64, 1 / 128, 1 / 256, 1 / 512):
            if lam > min(g.lam1, 0.5):
                continue
            for phi in (0.37, 2.1):
                wv = Wave(g, lam, phi)
                # nodes: graded near both tips
                xa, wa = grid(-0.5, -0.5 + 4 * lam, lam / 64)
                xb, wb = grid(-0.5 + 4 * lam, 0.5 - 4 * lam, lam / 8)
                xc, wc = grid(0.5 - 4 * lam, 0.5, lam / 64)
                x = np.concatenate([xa, xb, xc])
                w = np.concatenate([wa, wb, wc])
                o = np.argsort(x)
                x, w = x[o], w[o]
                xs, dT, m, mr = detT_u(wv, x)
                W0 = m * mr / g.M
                slice_ = math.pi * W0 * g.a2(xs)
                Zl = np.sum(w * slice_ * dT)
                Z0n = np.sum(w * slice_)
                tv = 0.5 * np.sum(w * slice_ * np.abs(dT / Zl - 1 / Z0n))
                bound = 4 * g.fmax * g.CE**2 / g.gamma * lam**4 * (1 + 2 * math.log(1 / (2 * lam)))
                dens = np.max(g.a2(xs) * dT / Zl) * Z0
                out.append(f"(a,b) q={q} lam=1/{round(1/lam)} phi={phi}: TV={tv:.3e} bound={bound:.3e}"
                           f" ok={tv <= bound}  TV/(lam^4 log(1/lam))={tv/(lam**4*math.log(1/lam)):.2f}"
                           f"  Z0 check {Z0n/Z0-1:.1e}  sup density*Z0={dens:.4f} (<=2)")
    # (c)
    for q in (0.25, 1.0):
        g = Geo(q)
        for phi in (0.37, 2.1):
            for lam in (1 / 32, 1 / 64, 1 / 128, 1 / 256):
                if lam > g.lam1:
                    continue
                wv = Wave(g, lam, phi)
                rows = []
                for (s, t) in ((-0.5, -0.5 + lam / 3), (-0.2, -0.2 + lam), (0.1, 0.1 + 3 * lam),
                               (-0.5, 0.2), (-0.25, 0.2), (-0.5, 0.5)):
                    Ll, L0 = lam_ratio(wv, s, t)
                    dev = abs(Ll - L0)
                    bd = 3 * g.fmax * g.CE**2 * lam**2 * Phi(t - s, lam)
                    rows.append(f"[{s:.3f},{t:.3f}] dev/bound={dev/bd:.3f} dev/lam^4={dev/lam**4:.3f}")
                out.append(f"(c) q={q} phi={phi} lam=1/{round(1/lam)}: " + "; ".join(rows))
    # (d)
    for q in (0.01, 0.25, 1.0):
        g = Geo(q)
        worst_a, worst_p = 0, 0
        for ell in (0.1, 0.3, 0.6, 1.0):
            s, t = -0.5, -0.5 + ell
            x, w = grid(s, t, 1 / 64)
            kern = g.a2(x) * g.m(s, x) * g.m(x, t) / g.m(s, t)
            L0 = np.sum(w * kern)
            worst_a = max(worst_a, abs(L0 / g.alpha(ell) - 1))
            L20 = np.sum(w * np.array([g.alpha(t - xi) for xi in x]) * kern)
            pf = L20 / (3 * g.alpha(ell) ** 2)
            z = mp.mpf(q) * ell**2
            C = lambda r: mp.hyp0f1(r + mp.mpf(1) / 2, -z / 4) / mp.factorial(2 * r - 1)
            p0 = float(C(1) * C(3) / (3 * C(2) ** 2))
            worst_p = max(worst_p, abs(pf - p0))
        out.append(f"(d) q={q}: max |Lambda_0/alpha - 1| = {worst_a:.1e}; max |pair fraction - p0| = {worst_p:.1e}")
    print("\n".join(out))


if __name__ == "__main__":
    main()
