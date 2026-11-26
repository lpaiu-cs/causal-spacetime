# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of gram_fourier_lemmas.md.
"""Referee check r1: Lemma E (canonical map) at the tips, with the stated constants.

Independent code. Deterministic. For q in {0.01, 0.25, 1}, several lambda <= lambda_1,
six phases, and tip distances s from 1e-10 to 0.5 on both sides, compute exactly
(up to quadrature) eps_L, eps_R, T = (W0 Q)^(-1/2), det T, ||T - I||, the sup of
|delta v| and |delta y| over the D_0 slice, and compare with the note's bounds:
  eps <= e(s), |det T - 1| <= 4 e^2, ||T - I|| <= 4 e,
  sup|delta v| <= 3 sig(1-sig)(epsL+epsR) <= 6 C_E lam^2/M <= 12 f_max C_E lam^2/M,
  sup|delta y| <= 4 e sqrt(2 f_max s).
Also: det T closed form, zeta preservation, Lemma E(iv) first order.
"""
import math
import numpy as np

GN, GW = np.polynomial.legendre.leggauss(24)


class Geo:
    def __init__(self, q):
        self.q = q
        self.k = math.sqrt(q)
        self.fmax = 1 / math.cos(self.k / 2) ** 2
        self.CE = 2 * self.k * (4 * self.fmax - 2)
        self.M = self.m(-0.5, 0.5)
        self.lam1 = 1 / (4 * self.CE)

    def f(self, w):
        return 1 / np.cos(self.k * w) ** 2

    def m(self, s, t):
        k = self.k
        return math.sin(k * (t - s)) / (k * math.cos(k * s) * math.cos(k * t))


class Wave:
    def __init__(self, g, lam, phi):
        self.g, self.lam, self.phi = g, lam, phi
        self.c = math.sqrt(1 + 4 * g.q * lam**2)
        self.pref = 2 * lam * g.k / self.c * complex(math.cos(phi), math.sin(phi))
        self.h = min(lam / 8, 1 / 128)
        n = math.ceil(1 / self.h)
        self.edges = np.linspace(-0.5, 0.5, n + 1)
        vals = np.array([self._gl(a, b) for a, b in zip(self.edges[:-1], self.edges[1:])])
        self.P = np.concatenate([[0], np.cumsum(vals)])

    def _gl(self, a, b):
        x = 0.5 * (b - a) * GN + 0.5 * (a + b)
        return 0.5 * (b - a) * np.sum(GW * self.g.f(x) * np.exp(1j * x / self.lam))

    def prim(self, u):
        j = min(int(np.searchsorted(self.edges, u, side="right")) - 1, len(self.edges) - 2)
        j = max(j, 0)
        return self.P[j] + self._gl(self.edges[j], u)

    def integral(self, s, t):
        if t - s <= self.h:
            return self._gl(s, t)
        return self.prim(t) - self.prim(s)

    def beta(self, s, t):
        return self.pref * self.integral(s, t)


def bm(b):
    return np.array([[b.real, b.imag], [b.imag, -b.real]])


def spow(A, p):
    w, V = np.linalg.eigh(A)
    return V @ np.diag(w**p) @ V.T


def frame(wv, u):
    g = wv.g
    m = g.m(-0.5, u)
    mr = g.m(u, 0.5)
    bL, bR = wv.beta(-0.5, u), wv.beta(u, 0.5)
    epsL, epsR = abs(bL) / m, abs(bR) / mr
    EL, ER = bm(bL) / m, bm(bR) / mr
    P = (np.eye(2) + EL) / (1 - epsL**2)
    S = (np.eye(2) + ER) / (1 - epsR**2)
    sig = m / g.M
    X = (1 - sig) * P + sig * S  # = W0 Q
    T = spow(X, -0.5)
    W0 = m * mr / g.M
    return dict(m=m, mr=mr, bL=bL, bR=bR, epsL=epsL, epsR=epsR, EL=EL, ER=ER,
                P=P, S=S, sig=sig, X=X, T=T, W0=W0)


def main():
    out = []
    worst = dict(eps=0, det=0, TI=0, dv1=0, dv2=0, dy=0, closed=0, zeta=0, dvid=0)
    for q in (0.01, 0.25, 1.0):
        g = Geo(q)
        lams = sorted({min(g.lam1, 0.5), min(g.lam1, 0.5) / 3, 1 / 64, 1 / 512})
        lams = [l for l in lams if l <= g.lam1]
        for lam in lams:
            for phi in (0.0, 0.37, 1.3, 2.1, 4.0, 5.5):
                wv = Wave(g, lam, phi)
                svals = np.unique(np.concatenate([np.logspace(-10, math.log10(0.5), 50),
                                                  lam * np.array([0.3, 1, 2, 2 * math.pi, 5])]))
                svals = svals[svals <= 0.5]
                for s in svals:
                    for side in (-1, 1):
                        u = -0.5 + s if side < 0 else 0.5 - s
                        if abs(u) >= 0.5:
                            continue
                        F = frame(wv, u)
                        e = g.CE * lam * min(1.0, lam / s)
                        worst["eps"] = max(worst["eps"], max(F["epsL"], F["epsR"]) / e)
                        T = F["T"]
                        dT = np.linalg.det(T)
                        a0 = (1 - F["sig"]) / (1 - F["epsL"] ** 2) + F["sig"] / (1 - F["epsR"] ** 2)
                        Ec = (1 - F["sig"]) * F["EL"] / (1 - F["epsL"] ** 2) + F["sig"] * F["ER"] / (1 - F["epsR"] ** 2)
                        En = np.linalg.norm(Ec, 2)
                        closed = (a0**2 - En**2) ** -0.5
                        worst["closed"] = max(worst["closed"], abs(dT / closed - 1))
                        worst["det"] = max(worst["det"], abs(dT - 1) / (4 * e**2))
                        worst["TI"] = max(worst["TI"], np.linalg.norm(T - np.eye(2), 2) / (4 * e))
                        # delta v: stable form A = T M_L^-1 T - I/m = -(1/M) X^-1/2 (S-P) X^-1/2
                        A = -(1 / g.M) * T @ (F["S"] - F["P"]) @ T
                        if s > 1e-3:
                            ML = F["m"] * np.eye(2) - bm(F["bL"])
                            A2 = T @ np.linalg.inv(ML) @ T - np.eye(2) / F["m"]
                            worst["dvid"] = max(worst["dvid"], np.linalg.norm(A - A2, 2) * F["W0"])
                        supdv = F["W0"] * np.linalg.norm(A, 2)  # sup over |z|^2 <= 2 W0 of |z^T A z|/2
                        b1 = 3 * F["sig"] * (1 - F["sig"]) * (F["epsL"] + F["epsR"])
                        worst["dv1"] = max(worst["dv1"], supdv / b1 if b1 > 0 else 0)
                        worst["dv2"] = max(worst["dv2"], supdv / (12 * g.fmax * g.CE * lam**2 / g.M))
                        supdy = np.linalg.norm(T - np.eye(2), 2) * math.sqrt(2 * F["W0"])
                        worst["dy"] = max(worst["dy"], supdy / (4 * e * math.sqrt(2 * g.fmax * s)))
                        # zeta preservation at a canonical point inside the slice
                        if s > 1e-6:
                            ang = 0.7 + s
                            r = 0.6 * math.sqrt(2 * F["W0"]) * 0.5
                            z = r * np.array([math.cos(ang), math.sin(ang)])
                            vt = 0.3 + z @ z / (2 * F["m"])
                            y = T @ z
                            ML = F["m"] * np.eye(2) - bm(F["bL"])
                            MR = F["mr"] * np.eye(2) - bm(F["bR"])
                            eta_l = vt - z @ z / (2 * F["m"])
                            v = eta_l + y @ np.linalg.solve(ML, y) / 2
                            zl = 1 - v - y @ np.linalg.solve(MR, y) / 2
                            z0 = 1 - vt - z @ z / (2 * F["mr"])
                            worst["zeta"] = max(worst["zeta"], abs(zl - z0))
        out.append(f"q={q}: lambda_1={g.lam1:.4f}, lams={['%.4g' % l for l in lams]}")
    out.append("max ratios (all must be <= 1 for the bound to hold):")
    out.append(f"  max(epsL,epsR)/e(s)                  = {worst['eps']:.4f}")
    out.append(f"  |det T - 1| / (4 e^2)                = {worst['det']:.4f}")
    out.append(f"  ||T - I|| / (4 e)                    = {worst['TI']:.4f}")
    out.append(f"  sup|dv| / (3 sig(1-sig)(eL+eR))      = {worst['dv1']:.4f}")
    out.append(f"  sup|dv| / (12 fmax C_E lam^2 / M)    = {worst['dv2']:.4f}")
    out.append(f"  sup|dy| / (4 e sqrt(2 fmax s))       = {worst['dy']:.4f}")
    out.append(f"identities: det T closed form rel err = {worst['closed']:.2e};"
               f" zeta preserved abs err = {worst['zeta']:.2e};"
               f" delta-v stable vs direct (x W0) = {worst['dvid']:.2e}")

    # Lemma E(iv): interior first order, s >= 0.1
    res_T, res_v = [], []
    g = Geo(1.0)
    for lam in (1 / 32, 1 / 64, 1 / 128, 1 / 256):
        wT = wv_ = 0
        for phi in (0.37, 2.1):
            wv = Wave(g, lam, phi)
            for u in np.linspace(-0.4, 0.4, 17):
                F = frame(wv, u)
                lin = -0.5 * ((1 - F["sig"]) * F["EL"] + F["sig"] * F["ER"])
                wT = max(wT, np.linalg.norm(F["T"] - np.eye(2) - lin, 2) / lam**4)
                A = -(1 / g.M) * F["T"] @ (F["S"] - F["P"]) @ F["T"]
                Alin = (F["EL"] - F["ER"]) / g.M
                wv_ = max(wv_, np.linalg.norm(A - Alin, 2) * F["W0"] / lam**4)
        res_T.append(wT)
        res_v.append(wv_)
    out.append("Lemma E(iv), q=1, |u|<=0.4: max ||T-I-lin||/lam^4 = "
               + ", ".join(f"{x:.3f}" for x in res_T)
               + "; max sup|dv - first order|/lam^4 = " + ", ".join(f"{x:.3f}" for x in res_v)
               + "  (lam = 1/32..1/256)")
    txt = "\n".join(out)
    print(txt)


if __name__ == "__main__":
    main()
