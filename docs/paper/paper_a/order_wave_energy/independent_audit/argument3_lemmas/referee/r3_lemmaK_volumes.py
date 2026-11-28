# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of gram_fourier_lemmas.md.
"""Referee check r3: Lemma K (exact interval volumes) for the wave, by direct quadrature.

(i)  vol(I(P1,P2)) by a 3-d quadrature in (u, y1, y2) of the v-length between the two
     slacks of Theorem 2 (no use of the closed-form h), times the y-chart density c a^2,
     against pi kappa^2 Lambda_lam(u1,u2).
(ii) ordered-pair measure: int_{R in I} vol(I(R,P2)) dR, with the inner v_R integral of
     kappa(R,P2)^2 done in closed form over the slice found from the two slacks, the
     (y1,y2) integral by a grid, and Lambda_lam(u_R, u2) by quadrature; against
     (pi^2 kappa^4 / 6) Lambda_2lam(u1,u2).
Wave in the y chart: M_lam(s,t) = m I - B(s,t), density c a^2.  q = 0.25, lam = 1/16.
"""
import math
import numpy as np

GN, GW = np.polynomial.legendre.leggauss(20)
q, lam, phi = 0.25, 1 / 16, 0.9
k = math.sqrt(q)
c = math.sqrt(1 + 4 * q * lam**2)
pref = 2 * lam * k / c * np.exp(1j * phi)


def f(w):
    return 1 / np.cos(k * w) ** 2


def m(s, t):
    return np.sin(k * (t - s)) / (k * np.cos(k * s) * np.cos(k * t))


def gl(a, b, n=1):
    e = np.linspace(a, b, n + 1)
    lo, hi = e[:-1, None], e[1:, None]
    return (0.5 * (hi - lo) * GN + 0.5 * (hi + lo)).ravel(), (0.5 * (hi - lo) * GW).ravel()


EDGES = np.linspace(-0.5, 0.5, int(math.ceil(1 / (lam / 8))) + 1)


def _panel_int(a, b):
    a = np.atleast_1d(a); b = np.atleast_1d(b)
    x = a[:, None] + (b - a)[:, None] * (GN[None, :] + 1) / 2
    return ((b - a)[:, None] * GW[None, :] / 2 * f(x) * np.exp(1j * x / lam)).sum(axis=1)


PTAB = np.concatenate([[0], np.cumsum(_panel_int(EDGES[:-1], EDGES[1:]))])


def prim(u):
    u = np.atleast_1d(np.asarray(u, float))
    j = np.clip(np.searchsorted(EDGES, u, side="right") - 1, 0, len(EDGES) - 2)
    return PTAB[j] + _panel_int(EDGES[j], u)


def Iosc(s, t):
    return complex(prim(t)[0] - prim(s)[0])


def detM(s, t):
    """det M_lam(s,t), vectorised in s or t."""
    b = pref * (prim(t) - prim(s))
    return m(s, t) ** 2 - np.abs(b) ** 2


def Lam(s, t, nn=None):
    n = nn or max(2, math.ceil((t - s) / (lam / 6)))
    x, w = gl(s, t, n)
    d = np.sqrt(detM(s, x) * detM(x, t) / detM(s, t))
    return float(np.sum(w * c * np.cos(k * x) ** 2 * d))


def Lam2(s, t):
    n = max(2, math.ceil((t - s) / (lam / 6)))
    x, w = gl(s, t, n)
    d = np.sqrt(detM(s, x) * detM(x, t) / detM(s, t))
    inner = np.array([Lam(xi, t) for xi in x])
    return float(np.sum(w * inner * c * np.cos(k * x) ** 2 * d))


def beta(s, t):
    return pref * Iosc(s, t)


def Mmat(s, t):
    b = beta(s, t)
    return m(s, t) * np.eye(2) - np.array([[b.real, b.imag], [b.imag, -b.real]])


def kap(P1, P2):
    (u1, v1, y1), (u2, v2, y2) = P1, P2
    dy = y2 - y1
    return v2 - v1 - dy @ np.linalg.solve(Mmat(u1, u2), dy) / 2


def slice_data(P1, P2, u, ny):
    (u1, v1, y1), (u2, v2, y2) = P1, P2
    A1 = np.linalg.inv(Mmat(u1, u))
    A2 = np.linalg.inv(Mmat(u, u2))
    # bounding box for the support: both slacks >= 0 needs |y - y1|^2 <= 2 dv eigmax(M1) etc.
    R = math.sqrt(2 * (v2 - v1) * max(np.linalg.eigvalsh(np.linalg.inv(A1)).max(), 1e-30))
    R = min(R, math.sqrt(2 * (v2 - v1) * np.linalg.eigvalsh(np.linalg.inv(A2)).max()))
    ctr = 0.5 * (y1 + y2)
    Rb = R + np.linalg.norm(y2 - y1)
    g = np.linspace(-Rb, Rb, ny)
    X1, X2 = np.meshgrid(g + ctr[0], g + ctr[1], indexing="ij")
    d1 = np.stack([X1 - y1[0], X2 - y1[1]], -1)
    d2 = np.stack([y2[0] - X1, y2[1] - X2], -1)
    lo = v1 + np.einsum("...i,ij,...j", d1, A1, d1) / 2
    hi = v2 - np.einsum("...i,ij,...j", d2, A2, d2) / 2
    return lo, hi, (g[1] - g[0]) ** 2


def main():
    out = []
    P1 = (-0.4, 0.05, np.array([0.01, -0.02]))
    P2 = (0.35, 0.95, np.array([0.06, 0.03]))
    K = kap(P1, P2)
    L = Lam(P1[0], P2[0])
    # (i)
    for ny in (401, 801):
        x, w = gl(P1[0], P2[0], 24)
        vol = 0.0
        for xi, wi in zip(x, w):
            lo, hi, dA = slice_data(P1, P2, xi, ny)
            vol += wi * c * math.cos(k * xi) ** 2 * np.sum(np.maximum(0, hi - lo)) * dA
        out.append(f"(i) ny={ny}: grid vol = {vol:.8f}, pi kappa^2 Lambda = {math.pi*K**2*L:.8f},"
                   f" rel diff = {vol/(math.pi*K**2*L)-1:.2e}")
    # (ii)
    L2 = Lam2(P1[0], P2[0])
    for ny in (301, 601):
        x, w = gl(P1[0], P2[0], 16)
        tot = 0.0
        for xi, wi in zip(x, w):
            lo, hi, dA = slice_data(P1, P2, xi, ny)
            hlen = np.maximum(0, hi - lo)
            # kappa(R,P2) = hi - v_R runs over [0, hlen]; int kappa^2 dv_R = hlen^3/3
            inner = np.sum(hlen**3 / 3) * dA
            tot += wi * c * math.cos(k * xi) ** 2 * math.pi * Lam(xi, P2[0], 24) * inner
        out.append(f"(ii) ny={ny}: pair measure = {tot:.10f}, pi^2 kappa^4 Lambda_2/6 = "
                   f"{math.pi**2*K**4*L2/6:.10f}, rel diff = {tot/(math.pi**2*K**4*L2/6)-1:.2e}")
    out.append(f"pair fraction Lambda_2/(3 Lambda^2) = {L2/(3*L**2):.6f} (flat value 0.1)")
    print("\n".join(out))


if __name__ == "__main__":
    main()
