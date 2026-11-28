# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of gram_fourier_lemmas.md.
"""Referee check r4: Proposition L (the signal (11)), independent computation.

Exact canonical map, exact wave slacks, exact perturbed anchor Gram (so the delta-ell term
is included), but calibration taken as exact (the O(lam^4 log) calibration error is the
separate Corollary K2 question, checked in r2). Computes
  E_num(omega) = int (a^2/Z0) w(u, v~, |z|^2) H(X) e^{-i omega u}  over canonical D_0
by Gauss quadrature in u, |z|^2 (v~ factorises: H does not depend on v~), and an exact
angular average (H is a quadratic polynomial in z_X; 8 angles). Compares with Step 5's
main term  E_pred = -(i k/c) lam^2 e^{i phi} int W f D_A e^{-i delta u} du.
Also: invariance of |E| under rotation and reflection of the anchor frame; the actual
value of |int W f D_A e^{-i delta u}| / int W f against the certified floor 0.0101678;
the normalisation against challenge (4).
"""
import math
import numpy as np

GN, GW = np.polynomial.legendre.leggauss(16)
q = 1.0
k = math.sqrt(q)
M = math.tan(k / 2) * 2 / k


def f(w):
    return 1 / np.cos(k * w) ** 2


def a2(w):
    return np.cos(k * w) ** 2


def mm(u):  # m(-1/2, u)
    return (np.tan(k * u) + math.tan(k / 2)) / k


def gl(a, b, n):
    e = np.linspace(a, b, n + 1)
    lo, hi = e[:-1, None], e[1:, None]
    return (0.5 * (hi - lo) * GN + 0.5 * (hi + lo)).ravel(), (0.5 * (hi - lo) * GW).ravel()


def bump(x):
    out = np.zeros_like(x)
    i = np.abs(x) < 1
    out[i] = np.exp(-1 / (1 - x[i] ** 2))
    return out


class Wave:
    def __init__(self, lam, phi):
        self.lam = lam
        self.c = math.sqrt(1 + 4 * q * lam**2)
        self.pref = 2 * lam * k / self.c * np.exp(1j * phi)
        self.E = np.linspace(-0.5, 0.5, int(math.ceil(8 / lam)) + 1)
        self.P = np.concatenate([[0], np.cumsum(self._pi(self.E[:-1], self.E[1:]))])

    def _pi(self, a, b):
        a, b = np.atleast_1d(a), np.atleast_1d(b)
        x = a[:, None] + (b - a)[:, None] * (GN[None, :] + 1) / 2
        return ((b - a)[:, None] * GW[None, :] / 2 * f(x) * np.exp(1j * x / self.lam)).sum(1)

    def prim(self, u):
        u = np.atleast_1d(np.asarray(u, float))
        j = np.clip(np.searchsorted(self.E, u, side="right") - 1, 0, len(self.E) - 2)
        return self.P[j] + self._pi(self.E[j], u)

    def bL(self, u):
        return self.pref * self.prim(u)

    def btot(self):
        return self.pref * self.prim(0.5)[0]


def bmats(b):
    """array of 2x2 matrices from complex array."""
    out = np.zeros(b.shape + (2, 2))
    out[..., 0, 0], out[..., 0, 1] = b.real, b.imag
    out[..., 1, 0], out[..., 1, 1] = b.imag, -b.real
    return out


def canon(wv, u):
    """T(u) and A(u) = T M_L^-1 T - I/m (so delta v = z^T A z / 2), vectorised in u."""
    u = np.atleast_1d(u)
    m = mm(u)
    mr = M - m
    bl = wv.bL(u)
    br = wv.btot() - bl
    eL, eR = np.abs(bl) / m, np.abs(br) / mr
    EL, ER = bmats(bl) / m[:, None, None], bmats(br) / mr[:, None, None]
    I = np.eye(2)
    P = (I + EL) / (1 - eL**2)[:, None, None]
    S = (I + ER) / (1 - eR**2)[:, None, None]
    sig = (m / M)[:, None, None]
    X = (1 - sig) * P + sig * S
    w, V = np.linalg.eigh(X)
    T = V @ (w[..., :, None] ** -0.5 * np.swapaxes(V, -1, -2))
    A = -(1 / M) * T @ (S - P) @ T
    return T, A, bl


def setup_anchors(wv, theta, reflect=False):
    d = 0.08
    R = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    F = np.diag([1, -1]) if reflect else np.eye(2)
    zs = [R @ F @ np.array(v) for v in ((d, 0.0), (0.0, d), (d, d))]
    us = np.array([-0.35, -0.25, -0.15])
    vts = np.array([0.15, 0.25, 0.35])
    T, A, bl = canon(wv, us)
    ys = [T[i] @ zs[i] for i in range(3)]
    vs = [vts[i] + zs[i] @ A[i] @ zs[i] / 2 for i in range(3)]
    return dict(us=us, vts=vts, zs=zs, ys=ys, vs=vs, bl=bl)


def kap_pair(wv, ui, vi, yi, uj, vj, yj):
    """exact wave slack, vectorised over j (uj array, vj array, yj (n,2))."""
    uj = np.atleast_1d(uj)
    dm = mm(uj) - mm(ui)
    b = wv.bL(uj) - wv.bL(ui)[0]
    Mt = dm[:, None, None] * np.eye(2) - bmats(b)
    dy = yj - yi
    sol = np.linalg.solve(Mt, dy[..., None])[..., 0]
    return vj - vi - 0.5 * np.sum(dy * sol, -1), dm


def DA(u, us, zs, ell):
    out = 0
    for i in range(3):
        zeta = zs[i][0] + 1j * zs[i][1]
        out = out - 0.5 * ell[i] * np.conj(zeta) ** 2 / (mm(u) - mm(us[i]))
    return out


def signal(lam, phi, theta, reflect=False, deltas=(0.0, 0.125, -0.125)):
    wv = Wave(lam, phi)
    an = setup_anchors(wv, theta, reflect)
    us, zs = an["us"], an["zs"]
    # perturbed anchor Gram (exact slacks, exact canonical radii)
    g = np.zeros((3, 3))
    for i in range(3):
        g[i, i] = zs[i] @ zs[i]
        for j in range(i + 1, 3):
            kl, dm = kap_pair(wv, us[i], an["vs"][i], an["ys"][i], np.array([us[j]]),
                              np.array([an["vs"][j]]), an["ys"][j][None, :])
            kc = an["vts"][j] - an["vts"][i] - np.sum((zs[j] - zs[i]) ** 2) / (2 * dm[0])
            g[i, j] = g[j, i] = zs[i] @ zs[j] + dm[0] * (kl[0] - kc)
    GA = g[:2, :2]
    l12 = -np.array([g[0, 2], g[1, 2]]) @ np.linalg.inv(GA)
    ell = np.array([l12[0], l12[1], 1.0])
    G0 = np.array([[zs[i] @ zs[j] for j in range(3)] for i in range(3)])
    l0 = -np.array([G0[0, 2], G0[1, 2]]) @ np.linalg.inv(G0[:2, :2])
    ell0 = np.array([l0[0], l0[1], 1.0])
    # quadrature
    u, wu = gl(0.1, 0.3, int(math.ceil(0.2 / (lam / 4))))
    r2n, r2w = gl(0.0025, 0.04, 2)
    vtn, vtw = gl(0.65, 0.75, 1)
    Wv = np.sum(vtw * bump((vtn - 0.7) / 0.05))
    T, A, bl = canon(wv, u)
    nang = 8
    Hbar = np.zeros(len(u))  # radial-weighted angular mean, times pi
    for r2, w2 in zip(r2n, r2w):
        wt = bump(np.array([(r2 - 0.02125) / 0.01875]))[0]
        if wt == 0:
            continue
        acc = np.zeros(len(u))
        for ang in np.arange(nang) * 2 * math.pi / nang:
            z = math.sqrt(r2) * np.array([math.cos(ang), math.sin(ang)])
            y = T @ z
            dvX = 0.5 * np.einsum("i,nij,j->n", z, A, z)
            H = np.zeros(len(u))
            for i in range(3):
                vX = 0.7 + dvX  # v~ cancels in kappa_lam - kappa_can, any value works
                kl, dm = kap_pair(wv, us[i], an["vs"][i], an["ys"][i], u, vX, y)
                kc = 0.7 - an["vts"][i] - np.sum((z - zs[i]) ** 2) / (2 * dm)
                H += ell[i] * (zs[i] @ z + dm * (kl - kc))
            acc += H
        Hbar += w2 * wt * math.pi * acc / nang
    Z0 = math.pi * (1 - 1 / math.tan(1.0)) / 2 / q  # pi alpha(1) at q=1
    Wfac = Wv / Z0
    W_u = a2(u) * bump((u - 0.2) / 0.1) * Wfac * np.sum(r2w * bump((r2n - 0.02125) / 0.01875)) * math.pi
    res = {}
    for dl in deltas:
        om = 1 / lam + dl
        Enum = np.sum(wu * a2(u) * bump((u - 0.2) / 0.1) * Wfac * Hbar * np.exp(-1j * om * u))
        Epred = -(1j * k / wv.c) * lam**2 * np.exp(1j * phi) * np.sum(
            wu * W_u * f(u) * DA(u, us, zs, ell0) * np.exp(-1j * dl * u))
        mod = abs(np.sum(wu * W_u * f(u) * DA(u, us, zs, ell0) * np.exp(-1j * dl * u))) / np.sum(wu * W_u * f(u))
        Ew = np.sum(wu * W_u)
        res[dl] = dict(Enum=Enum, Epred=Epred, rel=abs(Enum - Epred) / abs(Epred), mod=mod,
                       norm=Enum / (lam**2 * Ew))
    pred4 = -1j * k * np.exp(1j * phi) * np.sum(wu * bump((u - 0.2) / 0.1) * DA(u, us, zs, ell0)) / np.sum(
        wu * bump((u - 0.2) / 0.1) * a2(u))
    return res, pred4, np.max(np.abs(ell - ell0)) / lam**2


def main():
    out = []
    for lam in (1 / 64, 1 / 128, 1 / 256, 1 / 512):
        for (phi, theta, refl) in ((0.37, 0.0, False), (2.1, math.pi / 8, False), (0.37, 0.3, True)):
            res, pred4, dl_ell = signal(lam, phi, theta, refl)
            row = [f"lam=1/{round(1/lam)} phi={phi} theta={theta:.3f} refl={refl} |dell|/lam^2={dl_ell:.2f}:"]
            for dl, r in res.items():
                row.append(f"  delta={dl:+.3f}: |E_num|/lam^2={abs(r['Enum'])/lam**2:.6e}"
                           f" |E_pred|/lam^2={abs(r['Epred'])/lam**2:.6e} rel diff={r['rel']:.3e}"
                           f" |int W f D_A e|/int W f={r['mod']:.6f}")
            r0 = res[0.0]
            row.append(f"  challenge (4): E/(lam^2 E w) = {r0['norm']:.6f}, predicted limit = {pred4:.6f}")
            out.append("\n".join(row))
    print("\n".join(out))


if __name__ == "__main__":
    main()
