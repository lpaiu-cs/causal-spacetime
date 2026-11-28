# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of gram_fourier_lemmas.md.
"""Referee check r6: Lemma H (relation comparison), seeded test points (no order sampling).

Pairs of canonical points X, R in D_0 are placed near the background null cone (so that
mismatches occur), with many points near the tips. For each pair: physical points via the
exact inverse canonical map; wave relation of the physical pair (Theorem 2, y chart);
background relation of the canonical pair. On every mismatch check
  (1) |gt(X,R)| <= c0 lam^2 + sqrt2 (d(X) + d(R)),  d = exact (T,Z,y) displacement;
  (2) if |u_R - u_X| >= l0: |kappa_0| <= C (lam^2 + |dy_X| + |dy_R|); report the
      empirical constant max |kappa_0| / (lam^2 + |dy_X| + |dy_R|) for l0 = 0.07.
Also counts orientation flips (canonical X<R but physical background R<X), which the
'opposite sides of 0' step implicitly excludes.
"""
import math
import numpy as np

GN, GW = np.polynomial.legendre.leggauss(20)
q = 1.0
k = math.sqrt(q)
fmax = 1 / math.cos(k / 2) ** 2
CE = 2 * k * (4 * fmax - 2)
M = 2 * math.tan(k / 2) / k
c0 = 2 * math.sqrt(2) * CE


def f(w):
    return 1 / np.cos(k * w) ** 2


def mm(u):
    return (math.tan(k * u) + math.tan(k / 2)) / k


class Wave:
    def __init__(self, lam, phi):
        self.lam = lam
        c = math.sqrt(1 + 4 * q * lam**2)
        self.pref = 2 * lam * k / c * np.exp(1j * phi)
        self.E = np.linspace(-0.5, 0.5, int(math.ceil(8 / lam)) + 1)
        self.P = np.concatenate([[0], np.cumsum([self._pi(a, b) for a, b in zip(self.E[:-1], self.E[1:])])])

    def _pi(self, a, b):
        x = a + (b - a) * (GN + 1) / 2
        return (b - a) / 2 * np.sum(GW * f(x) * np.exp(1j * x / self.lam))

    def prim(self, u):
        j = min(max(int(np.searchsorted(self.E, u, side="right")) - 1, 0), len(self.E) - 2)
        return self.P[j] + self._pi(self.E[j], u)

    def beta(self, s, t):
        if t - s < 1e-3:
            return self.pref * self._pi(s, t)
        return self.pref * (self.prim(t) - self.prim(s))


def bm(b):
    return np.array([[b.real, b.imag], [b.imag, -b.real]])


def to_phys(wv, u, vt, z):
    m = mm(u)
    mr = M - m
    bl, br = wv.beta(-0.5, u), wv.beta(u, 0.5)
    ML = m * np.eye(2) - bm(bl)
    MR = mr * np.eye(2) - bm(br)
    W0 = m * mr / M
    X = W0 * (np.linalg.inv(ML) + np.linalg.inv(MR))
    w, V = np.linalg.eigh(X)
    T = V @ np.diag(w**-0.5) @ V.T
    y = T @ z
    eta = vt - z @ z / (2 * m)
    v = eta + y @ np.linalg.solve(ML, y) / 2
    return v, y


def wave_rel(wv, P1, P2):
    (u1, v1, y1), (u2, v2, y2) = P1, P2
    if u1 >= u2:
        return False
    Mt = (mm(u2) - mm(u1)) * np.eye(2) - bm(wv.beta(u1, u2))
    dy = y2 - y1
    return v2 - v1 - dy @ np.linalg.solve(Mt, dy) / 2 >= 0


def bg_rel(P1, P2):
    (u1, v1, y1), (u2, v2, y2) = P1, P2
    T1, T2 = (mm(u1) + v1) / math.sqrt(2), (mm(u2) + v2) / math.sqrt(2)
    Z1, Z2 = (mm(u1) - v1) / math.sqrt(2), (mm(u2) - v2) / math.sqrt(2)
    return T2 - T1 >= math.sqrt((Z2 - Z1) ** 2 + np.sum((y2 - y1) ** 2))


def gt(P1, P2):
    (u1, v1, y1), (u2, v2, y2) = P1, P2
    dT = (mm(u2) + v2 - mm(u1) - v1) / math.sqrt(2)
    dZ = (mm(u2) - v2 - mm(u1) + v1) / math.sqrt(2)
    return abs(dT) - math.sqrt(dZ**2 + np.sum((y2 - y1) ** 2))


def in_D0(u, vt, z):
    m = mm(u)
    return vt - z @ z / (2 * m) >= 0 and 1 - vt - z @ z / (2 * (M - m)) >= 0


def rand_point(rng, lam):
    r = rng.uniform()
    if r < 0.4:
        s = 10 ** rng.uniform(-7, math.log10(0.05))
        u = -0.5 + s if rng.uniform() < 0.5 else 0.5 - s
    else:
        u = rng.uniform(-0.45, 0.45)
    m = mm(u)
    W0 = m * (M - m) / M
    rad = math.sqrt(2 * W0) * math.sqrt(rng.uniform()) * 0.999
    ang = rng.uniform(0, 2 * math.pi)
    z = rad * np.array([math.cos(ang), math.sin(ang)])
    lo, hi = z @ z / (2 * m), 1 - z @ z / (2 * (M - m))
    vt = rng.uniform(lo, hi)
    return u, vt, z


def main():
    rng = np.random.default_rng(1126)
    out = []
    for lam in (1 / 32, 1 / 128):
        for phi in (0.37, 2.1):
            wv = Wave(lam, phi)
            n_mis = n_pairs = flips = 0
            worst1 = 0.0
            worst2 = 0.0
            tries = 0
            while n_pairs < 6000 and tries < 60000:
                tries += 1
                uX, vX, zX = rand_point(rng, lam)
                uR, _, zR = rand_point(rng, lam)
                if uR == uX:
                    continue
                if uR < uX:
                    uX, zX, uR, zR = uR, zR, uX, zX
                    m = mm(uX)
                    lo, hi = zX @ zX / (2 * m), 1 - zX @ zX / (2 * (M - m))
                    if lo > hi:
                        continue
                    vX = rng.uniform(lo, hi)
                dm = mm(uR) - mm(uX)
                # place R near the canonical cone of X (offset of order lam^2 .. lam^1.5)
                off = rng.choice([-1, 1]) * (10 ** rng.uniform(-1, 1)) * lam**2 * rng.uniform(0, 3)
                vR = vX + np.sum((zR - zX) ** 2) / (2 * dm) + off
                if not in_D0(uR, vR, zR) or not in_D0(uX, vX, zX):
                    continue
                n_pairs += 1
                pvX, yX = to_phys(wv, uX, vX, zX)
                pvR, yR = to_phys(wv, uR, vR, zR)
                Xc, Rc = (uX, vX, zX), (uR, vR, zR)
                Xp, Rp = (uX, pvX, yX), (uR, pvR, yR)
                w_rel = (wave_rel(wv, Xp, Rp), wave_rel(wv, Rp, Xp))
                b_rel = (bg_rel(Xc, Rc), bg_rel(Rc, Xc))
                if bg_rel(Xc, Rc) and bg_rel(Rp, Xp):
                    flips += 1
                if w_rel != b_rel:
                    n_mis += 1
                    dX = math.hypot(pvX - vX, np.linalg.norm(yX - zX))
                    dR = math.hypot(pvR - vR, np.linalg.norm(yR - zR))
                    worst1 = max(worst1, abs(gt(Xc, Rc)) / (c0 * lam**2 + math.sqrt(2) * (dX + dR)))
                    if uR - uX >= 0.07:
                        k0 = vR - vX - np.sum((zR - zX) ** 2) / (2 * dm)
                        worst2 = max(worst2, abs(k0) / (lam**2 + np.linalg.norm(yX - zX) + np.linalg.norm(yR - zR)))
            out.append(f"lam=1/{round(1/lam)} phi={phi}: pairs={n_pairs} mismatches={n_mis} orientation flips={flips};"
                       f" max |gt|/(c0 lam^2 + sqrt2(dX+dR)) = {worst1:.4f} (<=1);"
                       f" H(2) empirical C (l0=0.07) = {worst2:.3f}")
    print("\n".join(out))


if __name__ == "__main__":
    main()
