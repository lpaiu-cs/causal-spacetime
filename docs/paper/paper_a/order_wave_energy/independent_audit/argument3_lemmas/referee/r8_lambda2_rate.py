# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of gram_fourier_lemmas.md.
"""Referee check r8: the Lambda_2 comparison of Corollary K2(ii), rate in lambda.

dev2(lam) = Lambda_2lam(s,t)/c^2 - Lambda_20(s,t) for (s,t) = (-1/2, 0.2), q = 1, computed
as a difference integral (stable: r - 1 via log1p/expm1, background inner integral exact via
alpha), at lam = 1/32 ... 1/256. Fit dev2 = lam^4 (A log(1/lam) + B) to test whether the
order is lam^4 log (A finite, stable fit) rather than lam^4 log^2 or lam^3.
Also the same for Lambda (one integral) for comparison.
"""
import math
import numpy as np

GN, GW = np.polynomial.legendre.leggauss(12)
q = 1.0
k = 1.0
phi = 0.37


def f(w):
    return 1 / np.cos(k * w) ** 2


def a2(w):
    return np.cos(k * w) ** 2


def m(s, t):
    return np.sin(k * (t - s)) / (k * np.cos(k * s) * np.cos(k * t))


def alpha(l):
    x = k * np.asarray(l, float)
    out = np.where(x > 1e-4, (1 - x / np.tan(np.maximum(x, 1e-300))) / (2 * q), x**2 / (6 * q) * (1 + x**2 / 15))
    return out


def gl(a, b, h):
    n = max(1, math.ceil((b - a) / h))
    e = np.linspace(a, b, n + 1)
    lo, hi = e[:-1, None], e[1:, None]
    return (0.5 * (hi - lo) * GN + 0.5 * (hi + lo)).ravel(), (0.5 * (hi - lo) * GW).ravel()


class Wave:
    def __init__(self, lam):
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


def rho(wv, Ps, Pt, s, t):
    return np.abs(wv.pref * (Pt - Ps)) ** 2 / m(s, t) ** 2


def rm1(r1, r2, r0):
    return np.expm1(0.5 * (np.log1p(-r1) + np.log1p(-r2) - np.log1p(-r0)))


def dev_Lambda(wv, s, t, h):
    """Lambda_lam(s,t)/c - Lambda_0(s,t) as one difference integral."""
    x, w = gl(s, t, h)
    Px, Ps, Pt = wv.prim(x), wv.prim(s)[0], wv.prim(t)[0]
    K0 = a2(x) * m(s, x) * m(x, t) / m(s, t)
    d = rm1(rho(wv, Ps, Px, s, x), rho(wv, Px, Pt, x, t), rho(wv, Ps, Pt, s, t))
    return float(np.sum(w * K0 * d))


def dev_Lambda2(wv, s, t, h):
    x, w = gl(s, t, h)
    Px, Ps, Pt = wv.prim(x), wv.prim(s)[0], wv.prim(t)[0]
    K0 = a2(x) * m(s, x) * m(x, t) / m(s, t)
    d = rm1(rho(wv, Ps, Px, s, x), rho(wv, Px, Pt, x, t), rho(wv, Ps, Pt, s, t))
    D = np.array([dev_Lambda(wv, xi, t, h) if t - xi > 1e-12 else 0.0 for xi in x])
    return float(np.sum(w * (D * K0 * (1 + d) + alpha(t - x) * K0 * d)))


def main():
    out = []
    s, t = -0.5, 0.2
    rows = []
    for lam in (1 / 32, 1 / 64, 1 / 128, 1 / 256):
        wv = Wave(lam)
        h = lam / 4
        d1 = dev_Lambda(wv, s, t, h)
        d2 = dev_Lambda2(wv, s, t, h)
        d2b = dev_Lambda2(wv, s, t, h / 2) if lam >= 1 / 64 else float("nan")
        rows.append((lam, d1, d2))
        out.append(f"lam=1/{round(1/lam)}: dev1/lam^4={d1/lam**4:.5f}  dev2/lam^4={d2/lam**4:.5f}"
                   f"  (h/2 refinement dev2/lam^4={d2b/lam**4:.5f})"
                   f"  dev2/(lam^4(1+log))={abs(d2)/(lam**4*(1+math.log(1/lam))):.4f}")
    # fit dev/lam^4 = A L + B on consecutive pairs
    for idx, name in ((1, "Lambda"), (2, "Lambda_2")):
        Ls = [math.log(1 / r[0]) for r in rows]
        ys = [r[idx] / r[0] ** 4 for r in rows]
        As = [(ys[i + 1] - ys[i]) / (Ls[i + 1] - Ls[i]) for i in range(len(rows) - 1)]
        out.append(f"{name}: successive slopes d(dev/lam^4)/d log(1/lam) = " + ", ".join(f"{a:.4f}" for a in As)
                   + "  (constant slope => lam^4 log; growing slope => worse)")
    print("\n".join(out))


if __name__ == "__main__":
    main()
