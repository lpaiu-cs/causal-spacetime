# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written outside the repository by the fresh-context AI referee of one_point_strip_bound.md.
"""Referee r4.
(i) Step 1 inequalities with np.linalg.inv on seeded random symmetric trace-free E
    (inequality test only): A in [1/(m+e), 1/(m-e)], |A-B| <= e/(m(m-e)), and the
    sharper one-sided bound (B-A)^+ <= e/(m(m+e)).
(ii) Remark 4 numerically: E at phase phi equals R(phi/2) E(phase 0) R^T.
(iii) Cartesian brute force of the y-integral of the v-section lengths of D_lam,
    D_0, D_cap at u near the tips for q = 9.8 (pole region), against
    (h^2/2) int dtheta / S (theta by kink-aware GL).
"""
import math
import numpy as np

rng = np.random.default_rng(7)
worst_a = worst_ab = worst_sharp = 0.0
viol = 0
for _ in range(200000):
    m = rng.uniform(1e-3, 5)
    e = m * rng.uniform(0, 0.999)
    ang = rng.uniform(0, 2 * np.pi)
    E = e * np.array([[np.cos(ang), np.sin(ang)], [np.sin(ang), -np.cos(ang)]])
    th = rng.uniform(0, 2 * np.pi)
    n = np.array([np.cos(th), np.sin(th)])
    A = n @ np.linalg.inv(m * np.eye(2) + E) @ n
    B = 1 / m
    lo, hi = 1 / (m + e), 1 / (m - e)
    tol = 1e-12 * hi
    if not (lo - tol <= A <= hi + tol):
        viol += 1
    if abs(A - B) > e / (m * (m - e)) * (1 + 1e-9) + tol:
        viol += 1
    if max(B - A, 0) > e / (m * (m + e)) * (1 + 1e-9) + tol:
        viol += 1
    if e > 0:
        worst_ab = max(worst_ab, abs(A - B) / (e / (m * (m - e))))
        worst_sharp = max(worst_sharp, max(B - A, 0) / (e / (m * (m + e))))
print(f"(i) Step 1 random test: violations {viol}; max |A-B|/(e/(m(m-e))) = {worst_ab:.9f}; "
      f"max (B-A)^+/(e/(m(m+e))) = {worst_sharp:.9f}")


def T(x):
    return np.array([[np.cos(x), np.sin(x)], [np.sin(x), -np.cos(x)]])


def R(a):
    return np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])


worst = 0.0
for _ in range(2000):
    x, phi = rng.uniform(0, 2 * np.pi, 2)
    worst = max(worst, np.max(np.abs(R(phi / 2) @ T(x) @ R(phi / 2).T - T(x + phi))))
print(f"(ii) max |R(phi/2) T(x) R^T - T(x+phi)| = {worst:.2e}")

# (iii)
GX, GW = np.polynomial.legendre.leggauss(40)


def segG(lo, hi, s, lam, phi, n=4000):
    ed = np.linspace(lo, hi, n + 1)
    h = (ed[1:] - ed[:-1]) / 2
    md = (ed[1:] + ed[:-1]) / 2
    x = md[:, None] + h[:, None] * GX[None, :]
    fx = 1 / np.cos(s * x) ** 2
    return float(np.sum(h[:, None] * GW * fx)), complex(np.exp(1j * phi) * np.sum(h[:, None] * GW * fx * np.exp(1j * x / lam)))


def kinv(m, G, t):
    Mt = np.array([[m - t * G.real, -t * G.imag], [-t * G.imag, m + t * G.real]])
    return np.linalg.inv(Mt)


H = 1.0
for q, lam, phi, us in ((9.8, 1 / 8, 0.9, (-0.4995, 0.499, 0.4999)), (9.0, 1 / 32, 2.1, (-0.49, 0.0, 0.4999))):
    s = math.sqrt(q)
    t = 2 * lam * s / math.sqrt(1 + 4 * q * lam * lam)
    for u in us:
        m1, G1 = segG(-0.5, u, s, lam, phi)
        m2, G2 = segG(u, 0.5, s, lam, phi)
        K1, K2 = kinv(m1, G1, t), kinv(m2, G2, t)
        # polar: theta GL with many points over [0, pi) split in 64 pieces (kinks resolved by brute force)
        th_ed = np.linspace(0, np.pi, 4097)
        hh = (th_ed[1:] - th_ed[:-1]) / 2
        mm = (th_ed[1:] + th_ed[:-1]) / 2
        gx8, gw8 = np.polynomial.legendre.leggauss(8)
        th = (mm[:, None] + hh[:, None] * gx8).ravel()
        wt = (hh[:, None] * gw8).ravel()
        n = np.stack([np.cos(th), np.sin(th)])
        A1 = np.einsum("in,ij,jn->n", n, K1, n)
        A2 = np.einsum("in,ij,jn->n", n, K2, n)
        B1, B2 = 1 / m1, 1 / m2
        polar = {
            "lam": H * H * np.sum(wt / (A1 + A2)),
            "0": H * H * np.pi / (B1 + B2),
            "cap": H * H * np.sum(wt / (np.maximum(A1, B1) + np.maximum(A2, B2))),
        }
        # Cartesian: exact y2-integration is not used; dense midpoint grid
        rmax = 1.02 * math.sqrt(2 * H / min(np.linalg.eigvalsh(K1 + K2).min(), B1 + B2))
        ng = 3001
        yy = -rmax + (np.arange(ng) + 0.5) * 2 * rmax / ng
        y1, y2 = np.meshgrid(yy, yy, indexing="ij")
        cell = (2 * rmax / ng) ** 2

        def qf(K):
            return K[0, 0] * y1 * y1 + 2 * K[0, 1] * y1 * y2 + K[1, 1] * y2 * y2

        Q1, Q2 = qf(K1), qf(K2)
        r2 = y1 * y1 + y2 * y2
        lo_l, hi_l = Q1 / 2, H - Q2 / 2
        lo_0, hi_0 = r2 * B1 / 2, H - r2 * B2 / 2
        cart = {
            "lam": cell * np.sum(np.maximum(hi_l - lo_l, 0)),
            "0": cell * np.sum(np.maximum(hi_0 - lo_0, 0)),
            "cap": cell * np.sum(np.maximum(np.minimum(hi_l, hi_0) - np.maximum(lo_l, lo_0), 0)),
        }
        # differences (the quantities the TV is made of)
        dl_p, d0_p = polar["lam"] - polar["cap"], polar["0"] - polar["cap"]
        dl_c, d0_c = cart["lam"] - cart["cap"], cart["0"] - cart["cap"]
        rel = max(abs(cart[k] / polar[k] - 1) for k in polar)
        print(f"(iii) q={q} lam={lam} u={u}: max rel diff of volumes {rel:.2e}; "
              f"Delta_lam polar {dl_p:.6e} cart {dl_c:.6e}; Delta_0 polar {d0_p:.6e} cart {d0_c:.6e}; e/m = {t*abs(G1)/m1:.3f},{t*abs(G2)/m2:.3f}")
