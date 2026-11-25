# ruff: noqa: E501, B023
"""Referee-side deterministic checks for route_c/leading_coefficient.md.

No order is sampled. No Monte Carlo of the Section 6 claim. Only identities,
expansions and deterministic geometry (quadrature / root finding).

Copied into the repository from the referee's scratch folder. Two changes
only: the absolute path to the geometry module was made relative, and the
file-level lint exemption above was added.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent / "proof3b_wave_signal"
sys.path.insert(0, str(REPO))
import wavegeom as wg  # noqa: E402


def mu0(m1, v1, z1, m2, v2, z2):
    """Flat background margin in the (m, v, zeta) chart, m2 > m1."""
    return v2 - v1 - np.abs(z2 - z1) ** 2 / (2.0 * (m2 - m1))


# ---------------------------------------------------------------- A
def check_gauge_family():
    """Per-event shear field s(m) = 2 a0 m + a1 m^2 is pure gauge at first order.

    Map (m, v', z') -> (m, v' + Re(a1 conj z'^2)/2, z' + (a0 + a1 m) conj z').
    Claim: mu0(mapped) - mu0(orig) = (1/2) Re((s_Y - s_Z) conj(omega)^2) + O(a^2)
    for EVERY pair (not only null), omega = dz'/dm.
    """
    rng = np.random.default_rng(1)
    n = 5000
    m1 = rng.uniform(0, 0.5, n)
    m2 = m1 + rng.uniform(0.05, 0.5, n)
    v1, v2 = rng.uniform(0, 1, n), rng.uniform(0, 1, n)
    z1 = rng.normal(0, 0.3, n) + 1j * rng.normal(0, 0.3, n)
    z2 = rng.normal(0, 0.3, n) + 1j * rng.normal(0, 0.3, n)
    print("A. gauge family s(m) = 2 a0 m + a1 m^2 (first-order identity)")
    for eps in (1e-2, 1e-3, 1e-4):
        a0 = eps * (0.7 - 0.4j)
        a1 = eps * (-0.3 + 0.9j)

        def fwd(m, v, z):
            al = a0 + a1 * m
            return v + 0.5 * (a1 * np.conj(z) ** 2).real, z + al * np.conj(z)

        V1, Z1 = fwd(m1, v1, z1)
        V2, Z2 = fwd(m2, v2, z2)
        shift = mu0(m1, V1, Z1, m2, V2, Z2) - mu0(m1, v1, z1, m2, v2, z2)
        om = (z2 - z1) / (m2 - m1)
        s1 = 2 * a0 * m1 + a1 * m1**2
        s2 = 2 * a0 * m2 + a1 * m2**2
        pred = 0.5 * ((s1 - s2) * np.conj(om) ** 2).real
        err = np.max(np.abs(shift - pred))
        print(
            f"   |a|~{eps:.0e}: max|shift| {np.max(np.abs(shift)):.3e}, "
            f"max|shift - pred| {err:.3e}, err/eps^2 {err / eps**2:.3f}"
        )
    # Jacobian of z' -> z' + al conj z' as a real 2x2 map
    al = 0.01 * (0.3 + 0.8j)
    J = np.array([[1 + al.real, al.imag], [al.imag, 1 - al.real]])
    print(f"   det of transverse Jacobian = {np.linalg.det(J):.10f} = 1 - |al|^2 = {1 - abs(al) ** 2:.10f}")


# ---------------------------------------------------------------- B
def ray_root(A, B, C):
    """Positive root of A rho^2 + B rho + C = 0 (C < 0, A > 0)."""
    return (-B + np.sqrt(B * B - 4 * A * C)) / (2 * A)


def region_extent(c, w, r, b, theta, absorbed=True):
    """Radial extent from c along e^{i theta} of {u : |u-w|^2 < r^2 - E(u)}.

    absorbed=True : E(u) = (1/2) Re(b (conj u^2 - conj w^2))   (T_P residual)
    absorbed=False: E(u) = (1/2) Re(b conj w^2)                 (unabsorbed P shift)
    """
    e = np.exp(1j * theta)
    d = c - w
    # |d + rho e|^2 = |d|^2 + 2 rho Re(conj(d) e) + rho^2
    A = np.ones_like(theta)
    B = 2 * (np.conj(d) * e).real
    C = np.abs(d) ** 2 - r**2 + np.zeros_like(theta)
    if absorbed:
        # E(c + rho e) = (1/2)Re(b(conj c^2 - conj w^2)) + rho Re(b conj c conj e)
        #               + (1/2) rho^2 Re(b conj e^2)
        A = A + 0.5 * (b * np.conj(e) ** 2).real
        B = B + (b * np.conj(c) * np.conj(e)).real
        C = C + 0.5 * (b * (np.conj(c) ** 2 - np.conj(w) ** 2)).real
    else:
        C = C + 0.5 * (b * np.conj(w) ** 2).real
    return ray_root(A, B, C)


def area(c, ws, r, b, absorbed, nth=200000):
    th = (np.arange(nth) + 0.5) * 2 * np.pi / nth
    rho = np.min([region_extent(c, w, r, b, th, absorbed) for w in ws], axis=0)
    return 0.5 * np.sum(rho**2) * 2 * np.pi / nth


def check_sliver_first_order():
    """First-order (in P's shear b) change of sliver cross-section areas."""
    print("B. first-order area change of sliver cross-sections from P's own shear")
    w0 = 0.8 + 0.3j
    h = 1e-6
    print("   isolated disk: dA/db (max over phase), absorbed vs unabsorbed")
    for r in (0.2, 0.1, 0.05):
        best_abs, best_naive = 0.0, 0.0
        for psi in np.linspace(0, np.pi, 7):
            e = np.exp(1j * psi)
            da = (area(w0, [w0], r, h * e, True) - area(w0, [w0], r, -h * e, True)) / (2 * h)
            dn = (area(w0, [w0], r, h * e, False) - area(w0, [w0], r, -h * e, False)) / (2 * h)
            best_abs, best_naive = max(best_abs, abs(da)), max(best_naive, abs(dn))
        print(f"     r={r:5.3f}: absorbed {best_abs:.2e}   unabsorbed {best_naive:.3e} (= pi|w|^2/2 = {np.pi * abs(w0) ** 2 / 2:.3e})")
    print("   lens of two overlapping slivers (|w1-w2| = r): dA/db over phase")
    for r in (0.2, 0.1, 0.05, 0.025):
        w1, w2 = w0 - 0.5 * r * np.exp(0.4j), w0 + 0.5 * r * np.exp(0.4j)
        c = 0.5 * (w1 + w2)
        lens0 = area(c, [w1, w2], r, 0.0, True)
        best_abs, best_naive = 0.0, 0.0
        for psi in np.linspace(0, np.pi, 13):
            e = np.exp(1j * psi)
            da = (area(c, [w1, w2], r, h * e, True) - area(c, [w1, w2], r, -h * e, True)) / (2 * h)
            dn = (area(c, [w1, w2], r, h * e, False) - area(c, [w1, w2], r, -h * e, False)) / (2 * h)
            best_abs, best_naive = max(best_abs, abs(da)), max(best_naive, abs(dn))
        print(
            f"     r={r:6.4f}: lens area {lens0:.3e}; absorbed {best_abs:.3e}, unabsorbed {best_naive:.3e}; "
            f"ratio {best_abs / best_naive:.3e}; ratio/r^2 {best_abs / best_naive / r**2:.3f}"
        )


# ---------------------------------------------------------------- C
def check_boundary_channel():
    """Critical-pair sliver volume vs a shift of the diamond boundary.

    Flat (m, v, zeta) chart, B = origin, T = (1, 1, 0). In the T_P frame the
    past exit condition is mu0(B, Z) >= +delta and the future one
    mu0(W, T) >= -delta, delta = phi_P(w). V = pi eps^2 [int_0^Sp s(L+s)/L ds
    + int_0^Sf s(L+s)/L ds]. Report (dV/d delta)/V, and the same ratio for
    the unabsorbed shift of P's cone (2/eps), at eps = 1e-2.
    """
    print("C. boundary channel: d log V / d delta for a critical pair (flat diamond)")
    B = (0.0, 0.0, 0.0j)
    T = (1.0, 1.0, 0.0j)
    eps = 1e-2

    def exit_past(P, w, delta):
        lo, hi = 0.0, P[0] - 1e-12
        for _ in range(200):
            s = 0.5 * (lo + hi)
            Z = (P[0] - s, P[1] - s * abs(w) ** 2 / 2, P[2] - s * w)
            if mu0(*B, *Z) >= delta:
                lo = s
            else:
                hi = s
        return lo

    def exit_future(X, w, delta):
        lo, hi = 0.0, T[0] - X[0] - 1e-12
        for _ in range(200):
            s = 0.5 * (lo + hi)
            W = (X[0] + s, X[1] + s * abs(w) ** 2 / 2, X[2] + s * w)
            if mu0(*W, *T) >= -delta:
                lo = s
            else:
                hi = s
        return lo

    def vol(P, w, L, delta):
        X = (P[0] + L, P[1] + L * abs(w) ** 2 / 2 - eps, P[2] + L * w)
        sp, sf = exit_past(P, w, delta), exit_future(X, w, delta)
        F = lambda S: S**2 / 2 + S**3 / (3 * L)  # noqa: E731
        return np.pi * eps**2 * (F(sp) + F(sf)), sp, sf

    cases = [
        ((0.45, 0.50, 0.05 + 0.02j), 0.3 + 0.1j, 0.15),
        ((0.45, 0.50, 0.05 + 0.02j), -0.4 + 0.2j, 0.15),
        ((0.30, 0.40, -0.10 + 0.05j), 0.5 - 0.2j, 0.20),
        ((0.60, 0.55, 0.00 + 0.10j), 0.1 + 0.6j, 0.10),
    ]
    d = 1e-6
    for P, w, L in cases:
        v0, sp, sf = vol(P, w, L, 0.0)
        vp, _, _ = vol(P, w, L, d)
        vm, _, _ = vol(P, w, L, -d)
        dlog = (vp - vm) / (2 * d) / v0
        print(
            f"   P={P}, w={w}, L={L}: Sp={sp:.3f} Sf={sf:.3f}; dlogV/ddelta={dlog:+.3f}; "
            f"unabsorbed 2/eps={2 / eps:.0f}; ratio={dlog / (2 / eps):+.2e}"
        )


# ---------------------------------------------------------------- D
def check_transitivity():
    """Free per-event shears in the first-order rule need not give a partial order."""
    print("D. first-order rule with free per-event shears: transitivity")
    sig = 1e-3
    delta = 1e-3
    om1, om2 = 1.0 + 0j, 1.0 + delta
    L1 = L2 = 1.0
    sY, sZ, sW = 0.0, sig, 0.0

    def rel(p, q, sp, sq):
        om = (q[2] - p[2]) / (q[0] - p[0])
        return mu0(*p, *q) + 0.5 * ((sp - sq) * np.conj(om) ** 2).real

    a = 0.5 * sig + 1e-9
    b = -0.5 * (sig * np.conj(om2) ** 2).real + 1e-9
    Y = (0.0, 0.0, 0j)
    Z = (L1, L1 * abs(om1) ** 2 / 2 + a, L1 * om1)
    W = (L1 + L2, Z[1] + L2 * abs(om2) ** 2 / 2 + b, Z[2] + L2 * om2)
    print(
        f"   mu(Y,Z)={rel(Y, Z, sY, sZ):+.2e}  mu(Z,W)={rel(Z, W, sZ, sW):+.2e}  "
        f"mu(Y,W)={rel(Y, W, sY, sW):+.2e}   (Y<Z, Z<W but not Y<W if last < 0)"
    )


# ---------------------------------------------------------------- E
def check_short_pairs():
    """Eq. (1) for pairs with du ~ lambda: error vs |beta|^2 |omega|^2 / d."""
    print("E. eq. (1) at short separation du in [lam/4, 4 lam]")
    rng = np.random.default_rng(3)
    n = 4000
    for r in (64, 256, 1024):
        w = wg.Wave(1.0, r, 0.37)
        lam = 1.0 / r
        u1 = rng.uniform(-0.4, 0.3, n)
        u2 = u1 + lam * rng.uniform(0.25, 4.0, n)
        d = w.dm(u1, u2)
        z1 = rng.normal(0, 0.05, n) + 1j * rng.normal(0, 0.05, n)
        dz = d * (rng.normal(0, 1, n) + 1j * rng.normal(0, 1, n))
        v1 = np.zeros(n)
        v2 = np.abs(dz) ** 2 / (2 * d)
        exact = w.margin(u1, v1, z1, u2, v2, z1 + dz) - (v2 - v1 - np.abs(dz) ** 2 / (2 * d))
        om = dz / d
        b1 = w.beta(-0.5 * np.ones(n), u1)
        b2 = w.beta(-0.5 * np.ones(n), u2)
        first = 0.5 * ((b1 - b2) * np.conj(om) ** 2).real
        scale = np.abs(w.beta(u1, u2)) ** 2 * np.abs(om) ** 2 / d
        print(
            f"   r={r:5d}: max|err|/scale {np.max(np.abs(exact - first) / scale):.3f}; "
            f"max|err|/(lam^3 |om|^2) {np.max(np.abs(exact - first) / (lam**3 * np.abs(om) ** 2)):.3f}; "
            f"max|first|/(lam |om|^2 d) {np.max(np.abs(first) / (lam * np.abs(om) ** 2 * d)):.3f}"
        )


if __name__ == "__main__":
    check_gauge_family()
    check_sliver_first_order()
    check_boundary_channel()
    check_transitivity()
    check_short_pairs()
