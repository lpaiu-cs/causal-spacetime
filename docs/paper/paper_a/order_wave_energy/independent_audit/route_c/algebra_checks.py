"""Deterministic checks of the Route C derivation (leading_coefficient.md).

Generator coordinates only; no order is sampled and nothing here is a
measurement. Uses the referee-side geometry of ../proof3b_wave_signal.

    python algebra_checks.py
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

GEOMETRY = pathlib.Path(__file__).resolve().parent.parent / "proof3b_wave_signal"
sys.path.insert(0, str(GEOMETRY))
import wavegeom as wg  # noqa: E402

Q, PHI = 1.0, 0.37
SIZES = (64, 128, 256, 512, 1024)


def bg_margin(geometry, u1, v1, z1, u2, v2, z2):
    """Background margin of the pair (earlier 1, later 2)."""
    return v2 - v1 - np.abs(z2 - z1) ** 2 / (2.0 * geometry.dm(u1, u2))


def b_lead(w, u):
    """Leading oscillating part of the per-event shear beta(-1/2, u)."""
    lam = w.lam
    return -2j * w.k * lam**2 * w.f(u) * np.exp(1j * (u / lam + w.phi)) / w.cfac


def check_first_order():
    """Exact wave margin shift against phi_P - phi_X, same pairs at every r."""
    rng = np.random.default_rng(7)
    n = 4000
    u1 = rng.uniform(-0.4, 0.2, n)
    u2 = u1 + rng.uniform(0.1, 0.25, n)
    z1 = rng.normal(0, 0.05, n) + 1j * rng.normal(0, 0.05, n)
    dz = rng.normal(0, 0.08, n) + 1j * rng.normal(0, 0.08, n)
    v1 = rng.uniform(0, 0.2, n)
    print("1. max |exact shift - first order| / lam^4 (99th percentile)")
    for r in SIZES:
        w = wg.Wave(Q, r, PHI)
        d = w.dm(u1, u2)
        v2 = v1 + np.abs(dz) ** 2 / (2 * d)
        exact = w.margin(u1, v1, z1, u2, v2, z1 + dz) - bg_margin(
            w, u1, v1, z1, u2, v2, z1 + dz
        )
        om = dz / d
        b1 = w.beta(-0.5 * np.ones(n), u1)
        b2 = w.beta(-0.5 * np.ones(n), u2)
        first = 0.5 * ((b1 - b2) * np.conj(om) ** 2).real
        err = np.abs(exact - first) * r**4
        print(
            f"   r={r:5d}: max {err.max():8.2f}   p99 {np.quantile(err, 0.99):7.2f}"
            f"   (max |shift| / lam^2 {np.abs(exact).max() * r**2:.2f})"
        )


def check_leading_shear():
    """Differences of beta(-1/2, u) against the leading form b_lead."""
    print("2. relative error of beta differences against b_lead, and times r")
    u = np.linspace(-0.35, 0.35, 301)
    for r in SIZES:
        w = wg.Wave(Q, r, PHI)
        beta = w.beta(-0.5 * np.ones_like(u), u)
        exact = beta[1:] - beta[0]
        lead = b_lead(w, u[1:]) - b_lead(w, u[0])
        rel = np.max(np.abs(exact - lead)) / np.max(np.abs(lead))
        print(f"   r={r:5d}: {rel:.3e}   x r = {rel * r:.3f}")


def ugrid(a, b, r, per=12):
    """Gauss panels with about `per` nodes per wavelength, sin^2 end taper."""
    panels = int(np.ceil((b - a) * r / (2 * np.pi))) + 2
    x, wt = leggauss(per)
    edges = np.linspace(a, b, panels + 1)
    width = edges[1:] - edges[:-1]
    nodes = edges[:-1, None] + width[:, None] * (x + 1) / 2
    weights = width[:, None] * wt / 2 * np.sin(np.pi * (nodes - a) / (b - a)) ** 2
    return nodes.ravel(), weights.ravel()


def cone_rows(w, up, vp, zp, sigma, ux, wu, om, wo):
    """Weighted rows of the six-parameter fit for partners on P's exact cone.

    Unknowns: dv, dz_re, dz_im, dm of P, then Re s, Im s of its shear.
    Partners lie on the exact wave cone; the right-hand side is minus their
    background margin.
    """
    u = np.repeat(ux, om.size)
    weight = np.repeat(wu, om.size) * np.tile(wo, ux.size)
    o = np.tile(om, ux.size)
    lo, hi = (up * np.ones_like(u), u) if sigma > 0 else (u, up * np.ones_like(u))
    d = w.dm(lo, hi)
    a, b = w.inv_action(lo, hi)
    half = 0.5 * wg.sqf(a, b, o * d)
    if sigma > 0:
        mb = bg_margin(w, up, vp, zp, u, vp + half, zp + o * d)
    else:
        mb = bg_margin(w, u, vp - half, zp - o * d, up, vp, zp)
    design = sigma * np.stack(
        [
            -np.ones_like(u),
            o.real,
            o.imag,
            -0.5 * np.abs(o) ** 2,
            0.5 * (o.real**2 - o.imag**2),
            o.real * o.imag,
        ],
        axis=1,
    )
    root = np.sqrt(weight)
    return design * root[:, None], -mb * root


def check_cone_fit():
    """Six-parameter fit of one event's exact cone: fitted shear / b_lead(u_P)."""
    angles = 2 * np.pi * np.arange(16) / 16
    x, mw = leggauss(6)
    om = ((0.6 + 0.4 * x)[:, None] * np.exp(1j * angles)[None, :]).ravel()
    wo = np.repeat(mw, angles.size)
    events = ((-0.05, 0.03 + 0.01j), (0.08, -0.02 + 0.04j), (-0.12, 0.0j))
    print("3. fitted shear / b_lead(u_P) at three events; worst relative error")
    for r in SIZES:
        w = wg.Wave(Q, r, PHI)
        ratios, adjust = [], 0.0
        for up, zp in events:
            blocks = []
            for sigma in (+1, -1):
                a, b = (up + 0.08, up + 0.30) if sigma > 0 else (up - 0.30, up - 0.08)
                ux, wu = ugrid(a, b, r)
                blocks.append(cone_rows(w, up, 0.3, zp, sigma, ux, wu, om, wo))
            design = np.concatenate([blk[0] for blk in blocks])
            rhs = np.concatenate([blk[1] for blk in blocks])
            sol, *_ = np.linalg.lstsq(design, rhs, rcond=None)
            ratios.append((sol[4] + 1j * sol[5]) / b_lead(w, up))
            adjust = max(adjust, np.abs(sol[:4]).max() * r**2)
        worst = max(abs(z - 1) for z in ratios)
        shown = ", ".join(f"{z.real:.4f}{z.imag:+.4f}i" for z in ratios)
        print(
            f"   r={r:5d}: {shown}; worst {worst:.2e} (x r = {worst * r:.3f});"
            f" max |position adjustment| / lam^2 {adjust:.1e}"
        )


def check_sliver():
    """Minkowski identity behind the sector locality of critical pairs."""
    rng = np.random.default_rng(11)
    n = 1000
    big, small = rng.uniform(0.05, 1, n), rng.uniform(0.05, 1, n)
    t, eps = rng.uniform(0, 0.01, n), rng.uniform(0, 0.01, n)
    w = rng.normal(0, 1, n) + 1j * rng.normal(0, 1, n)
    u = rng.normal(0, 1, n) + 1j * rng.normal(0, 1, n)
    zm, zz, zv = -small, -small * u, -small * np.abs(u) ** 2 / 2 - t
    xm, xz, xv = big, big * w, big * np.abs(w) ** 2 / 2 - eps
    lhs = xv - zv - np.abs(xz - zz) ** 2 / (2 * (xm - zm))
    rhs = t - eps + big * small * np.abs(w - u) ** 2 / (2 * (big + small))
    print(f"4. sliver identity, max abs difference {np.max(np.abs(lhs - rhs)):.2e}")


if __name__ == "__main__":
    check_first_order()
    check_leading_shear()
    check_cone_fit()
    check_sliver()
