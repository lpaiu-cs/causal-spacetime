# ruff: noqa: E501, UP031
"""Check 01: the matrix primitive and the claimed expansion of B_lambda(u).

Claim (gram_fourier_challenge.md, Section 3):
  B_lambda(u) = -2 i sqrt(q) lam^2 f(u) e^{i(u/lam+phi)}
                + smooth endpoint term + O(lam^3),   f = a^-2,
in the complex trace-free representation, with Mt(s,t) = dm I - [B(t)-B(s)].

Deterministic.  ASCII only.  Exit code 1 on any FAIL.
"""

from __future__ import annotations

import math

import numpy as np
from wavegeom import Reporter, Wave, bump, gauss01, sapply, sqf

rep = Reporter()
Q, PHI = 1.0, 0.37

# ---------------------------------------------------------------- 1. primitive
try:
    import mpmath as mp

    mp.mp.dps = 30
    have_mp = True
except Exception:  # venv python has no mpmath
    have_mp = False


def reference_primitive(r, u):
    """Independent value of int_{-1/2}^u sec^2(k x) e^{i r x} dx."""
    k = math.sqrt(Q)
    if have_mp:
        pieces = max(8, int((u + 0.5) * r / 2) + 1)
        knots = [
            mp.mpf(-0.5) + (mp.mpf(u) + 0.5) * j / pieces for j in range(pieces + 1)
        ]
        val = mp.quad(lambda x: mp.exp(1j * r * x) / mp.cos(k * x) ** 2, knots)
        return complex(val)
    # fallback: a different rule (composite 6-point Gauss on 40 r panels)
    x, w = gauss01(6)
    n = int(40 * r)
    edges = np.linspace(-0.5, u, n + 1)
    nodes = edges[:-1, None] + (edges[1] - edges[0]) * x
    vals = np.exp(1j * r * nodes) / np.cos(k * nodes) ** 2
    return complex(np.sum(vals @ w) * (edges[1] - edges[0]))


worst = 0.0
for r in (64, 256):
    wave = Wave(Q, r, PHI)
    for u in (-0.3, 0.2, 0.5):
        worst = max(worst, abs(complex(wave.prim(u)) - reference_primitive(r, u)))
rep.check(
    "primitive table agrees with %s quadrature" % ("mpmath" if have_mp else "fallback"),
    worst < 1e-13,
    f"max abs err {worst:.2e}",
)

# ------------------------------------------ 2. complex representation helpers
rng = np.random.default_rng(20261124)
worst = 0.0
for _ in range(200):
    a, b1, b2, y1, y2 = rng.normal(size=5)
    mat = np.array([[a + b1, b2], [b2, a - b1]])
    y = np.array([y1, y2])
    zeta, b = y1 + 1j * y2, b1 + 1j * b2
    out = mat @ y
    worst = max(worst, abs(sapply(a, b, zeta) - (out[0] + 1j * out[1])))
    worst = max(worst, abs(sqf(a, b, zeta) - y @ mat @ y))
rep.check(
    "complex representation of symmetric 2x2 matrices",
    worst < 1e-12,
    f"max err {worst:.2e}",
)

# --------------------------- 3. sign convention of the wave action matrix
# gamma = a^2 R diag(e^{2b}, e^{-2b}) R^T, theta = (u/lam+phi)/2, sinh(2b)=2 lam k.
# Integrate gamma^-1 entrywise and compare with c (dm I - B).
r = 128
wave = Wave(Q, r, PHI)
lam, k = 1.0 / r, math.sqrt(Q)
bb = 0.5 * math.asinh(2 * lam * k)
s0, t0 = -0.27, 0.19
x, w = gauss01(12)
edges = np.linspace(s0, t0, 2000)
nodes = (edges[:-1, None] + (edges[1] - edges[0]) * x).ravel()
wts = np.tile(w * (edges[1] - edges[0]), len(edges) - 1)
th = 0.5 * (nodes / lam + PHI)
ct, st = np.cos(th), np.sin(th)
a2 = np.cos(k * nodes) ** 2
g11 = a2 * (math.exp(2 * bb) * ct**2 + math.exp(-2 * bb) * st**2)
g22 = a2 * (math.exp(2 * bb) * st**2 + math.exp(-2 * bb) * ct**2)
g12 = a2 * (math.exp(2 * bb) - math.exp(-2 * bb)) * ct * st
det = g11 * g22 - g12**2
minv = np.array(
    [
        [np.sum(wts * g22 / det), -np.sum(wts * g12 / det)],
        [-np.sum(wts * g12 / det), np.sum(wts * g11 / det)],
    ]
)
beta = complex(wave.beta(s0, t0))
dm = float(wave.dm(s0, t0))
mine = wave.cfac * np.array(
    [[dm - beta.real, -beta.imag], [-beta.imag, dm + beta.real]]
)
err = float(np.max(np.abs(minv - mine)))
rep.check(
    "int gamma^-1 = c (dm I - B) with B <-> beta (sign and helicity)",
    err < 1e-12,
    f"max err {err:.2e}",
)

# -------------------------------------------------- 4. the expansion of B
print("     r        u   |B-lead|/lam^3  bound3   |B-lead-next|/lam^4  bound4")
ok3 = ok4 = True
worst_ratio3 = 0.0
a_end = -0.5
for r in (64, 128, 256, 512, 1024, 2048):
    wave = Wave(Q, r, PHI)
    lam = 1.0 / r
    for u in (-0.2, 0.0, 0.2, 0.35):
        b = complex(wave.beta(a_end, u))
        eu = np.exp(1j * (u / lam + PHI))
        ea = np.exp(1j * (a_end / lam + PHI))
        lead = -2j * k * lam**2 * (wave.f(u) * eu - wave.f(a_end) * ea)
        nxt = 2 * k * lam**3 * (wave.fprime(u) * eu - wave.fprime(a_end) * ea)
        f2 = lambda z: 2 * k * k * (1 + 3 * math.tan(k * z) ** 2) / math.cos(k * z) ** 2  # noqa: E731
        bound3 = 1.1 * 2 * k * (abs(wave.fprime(u)) + abs(wave.fprime(a_end))) + 0.05
        bound4 = (
            1.1
            * (2 * k * (f2(u) + f2(a_end)) + 4 * Q * k * (wave.f(u) + wave.f(a_end)))
            + 0.05
        )
        r3 = abs(b - lead) / lam**3
        r4 = abs(b - lead - nxt) / lam**4
        ok3 &= r3 <= bound3
        ok4 &= r4 <= bound4
        worst_ratio3 = max(worst_ratio3, r3)
        if u in (0.2,):
            print(
                "%6d %8.2f %14.5f %8.3f %18.5f %8.3f" % (r, u, r3, bound3, r4, bound4)
            )
rep.check(
    "B = lead + O(lam^3): |B - lead| <= C lam^3 at all 24 (r,u)",
    ok3,
    f"max ratio {worst_ratio3:.4f}",
)
rep.check("B = lead + next + O(lam^4) with next = 2 sqrt(q) lam^3 [f' e^{i psi}]", ok4)

# ---------------------- 5. coefficient and phase of the local rotating term
x, w = gauss01(8)
prev = None
ok = True
for r in (128, 512, 2048):
    wave = Wave(Q, r, PHI)
    lam = 1.0 / r
    edges = np.linspace(0.1, 0.3, int(0.2 * r * 2) + 40)
    nodes = (edges[:-1, None] + (edges[1] - edges[0]) * x).ravel()
    wts = np.tile(w * (edges[1] - edges[0]), len(edges) - 1) * bump((nodes - 0.2) / 0.1)
    b = wave.beta(-0.5 * np.ones_like(nodes), nodes)
    coef = np.sum(wts * b * np.exp(-1j * (nodes / lam + PHI))) / (
        lam**2 * np.sum(wts * wave.f(nodes))
    )
    target = -2j * k
    print(
        "     r=%5d  local coefficient %.6f%+.6fi  (claimed %.6f%+.6fi)"
        % (r, coef.real, coef.imag, target.real, target.imag)
    )
    ok &= abs(coef - target) < 6 * lam + 2e-3
rep.check(
    "windowed local coefficient of B/lam^2 at e^{i(u/lam+phi)} is -2 i sqrt(q)", ok
)

rep.finish()
