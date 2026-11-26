# ruff: noqa: E501, UP031
"""Check 05: the intrinsic anchor neighbourhood of gram_fourier_noise.md
Section 6 and the rational guard certificate.

  1. brute-force scan of the outer neighbourhood (u_i, v_i within 0.002,
     Gram entries within 0.00002, 0 < q <= 1, 0.1 <= u <= 0.3): the
     sampled extremes of Im D_A and |D_A| must respect the certificate
     bounds (they are bounds, so sampled values must lie inside them);
  2. the text's step "|a_2x a_2y| = |G12| <= 0.00002" checked literally;
  3. the neighbourhood is stated through (u_i, v_i, R_i^2, g_ij) only:
     can a NON-timelike anchor pair (b_ij = 0) pass the Gram window?
     Exact rational bound valid for every q >= 0;
  4. background slack lower bounds against the certificate;
  5. INFO: order of magnitude of the physical probability of the
     quarter-tolerance neighbourhood (how slowly "probability tending to
     one" is reached with O(log n) candidates).

Deterministic (fixed seed).  ASCII only.  Exit code 1 on any FAIL.
"""

from __future__ import annotations

import math
from fractions import Fraction as F

import numpy as np
from wavegeom import Background, Reporter, gauss01

rep = Reporter()
CERT = {  # gram_fourier_evidence/anchor_guard_certificate.json
    "D_imaginary_lower": 0.01069797423346495,
    "D_absolute_upper": 0.042390706889430904,
    "grid_coefficient_lower": 0.010168090397347063,
    "anchor_pair_slack_lower": 0.028916666666666667,
    "active_probe_slack_lower": 0.1001972267437538,
}
TOL_G, TOL_UV = 2e-5, 2e-3
BASE_U = np.array([-0.35, -0.25, -0.15])
BASE_V = np.array([0.15, 0.25, 0.35])
TARGET = np.array([0.0064, 0.0064, 0.0, 0.0064, 0.0064])  # G11 G22 G12 G13 G23
rng = np.random.default_rng(20261124)


def sample(n):
    """Random points and all corners of the Gram box, both tolerances."""
    corners = np.array(
        [[2 * b - 1 for b in idx] for idx in np.ndindex(2, 2, 2, 2, 2)], dtype=float
    )
    sg = np.vstack((corners, rng.uniform(-1, 1, size=(n, 5))))
    suv = np.sign(rng.uniform(-1, 1, size=(sg.shape[0], 6)))  # extremes of u_i, v_i
    return (
        TARGET + TOL_G * sg,
        BASE_U + TOL_UV * suv[:, :3],
        BASE_V + TOL_UV * suv[:, 3:],
    )


def vectors(g):
    g11, g22, g12, g13, g23 = g.T
    a1 = np.sqrt(g11) + 0j
    a2x = g12 / a1.real
    a2 = a2x + 1j * np.sqrt(g22 - a2x**2)
    a3x = g13 / a1.real
    a3 = a3x + 1j * (g23 - a2x * a3x) / a2.imag
    return np.stack((a1, a2, a3), axis=1)


g, ua, va = sample(200000)
z = vectors(g)
keep = np.abs(np.abs(z[:, 2]) ** 2 - 0.0128) <= TOL_G
g, ua, va, z = g[keep], ua[keep], va[keep], z[keep]
print("admissible sampled anchor triples: %d" % len(g))
# ell from sum ell_i a_i = 0, ell_3 = 1 (2x2 real solve per sample)
det = z[:, 0].real * z[:, 1].imag - z[:, 1].real * z[:, 0].imag
ell1 = -(z[:, 2].real * z[:, 1].imag - z[:, 1].real * z[:, 2].imag) / det
ell2 = -(z[:, 0].real * z[:, 2].imag - z[:, 2].real * z[:, 0].imag) / det
ell = np.stack((ell1, ell2, np.ones_like(ell1)), axis=1)
resid = np.max(np.abs(np.sum(ell * z, axis=1)))
rep.check(
    "sum ell_i a_i = 0 on every sampled triple",
    resid < 1e-15,
    f"max residual {resid:.1e}",
)

uu = np.linspace(0.1, 0.3, 21)
im_min, abs_max, grid_min = np.inf, 0.0, np.inf
for q in (1e-4, 0.25, 0.5, 1.0):
    bg = Background(q)
    gaps = bg.m(uu)[None, :, None] - bg.m(ua)[:, None, :]
    d = -0.5 * np.sum(ell[:, None, :] * np.conj(z)[:, None, :] ** 2 / gaps, axis=2)
    im_min = min(im_min, float(d.imag.min()))
    abs_max = max(abs_max, float(np.abs(d).max()))
    grid_min = min(
        grid_min, float((d.imag.min(axis=1) - np.abs(d).max(axis=1) / 80).min())
    )
print(
    "sampled min Im D_A = {:.6f} (certificate lower bound {:.6f})".format(
        im_min, CERT["D_imaginary_lower"]
    )
)
print(
    "sampled max |D_A|  = {:.6f} (certificate upper bound {:.6f})".format(
        abs_max, CERT["D_absolute_upper"]
    )
)
rep.check(
    "certificate lower bound on Im D_A holds on the scan",
    im_min >= CERT["D_imaginary_lower"],
)
rep.check(
    "certificate upper bound on |D_A| holds on the scan",
    abs_max <= CERT["D_absolute_upper"],
)
rep.check(
    "frequency-offset coefficient bound holds on the scan",
    grid_min >= CERT["grid_coefficient_lower"],
    f"sampled min {grid_min:.6f}",
)

# 2. the literal step |a_2x a_2y| = |G12| <= 0.00002
prod = float(np.max(np.abs(z[:, 1].real * z[:, 1].imag)))
exact = TOL_G * math.sqrt((0.0064 + TOL_G) / (0.0064 - TOL_G))
print(
    f"max |a_2x a_2y| on the scan = {prod:.6e}; exact supremum {exact:.6e}; text claims <= {TOL_G:.1e}"
)
print(
    f"NOTE (MINOR) the identity |a_2x a_2y| = |G12| is not exact: the supremum exceeds 0.00002 by a factor {exact / TOL_G:.4f}; "
    f"effect on the Im D_A bound {1.02 * (exact - TOL_G) / 0.348:.1e}"
)
rep.check(
    "referee finding reproduced: sup |a_2x a_2y| exceeds 0.00002, by less than 0.4%",
    TOL_G < prod <= exact * (1 + 1e-12) and exact / TOL_G < 1.004,
)
rep.check(
    "corrected step: Im D_A lower bound still clears 0.01 with |a_2x a_2y| <= 2.0063e-5",
    CERT["D_imaginary_lower"] - 1.02 * (exact - TOL_G) / 0.348 > 0.01,
)

# 3. non-timelike pairs: b = 0 gives g = (R_i^2 + R_j^2)/2 - dm dv, dm >= du for q >= 0
r2 = [F(642, 100000), F(642, 100000), F(1282, 100000)]
target = {(0, 1): F(0), (0, 2): F(64, 10000), (1, 2): F(64, 10000)}
ok = True
for (i, j), t in target.items():
    du = F(j - i, 10) - F(4, 1000)
    dv = F(j - i, 10) - F(4, 1000)
    upper = (r2[i] + r2[j]) / 2 - du * dv
    ok &= upper < t - F(2, 100000)
    print(
        "pair (%d,%d): spacelike g_ij <= %.6f, window starts at %.6f"
        % (i + 1, j + 1, float(upper), float(t - F(2, 100000)))
    )
rep.check(
    "no non-timelike anchor pair can pass the Gram window (exact rational, all q >= 0)",
    ok,
)

# 4. background slacks on the scan
bg = Background(1.0)
ma = bg.m(ua)
pair = np.inf
for i, j in ((0, 1), (0, 2), (1, 2)):
    s = (
        va[:, j]
        - va[:, i]
        - np.abs(z[:, j] - z[:, i]) ** 2 / (2 * (ma[:, j] - ma[:, i]))
    )
    pair = min(pair, float(s.min()))
probe = np.inf
for i in range(3):
    s = 0.65 - va[:, i] - (0.2 + np.abs(z[:, i])) ** 2 / (2 * (bg.m(0.1) - ma[:, i]))
    probe = min(probe, float(s.min()))
print(
    "sampled min anchor-pair slack {:.5f} (certificate {:.5f}); min probe slack {:.5f} (certificate {:.5f})".format(
        pair, CERT["anchor_pair_slack_lower"], probe, CERT["active_probe_slack_lower"]
    )
)
rep.check(
    "certificate slack lower bounds hold on the scan (q = 1)",
    pair >= CERT["anchor_pair_slack_lower"]
    and probe >= CERT["active_probe_slack_lower"],
)

# 5. INFO: probability that a u-sorted random triple lies in the quarter-tolerance set
x, w = gauss01(64)
u = -0.5 + x
m = bg.m(u)
zbg = math.pi * float(np.sum(w * bg.a2(u) * m * (bg.M - m) / bg.M))
tg, tuv = TOL_G / 4, TOL_UV / 4
gq = TARGET + tg * rng.uniform(-1, 1, size=(400000, 5))
f33 = float(np.mean(np.abs(np.abs(vectors(gq)[:, 2]) ** 2 - 0.0128) <= tg))
vol_z = 2 * math.pi * (2 * tg) ** 5 / (2 * 0.0064 * 0.0064) * f33
vol_uv = (2 * tuv) ** 6
dens = float(np.prod(bg.a2(BASE_U))) / zbg**3
p0 = 6 * dens * vol_z * vol_uv
print(
    f"INFO one-point normalisation Z = {zbg:.4f}; fraction of the Gram box compatible with rank two = {f33:.3f}"
)
print(
    f"INFO probability of the quarter-tolerance neighbourhood per candidate triple ~ {p0:.1e}"
)
print(
    f"INFO with ceil(log m) disjoint triples, Pr(no admissible triple) ~ m^(-{p0:.1e}): tends to zero, but only for log m >> {1 / p0:.1e}"
)
rep.check(
    "the small neighbourhood has positive probability (limit statement is valid)",
    p0 > 0,
)
rep.finish()
