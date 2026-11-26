# ruff: noqa: E501, UP031
"""Check 04: independent population Fourier coefficient of the Gram residual.

Claims under test (gram_fourier_challenge.md Section 3, boxed (4);
gram_fourier_noise.md Section 6, (11)):

  S(1/lam) / (lam^2 E w) -> -i sqrt(q) e^{i phi} int w_u D / int w_u a^2,
  |E_lam[w H e^{-i omega u}]| >= c lam^2 for |omega - 1/lam| <= 1/8,
  for anchors in an intrinsic neighbourhood, any rotation or reflection.

Everything here is computed with the referee's own code (popsignal.py,
wavegeom.py); the author's source is not imported or re-run.  The stored
rows of the author (gram_fourier_evidence/rows.jsonl) are quoted as
constants and used only as a comparison target.  Weight and phase use
RECONSTRUCTED coordinates.

Tolerances: 1e-3 for reproducing the stored values (the author's v and
radius quadrature is coarser); 1% at r = 512 for agreement with the
limit, which the author's frozen row already showed (0.53%); these were
written down after an exploratory run of the referee code.

Deterministic.  ASCII only.  Exit code 1 on any FAIL.  About 75 s.
"""

from __future__ import annotations

import math

import numpy as np
from popsignal import (
    U_HI,
    U_LO,
    Probes,
    anchors_from_gram,
    composite,
    d_anchor,
    frozen_anchors,
    prediction,
)
from wavegeom import Reporter, bump

rep = Reporter()
Q, PHI = 1.0, 0.37
AUTHOR = {  # scaled_complex_signal rows of the author's frozen run
    (128, 0): complex(-0.022763969584423623, -0.06884304319051937),
    (256, 0): complex(0.015247779515825335, 0.005876091725943968),
    (512, 0): complex(0.017393017760465323, 0.00821974778981902),
    (512, 1): complex(0.018164225996713287, -0.006426814508791929),
}
CERT_GRID_LOWER = 0.010168090397347063  # anchor_guard_certificate.json


def rel(a, b):
    return abs(a - b) / abs(b)


def fmt(z):
    return f"{z.real:.8f}{z.imag:+.8f}i"


probes = {r: Probes(Q, r, PHI) for r in (128, 256, 512)}
frozen = frozen_anchors(0.0)
vals = {}
for r in (128, 256, 512):
    out, diag = probes[r].evaluate(frozen, [float(r)])
    vals[r] = out[0] * r * r
    print(
        "r=%d theta=0: S/(lam^2 Ew) = %s, limit %s, ell = %s"
        % (
            r,
            fmt(vals[r]),
            fmt(prediction(probes[r].wave, frozen)),
            np.array2string(diag["ell"], precision=8),
        )
    )

pred0 = prediction(probes[512].wave, frozen)
rep.check(
    "r=128 theta=0 reproduces the stored value of the author (a row 4.5x away from the limit)",
    rel(vals[128], AUTHOR[(128, 0)]) < 1e-3 and rel(vals[128], pred0) > 4,
    f"rel diff {rel(vals[128], AUTHOR[(128, 0)]):.2e}, distance from limit {rel(vals[128], pred0):.3f}",
)
rep.check(
    "r=256 theta=0 reproduces the stored value of the author",
    rel(vals[256], AUTHOR[(256, 0)]) < 1e-3,
    f"rel diff {rel(vals[256], AUTHOR[(256, 0)]):.2e}",
)
rep.check(
    "r=512 theta=0 reproduces the stored value of the author",
    rel(vals[512], AUTHOR[(512, 0)]) < 1e-3,
    f"rel diff {rel(vals[512], AUTHOR[(512, 0)]):.2e}",
)
rep.check(
    "r=512 theta=0 agrees with boxed limit (4): modulus, phase and sign",
    rel(vals[512], pred0) < 0.01,
    f"rel err {rel(vals[512], pred0):.5f}",
)
rep.check(
    "error against the limit falls by more than 10x from r=256 to r=512",
    rel(vals[256], pred0) > 10 * rel(vals[512], pred0),
    f"{rel(vals[256], pred0):.5f} -> {rel(vals[512], pred0):.5f}",
)

# rotation of the anchor frame: phase phi - 2 theta, same modulus
rot = frozen_anchors(math.pi / 8)
out, _ = probes[512].evaluate(rot, [512.0])
v_rot = out[0] * 512**2
p_rot = prediction(probes[512].wave, rot)
print(f"r=512 theta=pi/8: {fmt(v_rot)}, limit {fmt(p_rot)}")
rep.check(
    "r=512 theta=pi/8 reproduces the stored value of the author",
    rel(v_rot, AUTHOR[(512, 1)]) < 1e-3,
    f"rel diff {rel(v_rot, AUTHOR[(512, 1)]):.2e}",
)
rep.check(
    "rotation by theta multiplies the limit by exp(-2 i theta)",
    abs(p_rot / pred0 - np.exp(-2j * math.pi / 8)) < 1e-12 and rel(v_rot, p_rot) < 0.01,
    f"rel err {rel(v_rot, p_rot):.5f}",
)

# reflection: D -> conj(D); modulus preserved
ref = frozen_anchors(0.0, reflect=True)
out, _ = probes[512].evaluate(ref, [512.0])
v_ref = out[0] * 512**2
p_ref = prediction(probes[512].wave, ref)
print(f"r=512 reflected anchors: {fmt(v_ref)}, limit {fmt(p_ref)}")
rep.check(
    "reflected anchors: coefficient matches the limit with D replaced by conj(D)",
    rel(v_ref, p_ref) < 0.01,
    f"rel err {rel(v_ref, p_ref):.5f}",
)
rep.check(
    "reflection preserves the modulus",
    abs(abs(p_ref) - abs(pred0)) < 1e-12 * abs(pred0)
    and abs(abs(v_ref) / abs(vals[512]) - 1) < 0.01,
    "|S_ref|/|S| = %.5f" % (abs(v_ref) / abs(vals[512])),
)

# frequency offset |omega - 1/lam| <= 1/8 and the certified lower bound
u, w = composite(U_LO, U_HI, 200, 6)
w = w * bump((u - 0.2) / 0.1)
floor = math.sqrt(Q) * CERT_GRID_LOWER * np.sum(w) / np.sum(w * probes[512].wave.a2(u))
offs = (-0.125, 0.125)
out, _ = probes[512].evaluate(frozen, [512.0 + o for o in offs])
ok_pred, ok_floor, worst = True, True, 0.0
for o, val in zip(offs, out, strict=False):
    sc = val * 512**2
    pr = prediction(probes[512].wave, frozen, offset=o)
    worst = max(worst, rel(sc, pr))
    ok_pred &= rel(sc, pr) < 0.015
    ok_floor &= abs(sc) > floor
    print(
        f"r=512 omega=r{o:+.3f}: {fmt(sc)}, offset limit {fmt(pr)}, modulus {abs(sc):.6f}, certified floor {floor:.6f}"
    )
rep.check(
    "frequency offsets +-1/8: coefficient matches the offset limit",
    ok_pred,
    f"max rel err {worst:.5f}",
)
rep.check("(11): modulus stays above the certificate floor at offsets +-1/8", ok_floor)

# a worst-case corner of the outer neighbourhood, rotated and reflected
tol_g, tol_uv = 2e-5, 2e-3
du, dv = (-tol_uv, tol_uv, -tol_uv), (tol_uv, -tol_uv, tol_uv)
best = None
for signs in np.ndindex(2, 2, 2, 2, 2):
    s = [2 * x - 1 for x in signs]
    args = (
        0.0064 + s[0] * tol_g,
        0.0064 + s[1] * tol_g,
        s[2] * tol_g,
        0.0064 + s[3] * tol_g,
        0.0064 + s[4] * tol_g,
    )
    anc, g33 = anchors_from_gram(*args, du=du, dv=dv)
    if abs(g33 - 0.0128) > tol_g:
        continue
    im = float(np.min(d_anchor(probes[512].wave, anc, u).imag))
    if best is None or im < best[0]:
        best = (im, args)
print("worst admissible Gram corner: min Im D_A = {:.6f} at Gram args {}".format(*best))
corner, _ = anchors_from_gram(*best[1], du=du, dv=dv, angle=0.7, reflect=True)
out, diag = probes[512].evaluate(corner, [512.0, 512.125])
p_c = [prediction(probes[512].wave, corner, offset=o) for o in (0.0, 0.125)]
sc = out * 512**2
print(
    "corner anchors (rotated 0.7, reflected): {}, limit {}; min probe slack {:.4f}".format(
        fmt(sc[0]), fmt(p_c[0]), diag["min_slack"]
    )
)
rep.check(
    "neighbourhood corner, rotated and reflected: matches the general-anchor limit D_A",
    max(rel(sc[0], p_c[0]), rel(sc[1], p_c[1])) < 0.015,
    f"rel errs {rel(sc[0], p_c[0]):.5f} {rel(sc[1], p_c[1]):.5f}",
)
rep.check(
    "(11) at the neighbourhood corner: modulus above the certificate floor",
    min(abs(sc)) > floor,
    f"min modulus {min(abs(sc)):.6f}, floor {floor:.6f}",
)

# reconstruction: how much do the O(lam^4 log) calibration errors move the result?
geo = {}
for r in (256, 512):
    out, _ = probes[r].evaluate(frozen, [float(r)], reconstructed=False)
    geo[r] = out[0] * r * r
    dmax = float(np.max(np.abs(probes[r].u_hat - probes[r].u)))
    print(
        "r=%d: canonical-coordinate value %s (rel err to limit %.5f); calibration shift |full-canonical|/|limit| = %.5f; "
        "max|u_hat-u| = %.3e = %.0f lam^4 log(1/lam); max phase error %.2e rad"
        % (
            r,
            fmt(geo[r]),
            rel(geo[r], pred0),
            abs(vals[r] - geo[r]) / abs(pred0),
            dmax,
            dmax * r**4 / math.log(r),
            r * dmax,
        )
    )
rep.check(
    "canonical-coordinate (calibration-free) value agrees with the limit at r=512 to 0.2%",
    rel(geo[512], pred0) < 2e-3,
    f"rel err {rel(geo[512], pred0):.5f}",
)
rep.check(
    "calibration shift of the coefficient decays faster than lam^2 (more than 16x per halving of lam)",
    abs(vals[256] - geo[256]) > 16 * abs(vals[512] - geo[512]),
    f"{abs(vals[256] - geo[256]) / abs(pred0):.5f} -> {abs(vals[512] - geo[512]) / abs(pred0):.5f} of the limit",
)
# other parameters, not run by the author for this challenge: q = 1/4 and another phase
for q2, phi2, tol_full in ((0.25, PHI, 0.02), (1.0, 2.1, 0.01)):
    other = Probes(q2, 512, phi2)
    out, _ = other.evaluate(frozen, [512.0])
    can, _ = other.evaluate(frozen, [512.0], reconstructed=False)
    sc, sg = out[0] * 512**2, can[0] * 512**2
    pr = prediction(other.wave, frozen)
    print(
        f"q={q2:.2f} phi={phi2:.2f} r=512: {fmt(sc)}, canonical {fmt(sg)}, limit {fmt(pr)}"
    )
    rep.check(
        f"q={q2:.2f} phi={phi2:.2f} r=512: limit (4) holds (sqrt(q) and e^{{i phi}} dependence)",
        rel(sc, pr) < tol_full and rel(sg, pr) < 2e-3,
        f"rel err reconstructed {rel(sc, pr):.5f}, canonical {rel(sg, pr):.5f}",
    )

du256 = float(np.max(np.abs(probes[256].u_hat - probes[256].u)))
du512 = float(np.max(np.abs(probes[512].u_hat - probes[512].u)))
rep.check(
    "max |u_hat - u| scales like lam^4 log(1/lam) (ratio between 10 and 16)",
    10 < du256 / du512 < 16,
    "ratio %.2f" % (du256 / du512),
)
rep.check(
    "phase error omega |u_hat - u| is negligible at r=512 (< 1e-3 rad)",
    512 * du512 < 1e-3,
    "%.2e rad" % (512 * du512),
)
rep.finish()
