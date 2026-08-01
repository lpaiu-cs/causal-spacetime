"""P14 design check 2: the diamond volume is NOT flat, and by how much.

Review R1.2 (2026-08-01, in-session) refuted the first draft's control
C0 -- "the P12/P13 volume reading must return flat" -- by computing the
Alexandrov interval volume in the constant-A plane wave and finding it
0.4% to 7.3% ABOVE flat inside the conjugate-point bound. sqrt(-g) = 1
fixes the volume FORM; the diamond's boundary moves with the light
cones, and the interval volume moves with it. This script reproduces
that computation so the corrected Section 4.4 quotes a number its own
repository can regenerate.

Setup: A = w^2 constant, central-axis pair. On the axis H = 0, so the
pair's proper time is IDENTICAL in both arms: the comparison is at equal
clock reading by construction, and the volume still shifts -- the
deformation lives in the interval domain.

THE SIGN OF Dv, stated explicitly because getting it wrong picks
SPACELIKE anchors (review R7.2). With ds^2 = 2 du dv + ... and the
mostly-plus signature, an axis pair has

    tau^2 = -ds^2 = -2 Du Dv,

so tau = T at Du = T requires **Dv = -T/2**, negative. That is the same
convention the probe's predicate enforces -- future-related points have
DECREASING v -- and P2's fixed anchors must be built with it. The
quantity entering the volume below is the positive half-extent

    half_tau = -Dv = T/2,

and it appears squared, so the volume formula is unchanged; only a
caller reading this file for its anchor construction is affected, which
is exactly who P2 is.

Derivation (re-done from scratch for this check, then matched against
the reviewer's formula and the plane-wave world function of Harte &
Drivas, Phys. Rev. D 85, 124039):

  - transverse null-geodesic v-cost from (0, 0_perp) to (s, x_perp):
    (1/2) x_perp^T K(s) x_perp with K_x = w coth(w s) for x'' = +w^2 x
    (defocusing) and K_y = w cot(w s) for y'' = -w^2 y (focusing);
    flat limit K -> 1/s is the straight line.
  - a point at u-offset s must clear the cost from p and the cost to q,
    so the admissible v-extent at fixed (s, x_perp) is
    [half_tau*2 - q]_+ with q = (a_x x^2 + a_y y^2)/2 and
      a_x = w [coth(ws) + coth(w(T-s))]
      a_y = w [cot(ws)  + cot(w(T-s))]
  - the ellipse integral gives V = int_0^T pi (2*half_tau)^2/4 /
    sqrt(a_x a_y) ds, i.e. with tau = T and half_tau = T/2:
    V_A = (pi T^2 / 4) int_0^T ds / sqrt(a_x a_y).
  - flat limit, analytically: a -> T/(s(T-s)), integrand s(T-s)/T,
    integral T^2/6, V -> pi T^4 / 24 = the Alexandrov volume. That
    limit is asserted numerically below.

The two headline values are PINNED so this file and the design document
cannot drift apart (the C25/C28 lesson applied at draft stage). The
quartic-dominant growth of the ratio is asserted too: it is quadratic
in A, consistent with the small-diamond expansion's R and R_00 terms
(Roy, Sinha & Surya, arXiv:1212.0631) vanishing identically in vacuum.

Aborts on failure rather than recording it. numpy only.
"""

import numpy as np


def volume_ratio(wT: float, nodes: int = 6000) -> float:
    """V_A / V_0 for a central-axis pair at u-separation T, A = (wT/T)^2."""

    T, w = 1.0, wT
    xg, wg = np.polynomial.legendre.leggauss(nodes)
    s = 0.5 * T * (xg + 1.0)
    ww = 0.5 * T * wg
    a_x = w * (1.0 / np.tanh(w * s) + 1.0 / np.tanh(w * (T - s)))
    a_y = w * (np.cos(w * s) / np.sin(w * s)
               + np.cos(w * (T - s)) / np.sin(w * (T - s)))
    integral = float(np.sum(ww / np.sqrt(a_x * a_y)))
    return (np.pi * T ** 2 / 4.0 * integral) / (np.pi * T ** 4 / 24.0)


# the flat limit is recovered
flat = volume_ratio(1e-3)
assert abs(flat - 1.0) < 1e-9, f"flat limit broken: {flat}"

# the two values Section 4.4 quotes, both inside |Delta u| < pi/sqrt(A)
r1, r2 = volume_ratio(1.0), volume_ratio(2.0)
assert abs(r1 - 1.00400047) < 1e-7, f"wT=1: {r1:.8f}"
assert abs(r2 - 1.07300802) < 1e-7, f"wT=2: {r2:.8f}"

# monotone in wT, and quartic-dominant (quadratic in A): the ratio of
# excesses at wT=2 vs wT=1 sits at 16 for a pure quartic; higher terms
# push it above, and it must never fall below
excess_ratio = (r2 - 1.0) / (r1 - 1.0)
assert excess_ratio > 16.0, f"sub-quartic growth: {excess_ratio:.2f}"

# The leading coefficient, supplied by review R2 and confirmed here:
#     V_A / V_0 = 1 + (wT)^4 / 252 + O((wT)^8)
# This upgrades "quadratic in A" from a scaling observation to a closed
# number, so §4.4 quotes a coefficient rather than a trend. The window
# matters: above wT ~ 0.5 the (wT)^8 term is visible, and below wT ~ 0.1
# the ratio minus one falls to ~1e-9 where double-precision cancellation
# against 1 dominates -- the check is made where the physics, not the
# arithmetic, is the limit.
LEADING = 1.0 / 252.0
for wT in (0.5, 0.3, 0.2, 0.15):
    c = (volume_ratio(wT) - 1.0) / wT ** 4
    assert abs(c - LEADING) < 2.1e-6, f"wT={wT}: coefficient {c:.12f}"

# and the next term is (wT)^8, not (wT)^6: the residual over (wT)^4 is
# flat across that window, which it would not be at sixth order
resid = [((volume_ratio(wT) - 1.0) / wT ** 4 - LEADING) / wT ** 4
         for wT in (0.5, 0.3, 0.2, 0.15)]
assert max(resid) / min(resid) < 1.1, f"next term is not (wT)^8: {resid}"

print("flat-limit check (wT=1e-3):", f"{flat:.12f}")
print(f"wT=1: V_A/V_0 = {r1:.8f}   (pinned 1.00400047)")
print(f"wT=2: V_A/V_0 = {r2:.8f}   (pinned 1.07300802)")
print(f"excess ratio (quartic => 16): {excess_ratio:.2f}")
print(f"leading coefficient -> 1/252 = {LEADING:.12f}, "
      f"next term (wT)^8 with coefficient ~{np.mean(resid):.2e}")
print("=> the volume channel is NOT silent; C0 as drafted is dead: PASS")
