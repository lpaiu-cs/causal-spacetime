# Schwarzschild diamond-volume oracle: analytic reductions and open items

Status note (2026-08-10; Section 4 updated 2026-08-12 after the O3
campaign). This records what is DERIVED versus what remains open on
the path to a bounded-error diamond-volume oracle for the
Schwarzschild exterior — the prediction-anchored claim class that
Paper A's Section 6.7 deliberately does not use. Sections 1-3 are
mathematical statements about the continuum setup; the numerical
certification that closed Section 4's items lives in
`schwarzschild_volume_oracle_certification.md`, and Section 4 now
records its outcome.

## 1. [VERIFIED analytic] Staticity removes the time integral (4D -> 3D)

In the static exterior, an event x lies in J+(p) iff
t_x - t_p >= T(P, X), where T(P, X) is the minimal coordinate flight
time between the SPATIAL positions P and X — a function of space
alone. For the Alexandrov diamond of p = (t_p, P), q = (t_q, Q) with
Delta t = t_q - t_p, the t-slice at spatial point X is the interval
[t_p + T(P, X), t_q - T(X, Q)], of length [Delta t - T(P,X) - T(X,Q)]_+.
Hence, with the mass-independent volume element sqrt(-g) = r^2 sin(theta),

  V(p, q) = INT r^2 sin(theta) [ Delta t - T(P, X) - T(X, Q) ]_+ d^3X .

This holds for ANY anchor pair in the static exterior; the time
dimension never has to be discretized.

## 2. [VERIFIED analytic] Radially aligned anchors reduce to 2D

If P and Q lie on one radial ray, the whole configuration is invariant
under rotations about that ray. Taking it as the polar axis, T(P, X)
and T(X, Q) depend only on (r, theta), the phi integral contributes
2*pi, and

  V = 2*pi INT INT r^2 sin(theta) [ Delta t - T1(r, theta)
        - T2(r, theta) ]_+ dr d(theta) .

For general anchor placements only a reflection symmetry survives and
the integral stays 3D (still with the closed-form time factor of
Section 1) — which is why oracle anchors are chosen radially aligned.

## 3. Optical-metric structure of T and the Lipschitz skeleton

For a static metric ds^2 = -f dt^2 + h_ij dx^i dx^j, Fermat's
principle makes the spatial projections of null geodesics the
geodesics of the optical metric h_ij / f, and the minimal coordinate
flight time IS the optical distance. For Schwarzschild,

  dl_opt^2 = dr^2 / f^2 + (r^2 / f) (d(theta)^2 + sin^2(theta) d(phi)^2),
  f = 1 - 2M/r .

Consistency check: the radial optical distance is INT dr / f = the
tortoise-coordinate difference, exactly the S1 solver's exact radial
branch. As a Riemannian distance function, T is 1-Lipschitz in the
optical metric, giving explicit coordinate bounds
|dT/dr| <= 1/f(r) and |dT/d(theta)| <= r / sqrt(f); together with the
1-Lipschitz [ . ]_+ these bound the integrand's oscillation on any
(dr, d(theta)) cell — the skeleton of a deterministic cell-refinement
enclosure with O(h) convergence.

## 4. Former open items — closed (2026-08-12)

- [CERTIFIED] Numerical flight-time enclosure (PR #64). The
  Gauss-Legendre successive-refinement err was discarded as
  certification grounds. The per-call enclosure is now MPFR
  directed-rounding interval arithmetic end to end, with a certified
  composite-midpoint remainder (interval-evaluated second-derivative
  bounds) — see the certification document, L1-L2.
- [VERIFIED] Anchor-diamond containment (PR #64). The frozen anchors
  (12, 18, 8.5)M carry a containment certificate proven uniformly, in
  the optical-distance form anticipated here: margins
  1.5600 / 3.3326 / 2.9111 / 1.8923, box r in [11.3536, 18.6950] M,
  psi <= 0.817913, perihelion floor 10.417 M (L4).
- [SIZED, then EXECUTED] Grid and cost (PR #65-#67 + campaign). The
  measured neighbor price ladder (`docs/prereg/p14_oracle_price.json`)
  replaced the earlier estimates, and the frozen campaign is the
  final price point: 137,958 predicate calls, 3.72 h wall clock.

The oracle result. Executed once from the clean exact freeze checkout
`785148e` under the frozen caps (600,000 calls / 24 h / depth 18),
the frozen configuration's certified volume is

  V = [56.212737, 57.348019]   (outward binary64 endpoints),
  (V_hi - V_lo) / (V_hi + V_lo) = 0.009997 <= 0.01,

terminated `target-met` (max depth reached 12, no uncosted cells, the
L6d intersection never active — so the mode decomposition is exact).
Artifact: `docs/prereg/p14_o3_volume.json` (write-once); executed
freeze surface: `docs/prereg/p14_o3_executed_freeze_manifest.json`.

## 5. Relation to Paper A

Section 6.7's Schwarzschild confirmation is operationally anchored and
does not use this oracle. The note exists so that the manuscript's
"partially derived theoretical path" phrasing has a public referent:
Sections 1-2 are derived; Section 3 is a standard-result application
with explicit constants; Section 4's former open boundary is now
closed by the certification document and the executed campaign. Any
prediction-anchored USE of the certified volume (sprinkle counts vs
rho*V) is a separate, not-yet-designed stage; nothing in Paper A is
upgraded by this note.
