# Schwarzschild diamond-volume oracle: analytic reductions and open items

Status note (2026-08-10). This records what is DERIVED versus what
remains open on the path to a bounded-error diamond-volume oracle for
the Schwarzschild exterior — the prediction-anchored claim class that
Paper A's Section 6.7 deliberately does not use. Verdicts here are
mathematical statements about the continuum setup; nothing below is a
numerical certification.

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

## 4. Open items (the reasons this is not yet an oracle)

- [TO CERTIFY] Numerical flight-time enclosure. The S1 solver's
  reported err is the difference between successive Gauss-Legendre
  refinements — a stopping heuristic, not a proven upper bound on the
  true quadrature error. Until the per-call bound is certified (e.g.
  by derivative-bounded remainder terms or validated/interval
  quadrature), the pointwise integrand interval is not rigorous and
  the assembled volume interval is not a bounded-error statement.
- [TO VERIFY] Anchor-diamond containment. The chosen anchors' diamond
  must be PROVEN to lie inside the solver's validated domain (the
  exterior shell and angular cap); a containment criterion in terms of
  the optical distance to the domain boundary is the natural form.
- [TO SIZE] Grid and cost. The cell count needed for a target volume
  error has not been derived; per-anchor cost figures quoted so far
  (1e4-1e5 predicate calls, ~15 s - 2.5 min) are estimates, not a
  price table.

## 5. Relation to Paper A

Section 6.7's Schwarzschild confirmation is operationally anchored and
does not use this oracle. The note exists so that the manuscript's
"partially derived theoretical path" phrasing has a public referent:
Sections 1-2 are derived; Section 3 is a standard-result application
with explicit constants; Section 4 is the open boundary.
