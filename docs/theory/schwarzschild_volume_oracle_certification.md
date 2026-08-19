# Schwarzschild diamond-volume oracle: certification (L0-L6)

Status note (2026-08-11, PR-O1). This document carries the proofs that
back the certified flight-time contract implemented in
`experiments/oracle/certified_flight_time.py` (PR-O1) and states the
cell-enclosure and composition contract (L6) that the volume
integrator (PR-O2) must implement. It is the certification referent
the design review required before `[TO CERTIFY]` in
`schwarzschild_volume_oracle_note.md` can be closed; the note's status
tags flip only when the assembled oracle lands (PR-O3).

Everything below is stated for the S1 patch: geometric units G = c = 1,
mass M (M = 1 in the frozen configuration), exterior shell
r in [R_MIN, R_MAX] = [10, 20] M, pairwise center angle psi <= 2.
Throughout, f = 1 - 2M/r, the tortoise coordinate is
rho(r) = r + 2M ln(r/2M - 1), and the optical radius is
w(r) = r / sqrt(f).

## 0. Floating-point contract (what "certified" means)

A certified quantity is a closed interval [lo, hi] of MPFR numbers
with the guarantee that the true real value lies inside. Every lower
endpoint is computed under MPFR rounding mode RoundDown and every
upper endpoint under RoundUp, at 96-bit precision, through gmpy2.
The contract rests on exactly two assumptions, both documented
library guarantees rather than platform folklore:

1. MPFR returns correctly rounded results in the requested direction
   for every operation used (add, sub, mul, div, sqrt, log, sin, cos,
   const_pi). This is MPFR's defining contract and covers the
   transcendental functions that platform libm does not.
2. Python floats and integers embed into MPFR exactly (binary64 is a
   subset of 96-bit MPFR; embedding is checked by unit test, not
   assumed).

No certified path uses platform libm, bare `math.nextafter` steps, or
fixed inflation factors. Fixed-padding and summation-order-permutation
checks may exist as auxiliary regression diagnostics only; they are
never certification grounds. Fail-closed rule: any operation whose
precondition cannot be certified -- division by an interval containing
zero, sqrt or log reaching a non-positive endpoint, trigonometric
evaluation outside a certified monotone domain, or an undecidable
interval comparison at a point where the algorithm must branch --
raises `CertificationError`; no certified result is produced past it.

## L0. Reduction: causality -> optical distance -> a 2D
Cartan-Hadamard plane

**L0a (static causality).** In the static exterior, a future-directed
causal curve satisfies f dt^2 >= h_ij dx^i dx^j, i.e. dt >= dl_opt
with the optical metric dl_opt^2 = h_ij dx^i dx^j / f. Hence
q in J+(p) requires Delta t >= d_opt(P, Q); conversely the path that
follows a minimizing optical geodesic with dt = dl_opt is null, so
the condition is sharp: **causal order is Delta t >= T with
T = optical distance.** (This is the note's Section 1, restated so T
is unambiguously the metric distance, not merely "the direct arc".)

**L0b (spherical optical metric).** For Schwarzschild,
dl_opt^2 = dr^2 / f^2 + (r^2/f) dOmega^2
         = drho^2 + w(rho)^2 (dtheta^2 + sin^2 theta dphi^2),
with rho the tortoise coordinate and w = r/sqrt(f); the radial
optical distance is exactly the tortoise difference.

**L0c (meridian projection: 3D -> 2D).** Let P lie on the polar axis
and X be arbitrary; rotate so X sits in the meridian half-plane
phi = 0. The projection (rho, theta, phi) -> (rho, theta) fixes P
(a pole) and X and, by L0b, maps any path to one whose length drops
the nonnegative w^2 sin^2 theta dphi^2 term: the projection never
increases length. Therefore the distance from an axis point is
realized in the 2D metric

    g_2 = drho^2 + w(rho)^2 dpsi^2,   psi = polar angle,

and T(P, X) = d_{g_2}((rho_P, 0), (rho_X, psi_X)). Both oracle
anchors are on the axis (radial alignment), so every distance the
oracle needs is of this form. This is the promoted, separately
stated version of the note's Section 2.

**L0d (the plane is Cartan-Hadamard).** Extend g_2 to the full plane
(rho, psi) in R^2. It is smooth and complete (w > 0, rho-lines have
infinite length both ways), and its Gauss curvature (L3) satisfies
K < 0 everywhere on the exterior r > 2M. So the plane is a
Cartan-Hadamard surface: any two points are joined by a unique
geodesic, which is minimizing, there are no conjugate points and no
cut locus, and Hessian comparison holds globally. A minimizer from
(rho_P, 0) to (rho_X, psi_X) has Clairaut constant b = w^2 dpsi/dl
conserved, so psi is monotone along it and it stays in the strip
psi in [0, psi_X] subset [0, pi): the full-plane distance equals the
half-plane (physical) distance of L0c.

## L1. Eikonal identities and Lipschitz bounds

Let T(X) = d_{g_2}(P, X) with P on the axis. Away from P, T is smooth
(L0d) and:

- **Eikonal.** |grad T|_{g_2} = 1, i.e.
  (dT/drho)^2 + w^{-2} (dT/dpsi)^2 = 1.
- **Angular derivative is the impact parameter.**
  dT/dpsi = b, the Clairaut constant of the connecting geodesic,
  with 0 <= b and, since b = w sin(alpha) at every orbit point,
  b <= w at each endpoint: b <= min(w(r_P), w(r_X)).
- **Coordinate Lipschitz bounds.** |dT/drho| <= 1 and
  |dT/dpsi| <= w(r_X); in the original coordinates
  |dT/dr| <= 1/f and |dT/dtheta| <= r/sqrt(f).
- **Monotone optical radius.** w'(rho) = dw/drho = (1 - 3M/r)/sqrt(f)
  is positive for r > 3M and w attains its global exterior minimum
  w_glob = 3 sqrt(3) M at the photon sphere r = 3M (differentiate
  w^2 = r^3/(r - 2M)).
- **Angle correction rule.** For fixed endpoint radii, T as a
  function of the separation angle psi has |dT/dpsi| <= L with
  L = min(w(r_1), w(r_2)); hence for any computed geodesic arc that
  connects the radii at angle a,

      |T(psi) - T(a)| <= L |psi - a|.

  This is the certified conversion of residual angle mismatch into
  time error; it needs no monotonicity of the swept angle in b.

## L2. Orbit containment, family exhaustiveness, and positivity

The orbit family solved by S1 is, in u = 1/r,
R(u) = 1/b^2 - u^2 (1 - 2Mu), with turning point u_t solving
R(u_t) = 0 and the exact factorization R(u) = (u_t - u) Q(u),

    Q(u) = (u + u_t) - 2M (u^2 + u u_t + u_t^2).

**L2a (perihelion lower bound by flat comparison).** With
1/b^2 = u_t^2 (1 - 2M u_t),
R(u) = (u_t^2 - u^2) - 2M (u_t^3 - u^3) < R_flat(u) = u_t^2 - u^2
for 0 < u < u_t. Hence each side's swept angle
int du / sqrt(R) exceeds the flat value arccos(u_e / u_t) =
arccos(r_p / r_e). If a matched one-turn arc between endpoints with
r_e >= r_min had perihelion r_p < r_min cos(psi/2), its total swept
angle would exceed 2 arccos(cos(psi/2)) = psi -- contradiction. So
every matched direct arc has

    r_p >= r_min_endpoints * cos(psi / 2).

For patch inputs (r >= 10M, psi <= 2): r_p >= 10 cos(1) M = 5.403 M.
For the frozen anchor box of L4 (r >= 11.353M, psi <= 0.8180):
r_p >= 11.353 cos(0.409) M = 10.417 M -- **all oracle orbits stay
inside the S1 shell.** No-turn arcs are monotone in r, so their
minimum radius is an endpoint radius; the bound is trivial there.

**L2b (Q > 0 on all arcs).** Q is concave in u (Q'' = -4M), so on
[0, u_t] its minimum is at an endpoint: Q(0) = u_t (1 - 2M u_t) > 0
and Q(u_t) = 2 u_t (1 - 3M u_t) > 0 whenever u_t < 1/(3M). By L2a
every arc the oracle meets has u_t <= 1/(5.4M) << 1/(3M), so the
factored integrand 2/sqrt(Q) is smooth and bounded on every arc.
The implementation still certifies Q > 0 on every quadrature
subinterval at runtime (defense in depth); this lemma guarantees
those checks succeed on-patch.

**L2c (family exhaustiveness).** The unique connecting geodesic
(L0d) has at most one turning point in the region it visits: turning
points solve w = b, w is strictly monotone for r > 3M (L1), and L2a
keeps oracle orbits at r >= 10.4M > 3M. So the geodesic is either
monotone in r (S1's no-turn family, including the sub-critical
b < 3 sqrt(3) M case with no real turning point) or has exactly one
perihelion (S1's one-turn family), and its Clairaut constant lies in
the searched brackets: b = w(r_p) <= w(r_inner endpoint) = b_eq for
one-turn arcs, b in (0, b_eq) for no-turn arcs. S1's two families
cover the geodesic; the certified solver searches the same brackets.

**L2d (certified root brackets).** All root-finding is by bisection
on brackets whose signs are certified by interval evaluation:
R(0) = 1/b^2 > 0 and R(1/(3M)) = 1/b^2 - 1/(27 M^2) < 0 iff
b > 3 sqrt(3) M (for M = 0 the bracket top 2/b has
R = -3/b^2 < 0). A bisection step whose sign is not certifiable
stops the refinement and returns the current bracket -- sound, merely
wider. u_t, b, and every downstream quantity are intervals; nothing
is ever a point estimate.

**L2e (no-turn angle monotonicity: the Clairaut constant is
identifiable).** For a no-turn configuration the swept angle is
ang(b) = INT_{u_out}^{u_in} du / sqrt(R(u; b)) with integration
limits INDEPENDENT of b, and R(u; b) = 1/b^2 - u^2 (1 - 2Mu) is
strictly decreasing in b pointwise, so the integrand is strictly
increasing in b pointwise: **ang(b) is strictly increasing on the
no-turn family.** (This is pure pointwise comparison; no derivative
interchange is needed. It does NOT extend to the one-turn family,
whose limits move with b -- no monotonicity is claimed or used
there.) Consequence: a bisection bracket [b_lo, b_hi] whose
endpoint angles are CERTIFIED to straddle the requested dpsi
encloses the unique matched b*, so for no-turn calls the certified
b* enclosure is the final straddling bracket -- as tight as the
bisection converges -- rather than the family-wide hull [0, b_eq].
If the straddle cannot be certified (quadrature too coarse at an
endpoint), the contract falls back to [0, b_eq], sound and wider.

## L3. Curvature and the two-sided Hessian comparison

**L3a (closed-form curvature).** For a surface of revolution
g_2 = drho^2 + w^2 dpsi^2, K = -w''(rho)/w. With
w' = (1 - 3M/r)/sqrt(f) and w'' = (M/(r^2 sqrt(f)))(2 - 3M/r):

    K(r) = - 2M/r^3 + 3M^2/r^4  =  -(2M/r^3)(1 - 3M/2r),

negative for all r > 1.5M, hence on the whole exterior. |K| is
decreasing in r for r > 2M. The comparison in L3b runs along the
WHOLE minimizing geodesic from the anchor to the evaluation point,
not just through the evaluation cell -- a one-turn geodesic dips to
its perihelion, which can lie below the cell's own radial range. So
the certified curvature constant is

    kappa^2 = sup |K| = 2M/r_geo^3 - 3M^2/r_geo^4,

with r_geo a certified lower bound on the radius ALONG the geodesic:
by L2a, r_geo = min(r_anchor, r_cell_lo) cos(psi_max_cell / 2), and
on the frozen L4 box every geodesic satisfies r_geo >= 10.417M
(the frozen perihelion floor), giving kappa^2 <= 1.52e-3 / M^2,
kappa <= 0.0390 / M. Using the cell's own r_lo instead of r_geo is
NOT admissible.

**L3b (two-sided Hessian comparison).** Let d = d_{g_2}(P, .) on the
Cartan-Hadamard plane (no cut locus, L0d), and let
Pi = g_2 - dd (x) dd be the projection orthogonal to grad d -- the
direction with vanishing Hessian is grad d, the gradient of the
geodesic distance, not the coordinate direction d/drho. For
K in [-kappa^2, 0], Riccati comparison gives BOTH sides:

    (1/d) Pi  <=  Hess d  <=  kappa coth(kappa d) Pi.

The lower bound (from the curvature upper bound K <= 0) is what makes
the tangent-plane cell model two-sided; the upper bound (from
K >= -kappa^2) satisfies the elementary estimate
kappa coth(kappa d) <= 1/d + kappa^2 d / 3. In particular
0 <= Hess d <= lambda Pi <= lambda g_2 with
lambda = kappa coth(kappa d_min) on any region where d >= d_min > 0,
so |Hess d (v, v)| <= lambda |v|^2_{g_2} for every v.

**L3c (coordinate conversion).** The nonzero Christoffel symbols of
g_2 are Gamma^rho_{psi psi} = -w w' and
Gamma^psi_{rho psi} = w'/w. For a distance sum S = T_1 + T_2 (both
1-Lipschitz), the coordinate second partials obey

    |S_rho rho| <= lam,
    |S_rho psi| <= lam w + 2 w',
    |S_psi psi| <= lam w^2 + 2 w w',

with lam = sum_i kappa_i coth(kappa_i d_i,min), where each kappa_i
is the PATH-WISE constant of L3a for anchor i (r_geo, not the cell's
r_lo) and w, w', d_i,min are certified interval bounds over the
cell. (Each bound is |Hess| <= lam on the metric components plus
|Gamma . grad S| with |S_rho| <= 2, |S_psi| <= 2w.)

## L4. Anchor-diamond containment (frozen configuration)

The static reduction makes the spatial shadow of the Alexandrov
diamond of p = (t_p, P), q = (t_q, Q) exactly the **optical ellipse**
{X : T(P, X) + T(X, Q) <= Delta t}. Two closed-form 1-Lipschitz
functionals separate it from the domain boundary:

- rho(r) is 1-Lipschitz for g_2 (|grad rho| = |drho| / dl <= 1), so
  changing radius costs at least the tortoise difference;
- w_glob * psi with w_glob = 3 sqrt(3) M is 1-Lipschitz
  (|grad (w_glob psi)| = w_glob / w <= 1), so changing the polar
  angle costs at least w_glob times the angle change.

**Frozen anchors (design-review ruling): r_p = 12M, r_q = 18M on the
polar axis, Delta t = 8.5M, M = 1.** With
rho(10) = 12.7726, rho(12) = 15.2189, rho(18) = 22.1589,
rho(20) = 24.3944 (4 dp):

| check                    | closed-form lower bound            | margin vs 8.5 |
|--------------------------|------------------------------------|---------------|
| diamond nonempty         | T(P,Q) = rho(18) - rho(12) = 6.9400 | +1.5600      |
| exit inner shell r <= 10 | (rho12 - rho10) + (rho18 - rho10) = 11.8326 | +3.3326 |
| exit outer shell r >= 20 | (rho20 - rho12) + (rho20 - rho18) = 11.4110 | +2.9110 |
| exit polar cap psi >= 1  | 2 w_glob = 10.3923                 | +1.8923       |

Containment box (from the same functionals):
rho(X) in [(rho_P + rho_Q - Delta t)/2, (rho_P + rho_Q + Delta t)/2],
i.e. r in [11.3536, 18.6950], and psi <= Delta t / (2 w_glob)
= 0.817913. Consequences: every solver call has dpsi <= 0.8180 <= 2,
and by L2a every orbit has perihelion >= 10.417M -- strictly inside
the S1 shell.

**Entry gate.** The oracle recomputes this table as certified
intervals at entry; if any margin interval is not certainly positive,
it refuses to run (no numbers are returned). The gate is part of the
oracle contract, mirroring the S4/S5 verify-freeze pattern.

## L5. The solver's direct arc is the global minimizer
(non-circular)

By L0c-L0d the distance from an axis anchor is realized by the unique
geodesic of a Cartan-Hadamard plane, and by L2c that geodesic is in
the family the solver searches. Therefore, for every configuration
the oracle evaluates:

- the solver's matched arc at Clairaut constant b IS the unique
  minimizing geodesic for the endpoint pair it connects, so its
  certified arc length encloses T exactly at its own swept angle, and
- the L1 angle-correction rule transports that value to the requested
  angle with a certified Lipschitz constant.

No step of this argument assumes a bound on T inside the diamond, so
the circularity flagged in review is gone: uniqueness comes from
curvature (L0d), not from comparing T to Delta t.

*Remark (independent cross-check, review-supplied).* In the
alternative equatorial-cylinder reduction (both endpoints off-axis,
psi in S^1), winding paths are not excluded topologically and one
must compare lengths. On the L4 containment box a uniform direct-path
bound gives T_direct <= (rho_hi - rho_lo) + w(r_hi) psi_max
= 8.5 + 19.783 * 0.817913 = 24.681 M, while any winding path costs
at least w_glob (2 pi - psi_max) = 28.398 M. The direct class wins
with margin 3.72M -- consistent with, and independent of, the
projection argument above. (At Delta t = 9: 25.53 < 28.15 also
passes; the frozen configuration uses 8.5.)

## L6. Cell tangent-plane enclosure and composition (PR-O2 contract)

The volume is V = 2 pi INT INT r^2 sin theta
[Delta t - S]_+ dr dtheta with S(X) = T_1 + T_2. The integrator must
implement exactly the following enclosure; nothing weaker counts as
certified.

**L6a (per-cell enclosure).** For a coordinate cell
C = [rho_0, rho_1] x [psi_0, psi_1] with center c whose optical
distances to both anchors are certified >= d_i,min > 0:

1. Certified center values: S(c) in [S_lo, S_hi] via the L1-L5
   flight-time contract (two calls).
2. Certified gradient at c via the eikonal (L1):
   dT_i/dpsi = b_i, where b_i is the certified b* hull the
   flight-time contract returns -- one-turn calls certify
   b* in [w(r_floor), b_eq] (L2a), while no-turn and
   equal-perihelion-corridor calls certify only b* in [0, b_eq];
   the tangent model must hull the gradient over the WHOLE returned
   interval. For no-turn calls, L2e upgrades the enclosure to the
   certified straddling bisection bracket (tight); the family-wide
   hull remains the fallback whenever the straddle is not certified.
   dT_i/drho = sigma_i sqrt(1 - b_i^2 / w^2), with the sign sigma_i
   determined by the arc geometry when certifiable (no-turn: the
   sign of r_X - r_anchor, since r is monotone along the arc;
   one-turn: +1 at the non-anchor endpoint, which sits past the
   perihelion on the outgoing side) and otherwise hulled over
   {-1, +1} (sound, wider).
3. Remainder: for x in C, the straight coordinate segment from c
   stays in C, so by L3c

       |S(x) - S(c) - grad S(c) . (x - c)|
         <= 1/2 ( M_rr drho^2 + 2 M_rp drho dpsi + M_pp dpsi^2 ),

   with (M_rr, M_rp, M_pp) the certified L3c bounds over C and
   (drho, dpsi) the half-widths. Call this remainder interval E_C.
4. The hinge is 1-Lipschitz and monotone:
   [Delta t - S]_+ on C is enclosed by applying the hinge to the
   affine model +- E_C endpoints.
5. The cell integral carries the COORDINATE JACOBIAN of the
   (rho, psi) chart: with dr = f drho,

       V = 2 pi INT INT r^2 sin psi dr dpsi
         = 2 pi INT INT r^2 f(r) sin psi drho dpsi,

   so the per-cell weight is r^2 f sin psi = r (r - 2M) sin psi --
   strictly increasing in r (hence in rho) for r > M and monotone in
   psi on the cap, so its exact per-cell interval sits at the cell
   corners. The cell integral encloses
   INT INT (r^2 f sin psi) (affine hinge) drho dpsi either
   (i) analytically -- the integral of an affine function's positive
   part over the (rho, psi) box with the weight enclosed by its
   exact per-cell interval -- or (ii) by an interval subcell sum
   that COVERS the cell (directed Riemann enclosure on every
   subcell, no point sampling). Dropping the f factor, or applying
   the affine-over-a-box rule in (r, psi) where the model is affine
   in rho, is NOT admissible. Finite pointwise sampling of the
   model is not admissible as a bound.

**L6b (anchor neighborhoods).** Near an anchor the tangent model's
curvature constant 1/d blows up, so cells within a declared optical
distance d_switch of either anchor are handled by one of two
admissible certified treatments:

1. (implemented) the FIRST-ORDER mode: center flight-time calls
   plus the L1 coordinate Lipschitz bounds |dS/drho| <= 2,
   |dS/dpsi| <= 2 w -- valid at every distance because no curvature
   constant enters; anchor neighborhoods then need no excision at
   all, and their enclosures tighten under the same adaptive
   refinement as every other cell (the small sin psi weight near the
   axis does the rest);
2. (admissible alternative) explicit excision of an optical ball of
   radius delta_b with its contribution ADDED as [0, U_ball],
   U_ball = 2 pi * Delta t * r_hi^2 * (1 - cos(delta_b / w_lb))
   * 2 delta_b f(r_hi), from the ball's coordinate bounding box
   (|rho - rho_anchor| <= delta_b, psi <= delta_b / w_lb, r-extent
   <= 2 delta_b f(r_hi) since dr/drho = f).

Either way, nothing is ever dropped as "negligible": option 1 keeps
the cells in the certified sum, option 2 keeps the ball term in it.

**L6b' (shell-local angular cost).** The pruning and anchor-distance
functionals may use, in place of the global w_glob = 3 sqrt(3) M,
the sharper certified bound

    d >= max( |Delta rho|,
              min( w(R_MIN) dpsi,
                   (rho_a - rho(R_MIN)) + (rho_X - rho(R_MIN)) ) ):

any path either stays in the shell r >= R_MIN, paying at least
w(R_MIN) per unit angle (w monotone above the photon sphere), or
leaves it, paying the tortoise exit cost from both endpoints. For
M = 0 the constant is 2 R_MIN / pi instead, from the chord bound
chord >= 2 sqrt(r1 r2) sin(psi/2) >= (2 R_MIN / pi) psi -- flat
space has NO photon-sphere floor and a naive "w >= R_MIN" would be
unsound.

**L6c (pruning tier).** A cell where the closed-form L4 lower bounds
already certify T_1 + T_2 > Delta t contributes exactly [0, 0]
without any solver call. Pruning uses only certified lower bounds, so
it can never discard a cell that intersects the diamond.

**L6d (composition and refinement).** The global interval is the
directed-rounding interval sum of all cell intervals plus all
[0, U_ball] terms -- a finite sum, no statistical step. Under
refinement the new global interval is INTERSECTED with the previous
one; an empty intersection is a certification failure and hard-stops
the run (it means an enclosure was wrong, not that precision
improved). Monte Carlo cross-checks are diagnostic only: a
non-intersection raises an investigation flag, never a pass/fail
verdict, and no probabilistic statement enters the certified
interval.

**Target (frozen).** The oracle's published goal is

    (V_hi - V_lo) / (V_hi + V_lo) <= 0.01,  with V_lo > 0 certified
    as a precondition of the ratio.

Failure to reach it under the cost caps (solver-call cap, wall-clock
cap, refinement-depth cap) returns the certified interval with status
`target-not-met`; the interval itself remains certified regardless.

**Price-measurement discipline (the [TO SIZE] deliverable).** The
cost of reaching each ratio is measured on a NEIGHBOR anchor
configuration, not on the frozen one: the oracle is deterministic,
but the freeze discipline exists so that nothing about the frozen
run is chosen after seeing its number, so the frozen configuration's
volume stays unobserved until the PR-O3 execution. The ladder is ONE
adaptive pass that records the cost at each first crossing (a
restarted run per rung discards every earlier rung's refinement and
measures the same numbers at several times the price), it streams
its trace so a slow run is distinguishable from a hung one, and it
publishes the full (calls, cells, ratio, wall) curve so the
convergence exponent is auditable.

The convergence is **floor-limited**, and the analysis has to say
so. Every cell strictly inside the diamond keeps its center flight
time's uncertainty however small the cell becomes, so cell
refinement converges not to zero but to

    W_floor(n_sub) ~ width(T1 + T2) * 2 pi INT INT r^2 sin psi
                      dr dpsi   over the diamond's support,

which only a larger `n_sub` reduces. The artifact reports this floor
per quadrature setting and extrapolates unreached targets with the
floor-aware model `ratio = floor + C * calls^slope`; a plain power
law would contradict the floor diagnostic printed beside it. A
target at or below the floor is reported as unreachable at that
`n_sub`, never as a call count. The floor is measured STANDALONE
(support integral on a fixed grid with the fast non-certified
solver, widths from the certified solver): a per-cell tally
accumulated inside the integrator also counts cells OUTSIDE the
diamond, whose true contribution width is zero, and overstates the
floor by orders of magnitude -- an error that inverted the verdict
on which lever binds before it was caught.

Extrapolations live in fields that cannot be mistaken for measured
crossings, and the artifact carries a per-`n_sub` plan for the
FROZEN target so PR-O3 has configuration GUIDANCE instead of a
number to guess.

**Epistemic grading of the plan (review ruling on PR #65).** The
plan rows are model-based, and the wording is bounded accordingly:

- `n_sub = 16` is the current RECOMMENDED PLAN, not a frozen
  parameter;
- "a larger `n_sub` is counterproductive" is a MODEL-BASED
  EXPECTATION -- the `n_sub = 32/64` floor rows combine flight-time
  widths at a few probe points with a linear cost model, they are
  not full-integrator runs;
- the frozen target is ESTIMATED reachable; the extrapolated call
  count is reported with its fit-window range (the spread across
  windows is the model's honest error bar) and is NEVER frozen as
  an execution cap;
- the allowed bottleneck statement is that cell refinement is the
  LIKELY binding lever and the first-order mode a CANDIDATE carrier
  of the remaining width -- until the mode-width diagnostic
  (`o2_mode_width.py`, run on the NEIGHBOR configuration only, with
  raw/certified width decomposition and intersection accounting)
  has been executed and committed. That diagnostic is a
  PREREQUISITE for any O3 budget or algorithm decision; running the
  frozen configuration first and deciding afterwards would collapse
  the freeze boundary.

**MC is never a CI gate.** The Monte Carlo cross-check's overlap
with the certified interval is computed and RECORDED by the
runners, but no test may pass or fail on it: a probabilistic
comparison can reject a perfectly correct oracle, so gating CI on
overlap would contradict the "never a verdict" contract above. If a
probabilistic gate is ever wanted, it must be redesigned with
frozen coverage and failure-rate accounting -- as a separate,
explicitly approved change.

## Relation to the note and to Paper A

This document is PR-O1 of the oracle arc: L0-L5 are implemented and
tested by `experiments/oracle/`; L6 is the frozen contract for the
PR-O2 integrator; the note's Section 4 status tags flipped in PR-O3
when the assembled bounded-error volume landed
(`docs/prereg/p14_o3_volume.json`). The certified volume has since
been consumed once, by the auxiliary O4b direct-MC instrument audit
(CONCORDANT at the single frozen configuration, verdict recovered
from the completed run's preserved statistics —
`docs/prereg/p14_o4b_results.json`); per the Section 6.7 claim
boundary that audit upgrades nothing in Paper A. The
prediction-anchored Poisson-count stage is separate from that audit and
has since been executed across the S6 mass ladder — CONCORDANT on each
of the three preregistered rungs, per rung and with no joint verdict.

## S6 addendum: mass-general instantiation (the ladder rungs)

The lemmas above were STATED for the frozen configuration M = 1; this
addendum records which of them are mass-generic, under exactly what
conditions, and how a new rung re-certifies them. Nothing here
weakens the M = 1 statements, and nothing is inherited silently: a
rung exists only after every row of its own certified table passes.

**Ladder convention (PI ruling).** The shell [10, 20], the cap
psi <= 2 and the anchors (12, 18) stay FIXED in absolute
coordinates; only M changes, so a rung is never an isometric copy of
another. The time window is not free either:

    dt(M) = 8.5 * T_min(M) / T_min(1),     KAPPA = 8.5 / T_min(1)
                                                 = 1.224782580680992,

the central rung's dimensionless slack carried to every rung (the
ratio form reproduces dt(1) = 8.5 exactly in binary64). The
pre-frozen dimensionless curvature indicator is the compactness
MU = 2M/r_c at the anchor midpoint r_c = 15.

**Mass-generic lemma conditions.** For a rung at mass M on this
shell (M = 0 is out of scope -- its angle-cost constant is the
different formula 2 R_MIN / pi):

| lemma | condition at mass M | why |
|---|---|---|
| exterior | 2M < R_MIN | horizon below the shell |
| L3 (K < 0) | 1.5M < R_MIN | K = -(2M/r^3)(1 - 3M/2r) |
| L1/L2c (w monotone) | 3M < R_MIN | photon sphere below the shell |
| L2b (Q > 0) | perihelion floor > 3M | u_t <= 1/r_p < 1/(3M) |
| L2a (patch orbits) | perihelion floor vs R_MIN cos(1) | certified per rung |
| L4 (containment) | all four margins certainly positive | the entry gate, re-run per rung |
| L5 cross-check | w_glob (2pi - psi_max) > dt + w(r_hi) psi_max | winding vs direct, re-verified per rung |

`experiments/oracle/s6_rungs.py` evaluates every row as a certified
interval comparison through the SAME parametric entry gate the
runners use (`containment_certificate(12, 18, dt(M), m=M)`), derives
the box and SCALE through `o4_sizing`'s exact path (o4_sizing itself
is unchanged -- it is the O4b recovery's frozen decision surface and
remains the M = 1 instance, reproduced bit-exactly by the parametric
path), and pins the frozen per-rung constants as import-checked
literals.

**Certified ladder tables (frozen).** All rows PASS for every rung;
margins are the certified lower bounds:

| M | mu | dt | box r | psi_max | SCALE | binding margins |
|---|---|---|---|---|---|---|
| 1.0 | 0.1333 | 8.5 (exact) | [11.3536, 18.6950] | 0.817913 | 28546.5704734038 | nonempty +1.5600, polar cap +1.8923, winding +3.72 |
| 1.4 | 0.1867 | 9.070565190742672 | [11.3672, 18.7053] | 0.623440 | 18141.08004658374 | nonempty +1.6647, polar cap +5.4787, winding +19.46 |
| 1.8 | 0.2400 | 9.725248174609407 | [11.3825, 18.7174] | 0.519896 | 13679.093767152488 | nonempty +1.7849, polar cap +8.9809, winding +33.35 |
| 3.0 | 0.4000 | 12.442423039673733 | [11.4429, 18.7691] | 0.399091 | 10472.024224793164 | nonempty +2.2835, polar cap +18.7345, winding +70.20 |

**Where the ladder stops, and why (the deep end).** The rungs are not
bounded by taste but by the table above. Sweeping M through the same
certified entry gate, every row passes exactly on

    M in [0.92, 3.33],   i.e.   MU in [0.1227, 0.4440],

bounded BELOW by the L5 winding cross-check and ABOVE by L1/L2c: the
photon sphere 3M reaches the shell floor R_MIN = 10 at M = 10/3, and
`w_monotone` refuses every rung past it. Inside that window the
binding margin shrinks monotonically -- w_monotone clears by 4.6 at
M = 1.8, 1.0 at M = 3.0, 0.4 at M = 3.2 and 0.1 at M = 3.3 -- while
every OTHER margin (L4 containment, L5 winding) GROWS with depth. The
strong-curvature rung is therefore set at M = 3.0: it keeps a full
tenth of R_MIN on the one condition that binds, where 3.2 or 3.3 would
buy at most 10% more MU in exchange for standing on the cliff. A
contract test pins both the ceiling and that margin, so no later rung
can be parked past the photon sphere by editing a literal.

The instrument paths are mass-parametric end to end: the wrapper
contract (`o4b_g3a`) accepts the rung mass explicitly and the S6
contract tests reject any path that silently falls back to M = 1.
Per-rung wrapper-preflight call counts, oracle targets and
count-stage sizings are frozen by each rung's own freeze, from that
rung's certified oracle endpoints -- never copied from another rung.
