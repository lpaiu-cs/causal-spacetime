# Window-response study: G1 theory and validation protocol

Date: 2026-09-21. Parent design: `supplementary_research_plan_2026-09-21.md`.
The user has authorized execution of that gated plan. No further routine phase
approval is needed; its stopping conditions and CPU caps still apply.
This document starts before the numerical validation and new poset campaigns.

## 1. Model and exact reduction

Use `g_A = 2 du dv + A(x²-y²) du² + dx² + dy²`. In this note the box has
`Du=Dv=1`, transverse full widths `a,b`, and its transverse centre is zero.
The first conjugate separation is excluded: `sqrt(abs(A)) < pi`.
Negative A means exchange x and y, not a negative oscillator frequency.
Here A equals the design's dimensionless lambda because Du=1.

Sort an independent pair by u. Its separation s has density `2(1-s)` on (0,1).
Conditional on the other coordinates, integrating the two uniform v coordinates
gives the exact causal-pair probability `G(C_A)` where

```text
G(C) = 0                         C >= 1
       (1-C)^2/2                 0 <= C <= 1
       1-(1+C)^2/2              -1 <= C <= 0
       1                         C <= -1.
```

G is continuously differentiable and its derivative is Lipschitz. For a point
pair write `d_x=x_q-x_p`, `m_x=(x_q+x_p)/2`, similarly for y. The exact action is

```text
C_A = w tanh(ws/2) m_x² + (w/4)coth(ws/2) d_x²
    - w tan(ws/2) m_y²  + (w/4)cot(ws/2) d_y²,   w=sqrt(A), A>0.
```

This is algebraically the same geodesic action as the existing certified P14
predicate. It is a different evaluation form, not an independent physical model.
Independent validation below uses numerical shooting and action integration.

The exact mean is `p_R(A)=E[G(C_A)]`, independent of N at fixed N>=2.
This is a five-dimensional integral: s, |d_x|, |d_y|, and two midpoint coordinates.
The difference density is `2(a-|d_x|)/a²`; conditional on |d_x| the midpoint is
uniform on `[-(a-|d_x|)/2, +(a-|d_x|)/2]`, and likewise for y.

## 2. Taylor coefficients without fitting sprinklings

Expanding the action gives `C_A=C0+A C1+A² C2+O(A³)` with

```text
C0 = (d_x²+d_y²)/(2s)
C1 = s [(m_x²-m_y²)/2 + (d_x²-d_y²)/24]
C2 = -s³ [(m_x²+m_y²)/24 + (d_x²+d_y²)/1440].
```

On every compact A neighbourhood inside the conjugate bound the action
derivatives are bounded uniformly in the box, including s approaching zero:
the A-dependent terms start at s, whereas only the A-independent C0 has 1/s.
G' is bounded and Lipschitz. Away from the zero-measure surfaces C0=0 and C0=1,
the second derivative of G exists. The Lipschitz bound controls second difference
quotients, so dominated convergence gives a twice differentiable integrated mean:

```text
p_R(A) = p0 + c1 A + c2 A² + o(A²)
c1 = E[-(1-C0) C1 1{C0<1}]
c2 = E[(-(1-C0) C2 + C1²/2) 1{C0<1}].
```

This establishes a local expansion, not a finite-A error certificate. In
particular, evenness does not license an O(A^4) remainder: G'' changes at C=0.
Possible nonanalytic higher terms remain unadjudicated. Finite-A predictions
will use the full exact integral, not extend this truncated expansion to A=2.

For the three boxes `min(a,b)>sqrt(2)`. Thus `C0<1` is a complete displacement
disk `d_x²+d_y²<2s` throughout 0<s<1. There is no clipping by the transverse
difference support in the coefficient integrals. Midpoint moments are

```text
E[m_x² | |d_x|] = (a-|d_x|)²/12
E[m_x^4 | |d_x|] = (a-|d_x|)^4/80.
```

After taking these moments, all terms reduce to integrals of
`s^k |d_x|^i |d_y|^j` over the disk with density `2(1-s)`. With
`h=(i+j+2)/2` and `t=k+h`, that integral equals

```text
4 * B((i+1)/2,(j+1)/2) * 2^h / [(i+j+2)(t+1)(t+2)].
```

Multiply by the displacement weight `(a-|d_x|)(b-|d_y|)/(a²b²)` before
integrating. Exact symbolic polynomial integration produced the expressions
below; independent Gauss integration of the unreduced expectation is a check.

| Region | p0 | c1 | c2 |
|---|---|---|---|
| R+, (a,b)=(2,6) | `pi/108 + 1/5184 - 128sqrt(2)/33075` | `pi/54 + 1/540 - 352sqrt(2)/25515` | `1613pi/48600 + 16733/2177280 - 1316008sqrt(2)/30405375` |
| R0, (sqrt(12),sqrt(12)) | `pi/108 + 1/5184 - 64sqrt(6)/33075` | 0 | `257pi/48600 + 73/435456 - 137204sqrt(6)/30405375` |
| R-, (6,2) | same as R+ | negative of R+ | same as R+ |

Numerical values from these expressions (not fitted simulation results):

| Region | p0 | c1 | c2 |
|---|---:|---:|---:|
| R+ | 0.0238087263042011 | +0.0405192770866489 | 0.0507424289511762 |
| R0 | 0.0245419687061855 | 0 | 0.00572728721822277 |
| R- | 0.0238087263042011 | -0.0405192770866489 | 0.0507424289511762 |

The axis-exchange identity is exact: `p_(a,b)(A)=p_(b,a)(-A)`.
The square response is even; both anisotropic linear coefficients are nonzero.
The second coefficient is positive because C2<0 and both terms of its integrand
are nonnegative, strictly positive on a set of positive measure.

These explicit coefficients and their explanation are additional to P3-E's
empirical exponent. G0 can retain this narrow contribution for validation.
This is not a proof of literature priority or of publication significance.

## 3. Fixed-N noise and the information supplied

For a single poset let `f_N = sum_{i<j} h_ij / binom(N,2)`. Under iid sampling,

```text
zeta = E[h12 h13] - p² = Var(E[h12|X1])
Var(f_N) = 2p(1-p)/[N(N-1)] + 4(N-2)zeta/[N(N-1)].
```

This identity is standard. Neither millions of pair evaluations nor shared
points across A create independent poset replications. We compute p and zeta
with separate two/three-point integrals; simulation validates the resulting
distribution approximation. Single-poset thresholds know region and N, not A.

For integration, condition on each point's u,x,y coordinates. Its v coordinate
is uniform on an interval. Given the anchor v, each relation's conditional
probability is a clipped affine function of that v. Their product is piecewise
quadratic. Integrate it exactly using its breakpoints and Simpson's rule on
each piece. This removes discontinuous pair indicators from the nine-dimensional
triple integral. Independent v values are integrated, not shared between points.

## 4. The separate diamond population

Tips are `(0,+1/2,0,0)` and `(1,-1/2,0,0)`. At 0<u<1 define

```text
Bx = (w/2)[coth(wu)+coth(w(1-u))]
By = (w/2)[cot(wu)+cot(w(1-u))]
v_low = -1/2 + action((u,x,y),(1,0,0))
v_high = +1/2 - action((0,0,0),(u,x,y)).
```

The allowed ellipse has `Bx*x²+By*y²<1`, with v-width
`1-Bx*x²-By*y²`. Its slice volume is `pi/(2sqrt(Bx By))`.
Both B are positive before the first conjugate point, even where an individual
cotangent changes sign. The marginal u density is proportional to `1/sqrt(Bx By)`.
Conditional on u, draw angle uniformly and radius squared as Beta(1,2), then
scale the axes by `1/sqrt(Bx), 1/sqrt(By)`; draw v uniformly between its limits.
At A=0, `Bx=By=1/[2u(1-u)]` and volume is pi/6 for tip proper time sqrt(2).

The u density is maximal at 1/2: its reciprocal-square B product can be written
as the product of `sinh(wu)sinh(w(1-u))` and `sin(wu)sin(w(1-u))`, each maximized
there by concavity of log(sinh) and log(sin). Thus simple uniform-u rejection
sampling is available for actual posets without numerical inverse-CDF bias.
Quadrature uses a numerical inverse CDF, which is checked at two resolutions;
that discretization error is recorded separately from scramble uncertainty.

This geometrically defined population is not discovered from an arbitrary
unlabelled poset. No coordinate eligibility filter is silently added to f_N.

## 5. G1 checks fixed before running them

Development streams use namespace `[20260921, 6101, ...]`. Later quadrature,
pilot and confirmation use distinct namespaces. G1 creates no holdout data.

- Coefficients: closed expressions versus direct Gauss integration at orders
  24 and 40; symmetry and sign identities; report maximum numerical differences.
- Predicate: 256 predeclared random pairs per A in {0, 1/64, 1/4, 1, 2};
  compare half-angle evaluation with existing certified margins. Add exact
  near-cone cases, coincident u, and x/y reversal. Undecided pairs escalate.
- Independent geometry: numerical ODE shooting plus action integration for
  24 pairs per A in {0, 1/4, 1, 2}. State integration tolerances and residuals.
- Conditional-v integrator: compare with independent one-dimensional adaptive
  integration on deterministic intervals, including opposite temporal directions.
- Diamond: check tip membership with the existing predicate, flat volume pi/6,
  and moments against independent cross-section integration. Use 4,096 sampled
  points per A in {0,1,2}; these are sampler diagnostics, not test posets.
- Positive control and end-to-end timing belong to the bounded pilot only after
  the above checks. No new formal physical verdict is produced by G1 checks.

Use existing NumPy and add SciPy to the experiment extra for standard scrambled
Sobol integration, ODE/quadrature and probability distributions. No ML model,
GPU, external service, new orchestration framework, or paid computation is needed.
