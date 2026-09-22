# Bounded rotating-polarization vacuum generator

2026-09-22. This generator prepares a genuine varying-profile test of the
unchanged order-only Weyl estimator. It is not a reconstruction result.
The general Jacobi-propagator/action framework is established; see
[Harte and Drivas, *Caustics and wave propagation in curved spacetimes*](https://arxiv.org/pdf/1202.0540),
sections III--IV. The implementation specializes that framework to a
homogeneous rotating wave. No novelty is claimed for the generator.

## Geometry, domain and exact propagator

Use `g=2du dv+dx^T dx+x^T H(u)x du^2`, with

```text
H(u) = R(w(u-1/2)) D R(w(u-1/2))^T,
D = diag(A,-A), J = [[0,-1],[1,0]], R(theta)=exp(theta J),
0 <= A <= .75, |w| <= .5, 0 <= u <= 1.
```

The trace is zero, so this is Ricci-flat. When `A*w != 0`, H changes with
u in the parallel transverse Brinkmann frame. The midpoint has H=D,
so `C_uiuj=-D_ij` there, independent of w. The tips
`p=(0,.5,0,0), q=(1,-.5,0,0)` have proper separation `sqrt(2)`;
the dimensionless full tidal norm at their midpoint is `2A`.

Store transverse coordinates `y=R^T x` and momentum `P=R^T dx/du`.
Then

```text
d/du [y;P] = K [y;P],    K=[[-wJ,I],[D,-wJ]],
F(L) = exp(KL) = [[a(L),b(L)],[c(L),d(L)]].
```

This is a canonical change of transverse variables. K is Hamiltonian and
F is symplectic. The production code evaluates its blocks as degree-12
Taylor polynomials in L^2, with `r_+=A-w^2`, `r_-=-A-w^2`. The divided
difference coefficient uses `sum(r_+^j*r_-^(n-j), j=0..n)`, avoiding both
division by A at the flat limit and subtraction of nearly equal cosines.
It returns `b(L)/L` without dividing a small computed matrix by L.
In the bounded domain, `|r_+|,|r_-|<=1`; the first omitted terms give
absolute analytic truncation tails below `1e-24` for each block entry.

There are no conjugate endpoints in this domain. For a zero-endpoint
variation z on an interval of length L<=1, the transverse second variation
is bounded below by `(1-A*L^2/pi^2)*integral|z'|^2`, using the Dirichlet
Poincare inequality. This is strictly positive at A<=.75. The unique
stationary transverse path therefore minimizes the action. Caustics,
singular endpoints and periods outside this domain are not implemented.

## Action and interval volume

For L>0 and endpoints y0,y1,

```text
S = 1/2 [y0^T b^-1 a y0 - 2 y0^T b^-1 y1 + y1^T d b^-1 y1].
```

The symplectic identities make both quadratic matrices symmetric. A
future timelike relation has `v0-v1-S>0`. Negative-u queries swap the
endpoints before evaluating S; treating the two transverse endpoint orders
as interchangeable would be wrong for this rotating frame.

At an intermediate u in a length-L interval, the Hessian of the sum of
the two endpoint actions is
`Q=d(u)b(u)^-1 + b(L-u)^-1 a(L-u)`. Block composition gives
`b(L)=b(L-u) Q b(u)` and hence its determinant. If `h=v0-v1-S`, integrating
the remaining v width over the transverse ellipse gives

```text
V = pi L^2 h^2 integral_0^1 t(1-t)
       sqrt(det Bs(Lt) det Bs(L(1-t)) / det Bs(L)) dt,
Bs(L)=b(L)/L.
```

The formula is exact in this domain. Its final smooth one-dimensional
integral is evaluated by Gauss quadrature, with 32 versus 64 nodes checked
on every oracle evaluation. It is not replaced by the flat volume formula.

## Physical-volume sampling in the actual diamond

The volume form is `du dv d^2y` because the transverse rotation has unit
Jacobian. Equality of volume forms alone does not specify the diamond.
For the fixed tips, put `Q=left+right`, where
`left=d(u)b(u)^-1`, `right=b(1-u)^-1 a(1-u)`. Its cross section has
width `1-y^T Q y/2`, and its integrated volume density is
`f(u)=pi/sqrt(det Q)`.

A one-endpoint inequality
`integral|x|^2 <= L^2/2 integral|x'|^2` bounds each endpoint action below by
`(1-A/2)|y|^2/(2L)`. Thus

```text
Q >= (1-A/2)/(u(1-u)) I,
f(u) <= pi*u*(1-u)/(1-A/2) <= pi/[4(1-A/2)].
```

Reject a uniform u against this analytic ceiling. Given u, draw
`r^2=1-sqrt(1-U)` (implemented by its cancellation-free equivalent), an
independent uniform angle and a uniform v across the allowed width;
map the transverse disk by `sqrt(2) Q^-1/2`. The radial law is Beta(1,2),
not uniform radius or uniform disk area. The resulting density is uniform
with respect to physical four-volume in the curved diamond.

## Oracle normal coordinates and numerical checks

The inverse exponential map is used only after order-only fitting returns.
Let `U=u-1/2`, `X=Bs(U)^-1 y`, and `y(s)=s Bs(Us) X`. Then

```text
V = v + U A integral_0^1 [yx(s)^2-yy(s)^2
                          -2 w U (1-s) yx(s) yy(s)] ds,
normal = ((U-V)/2, X/sqrt(2), (U+V)/2).
```

The term containing w comes from H'(u) in the v-geodesic equation. Dropping
it would bias the evaluator's reference frame.

Before random data, `rotating_plane_wave.py` checks A=0,.25,.75;
w=-.5,0,.25,.5; positive, zero, negative and very short propagation lengths.
Observed maximum errors were:

- closed propagator versus SciPy expm: `2.22e-16`;
- normal-map round trip via independent original-Brinkmann geodesic ODE:
  `1.11e-15`;
- action versus an independently integrated path action: `2.47e-16`.

The check also covers symplecticity, time-reversed queries, the physical
slice density, rejection bound, quadrature refinement, and agreement with
the earlier constant-wave action, volume and normal map when w=0. A flat
slice attains the bound, so its numerical comparison allows 64 machine
epsilons; the measured discrepancy was `4.44e-16`.

The floating action guard assumes per-entry block error <=1e-12 and
propagates inverse and quadratic-form error conservatively. Unresolved
causal margins stop the generator; they are not dropped or assigned a sign.
This is a floating-arithmetic error model plus cross-checks, **not a
directed-rounding interval certificate**. The sampling arithmetic likewise
is not an exact-real random generator. These limitations remain explicit.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/rotating_plane_wave.py
```
