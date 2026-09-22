# Curvature-squared correction without new fitted components

2026-09-22. Development in an isolated checkout while the frozen
fifth-degree coverage campaign continues unchanged. This uses the existing
ten Weyl coefficients; the number of nuisance gradient components remains
24. It is **not a complete sixth-degree model**: second covariant curvature
derivatives remain omitted. The 2% bias allowance is not enlarged.

## Two terms already fixed by the fitted curvature

Work in dimensionless normal coordinates, with tips (-1/2,0) and (1/2,0).
For each pair x,y let eta=diag(-1,1,1,1), D=x-y, S2=D^T eta D and define

```text
u_a = C(x,y,x,e_a),     v_a = C(x,y,y,e_a).
S6_CC = (4 u.eta.v - u.eta.u - v.eta.v)/45.
```

This follows by collecting the curvature-product terms in
[Brewin, corrected v3, equation 11.21d](https://arxiv.org/html/0903.2087).
An independent action derivation already appears in the repository's
[vacuum generalization, G5](../window_response/vacuum_generalization_derivation.md).
There P=u and W=v-u, so `(2P^2+2P.W-W^2)/45` gives the displayed expression.
All contractions use the Lorentz metric, not a Euclidean norm.

The second correction comes from replacing squared proper times by volume
clocks. The repository's independently checked four-dimensional vacuum
volume expansion, [G2](../window_response/vacuum_generalization_derivation.md), is

```text
V = (pi/24) tau^4 [1 + tau^4 (17 E_U^2 + 3 B_U^2)/10080 + higher orders].
```

Use that derived coefficient, not the inconsistent printed coefficient in
Wang's equation 79; the existing discrepancy record remains unchanged.
Define `K=C_abcd C^abcd/8 = E_U^2-B_U^2` and
`E_D,ab=C_acbd D^c D^d`. Then

```text
Q(D) = (20 E_D,ab E_D^ab - 3 S2^2 K)/10080
     = tau^4 (17 E_U^2 + 3 B_U^2)/10080,   U=D/tau.
```

The polynomial expression has no division by tau. It is also well defined
during a local optimizer step. With `D0=(1,0,0,0)`, expanding BOTH interval
and outer-reference volumes gives the measured-clock model

```text
-sqrt(V_xy/V_outer) = F5(x,y,C,grad C)
                     + S6_CC + S2 [Q(D)-Q(D0)]/2
                     + omitted second-derivative and higher terms.
```

Replacing true pair tangents by D and midpoint curvature by the central C
is sufficient at the retained curvature-squared order. Corrections to
those replacements enter later. The omitted second-derivative terms are
of the SAME degree as S6_CC, so this model does not improve the generic
asymptotic bias order proved conditionally in the earlier derivation. Its
candidate benefit is reducing the known part of the finite-window bias
without adding unknown coefficients or imposing a wave type/profile.

## Checks and diagnostic fixed before execution

The production helper uses only inferred coordinates and fitted Weyl
coefficients. Its analytic derivatives cover all coordinates and all ten
Weyl coefficients; gradient-coefficient derivatives of this correction are
zero. A generic electric-plus-magnetic tensor complex-step check, pair
exchange symmetry and the flat zero case are executable in the module.

Check the physical coefficient and sign on the fixed five-point constant
wave configuration used by the earlier gradient check. Use exactly the
scale ladder `[1,.5,.25,.125]` with A=.75*s^2 and w=0. Compare both exact
proper-time squares and exact volume clocks before/after the correction.
Since this locally symmetric geometry has vanishing covariant curvature
derivatives, the corrected residual should have the next higher order.
Preserve all four scales; do not change the ladder to hide a discrepancy.

For failure diagnosis, use exactly the three already completed coverage
cases **27,29,49**, identified after their original results. They all have
A=.75 and w=-.5. Keep each case's landmarks, clock mask and covariance.
Compare the saved fifth-degree continuum fit with: exact proper-time input,
exact fifth-degree polynomial input, true-C subtraction of the world term,
true-C subtraction of the volume term, and true-C subtraction of both.
These are explicitly ORACLE diagnostics, not operational recovery.

Separately fit the new correction as a function of the UNKNOWN fitted C,
first to continuum volumes and then to the actual saved counts. The latter
uses only the order archive. These are post-result development comparisons
on selected failures, not independent validation. Record all failures,
errors, fit ranks and covariance changes. No fifth-degree coverage row or
verdict is changed by this diagnostic. Any adoption requires fresh data.
