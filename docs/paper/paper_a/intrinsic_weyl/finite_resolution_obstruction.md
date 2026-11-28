# Arbitrarily rapid vacuum waves prevent uniform finite-sample recovery

Development derivation, 2026-09-22. This is a scope obstruction, not a
replacement for the ongoing positive reconstruction and coverage work.
No priority claim is made. The general comparison of finite causal-order
probability laws is established in [Bombelli (2000), sections II--III](https://arxiv.org/pdf/gr-qc/0002053).
General plane-wave geometry and Jacobi propagators are reviewed by
[Harte and Drivas (2012)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.85.124039).
The explicit action/coupling bound below is derived here.

## Statement

Fix any finite event budget N. Without a common regularity or frequency
restriction on a vacuum wave profile, no estimator from a uniformly sampled
finite causal order can recover midpoint observer-tidal curvature with a
uniform small error probability over all smooth profiles. This remains true
in regular preconjugate diamonds with bounded curvature amplitude. Singular
points and complicated topology are not needed for this obstruction.

It does NOT rule out consistency for each fixed smooth profile as N grows,
or finite recovery on a specified regularity class. It rules out a uniform
claim over arbitrarily rapid profiles at any fixed N. The current constant-
wave campaign therefore cannot establish the unrestricted claim.

## A smooth oscillating family with an unchanged target

Use the repository sign convention

```text
g_w = 2 du dv + dx^2 + dy^2 + X^T H_w(u) X du^2,
H_w(u) = a cos(w(u-1/2)) diag(1,-1),   0 < a <= 1/4,
p = (0,1/2,0,0),  q = (1,-1/2,0,0).
```

The comparison metric g_0 is flat, with H identically zero; it is NOT the
w=0 member of this oscillating family. Each wave is smooth and Ricci-flat,
with `C_uiuj=-H_ij`. The axis geodesic between the tips has proper time
`T=sqrt(2)` for every w. At its midpoint H=a*diag(1,-1), so

```text
theta(g) = T^2 sqrt(||E||_F^2 + ||B||_F^2),
theta(flat) = 0,       theta(g_w) = 2a.
```

The observer is the tip geodesic's unit tangent. The norm does not depend
on a spatial-axis convention. This fixed scalar separation is already enough
to obstruct a uniform full-tensor recovery claim.

## 1. Uniform control of the exact causal action

For transverse endpoints bounded by R, and `0 < L = u_q-u_p <= 1`, define

```text
S_H = min_X (1/2) integral ( |X'|^2 + X^T H X ) du.
```

Future timelike relatedness is `v_p-v_q > S_H`. Dirichlet coercivity follows
from `pi^2/L^2-a > 0`, so there is a unique minimizing transverse path, with
`X''=HX`. This is an exact quadratic-action statement, not a small-diamond
curvature expansion.

Write X=l+y, where l is the linear interpolant and y vanishes at the two
endpoints. The Dirichlet inverse of the second derivative has sup-norm at
most `L^2/8`. Hence, with `K=(1-a/8)^(-1)`,

```text
||X||_infinity <= K R.
```

Since `y''=HX`, its derivative can be written as the difference of the two
weighted integrals on either side of u. Their total absolute weight is at
most L/2. Thus

```text
||y'||_infinity <= a K R L/2,
integral |X'| du <= 2R + a K R L^2/2 <= R(2+aK/2).
```

A primitive of the oscillating matrix H is
`F=(a/w) sin(w(u-1/2)) diag(1,-1)`, with norm at most a/w. Integration by
parts along the minimizing path gives

```text
|(1/2) integral X^T H X du|
 <= (a/w) [K^2 R^2 + K R * integral |X'| du]
 <= (a R^2/w) [K^2(1+a/2)+2K].
```

For the flat linear path the analogous upper bound is `3aR^2/w`.
Evaluating each quadratic functional on the other's minimizing path therefore
gives the two-sided comparison

```text
|S_H - S_flat| <= C(a) R^2/w,
C(a) = a [K^2(1+a/2)+2K].
```

The bound is uniform even when L is small. For a=1/4, K=32/31 and
`C(a)=784/961 < 0.816`. These are deliberately loose constants.

## 2. The diamond sampling measures also converge

It is insufficient to note that `sqrt(-g)=1`: the sampling domains differ.
Here that difference is bounded explicitly.

For a tip-to-point path, one endpoint of X is zero. The one-endpoint
Poincare inequality gives

```text
S_H >= (1-4a/pi^2) |X_endpoint|^2/(2L) = kappa |X_endpoint|^2/(2L).
```

Consequently every point in either diamond has `|X|^2 <= 1/(2*kappa) < 1`
and `-1/2 < v < 1/2`. Both diamonds lie inside the common bounded cylinder
`B=(0,1) x (-1/2,1/2) x { |X|<1 }`. We may use R=1 above.

For each fixed (u,X), the permitted v interval has endpoints given by the
two tip actions. Each changes by at most C/w. The symmetric-difference
four-volume therefore obeys

```text
Delta = Vol(D_H symmetric_difference D_flat) <= 2*pi*C/w.
```

The flat diamond has `V_flat=pi/6`. For the normalized uniform physical-volume
sampling laws mu_H and mu_flat, the exact overlap formula for uniform
densities yields

```text
TV(mu_H,mu_flat) <= Delta/V_flat <= 12*C/w.
```

Thus this argument includes both order changes AND changed sampling domains.

## 3. Coupling complete finite causal orders

Couple N labelled iid samples point by point using maximal couplings. The
probability that at least one coordinate differs is at most `12*N*C/w`.
On shared points, an internal pair can change its causal relation only when
its flat v-margin is within C/w of zero.

The flat point density is `1/V_flat`. At fixed u and transverse coordinates,
the area of this strip in the unit (v_i,v_j) square is at most `2*C/w`.
Integrating the other coordinates over the bounding cylinder gives

```text
P_flat(a given pair can flip) <= 2*pi^2*C/(V_flat^2*w) = 72*C/w.
```

Tip relations are already identical when points belong to both diamonds.
Union bounds and the coupling inequality give, for the complete labelled
order laws P_H,N and P_flat,N,

```text
TV(P_H,N, P_flat,N)
 <= min(1, (C/w) [12N + 72*binom(N,2)]) -> 0 as w -> infinity.
```

Removing labels, selecting landmarks, splitting event roles, counting
intervals or making any other data-dependent summary cannot increase total
variation. The obstruction therefore applies to the present pipeline as
well as to an estimator that uses the entire sampled order.

For a fixed known Poisson intensity rho the same conclusion follows by
coupling a common Poisson process in B and restricting it to the two
diamonds. With `nu=rho*V_flat`, an explicit bound is
`TV <= min(1, (C/w)(12*nu+36*nu^2))`. Knowing the density or tip time does
not remove this finite-resolution obstruction.

## 4. Statistical consequence

The target values are separated by 2a while their finite-data laws become
arbitrarily close. The success intervals `|theta_hat|<a` and
`|theta_hat-2a|<a` are disjoint. For any estimator, the elementary two-point
testing bound therefore implies

```text
max(P_flat(|theta_hat|>=a), P_H(|theta_hat-2a|>=a))
 >= (1-TV(P_H,N,P_flat,N))/2.
```

For every finite N, choose w large enough to make the right side arbitrarily
close to 1/2. Uniform high-probability recovery over all smooth oscillating
profiles is therefore impossible. A uniformly valid confidence set also
cannot be uniformly narrow enough to separate these two target values.

The flat member is not essential. Take two positive amplitudes, for example
a=1/4 and b=1/8, at the same growing frequency. Both order laws converge to
the flat law, hence to one another by the triangle inequality, while their
tidal norms remain 1/2 and 1/4. Their relative +/-25% accuracy bands are
disjoint. The obstruction therefore also applies when both candidate
spacetimes have nonzero oscillating Weyl curvature.

This is an existence obstruction, NOT a useful numerical cutoff or a tight
sample-complexity rate: the N-squared union bound is intentionally crude.
Its practical conclusion is that a declared modulus of curvature variation,
frequency class or comparable local regularity control is essential for a
uniform pointwise accuracy/cost claim. More test profiles alone cannot
replace that assumption. Positive finite-order reconstruction under an
explicit class and a tested uncertainty procedure remains the active task.
