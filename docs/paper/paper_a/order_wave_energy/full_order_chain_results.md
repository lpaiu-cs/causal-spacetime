# Resolution-law progress and the fixed-chain bottleneck

2026-09-24. **The sharp full-order resolution law remains unresolved.**
This record establishes a concrete sufficient side and explains why merely
adding a fixed list of chain counts is not a promising way to close the gap.
It is a same-author derivation and numerical check, not an independent audit
or a claim of priority.

## Frozen run and what it checks

The [derivation/design](full_order_chain_asymptotics.md) and diagnostic source
were frozen at `e5cfddfdd7bab2d744756fc00cef05ed4b219f6f`. The run started from
a clean worktree. Raw rows have SHA256
`33aa0c4531b92a4f61e1533feb23334ac04fd7b2695d3cb8d67c6b3cec7b1fdc`.
The [portable evidence](resolution_evidence/manifest.json) includes 18 primary
conditions and 2 refinements, 20 oscillatory-primitive ODE solves in 11.10 s
on CPU. No finite causal set was sampled by this diagnostic.

The earlier independent Rosen ODE agreed with the stable correction kernel
and pair probability in the executable controls. The weak-q limit also
agreed with the earlier transfer expression. Every primary pair change was
negative, as predicted. At q=1 the two refinements changed the pair difference
by relative amounts `9.77e-15` (r=128) and `2.80e-13` (r=1024).
Those are numerical consistency differences, not rigorous error enclosures.

For `Delta p=p_lambda-p_0`, the theorem predicts
`Delta p=-A(q) r^-4 log r+O(r^-4)`. The deliberately retained finite-r
coefficient and remainder were:

| q | Predicted A(q) | -Delta p r^4/log r over the grid | Delta p r^4+A(q)log r over the grid |
|---|---:|---:|---:|
| 0.25 | 3.562421 | 1.7104 to 2.6331 | 6.2280 to 6.4417 |
| 1 | 13.776152 | 6.6136 to 10.1947 | 23.7503 to 24.8246 |
| 4 | 45.625205 | 20.5054 to 33.6356 | 74.7951 to 87.0586 |

The normalized coefficient is **not close enough to its limit to claim
percent-level asymptotic accuracy**. Its slow `1/log r` correction was predicted
before the run, and the bounded-scale remainder is compatible with it on
the inspected grid. In particular q=4 is not monotone between r=512 and 1024;
endpoint oscillations are part of the allowed remainder. No fitted exponent
is used in the proof.

Verification/reproduction of the saved summary:

```powershell
E:/lab/causal-spacetime/.venv/Scripts/python.exe experiments/positive_control/order_wave_resolution.py --summarize docs/paper/paper_a/order_wave_energy/resolution_evidence
```

The verifier checks raw hashes, all source/design hashes at the frozen commit,
its ancestry, and the complete specified grid before summarizing.

## The two-sided bracket, with its gap intact

For the exact circular family at fixed known `0<q<pi^2`, with the fixed-N
physical-volume protocol and observer clock of the design:

| Region | Proven/derived conclusion within the stated model |
|---|---|
| `N lambda -> 0` | Full-order TV tends to zero by the earlier coupling argument. |
| `N lambda^8 log^2(1/lambda) -> infinity` | Full-order TV tends to one; a pair-count statistic suffices. A wavelength-free test has fixed false-positive bound 0.025 and power tending to one. |
| Between those regions | Undetermined. No sharp event-count exponent or attaining full-order statistic has been found. |

The sufficient side no longer uses a weak-q-first limit: the analytic
coefficient `A(q)>0` is proved at fixed q. The executable order statistic
uses only a comparable-pair count, N, and the declared q calibration. It
does not consume coordinates or wavelength. Its finite-N error bound comes
from a bounded U-statistic inequality, not from the deterministic ODE run.

These statements do not address an unknown-q composite comparison,
arbitrary wave profiles, caustics, an unknown observation window, or a
known-density Poisson count. Neither side has had an independent proof audit.

## Why a fixed menu of longer chains does not supply the missing lever

The following is a post-run mathematical deduction, not a new numerical
experiment. Fix any finite m and collect the fractions of k-chains for
`k=2,...,m`. Their population shifts all have order `lambda^4 log(1/lambda)`
by equation (6) of the derivation. A possible loophole would be a linear
combination whose leading sampling variance vanished. It does not exist
for these chain fractions in the effective background.

For one interior event X, write its past and future diamond volumes relative
to the large diamond as p(X), f(X). Condition on X being one of the k sampled
events. The probability that all k form a chain is

\[
h_k(X)=\sum_{j=0}^{k-1}\binom{k-1}{j}
  p(X)^j f(X)^{k-1-j}\,
  P_j(p_{\rm tip},X)P_{k-1-j}(X,q_{\rm tip}),
\tag{1}
\]

where `P_0=P_1=1`. The subdiamond probabilities depend on their u endpoints,
not on their transverse location or positive longitudinal slack, by the
plane-wave hierarchy. At fixed interior u, use the past and future slacks
eta,zeta>0, with `eta+zeta<h`. Then
`p=alpha(u) eta^2`, `f=beta(u) zeta^2` with positive coefficients.
The physical-volume distribution has a positive density on this open
slack triangle after integrating the transverse angle.

Thus h_k is a nonzero homogeneous polynomial of degree `2(k-1)` in eta,zeta,
with positive coefficients. If a finite linear combination of the centered
first projections `h_k-P_k` were constant almost everywhere, that polynomial
identity would hold on an open set. Its highest-degree term forces the
coefficient of the largest k to vanish; induction forces every coefficient
to vanish. Therefore their covariance matrix is strictly positive definite.

The standard multivariate Hoeffding decomposition consequently has
nondegenerate `N^-1/2` fluctuations in every nonzero fixed linear direction.
This remains true along lambda_N->0: first-projection covariances are
combinations of fixed finite-order probabilities, which converge by the
coupling bound, and bounded kernels give the triangular-array CLT and a
uniform `O(N^-2)` remainder variance. If
`sqrt(N) lambda^4 log(1/lambda) -> c`, then

\[
\sqrt N(\widehat{\boldsymbol P}-\boldsymbol P_0)
\Rightarrow\mathcal N(c\,\boldsymbol b(q),\Sigma_0(q)),
\qquad\Sigma_0(q)>0.
\tag{2}
\]

Here b is the explicit coefficient vector from the chain asymptotics and
`Sigma_0[j,k]=jk Cov(h_j,h_k)`.

Hence a fixed number of chain fractions can change constants but cannot
improve the power-law event price of **regular Gaussian/Wald tests or fixed
linear-threshold tests**. At c=0 such tests have the same limiting rejection
probability under the two models; at c tending to infinity the pair statistic
already separates them. This is not a total-variation lower bound for the
exact discrete count vector: weak convergence does not exclude every
irregular test of fine lattice probabilities. It also says nothing about
chain lengths increasing with N, other motifs, or the full order.

The practical research consequence is narrow and useful: adding P3, P4,
and a few more ordinary chain features is not a demonstrated route to the
sharp boundary. The unresolved target must use more of the order or improve
the full-order impossibility argument. This partial result is not being
promoted to completion or automatically recast as a small publication.

Subsequent progress: the [conditional null-coordinate test](null_coordinate_results.md)
uses a growing collection of event-specific statistics and improves the
sufficient exponent to 4, with a logarithmic loss, for known `0<q<=1`.
It does not contradict the fixed-vector limitation above. The full-order
boundary remains unresolved.
