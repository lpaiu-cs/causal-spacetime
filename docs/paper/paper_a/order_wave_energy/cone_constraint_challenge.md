# Joint cone constraints: test the coordinate-adjustment obstruction first

2026-09-24. Pre-execution diagnostic, not a replacement objective. The target
is still the sharp finite-order law and an order-only attaining observable.
The current bracket is necessity below exponent 4/3 and sufficiency at 3.

The proposed next mechanism is to combine many cone constraints so that
coordinate errors can be eliminated jointly. Before constructing a test,
check whether the leading wave perturbation survives **arbitrary simultaneous
infinitesimal adjustments of all event coordinates**. If it does not, a
linear coordinate-elimination statistic cannot exploit that leading term.
This diagnostic deliberately uses generator coordinates. Its outputs are
not an observer and cannot establish a better sufficient region.

## Mathematical quantity

Use background Minkowski coordinates `x=(T,Z,y1,y2)`, with
`T=(m+v)/sqrt(2)`, `Z=(m-v)/sqrt(2)`. For an oriented pair put
`F=Delta T-|Delta z|`, `n=Delta z/|Delta z|`. Its differential is the row

\[
(A\xi)_{ij}=(-1,n)\cdot\xi_i+(1,-n)\cdot\xi_j.
\]

Let `E_ij=Delta m I-M_lambda(ui,uj)` be the trace-free primitive correction
already used by the exact circular generator. The first-order change in F is

\[
b_{ij}=-\frac{\Delta y^T E_{ij}\Delta y}
 {\Delta m(\Delta T+|\Delta z|)}.
\]

The exact change in the same normalized null action is
`-Delta y^T M_lambda^-1 E_ij Delta y/(Delta T+|Delta z|)`.
On macroscopic longitudinal segments b=O(lambda^2), and their difference
is O(lambda^4). We fit **all four coordinate changes at every event**, not
only time shifts, and retain `r=b-A argmin_xi ||b-Axi||_2`.
At an exact minimizer, `A^T r=0`; this is the linear combination that cancels
coordinate adjustments. A nonzero residual is only a local linear-algebra
obstruction. It does not exclude remote nonlinear embeddings, recover the
coordinates from an order, or provide a sampling distribution for r.

## Frozen campaign

- Fixed q=1, phase 0.37, N in `{128,256,512}`, seeds `{118,819}`, and
  r=1/lambda in `{32,128,512}`: 18 rows, all retained.
- Generate each (N,seed) once from the **background** physical-volume law.
  In the m,v,y chart this has density proportional to a(u(m))^4 on the
  fixed Minkowski diamond. Reuse the same coordinates for all three rungs.
  This conditional geometric diagnostic is not a wave-sprinkling campaign.
- Select oriented pairs with `|F|<=L*4 log(N)/N`, where L=sqrt(2M), and
  `Delta m>=M/10`, `Delta v>=1/10`. These fixed macroscopic restrictions
  avoid a short-segment singularity. Record edge count, degrees, and width.
  The selection uses coordinates and is not claimed order-measurable.
- Solve the sparse least-squares problem with LSMR, `atol=btol=1e-11`,
  `conlim=1e12`, `maxiter=20000`. Record the scaled forcing RMS, scaled
  residual RMS, residual/forcing ratio, stationarity residual, iterations,
  stopping code, and norm of the fitted coordinate adjustment.
- Before interpreting each graph, feed `A xi_control` from a fixed random
  coordinate adjustment. Require relative residual below `1e-7` and a
  successful LSMR stop in `{1,2,4,5}`. A failure marks that graph's diagnostic
  numerically unresolved; do not silently increase the iteration cap.
- At N=128,seed=118,r=512 compare against dense SVD least squares with
  singular-value cutoff `1e-12`. At N=512,seed=819,r=512 rerun only the
  linear solve with tolerances `1e-13` and cap 40000. Record differences
  without treating a numerical agreement as independent researcher review.
- Self-check the cone Jacobian against central finite differences on a
  fixed coordinate perturbation; check zero wave forcing at q=0 through
  direct zero matrices, translation null vectors, and the control projection.

No finite-sample accuracy or power gate is present. Residuals and their
frequency dependence are diagnostic outputs, not confidence statements about
the asymptotic exponent. Do not relabel this as a full-order success even
if all residuals remain large. A usable order-only statistic would still
need latent-coordinate elimination and a valid likelihood or confidence
construction, and a sharp boundary still needs matching impossibility.

## Why existing reconstruction results do not settle this

[Dani et al. (2022)](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2022.48)
give quantitative reconstruction results for Euclidean random geometric
graphs, using graph distances and common-neighbor estimates. Those results
are not a theorem of 1/N coordinate recovery from four-dimensional causal
orders. We will not assume such a recovery rate to manufacture the missing
upper bound. The earlier null-incidence and causal-set embedding references
are discussed in [the current results](local_null_results.md).
