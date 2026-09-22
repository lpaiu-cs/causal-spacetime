# Frozen blind fifth-degree vacuum reconstruction pilot

2026-09-22. Commit code and this protocol before new random orders. This
prospective development pilot tests the operational version of the
[gradient diagnosis](gradient_model_diagnosis.md). It is not an accuracy
or coverage certification, and none of the earlier failures is re-scored.

## Model and blind boundary

The world-function model is

```text
S_ij = (x_i-x_j)^T eta (x_i-x_j) - C(x_i,x_j,x_i,x_j)/3
       - (nabla_e C)(x_i,x_j,x_i,x_j)(x_i+x_j)^e/12.
```

Jointly fit all unknown landmark coordinates (same tip/axis gauge), all ten
Weyl coefficients, and all **24** vacuum-compatible derivative coefficients.
The derivative space is the null space of the 16 independent differential
Bianchi constraints on 40 derivative components. No plane-wave algebraic
type, propagation direction, polarization, amplitude or frequency is imposed.
No oracle correction or true derivative enters the fit.

`fit_order_archive(path, model_order=5)` reads exactly the same six order
fields as before. The fifth-degree model is a fixed algorithm choice, not
case-specific metadata. The connected selector, correlated shared-membership
clock covariance, four flat initialization starts, optimizer tolerance,
iteration budget and all order counts remain unchanged. The enlarged
curved solve starts with zero curvature and gradient after the best flat fit.
Analytic Jacobians include the coordinate dependence of the gradient term.

Project curvature columns off BOTH coordinate and gradient nuisance columns.
Require full nuisance rank and ten projected Weyl components. Compute the
ten-component covariance by marginalizing the full joint Jacobian; do not
mistake the final ten gradient parameters for the Weyl coefficients. The
existing internal-frame/reflection confidence ball uses this marginal
covariance. Its coverage remains nominal until separately tested.

## Fixed data and paired reference

Use the unchanged exact rotating-wave generator and the ordered cases

```text
(A,w) = (0,0), (0,.5), (.25,0), (.75,0),
        (.25,.5), (.75,.5), (.25,-.5), (.75,-.5).
```

Case i has seed **2026092800+i**, with independent selection/measurement
streams `[seed,0]` and `[seed,1]`. Each uses 4096 selection events and
262144 measurement events, chunks of 4096, plus the supplied tips. The
maximum event budget is **2,129,920 plus tips**. No added events, seed
replacement, changed frequency or extra case after outcomes.

For every generated order, also fit the unchanged quartic procedure to the
SAME archive as a reference. It is another reconstruction of that case,
not a flat control. Evaluate the two models separately; failure of the
fifth-degree fit must not prevent the quartic reference from being attempted.
Preserve each failure in its assigned method. No post-result optimizer or
cutoff change is allowed within this run.

## Evaluation

The target remains all ten signed dimensionless Weyl components at the
midpoint in the order-defined internal frame, modulo residual reflection.
Only after a pure fit returns does its evaluator open the separate truth
archive. Report actual finite-order error, nominal radius/covariance,
optimizer, rank and boundary diagnostics, and the full fitted nuisance
gradient. No precision claim about the gradient itself is made.

For each model refit exact continuum volume clocks on the SAME selected
landmarks, selected pair mask and covariance weights. Report its full-tensor
continuum bias, with 32/64 quadrature refinement required below 1e-10.
Record the unchanged relative 2% bias condition (1e-10 numerical flat zero).
Do not relax that condition if the extra term fails to repair it.
Report actual finite errors and confidence widths alongside continuum bias;
improved continuum behavior alone does not meet the reconstruction objective.

This small stratified pilot has no powered accuracy or coverage gate. Any
later accuracy claim needs an explicit remainder/regularity scope, confidence
containment and 90% power sizing. The fifth-degree approximation still omits
second derivatives and curvature-squared terms; it is not an exact formula
for arbitrary oscillation frequencies, caustics or singularities.

Before data, the deterministic inverse self-check compares every coefficient
Jacobian and representative coordinate derivatives against finite differences,
checks the zero-gradient reduction, and recovers zero Weyl/gradient from an
independently constructed exact-flat interval matrix with full nuisance rank.
Replay an earlier archive in default quartic mode to verify compatibility.
Record reachable pre-data SHA, source hashes, all artifacts and failures.
After completion replay a fifth-degree archive through the pure estimator
and independently verify hashes, counts, ranks, covariance and comparisons.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_rotating_pilot.py --model-order 5 --data-seed 2026092800 --protocol jet_pilot_protocol.md --output outputs/intrinsic_weyl_jet_pilot_v1
```
