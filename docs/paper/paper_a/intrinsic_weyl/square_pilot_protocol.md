# Frozen fresh pilot for the curvature-squared correction

2026-09-22. Commit this protocol and code before drawing the eight new orders.
This is a development comparison of the [derived correction](square_correction.md),
motivated by three selected failures of the fifth-degree coverage campaign.
Their outcomes are not reused as validation and that campaign is unchanged.
This pilot has no accuracy or coverage certification gate.

## Estimator and comparator

Fit the fifth-degree vacuum jet plus the two fixed curvature-square terms
in squared distance and interval/reference volume clocks. These terms depend
only on the already fitted ten Weyl coefficients; there are still 24 gradient
coefficients and the same unknown coordinates. No second-curvature-derivative
parameters or wave-profile constraints are added. The generic remainder and
conditional asymptotic bias order remain unchanged.

The corrected estimator receives only the same six-field order archive.
The connected selector, four initial flat starts, clock mask, covariance
construction, optimizer, tolerances and confidence-ball formula remain fixed.
The covariance uses the corrected analytic Jacobian and marginalizes both
coordinate and gradient nuisance columns. Curvature-square terms have zero
derivative at C=0, preserving the flat-limit information matrix.

Fit the unchanged fifth-degree estimator to the SAME archive as a paired
reference, irrespective of whether the corrected fit succeeds. This comparator
has the same parameters and data, not a paired flat control. Keep method-specific
failures and all full fits. Do not tune either method during this run.

## Data fixed before generation

Use the unchanged exact rotating-wave generator and ordered cases

```text
(A,w) = (0,0), (0,.5), (.25,0), (.75,0),
        (.25,.5), (.75,.5), (.25,-.5), (.75,-.5).
```

Case i uses seed **2026100500+i**; selection and measurement use independent
streams `[seed,0]` and `[seed,1]`. Each uses 4096 selection events and 262144
measurement events in chunks of 4096, plus the supplied tips. The fixed budget
is **2,129,920 events plus tips**. No replacement, extended count, extra case,
changed frequency or new initialization after seeing outcomes.

## Evaluation and checks

Open each separate truth archive only after its order-only fit returns.
Report ten-component tensor error in the order-defined internal frame modulo
reflection, confidence radius, largest covariance eigenvalue, ranks and fit
diagnostics for BOTH methods. For each, refit exact continuum volume clocks
on the same landmarks, clock mask and empirical weights. Require quadrature
refinement below 1e-10. Preserve the relative 2% bias allowance and numerical
flat-zero tolerance 1e-10 without relaxation.

Report every paired bias, finite-error and variance change; do not infer
finite accuracy from continuum improvement. Eight stratified cases cannot
validate population coverage. Any later accuracy campaign still requires
fresh calibration, a stated regularity/remainder scope, 95% two-sided confidence
containment and 90% power sizing before its data. Larger event counts alone
cannot remove remaining second-derivative bias.

Before generation run the existing generator and fifth-degree inverse checks,
the curvature-square full complex-step derivative/symmetry/flat checks and
the fixed constant-wave physical scaling check. Replay an earlier uncorrected
fifth-degree archive to verify the default path. Preserve clean generating SHA,
source hashes, all six-field archives, separate truth and every full fit/failure.
After completion audit hashes and replay at least one corrected case through
the pure estimator without generator/truth imports.

The original 256-case coverage uses an earlier checkout while this experiment
runs from an isolated worktree. Do not merge estimator changes into that active
checkout until its frozen run has finished.

```powershell
rtk proxy E:/lab/causal-spacetime/.venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_rotating_pilot.py --model-order 5 --square-correction --data-seed 2026100500 --protocol square_pilot_protocol.md --output E:/lab/causal-spacetime/outputs/intrinsic_weyl_square_pilot_v1
```
