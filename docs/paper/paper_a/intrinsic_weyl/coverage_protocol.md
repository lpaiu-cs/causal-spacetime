# Frozen independent end-to-end coverage validation

2026-09-22. Commit code and this protocol before generating any of the new
validation orders below. The preceding goal turn made progress: it repaired
central retention and measured all-component finite and continuum errors.
Those results are development data. No existing point error is re-scored as
an accuracy success.

## Estimand and internal frame

Recover the ten components of dimensionless Weyl curvature in the outer-tip
observer's normal frame, modulo an arbitrary spatial orientation. The first
selected landmark vector defines the positive x axis; the second defines
the positive y half-plane. The remaining z reflection is explicitly
quotiented. Do NOT orient z with the third landmark: one existing pilot has
its third z coordinate about 0.00064, making that choice unstable. This is
an axis convention, not a change to landmark selection or measured clocks.

Normalize the estimator's x/y signs, transforming both the coefficients and
their covariance. For a spatial diagonal transformation (sx,sy,sz), the
electric coefficient signs are `(1,1,sx*sy,sx*sz,sy*sz)`; magnetic signs also
include `sx*sy*sz`. A full four-index tensor transformation independently
checks this parity rule. The distance used in evaluation is

```text
d(c,t) = min( ||c-t||_G, ||c-Jt||_G ),
||c||_G^2 = ||E(c)||_F^2 + ||B(c)||_F^2,
J = the coefficient transform for diag(1,1,-1).
```

This compares all ten signed components and relative E/B directions. It
does not collapse the target to separate spectra. The z-reflected tensors
are representations of the same orientation-free internal tensor.

## Nominal confidence construction

Keep the connected selector, joint quartic fitter, covariance eigenvalue
cutoff, and optimizer unchanged. Let Sigma be its nominal full-component
covariance, G the tidal norm matrix and lambda the largest eigenvalue of
`G^(1/2) Sigma G^(1/2)`. The statistical radius is
`r_stat = sqrt(chi2_10(0.95) * lambda)`.

The bound is a 95% simultaneous Gaussian ball in ten components, before
quotienting the residual reflection. It is not ten uncorrected marginal
intervals. The covariance is estimated, and the model is nonlinear, so
Gaussian coverage is a hypothesis to validate, not an established fact.

Freeze a relative asymptotic model-bias allowance **delta=0.02**. This was
chosen after the previous continuum diagnostics; it is NOT a proved uniform
bound on arbitrary spacetimes. For an assumed bias no greater than
`delta*||true||_G`, the reported radius is

```text
r = (r_stat + delta*||estimate||_G)/(1-delta).
```

The triangle inequality gives this expression without feeding the true norm
to the estimator: on the statistical event,
`d(estimate,true) <= r_stat + delta*||true||`, and
`||true|| <= ||estimate|| + d(estimate,true)`.

The evaluator checks the bias assumption on EACH selected configuration by
fitting continuum interval volumes with that case's actual clock mask and
covariance weighting. It uses a true normal-frame construction fixed by the
same first two landmarks, then the residual-reflection distance. A violated
bias condition is an end-to-end validation failure, even if the ball covers
the truth. For the flat population, numerical zero means absolute bias at
most `1e-10`; no relative flat error is defined. Continuum quadrature must
stabilize to `1e-10` relative between orders 32 and 64.

## Independent population, budget and primary gate

- Exactly **256** cases. Amplitudes are iid choices from `(0,0.25,0.75)`
  with equal probabilities, using `default_rng(2026092250)`. Case i uses
  data seed `2026092300+i`, with independent selection/measurement streams
  `[seed,0]` and `[seed,1]`. None of these seeds was used in the pilots.
- Every case uses 4096 selection events and 262144 measurement events.
  The event budget is 68,157,440 plus the supplied tips. No sample increase,
  seed replacement, early stopping, or excluded failure after results.
- Preserve a separate evaluator-only truth archive. The estimator loads
  ONLY the order archive. Its confidence function uses only the returned
  coefficients, inferred positions and covariance. Neither function reads
  a profile, true coordinates, amplitude, or a paired flat observation.
- A case counts as valid-and-covered only if selection, generation, fitting,
  rank and boundary checks, continuum audit and bias condition succeed,
  and `d(estimate,true) <= r`. Failures remain in the denominator.
- Compute the **95% two-sided exact Clopper--Pearson interval** for the
  valid-and-covered probability of this frozen iid mixture. The primary
  gate requires its LOWER endpoint to be at least **0.95**. With n=256,
  this requires **251 or more** valid-and-covered cases.
- Exact power is **0.95471764 at true end-to-end probability 0.99**.
  At 0.98 the power is only **0.59470691**. The optimistic planning
  alternative is stated explicitly because the simultaneous norm ball is
  conservative under Gaussian sampling. Ninety-percent power is NOT
  claimed at 0.98 or at the 0.95 boundary. This validation can fail.
- Save per-case full outputs, order/truth hashes, runtimes, failures,
  source hashes and a reachable pre-data commit. Do not rewrite pilot data.
  No pure coverage claim is made by discarding invalid models or fit failures.

An outcome from this mixture is not a separate coverage guarantee for every
amplitude or every possible landmark arrangement. Preserve stratum counts
and outcomes descriptively. The study validates coverage/availability at
its stated event count; it does not certify 25% accuracy or varying profiles.

## Price of the later accuracy study: not authorization to run it yet

The six pilot covariances have `lambda*N` at most about 16659.4. A provisional
design envelope `lambda*N <= 0.08*262144 = 20971.52` adds about 26% over that
observed maximum. It is a planning assumption, not a uniform variance bound.
It must be checked on any later accuracy data, without raising its cap.

For a six-case simultaneous accuracy study, allocate confidence error
0.05/6 and power failure 0.1/6 per case. A sufficient Gaussian price for the
reported ball to lie wholly within a relative 25% target ball, allowing the
same relative bias delta, is

```text
N >= V * (sqrt(chi2_10(1-.05/6)) + sqrt(chi2_10(1-.1/6)))^2
         / ( ((1-delta)*.25 - 2*delta) * ||true||_G )^2.
```

The TWO appearances of bias are intentional: it affects the center's error
and the confidence-set radius. This is a sufficient 90% JOINT power bound
under the stated Gaussian/variance/bias assumptions, not a break-even price.
For the plane-wave norms `2A`, it predicts approximately **181.3 million**
measurement events per weak-wave case (A=0.25) and **20.15 million** per
strong-wave case (A=0.75). These prices concern ALL components and a joint
interval-containment gate; they are not comparable to the earlier one-
projection price without those qualifications.

Do not launch that larger campaign from this protocol. First complete this
coverage validation, inspect actual cost and failure modes, then freeze its
exact population, accuracy targets including the flat case, budget, bias
checks and independent seeds in a separate protocol. The full goal also
still requires finite-order validation beyond constant wave profiles.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_coverage.py --check
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_coverage.py
```
