# Noise-cost decision and covariance weighting

## Completed first diagnostic

`noise_cost_diagnostic.json` was produced from clean commit `5dc22a9` with
the prospective protocol in that commit's README. It has four independent
262144-event replicates. This is a local flat-limit cost projection, not a
successful recovery experiment or a sample-complexity lower bound.

| Scale | Predicted relative tensor RMS at N=1e6 | Optimistic event count for one projection |
|---|---:|---:|
| 0.04 | 221--259 | 5.32e11--6.86e11 |
| 0.02 | 886--1037 | 8.51e12--1.10e13 |

Ranges are the four observed replicate values, **not confidence intervals**.
The projection price assumes zero bias and normal errors, 95% two-sided
intervals contained in a relative +/-25% band, and 90% probability of that
containment. A full-tensor gate and truncation bias would require more work.
The small continuum error has therefore not established an affordable method.

Two of four actual noisy flat fits failed the reference Gram signature check.
The other two returned large spurious tensors. Minimum interval counts were
only 1--3. Those fits are retained as failures of this development design;
they must not be filtered into a favorable accuracy statement. The asymptotic
linearization is not a calibrated error model at these counts.

The point of the next diagnostic is narrow: determine whether the very large
cost is primarily an avoidable consequence of ordinary least squares ignoring
correlated residual noise. No parameter of the first diagnostic is changed.

## Prospective second diagnostic

Freeze before its validation seeds are used:

1. Keep the SAME 33-event configuration, intervals, normalization, and linear
   curvature response. Reuse the first diagnostic's four seeds as development
   data to estimate the full 224-by-224 residual covariance. Their reuse is
   explicit; they are not independent validation.
2. Use generalized least squares. Standardize the residual covariance by its
   diagonal, eigendecompose it, and keep eigenvalues greater than `1e-9` times
   the largest. Record the spectrum and retained rank. Require rank ten for
   the whitened curvature design; otherwise report failure without retuning.
3. Fix the weights before validation. Validate the projected noise cost on
   four NEW independent seeds `2026092210..2026092213`, each with `262144`
   events in chunks of `16384`. No extension or seed removal is permitted.
4. Compare ordinary and generalized least squares on each validation sample,
   using the same influence covariance, scales, and power price as before.
   The comparison is a development noise projection, not a preregistered
   accuracy verdict. Report all four price ratios and both absolute prices.
5. Check that both linear estimators return the identity response to the ten
   basis signals and annihilate common count scaling. Bind protocol, sources
   and result to the clean pre-run commit. Do not overwrite the first result.

This diagnostic still knows the limiting landmark configuration when deriving
the influence matrices. Thus even a large improvement would not demonstrate
automatic landmark selection, finite curved recovery, valid confidence
coverage, or an immediately deployable estimator. The estimator of the active
goal must eventually infer its geometry and covariance from order data alone.
The calculation below is a cost/design decision before spending resources on
that implementation.
