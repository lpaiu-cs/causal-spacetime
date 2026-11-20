# Calibration of the nominal covariance at high event counts

Frozen before any data of this stage were generated. Runner:
`experiments/positive_control/intrinsic_weyl_calibration.py`.

## 1. Question

The estimator reports the local linearized covariance `S` of the ten Weyl
coefficients under the fitted jet, built from the shared-membership covariance
of the interval clocks. Its coverage was tested at `N = 2.6e5` on one
population (509 of 512). At larger `N` the evidence is thin:

- the benchmark and the Schwarzschild test had two flat controls each, whose
  detection statistics averaged 11.2 and 10.2 over eight values each;
- one engineering diagnostic (Schwarzschild, `T = 4M`, sixteen blocks of
  `2^20` events with one landmark selection) gave a mean truth statistic of
  13.4 against the nominal 10, a hint of variance understatement by about a
  third at `N ≈ 1e6`.

This stage asks whether `S` is calibrated at `N = 2^20` and `N = 2^24`, for the
published fifth-degree jet and for the sixth-degree jet added after the
benchmark (commit `b5e14e1`, which adds the 42 vacuum second curvature
derivatives and removes the half-cycle continuum bias at a large variance cost;
see `second_derivative_results.md` once written).

## 2. Design

- **Geometry.** Flat spacetime from the plane-wave generator with zero
  profile, exactly as the benchmark controls. The true tensor is zero, so the
  truth statistic equals the detection statistic `Q = ĉᵀ S⁻¹ ĉ` and needs no
  reflection convention, and no jet has model bias.
- **Replicates.** `R = 64` independent replicates with seeds
  `2026112100 + r`, `r = 0..63`. Each draws its own landmark selection from
  4096 events with the unchanged nested design, and one measurement stream
  whose prefixes of `2^20` and `2^24` events are the two levels.
- **Fits.** Order-only pure fits with the fifth-degree jet plus
  curvature-square terms (as published) and with the sixth-degree jet plus
  curvature-square terms, from the same clocks. Exact integer second moments
  on GPU, as in the benchmark worker. The pure fits are written before any
  generator module is imported.
- **Secondary population (descriptive, not gated).** 16 replicates of the
  constant plane wave with the benchmark crest amplitude `A0 = π² h0 / 2`,
  `h0 = 0.01`, seeds `2026113100 + r`. Its truth statistic uses the reflection
  candidate with the smaller Mahalanobis distance, which biases it slightly
  downward; it is reported, not judged.

## 3. Statistic, interval and verdict

For each jet and level, `κ̂ = Σ_r Q_r / (10 R)`. If `Q_r = κ χ²₁₀`
independently, `Σ Q_r ~ κ χ²_{10R}`, and the 95% two-sided interval for `κ` is
`[Σ Q / χ²_{10R}(0.975), Σ Q / χ²_{10R}(0.025)]` (house default: 95%,
two-sided). For `R = 64` this is `[0.8989 κ̂, 1.1193 κ̂]`; under the model its
coverage is 94.8% in 20,000 simulations at `κ = 1.3`.

Two outputs per jet and level, both fixed now:

1. **Tolerance verdict.** `CALIBRATED` if the interval lies inside
   `[0.8, 1.25]`; otherwise `NOT_CALIBRATED`. The tolerance is on the variance
   scale, i.e. radii correct to about ±11%.
2. **Direction.** `UNDERSTATED` if the interval lies above 1, `OVERSTATED` if
   below 1, otherwise `CONSISTENT_WITH_ONE`.

Consequence, fixed now: for any jet and level that is `NOT_CALIBRATED`, the
manuscript scales its nominal radii at comparable `N` by the square root of the
interval's upper bound and says so. A `CALIBRATED` but `UNDERSTATED` outcome is
reported with `κ̂` and its interval.

Descriptive companions, not verdicts: a percentile bootstrap interval for the
mean of `Q` (10,000 resamples, generator seed 1); the number of `Q > q*`, with
`q* = χ²₁₀(0.999) = 29.5883`; the Kolmogorov–Smirnov distance of `Q/κ̂` to
`χ²₁₀`.

## 4. Power (before data)

Under the scaled-χ² model with `R = 64`, the probability of `CALIBRATED` is
0.958 at `κ = 1`, 0.614 at `κ = 1.1`, 0.105 at `κ = 1.2` and 0.001 at
`κ = 1.34`. The repository's convention of 90% power is met at `κ = 1`; a
variance understatement of the size hinted at would not be mistaken for
calibration. `R = 48` would give only 0.865 at `κ = 1`.

## 5. Predictions (before data)

- Fifth-degree jet: `κ̂ ∈ [0.95, 1.35]` at `2^20`, `κ̂ ∈ [0.95, 1.20]` at
  `2^24`. If the engineering hint is real, `NOT_CALIBRATED` and `UNDERSTATED`
  at `2^20`.
- Sixth-degree jet: at least as large as the fifth-degree `κ̂` at each level,
  because its 42 extra nuisance parameters make the linearization worse at low
  signal; `κ̂ ∈ [0.95, 1.5]`.
- Exceedances of `q*`: zero to two per jet and level out of 64.

## 6. Execution

From a clean worktree after this protocol and the runner are committed. Output
in `outputs/intrinsic_weyl_calibration_v1` (archives kept, about 50 GB). Six
worker processes with two BLAS threads each, because one fit at `2^24` peaks
near 13 GiB of memory. Nothing is rerun, reseeded or
extended; a failed replicate counts as a failure of its row and is reported,
and the verdict is computed over completed replicates only if at least 60 of
64 complete, otherwise `INCOMPLETE`.

## 7. Scope

Calibration of the nominal covariance in flat spacetime at two event counts,
for one landmark design family and two jets. It does not establish coverage
for curved fields, where model bias adds to the variance, nor at event counts
outside the two levels.
