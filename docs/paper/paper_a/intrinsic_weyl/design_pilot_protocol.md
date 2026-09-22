# Frozen paired comparison of central landmark spread

2026-09-22. Freeze the source and this protocol at a clean, reachable commit
before generating new data. This is a prospective development comparison of
the main noise lever, not an accuracy or coverage certification. The preceding
256-case coverage result, eight-case rotating pilot and their failures remain
unchanged. No old sample is relabelled as an independent validation sample.

## One design change and a common physical experiment

The previous completed fits had magnetic tidal variance about five times
electric variance at the median. In the three-band geometry, the magnetic
quartic against a central vector z contains a factor quadratic in |z|.
The competing effect is that wider centers may have fewer, smaller intervals
to the outer bands. Therefore a better signal-to-noise ratio is a hypothesis,
not a consequence of the larger radius alone.

Compare two order-only selectors:

- `connected`: the original algorithm, unchanged from `e057d10`.
- `wide`: first run `connected` and keep its past and future landmark IDs,
  in the same order. Expand only the central candidate radius cap from .20
  to .30. Keep `|t|<=.075`, radius>=.075, and at least two selection events
  in each endpoint interval. A candidate needs at least four qualified
  clocks to each retained outer side. Rank by
  `minimum_side_degree * radius_squared`, with stable input-ID tie handling;
  retain enough central points to bring the total to at most 64 interior
  landmarks. Require at least eight centers and degree at least eight in
  the retained graph. Failure is reported without changing the outer set,
  refilling a selection sample, or changing a cutoff.

Here t and radius are the existing proxies calculated from interval
cardinalities, not coordinates supplied by the generator. The outer tips,
first two frame vectors, clock threshold, ten-component quartic model,
optimizer, covariance calculation and nominal confidence construction
remain identical. Neither fit assumes a wave type, amplitude or profile.

Generate one independently labelled selection order per case. Select both
sets before drawing the measurement stream. Query the union of landmark
IDs against the SAME measurement events; save a separate order-only archive
for each design with only the six existing fields. Shared landmark queries
must be byte-identical across the two archives. True coordinates and profile
parameters are written to separate evaluator archives and opened only after
that arm's pure estimator returns. A flat case never supplies an input or
normalization to a curved case.

## Fixed population and budget

The ordered `(amplitude, rotation_rate)` cases are

```text
(0,0), (0,.5), (.25,0), (.75,0),
(.25,.5), (.75,.5), (.25,-.5), (.75,-.5).
```

Use data seed `2026092700+i`, i=0..7, with selection/measurement streams
`default_rng([seed,0])` and `default_rng([seed,1])`. These are new seeds.
Each case has **4096 selection events and 262144 measurement events**, in
measurement chunks of 4096, shared by the two arms. The maximum distinct
sampled-event count is **2,129,920 plus tips**, not twice that number.
If both selectors fail, record that case and generate no unused measurement
stream. If only the wider selector fails, still run the baseline arm.

The cases are deliberately small and stratified. There is no powered
population comparison or accuracy gate in this pilot. No seed replacement,
additional case, event-count increase, parameter change or favorable subset
may be used to rescue its outcome.

## Readouts, bias and failures

For each arm preserve the complete finite-order fit, nominal covariance,
confidence radius, optimizer/boundary/rank diagnostics, and the same-landmark
continuum fit from the existing rotating evaluator. Use the full signed
Weyl norm and shared internal frame modulo residual reflection. The two
arms' true canonical coefficient vectors must agree within 1e-10; otherwise
stop because the paired comparison has changed its target.

Report, per case and descriptively across returned pairs:

- magnetic and electric tidal variance traces, using the Frobenius Gram
  metric, and the largest full-tidal covariance eigenvalue;
- ratios `wide/connected` of those quantities and actual full-tensor error;
- full-tidal finite-order error and continuum model bias for each arm;
- whether each continuum bias satisfies the unchanged relative 2% allowance
  (1e-10 absolute numerical-zero tolerance at A=0), and whether its nominal
  ball contains the true tensor.

A small point error is not accuracy support. A nominal ball can contain
the truth while its bias assumption is violated. The four varying-profile
failures of the prior pilot are not repaired by this new run, and its 2%
allowance is not inflated here. Preserve every empty-clock, selection,
ambiguous-causality, fitting, rank or boundary failure in its assigned arm.
Continue the fixed remaining cases after recorded case/arm failures; no
failure is silently excluded when describing practical availability.

If variance decreases but continuum bias grows, record that tradeoff.
If the wider arrangement loses rank or increases variance, reject the
radius-only mechanism for this tested implementation. A favorable proxy
covariance must also be distinguished from actual smaller finite errors;
eight paired cases cannot establish ensemble coverage or a general rate.
Only after this comparison may a new precision design be frozen, with a
regularity/window bias condition, 95% simultaneous uncertainty and explicit
90% power sizing. No existing variance cap is treated as validated.

## Pre-data and post-data checks

A deterministic exact-flat configuration checks that the outer landmark
IDs/frame are preserved, the chosen central radius really increases,
the projected curvature rank remains ten and the fitted flat tensor is
numerically zero. Replay a completed old archive to check unchanged
`connected` behavior before the new data. The existing rotating generator
has independent geodesic, action and volume checks; it is unchanged here.

Record all source hashes, runtime versions, per-case source seed, artifacts
and hashes, failures and timings. After completion independently check
shared order/query bytes and frame IDs, recalculate comparison quantities,
and replay a `wide` archive through the pure estimator with no truth input.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_design_pilot.py --check
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_design_pilot.py
```
