# Blind curvature-and-gradient pilot: bias reduced, precision unfinished

2026-09-22. The [protocol](jet_pilot_protocol.md) and enlarged inverse model
were committed before fresh data at
`8e56607c1d9a25d4b598a78f4e909051d39c18c0`. All eight fixed cases completed.
Each model inferred coordinates and Weyl curvature from the same saved
order archive. The fifth-degree fit additionally inferred all 24 allowed
vacuum derivative coefficients; no true coordinate, profile, derivative,
polarization, amplitude or paired flat observation was supplied to it.

**The four varying-profile continuum biases fell from 2.14--4.31% in the
paired quartic fits to 1.19--1.83% in the blind fifth-degree fits.** All four
now met the unchanged 2% condition. This is a genuine inverse-model result,
unlike the earlier oracle subtraction. However, these bias numbers come
from the same-landmark continuum audit; they are not the errors of the
finite causal orders. The finite errors remain too large for a 25% accuracy
claim, and nominal confidence coverage has not been independently calibrated
for the enlarged model.

## Complete results

The two methods use identical landmarks, clock masks, physical events,
reference frame and covariance weights. Percent errors use the combined
electric/magnetic norm and compare all signed components modulo reflection.

| Case | A | w | Continuum bias, quartic -> fifth | Finite error, quartic -> fifth | Fifth-degree absolute confidence radius |
|---|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | numerical zero -> zero | .71958 -> .84155 (absolute) | 1.08314 |
| 1 | 0 | .5 | numerical zero -> zero | .69448 -> .61989 (absolute) | 1.23582 |
| 2 | .25 | 0 | .2119% -> .2572% | 100.41% -> 101.30% | 1.27254 |
| 3 | .75 | 0 | 1.0651% -> 1.0060% | 25.77% -> 27.72% | 1.10141 |
| 4 | .25 | .5 | 2.1436% -> 1.5014% | 88.60% -> 105.99% | .95153 |
| 5 | .75 | .5 | 4.3089% -> 1.5787% | 27.25% -> 33.11% | 1.30389 |
| 6 | .25 | -.5 | 2.4774% -> 1.1869% | 79.25% -> 94.69% | 1.06154 |
| 7 | .75 | -.5 | 3.6022% -> 1.8331% | 38.86% -> 36.50% | 1.18481 |

There were no selection, empty-clock, optimizer, boundary or rank failures.
Both the finite and continuum fifth-degree fits had nuisance rank 277
(253 coordinate parameters plus 24 gradients) and projected Weyl rank ten.
All reported nominal balls contained their targets, but eight cases are
not a coverage study. None of these broad balls certifies 25% accuracy.

The extra nuisance parameters raise variance: all eight largest-eigenvalue
ratios are above one. Actual finite error improves in only two of eight
paired cases. This is a bias--variance tradeoff, not a demonstrated
finite-sample precision gain. No radius, noise threshold, frequency, seed
or event count was changed after seeing outcomes. Earlier pilot failures
and the quartic-only coverage verdict remain unchanged.

## Independent checks

Before data, analytic derivatives of all 34 curvature/gradient parameters
and representative coordinate parameters agreed with finite differences
to at most `3.41e-12`. The zero-gradient model reduced to the original
quartic formula. An independently constructed exact-flat interval matrix
returned zero curvature and gradient with full nuisance rank. A complete
old quartic archive replay was exactly equal to its earlier saved fit.

After data, all **32 artifact hashes**, **10 source hashes**, seeds, archive
field schemas, selected IDs, shared order hashes, recorded confidence
calculations and paired evaluations were checked. The pure fifth-degree
replay of case 7 matched curvature, gradient coefficients, inferred
coordinates and covariance exactly (maximum difference **0.0**). No generator
or truth-evaluator module was imported in that replay.

For the same case the ten-component marginal covariance was independently
reconstructed by projecting off position AND gradient nuisance columns and
inverting the projected information. It agreed with the full-Jacobian
pseudoinverse covariance to relative maximum difference **4.96e-14**.
This checks that the reported covariance belongs to the ten Weyl components,
not to the final ten nuisance-gradient parameters. Volume quadrature
refinements were below `3e-15`.

Raw data and full fits remain in `outputs/intrinsic_weyl_jet_pilot_v1`.
`jet_pilot_summary.json` is the complete summary normalized to LF. Raw hash:
`a079dc217c44f96a6f3a32ec5c3d40d52af5795ea30f8c8c0b90140ae3939d56`.
Committed LF hash:
`d9d97d56f3f13968d2a31d75fdadcc4c0a02bf015c3a42e45301bf09d821a439`.
The run used Python 3.14.0, NumPy 2.5.3, SciPy 1.18.1, and 2,129,920 sampled
events plus tips. Sum of per-case runtimes was 254.05 seconds.

## Remaining requirement

The code now actually fits unknown coordinates, C and its vacuum-compatible
first derivatives from finite order data. This repairs a material source
of varying-profile model bias in the tested independent configurations.
The full goal remains open: actual full-tensor precision, uncertainty
coverage for this enlarged model, justified regularity/remainder scope and
a powered accuracy demonstration have not yet been achieved. The remaining
work must address sampling information/cost and these validation requirements;
adding another generator formula alone would not do so.
