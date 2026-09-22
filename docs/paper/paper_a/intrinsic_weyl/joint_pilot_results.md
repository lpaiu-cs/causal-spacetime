# Actual finite-order reconstruction: first pilot

2026-09-22. The full objective remains open. This pilot supplies a previously
missing piece: an executable estimator whose landmarks, observer frame,
geometry and weighting are inferred from causal-order data. Its output is a
quantitative signed tidal spectrum, not a classification trained on known
wave amplitudes. It is still a development result, with no accuracy verdict.

## Provenance and blind input

The prospective [protocol](joint_pilot_protocol.md) and implementation were
committed at `0bb682d` before sampling. All six fixed cases ran from that clean
commit; none was dropped or replaced. Each used 1024 independent selection
events and 262144 measurement events. The selected interior landmark counts
were 58, 46, 64, 64, 64, 58. Their minimum observed clock counts were respectively
99, 131, 69, 39, 104, 37. No ambiguous causal query needed escalation.

[The complete result](joint_pilot_result.json) preserves every fit, covariance,
optimizer diagnostic and archive hash. The packed order archives remain in
`outputs/intrinsic_weyl_joint_pilot_v1/case_0.npz` through `case_5.npz` and are
reproducible from the bound seeds and generating commit.

Each archive contains exactly `selection_order`, `selected_ids`, `after`,
`before`, `measurement_n`. It contains no coordinates, profile, generator
amplitude, physical density, true tensor, or paired flat sample. The reader
recomputes landmark selection from the saved selection order and verifies IDs.

After the run, the reader was moved into `intrinsic_weyl_joint.py` to make the
boundary independently executable. The estimator imports only standard-library
modules, NumPy, SciPy's optimizer, and the algebraic Weyl basis. It does not
import the generator. In a fresh process, replaying case 4 matched the original
ten coefficients to `1e-11` and exactly reproduced its selected IDs. This
refactoring does not replace the original frozen implementation or result.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_joint.py outputs/intrinsic_weyl_joint_pilot_v1/case_4.npz
```

## Observed finite-order estimates

The target is the electric Weyl spectrum in squared-tip-proper-time units,
relative to the timelike geodesic joining the diamond tips. For these generator
populations it is `(-A,0,A)`. The estimator fits a general ten-component Weyl
tensor; it is not supplied with this spectral form.

| A, evaluator only | Repetition | Estimated electric eigenvalues |
|---|---:|---|
| 0 | 0 | (-0.15002, 0.03712, 0.11291) |
| 0 | 1 | (-0.15194, -0.05557, 0.20751) |
| 0.25 | 0 | (-0.16300, -0.13303, 0.29603) |
| 0.25 | 1 | (-0.23596, -0.10948, 0.34545) |
| 0.75 | 0 | (-0.64051, -0.02833, 0.66884) |
| 0.75 | 1 | (-0.72811, -0.05101, 0.77912) |

Every optimizer reported success, every coordinate Jacobian had full rank,
and the ten curvature columns retained rank ten after removing coordinate
directions. No inferred time coordinate exceeded the tip interval. These are
mechanical checks, not statements that all curvature components are accurate.
In particular, the magnetic spectra for the two A=0.75 cases were
`(-1.30184, 0.12765, 1.17420)` and `(-1.05794, -0.40140, 1.45934)`, against
the same expected spectrum `(-0.75,0,0.75)`. Their errors are substantially
larger than the electric errors. The pilot has **not** demonstrated recovery
of the full Weyl tensor.
The six complete generation-and-fit times were about 1.8--10.8 seconds each
on this machine; no sample-size extrapolation is used in the table.

The first two spectra also show why zero is not recovered exactly from a
finite flat order. Eigenvalue error bars based on nondegenerate perturbation
theory are invalid at the triply degenerate flat spectrum. The nominal full-
tensor Wald statistics were 8.37 and 9.27 for the flat cases, versus the
chi-square(10) 95% reference value 18.31. This is a diagnostic under the local
normal model, not an independently validated type-I-error claim.

## Bias and uncertainty audit

The separate [continuum audit](joint_pilot_bias_audit.json) uses generator truth
to evaluate mean interval volumes on the SAME selected landmarks. It retains
the original data-derived covariance weights and repeats the fit with those
continuum clocks. This evaluator-side diagnostic is deliberately labelled as
an oracle calculation; it is not counted as another blind recovery or an
independent confirmation. Two quadrature orders agreed to about `3e-15` on the
used volume entries; this is a refinement check, not a certified enclosure.

The largest absolute electric-eigenvalue bias on the four curved configurations
was `0.002483`, occurring at A=0.75. The two A=0.75 continuum spectra were

```text
(-0.747919, -0.001912, 0.749831)
(-0.749191,  0.001674, 0.747517).
```

In comparison, nominal finite-sample standard errors of their three eigenvalues
were approximately 0.061--0.069. Sampling noise, rather than the measured
continuum truncation bias on these cases, is the larger remaining error.
This small bias must not be promoted to a uniform bound over other landmark
sets, profiles, or spacetimes.

Magnetic continuum errors are larger: the maximum absolute eigenvalue bias
in these four curved configurations is about 0.0201. That is still much
smaller than the magnetic errors in the finite pilot, but the electric bias
number must not be substituted for a full-tensor bias bound.

For example, the A=0.75 positive-eigenvalue marginal normal intervals were
approximately `[0.5480, 0.7897]` and `[0.6592, 0.8991]`. These are nominal 95%
two-sided delta-method intervals, not a simultaneous band or validated
coverage result. No frozen accuracy gate was applied to the pilot, and we
must not invent one after seeing these intervals.

## What remains before the full objective is established

A new design obstruction explains why increasing sample size is not the only
next step. Write two events as `x=(t,r)` and `y=(s,z)`, and let `B` be a
symmetric tracefree magnetic tensor. With `F(r)=r cross (B r)`, the purely
magnetic quartic is exactly

```text
C_B(x,y,x,y) = 2 [t r dot F(z) + s z dot F(r)].
```

If all interior landmarks lie on exactly two time slices `t` and `s`, then
the spatial shifts `delta r = s F(r)/3` on the first slice and
`delta z = t F(z)/3` on the second reproduce the curvature correction
`-C_B(x,y,x,y)/3` in every cross-slice squared interval. They preserve the
tip clocks because `r dot F(r)=0`. Within each slice there are no timelike
pair clocks. Thus **all five magnetic components are absorbed by unknown
position changes at first order**. This is not a claim of indistinguishability
at every higher order in curvature.

The symbolic identity is checked exactly in
`check_two_slice_obstruction.py`. Its rational configurations have projected
Weyl rank five on two slices and ten on three slices, as recorded in
[the certificate](two_slice_obstruction.json). Rational ranks are certified by
exact finite-field lower bounds meeting analytic rotation/null-mode upper
bounds; no floating rank tolerance is involved. A slow unrestricted rational
elimination was stopped and replaced by this equivalent bounded certificate.

The pilot uses time bands, not exact slices, so this does not make its rank-ten
result contradictory. It identifies a nearby degeneracy and a concrete design
lever: include a third time band to constrain the magnetic deformation. The
three-slice rank certificate alone does not guarantee a lower finite-sample
error; that comparison still needs independent sampled orders.

1. Freeze and execute an independently powered accuracy study. Interval-based
   verdicts must use the selected coverage and sidedness consistently in
   sizing and evaluation. Validate uncertainty rather than relying on point
   estimates or successful optimization; retain failed selections and fits.
2. Test genuinely varying wave profiles and more general vacuum curvature.
   The fitter accepts all ten components, but these six generator populations
   are constant plane waves. Rank ten does not prove finite-sample accuracy
   for the magnetic components or all vacuum algebraic types.
3. Treat locality/scale choice and truncation bias as part of the instrument.
   The global/singular counterexamples from earlier work remain in force.
   Unknown nonuniform sampling and nonvacuum Ricci contamination remain outside
   the present estimator's justified domain.

The previous two-billion-event cost projection and this 262144-event table
refer to different designs AND different statements: one was a powered
single-projection accuracy price, the other contains pilot point estimates.
They do not establish a several-thousand-fold efficiency improvement.
The actual advance is that the external landmark geometry and flat-control
inputs have now been removed in an executed finite-order reconstruction.
