# Shared-reference differences: evidence and remaining resolution gap

2026-09-24. The [analytic construction](local_null_challenge.md) gives the
proposed sufficient condition `N lambda^3/log N -> infinity`, for known
`0<q<=1` and this exact circular vacuum family. It uses correlated moment
increments and all causal probe pairs. It does **not** identify the sharp
boundary: the impossibility result still covers only each fixed exponent
below 4/3. Both proofs await independent researcher audit.

## Frozen geometry checks

Source and design were committed at
`7dec69a787164b3f9b0c9ee70754879f6b15fc3d` before the clean-tree run.
All 12 constructed pairs had positive exact wave causal margin. Ten had
positive population-coordinate decrease; the two r=64 failures remain in
the record. These are constructed population witnesses, **not** 10/12
classification accuracy or a random-sampling success probability.

| r=1/lambda | q=0.25: decrease/lambda^2 | q=1: decrease/lambda^2 |
|---:|---:|---:|
| 64 | -0.271854 | -0.179415 |
| 128 | 0.016543 | 0.092361 |
| 256 | 0.075523 | 0.147062 |
| 512 | 0.090478 | 0.160516 |
| 1024 | 0.092035 | 0.162671 |
| 2048 | 0.091750 | 0.162087 |

The q=1,r=2048 population decrease is `3.86446596945e-8`, and the wave
causal margin is `3.86450424167e-8`. Its three nonnegative unconditional
moment increments divided by lambda are `1.035476, 0.386949, 0.021323`.
Tighter integration changed the stored moments and decrease by zero at
double precision. The predeclared distinct five-state Rosen ODE and direct
action quadrature gave `3.86446692424e-8`: absolute difference `9.55e-15`,
relative difference `2.47e-7`, below the frozen 1% criterion. The check used
the same author and is neither independent researcher review nor a rigorous
numerical enclosure. Its reused helper labels its output "post-run", but
this geometry check was specified in the frozen design before the run.

Raw SHA256: `108d77a2ccf2fa720223b7e86a7cbcc8ad1a916107c945fb29aedc9a84377c83`.
The primary CPU run took 17.60 seconds; this excludes the separate
five-state ODE cross-check.

## The original synthetic decision control failed

The frozen q=1,r=2048 control supplied rounded population moments with a
**synthetic** reference count `10^30`. It returned the difference interval
`[-8.70523402561e-8, +9.76304409615e-9]`. Since the upper endpoint was
positive, the observer did **not** reject in this control. No finite sample
of this size was generated, and no detection-power result follows.

A subsequent arithmetic inspection found that applying the inverse
calibration once over the whole moment rectangle inflated the Hessian
absolute bound to about `1.50e12` in its largest entry. The calibration has
cancellation among interval-valued factors; enclosing every factor
separately remains safe but loses precision. This diagnoses conservative
interval dependency, not failure of the negative population signal.

## Pre-execution plan for one post-run enclosure check

This revision is explicitly **post-run**. Keep the raw geometry rows,
original failed interval, synthetic count, moments, confidence construction,
rounding, and 128-bit precision unchanged. Split the allowed `c=a/p^2`
interval into exactly eight equal bands. On each band bound the same first
and second derivatives, using the exact identities `c_p=-2c/p` and
`c_pp=6c/p^2`; then take the hull of the eight entrywise enclosures.
Every true moment point lies in a band, so coverage is preserved.

Run the existing independent finite-difference derivative self-check.
Commit this revision and plan, then run `--recheck-enclosure` once from a
clean tree on the same saved q=1,r=2048 moments. Save its own stamp, source
hashes, original raw hash, and outcome in a separate JSON. The criterion is
that the new difference interval's upper endpoint is negative. If it fails,
retain the failure without changing the sample count, band count, or pair.
Even a pass only validates this algebraic decision path for synthetic
moments; it does not retrospectively convert the frozen control to a pass.

The unresolved scientific work remains the optimal full-order boundary and
an observable attaining it, with a proof audit and finite-sample validation.
