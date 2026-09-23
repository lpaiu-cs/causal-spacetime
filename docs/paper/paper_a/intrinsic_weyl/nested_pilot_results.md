# Nested landmarks: precision gain with an unresolved coverage warning

The [frozen eight-pair pilot](nested_pilot_protocol.md), generated at
**c59cdb3edd359d4c20d2ffddb35d26210e7f810f**, completed both arms on the
same 4096 selection and 262144 measurement events per case. The candidate
retained every baseline landmark and used **122--131** interior landmarks,
versus 64 in the baseline. Both use the identical corrected fifth-degree
model, with ten free Weyl and 24 free gradient components.

Actual full-tensor error decreased in **all eight pairs**, including all six
curved cases. The median paired error ratio is **0.5628**: a 43.7% reduction.
The median largest-variance ratio is **0.1837** and the median magnetic-
variance ratio is **0.2036**. This is a substantial gain in these new data,
unlike the earlier wider-center or curvature-square-only comparisons.
These are descriptive pilot findings, not population precision guarantees.

| A | w | Baseline relative error | Nested relative error | Nested/baseline largest variance | Nested continuum bias |
|---:|---:|---:|---:|---:|---:|
| .25 | 0 | 96.60% | 44.00% | 0.2635 | 0.0041% |
| .75 | 0 | 29.70% | 7.23% | 0.1632 | 0.0425% |
| .25 | .5 | 119.30% | 74.25% | 0.1885 | 1.2362% |
| .75 | .5 | 38.70% | 18.20% | 0.2157 | 1.2244% |
| .25 | -.5 | 94.47% | 88.35% | 0.1440 | 1.2801% |
| .75 | -.5 | 26.75% | 22.88% | 0.2096 | 1.2379% |

The flat cases' absolute errors decrease from 0.67287 to 0.33864 and from
0.35830 to 0.24102. Relative error at zero remains undefined. All errors
include all ten tidal components in the shared internal frame modulo the
remaining reflection, not just eigenvalues or a fitted wave amplitude.

## Coverage warning and limit of the result

Both arms satisfy the unchanged 2% continuum-bias condition in every case.
The baseline nominal balls contain truth in eight cases, but the nested
balls contain truth in **seven of eight**. Case **6** (A=.25,w=-.5) has
absolute error **0.44176**, exceeding radius **0.41321**. Its point error
improves by only 6.5% even though its nominal maximum variance falls by 85.6%.
Keep this noncoverage; do not enlarge its radius or remove it after the fact.
Eight cases do not establish or refute 95% population coverage.

The three strong cases have point errors below 25%, but their confidence
radii are 0.459--0.538, exceeding the 0.375 relative-accuracy target radius
at truth norm 1.5. Thus **no 25% confidence-contained accuracy claim is made**.
Weak-case errors remain large. The objective is not completed by this pilot.

The data support retaining more relational information as a precision lever.
They do not yet support the new uncertainty construction's reliability.
The next step is independent coverage validation of this unchanged candidate,
including zero-clock, bias and fit failures. Any later accuracy study still
needs explicit regularity/bias scope, independent high-count calibration and
90% power sizing before its data. No unrestricted-frequency or global result
follows, and the old fifth-degree failure remains unchanged.

## Conditional price of a later accuracy test

Using only these pilot covariances, reuse the sufficient Gaussian planning
rule in [the earlier coverage protocol](coverage_protocol.md#price-of-the-later-accuracy-study-not-authorization-to-run-it-yet).
For six nonzero cases, allocate confidence error .05/6 and power failure .1/6
per case. With delta=.02 and V=N*lambda_max, the sufficient measurement count
for 25% confidence containment and at least 90% joint power is

```text
N >= V * (sqrt(chi2_10(1-.05/6)) + sqrt(chi2_10(1-.1/6)))^2
         / (.205 * ||true||_G)^2,
.205 = (1-delta)*.25 - 2*delta.
```

This includes center error and interval radius; it is not the smaller
break-even count. Applying it separately to each observed configuration gives:

| Case | A | w | Pilot V | Conditional measurement count | Nested/baseline count |
|---:|---:|---:|---:|---:|---:|
| 2 | .25 | 0 | 4258.32 | 36,819,113 | .2635 |
| 3 | .75 | 0 | 2510.81 | 2,412,159 | .1632 |
| 4 | .25 | .5 | 4370.47 | 37,788,829 | .1885 |
| 5 | .75 | .5 | 3537.42 | 3,398,436 | .2157 |
| 6 | .25 | -.5 | 2160.95 | 18,684,391 | .1440 |
| 7 | .75 | -.5 | 3213.10 | 3,086,856 | .2096 |

The total is **102,189,784 measurement events**, excluding selection and any
flat-case accuracy target. These are conditional extrapolations, not a frozen
future budget or a verified power result. They assume valid Gaussian
uncertainty, the 2% bias condition and stable V at higher counts; none is
established for new high-count cases or new landmark configurations. In
particular, case 6's observed noncoverage remains unresolved. The ten-component
confidence ball would use the stricter .05/6 error allocation in that later
study. No high-count accuracy run has been started from this calculation.

The apparent saving concerns sampled events. It does not establish lower CPU
cost: the nested estimator uses about four times as many interval clocks, and
the pure-fit replay below is substantially slower. High-count calibration and
measured computational feasibility are still necessary; weak cases cannot be
dropped to turn the cheaper strong-case counts into the full result.

## Resource cost and audit

The total sample is **2,129,920 events plus 16 supplied tips**; events were
shared, not doubled. The nested fits use **3727--4185** interval clocks versus
**948--1089**. Querying both directional membership bits requires roughly
65.0--69.7 million landmark/event predicates per nested case versus 34.6
million for the baseline; these counts are not CPU action-evaluation counts.
The actual paired generator queries only the union once. Compressed order
archives are approximately 4.24--4.49 MB versus 2.58--2.60 MB.

Computational cost also increases. One sequential pure-estimator replay of
the same case 5 took **4.37 seconds** for the baseline and **60.03 seconds**
for the nested fit on this host. No concurrent research job was running in
that replay, but it is one descriptive timing, not a controlled benchmark.
The sum of pilot per-pair times, including generation, both finite fits and
both continuum evaluations under variable load, is **1116.27 seconds**.

All 32 fits have Weyl rank ten, full nuisance rank for their actual landmark
count, successful optimizers and no boundary flags. The audit verified all
11 historical source hashes, all 48 artifact hashes, all eight shared-order
and shared-measurement bit comparisons, preserved frame prefixes, and all
saved error/variance ratios. Separate pure-estimator replays of **both**
case-5 arms match every saved fit field exactly without generator/profile/
truth imports. The additional coordinates are inferred by each estimator.

[Recorded results and audit](nested_pilot_summary.json) retain every case,
including case 6. Raw files are in `outputs/intrinsic_weyl_nested_pilot_v1`;
the raw summary SHA256 is
`981b88fee30c04ec5e3ca223ffbe5f60d182f60968c4f7ec196a4c36e480ecb7`.
The pure-replay timing record is `pure_replay_audit.json` in that folder.

## Exact reuse of repeated membership patterns

A subsequent engineering check, committed at
**63e1b6b0428905cc6b7532b2786f0eb6bef5bb66**, removes repeated arithmetic
without changing the measured moments. If z_e is the vector of interval
indicators for event e, and identical after/before landmark signatures have
representative z_r and multiplicity w_r, then

```text
sum_e z_e z_e^T = sum_r w_r z_r z_r^T.
```

The denominator remains the original event count N, not the number of distinct
signatures. The implementation multiplies by integer w_r, without square
roots; products and partial sums remain exactly representable in float64 for
N below 2^53. Interval counts, shared-membership covariance, selection, fitted
model and confidence construction are unchanged. No coordinates or generator
information are used. The existing inverse self-check now verifies this
identity against direct event sums, including repeated/singleton patterns and
landmark counts not divisible by a byte boundary.

Old case 5 has **107,961** distinct signatures among 262,144 events. The
prototype reproduces every covariance entry exactly. The integrated path then
reproduces **every saved fit and confidence field exactly in old cases 5 and
6**, without importing generator or truth modules. In particular, case 6
still misses truth: error 0.44176075 exceeds the unchanged radius 0.41320850.
Computational reuse does not repair this statistical failure.

During the simultaneous coverage campaign, one profiled baseline replay of
case 5 took 300.71 seconds, including 96.43 seconds for reading clocks and
constructing covariance. The grouped-moment prototype took 32.12 seconds;
integrated replays took 179.08 seconds for case 5 and 117.39 seconds for case
6. These unequal workloads under variable concurrent load are descriptive
timings, **not a controlled speedup estimate or a high-count feasibility
result**. The grouping requires additional memory linear in N times the
number of landmarks, plus sorting temporaries; its benefit depends on repeated
signatures. It does not remove fitting or dense covariance factorization.

## Fixed thread-setting replay before high-count work

2026-09-23, before calculation: replay the existing nested cases 0--7 in
separate processes with `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`.
Bracket these eight replays with two default-environment replays of case 5.
Each process receives only the same order archive and keeps the nonlinear
fifth-degree-plus-square estimator. No new events, landmarks, model terms,
selection changes or confidence changes are permitted. The ongoing 512-case
campaign keeps its original sources, processes and unset thread variables.

Check every returned fit and confidence field against the saved pilot at
absolute and relative tolerance 1e-8, with identical discrete fields and
rank/success flags. Record optimizer function-evaluation counts separately:
they describe computation, not the estimator output. Retain failures and all
numeric discrepancies, and recompute the original ball-containment flag after
the pure process returns. Case 6's missed ball may not disappear from the
record. Store wall time, process CPU time and clock-building/fitting profile
times, together with source and input hashes. Preserve the complete outputs.

This is a numerical-portability and execution-cost check at N=262144, not
high-count calibration or an accuracy test. Concurrent campaign load is
variable; two bracketing timings cannot establish a controlled or universal
speedup. A successful replay could justify considering a thread limit for
a later, separately frozen run, not changing the active run. OpenBLAS checks
these variables on startup; see its [runtime documentation](https://www.openmathlib.org/OpenBLAS/docs/runtime_variables/).

### Completed replay and its limits

The plan was committed as `75b9356`; the ten replays ran from clean commit
`85d795a12a4e283ec78e35a6356596c40b0183e7`. Five source files and all replay
inputs were bound by hash. The [compact record](thread_replay_result.json)
preserves all outcomes and hashes the complete per-field comparisons in
`outputs/intrinsic_weyl_thread_replay_v1`. Its raw summary SHA256 is
`46b29b213d4acdcbb1525c6a851d473778e84ed27b0c54993ac29cd023a87522`.

Both default case-5 replays match every original fit and confidence field
exactly. Under the one-thread setting, **six of eight** satisfy the frozen
direct field comparison. Cases **1 and 5 fail it**, and the recorded aggregate
`all_single_replays_match` remains false. Their differences are spatial-frame
signs, not a changed curvature solution:

- Case 1 is related to the original by `(x,y,z) -> (-x,-y,-z)`.
- Case 5 is related by `(x,y,z) -> (x,-y,-z)`.

A separate **post-result** audit chose signs using inferred landmark
coordinates only and transformed Weyl coefficients, the full covariant
gradient, and covariance by their tensor rules. It did not choose a frame
using truth or change the direct comparison rule. All eight fits and balls
then agree within the original 1e-8 tolerance. Across them, the largest
absolute discrepancies are 4.63e-15 in coordinates, 1.65e-13 in Weyl
coefficients, 4.59e-12 in gradient coefficients, and 6.17e-16 in covariance.
The all-pair model values are invariant under the sign transformations to
5.56e-17. All original containment flags are preserved, including the
**case-6 missed ball**. This uses the reflection quotient already defined in
[the original estimand](coverage_protocol.md#estimand-and-internal-frame);
it is not a new way to rescue a statistical failure. Direct byte/numeric
reproducibility and equivalence of the defined physical output are distinct
claims, and only the latter holds across all eight one-thread replays.

For the same case-5 input, the measured profiled costs were:

| Setting | Wall seconds | Process CPU seconds | Clock building | Nonlinear fit |
|---|---:|---:|---:|---:|
| Default, before | 164.11 | 417.44 | 32.52 | 131.55 |
| One thread | 160.47 | 129.20 | 46.32 | 114.12 |
| Default, after | 303.55 | 615.42 | 65.82 | 237.67 |

CPU time is aggregated over process threads; it is not elapsed time. The
one-thread run used 69--79% less CPU time than the two bracketing default
runs on this loaded host. Wall times and stage costs varied substantially,
and the first default wall time is close to the one-thread time. These
observations do **not** establish a controlled wall-time speedup, a best
worker/thread allocation, or a high-count event price. They support testing
a thread limit in future, separately frozen work without the affine bias
introduced by the rejected count-map shortcut. The active 512-case campaign
was not reconfigured or re-scored. No new events or improved physical
accuracy claim resulted from this engineering replay.

The [recorded audit](nested_pilot_summary.json) includes source and artifact
hashes under `engineering_clock_reuse`. Raw engineering files are retained at
`outputs/weyl_clock_reuse` in the isolated `weyl-curvature-square` worktree.
The active 512-case coverage run continues on its original **f001de6** code;
this implementation is held in the isolated worktree until that run ends.
No independent scientific cases or verdicts were added by the engineering
check.
