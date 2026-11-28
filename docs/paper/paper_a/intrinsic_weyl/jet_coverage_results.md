# Independent fifth-degree coverage: NOT_SUPPORTED

The [frozen 256-case campaign](jet_coverage_protocol.md), generated from
**567d373c6ccf8a15d00debf19d45845065999f78**, completed without changing
its algorithm, population, event counts or failure rules. Its result is

```text
valid and covered: 246 / 256
exact 95% two-sided Clopper-Pearson interval:
    [0.9293377140998706, 0.9811117242805399]
frozen gate (lower bound >= .95): NOT_SUPPORTED
```

The observed proportion is 96.09%, but its lower confidence bound is below
95%. The campaign therefore does not establish the preregistered reliability
claim. It does not prove that the population probability is below 95% either.

All **255 successfully fitted** cases have nominal balls containing truth.
That alone is not the frozen gate: **nine** cases violate the assumed 2%
continuum-bias condition, and **one** cannot be fitted because a selected
interval has zero independent measurement count. All ten remain failures.
The broad nominal balls cannot validate a false bias assumption by happening
to cover the truth. No later correction re-scores this campaign.

## Outcomes and failure mechanism

| A | w | Cases | Fitted | Valid and covered | Median finite error among fitted cases |
|---:|---:|---:|---:|---:|---:|
| 0 | -.5 | 28 | 28 | 28 | 0.47624 absolute |
| 0 | 0 | 32 | 32 | 32 | 0.51845 absolute |
| 0 | .5 | 29 | 29 | 29 | 0.46514 absolute |
| .25 | -.5 | 22 | 22 | 22 | 100.66% relative |
| .25 | 0 | 26 | 25 | 25 | 107.15% relative |
| .25 | .5 | 33 | 33 | 33 | 101.55% relative |
| .75 | -.5 | 33 | 33 | 25 | 31.35% relative |
| .75 | 0 | 28 | 28 | 28 | 35.08% relative |
| .75 | .5 | 25 | 25 | 24 | 35.54% relative |

These strata are descriptive, not independently powered coverage claims.
The median relative full-tidal errors pooled by amplitude are **101.51%**
(80 fitted weak cases) and **35.08%** (86 strong cases). Relative error at
zero curvature is undefined. This is not a 25% accuracy result.

Bias failures are cases **27,29,49,124,137,143,146,244,251**. All have A=.75;
case 124 has w=.5 and the other eight w=-.5. Their relative continuum biases
range from 2.0344% to 2.6699%. No directional asymmetry claim follows from
this post-result list. Case **150**, A=.25,w=0, has the zero-count failure.
The retained mask was chosen before seeing the measurement sample, so that
clock is not removed and the case is not replaced.

All 510 returned finite/continuum fits have projected Weyl rank ten, full
nuisance rank 277, successful optimizers and no boundary flags. The bias
violations persist when exact continuum clocks replace noisy counts; they
cannot be solved merely by adding more measurement events at this fixed
window. The maximum observed largest tidal variance is **0.1284773** at
N=262144. The earlier eight-case maximum was not a validated variance cap.

The [subsequent curvature-square pilot](square_pilot_results.md) studies a
new estimator on fresh data. Its bias improvement and the separate
[nested-landmark pilot](nested_pilot_protocol.md) are development evidence;
neither inherits coverage validation from this failed campaign. Accuracy,
regularity/remainder control and high-count calibration remain open.

## Verification and provenance

The actual process terminated with exit code zero after all 256 assigned
cases. Its 12 frozen source files remained unchanged throughout the run.
The iid population was independently regenerated with seed 2026092950 and
matched exactly; case seeds remain 2026100100+i. All 767 input/full-fit hashes
match, all archives have only the allowed six order fields, and all per-case
records agree with the final summary. The confidence/bias fields and exact
binomial gate were independently recomputed.

A pure-estimator replay of unfavorable case 124 reproduced every saved fit
field exactly without loading generator, profile or truth modules. Historical
source hashes were verified against the reachable generating commit.
[Recorded summary and audit](jet_coverage_summary.json) preserve every case,
stratum, failure and verification. Raw archives are in
`outputs/intrinsic_weyl_jet_coverage_v1`; the raw summary SHA256 is
`b8027f48b29946600c8dee070cb31900d027b8fd07759751955d48405d09fc55`.

The fixed sample count is **68,157,440 events plus 512 supplied tips**.
The sum of per-case wall times under four-worker execution and changing
machine load is 18,498.73 seconds; it is not the campaign's elapsed time.
There was no early stop, replacement or event-count extension.
