# Frozen independent coverage check for the blind fifth-degree model

2026-09-22. Freeze this protocol and its runner before fresh data. The
[eight-case pilot](jet_pilot_results.md) repaired varying-profile continuum
bias on its four nonconstant cases, while actual finite-order precision
remained poor. Its eight broad confidence balls do not validate coverage.
This experiment tests the unchanged enlarged model before an expensive
precision campaign. It does not add a new inverse model or re-score a pilot.

## Estimator, target and confidence construction

Use the same connected selector, 64-landmark cap, fifth-degree joint fit,
24 free vacuum gradient components, four flat initialization starts,
600-evaluation limit, tolerances and covariance cutoff. The estimator
receives only the six-field order archive; coordinates, amplitude, rotation
rate and separate truth archive remain unavailable until its fit returns.
No paired flat observation is supplied. The target is all ten signed
dimensionless Weyl components in the tip/landmark frame, modulo z reflection.

Keep the existing 95% simultaneous ten-component Gaussian norm ball,
including the relative bias allowance **delta=0.02**. If lambda is the
largest eigenvalue of the tidal-norm covariance and c is the estimate,

```text
r_stat = sqrt(chi2_10(.95) * lambda)
r = (r_stat + .02 * ||c||_G) / .98.
```

The covariance includes both unknown coordinate and gradient nuisance
parameters. Full nuisance rank, Weyl rank ten, optimizer success and no
boundary excess are required. This is an estimated-covariance Gaussian
construction whose performance is being tested, not presumed.

Refit exact continuum-volume clocks on each case's SAME landmarks, clock
mask and weights. The generator/evaluator checks 32/64 quadrature agreement
below 1e-10 relative. A relative continuum bias above 2% is a validation
failure even when the broad ball happens to cover the target. For A=0,
absolute continuum bias must be at most 1e-10. This does not prove a uniform
bias bound and does not use truth to tune or correct the finite estimate.

## Independent population and fixed budget

Exactly **256 iid cases**. Initialize `default_rng(2026092950)`, draw 256
amplitudes uniformly from `[0,.25,.75]`, then independently draw 256 rotation
rates uniformly from `[-.5,0,.5]`. This defines a nine-cell product mixture;
the three flat cells have identical physical laws. Do not choose cases to
balance observed outcomes. Preserve every stratum descriptively.

Case i uses seed **2026100100+i**, with independent selection/measurement
streams `[seed,0]` and `[seed,1]`. These seeds have not been used by the
development pilots. Each case receives 4096 selection events and 262144
measurement events in chunks of 4096. The total cap is **68,157,440 events
plus tips**. No sample increase, new seed, excluded failure, interim stopping
rule, parameter adjustment or replacement case. The old eight pilot cases
are excluded from the validation denominator.

Use the existing exact rotating-wave generator and evaluator. Run four
independent worker processes, with case-local random streams and files.
Worker scheduling must not affect the samples. Store and report all 256
results in case order. A generator or fit failure remains a failed case;
a process or storage failure that prevents completion makes the run
incomplete, never a success. Do not silently relaunch or replace it.

## Primary gate and power

A case is valid-and-covered only if generation, selection, fitting, ranks,
boundary checks, continuum audit and the frozen bias condition pass, and
the full-tidal error modulo reflection is no greater than r. Include ALL
256 cases, including failures, in the denominator.

Compute the **95% two-sided exact Clopper--Pearson interval** for the
mixture's end-to-end valid-and-covered probability. The gate requires the
lower endpoint to be at least **.95**, hence **at least 251/256 successes**.
This is the same construction and threshold used for the earlier quartic
coverage validation, now applied to an independent population and model.

The exact power at true probability .99 is **.95471764**, exceeding the
90% planning requirement. At .98 it is only **.59470691**; the planning
alternative .99 is optimistic and explicitly unvalidated for this model.
This gate can fail. No separate per-stratum or uniform-in-frequency
coverage guarantee follows from a mixture-level success.

Preserve source hashes and a reachable clean pre-data commit; all order,
truth, full-fit and failure archives; observed bias, tensor error, radius,
variance and runtime. The executable gate self-check checks 250 vs 251
successes and rejects bias failures despite nominal coverage. Before data,
replay an old case through one spawned worker and compare all six order
fields and the complete fifth-degree fit/evaluation to the saved pilot.
After completion verify hashes, population, counts, gate calculation and
at least one pure order-only estimator replay independently.

## Why coverage comes before a precision campaign

The eight-case pilot's largest observed `lambda*N` is **22120.75214**;
this is a development observation, not a variance envelope. Applying the
existing sufficient joint-confidence/power calculation with that maximum
gives, for ONE case and a 25% full-tensor containment target,

```text
N >= V * (sqrt(chi2_10(.95)) + sqrt(chi2_10(.90)))^2
         / (((1-.02)*.25 - 2*.02) * ||true||_G)^2.
```

The conditional prices are approximately **16.03 million** measurement
events at norm 1.5 and **144.25 million** at norm .5. For four simultaneous
cases, replacing the tail probabilities by .05/4 and .10/4 raises these
to approximately **20.13 million** and **181.14 million** per case.
These are Gaussian, bias- and variance-conditional planning calculations,
not validated prices or authorization to launch that accuracy study.

A timing replay of the unchanged old case 7 took about 16.93 seconds to
regenerate identical order fields and 4.45 seconds for the pure estimator.
Hardware load changes timing, so extrapolations are approximate. The
coverage study retains the smaller event count and first tests whether the
enlarged inverse and its bias/uncertainty construction survive fresh cases.
Even success leaves high-count calibration, actual interval-contained
accuracy and regularity/remainder scope open. No claim about arbitrary
oscillations, singularities or global topology is made by this experiment.

```powershell
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_jet_coverage.py --check
rtk proxy .venv/Scripts/python.exe experiments/positive_control/intrinsic_weyl_jet_coverage.py
```
