# Fixed high-count diagnostic on eight retained landmark configurations

2026-09-23. Commit this protocol and its runner before new measurement data.
This is the next check of the event-count extrapolation in
[the nested pilot](nested_pilot_results.md), not a powered accuracy or coverage
campaign. The active 512-case campaign and all its frozen rules stay unchanged.

## Fixed comparison

Retain **all eight** nested-pilot configurations, in original case order,
including both flat cases and the case-6 missed confidence ball. Reuse their
4096-event selection orders and exactly their nested selected IDs. Do not
regenerate selection, change the clock mask, choose easier shapes, or supply
coordinates to the estimator. The generator may read the old selected points
and profiles solely to generate new relational observations.

For case i draw a new measurement stream with `default_rng([2026103000+i,1])`.
None of the old measurement events is reused. Generate in chunks of 4096 and
save the nested prefixes at **262144, 1048576 and 4194304** events. The three
levels share prefixes; they are not independent replicates. The maximum new
sample is **33554432 events**, with no new selection events or supplied tips.
Every assigned case and level remains in the record. If a guarded causal
comparison fails, retain already completed prefixes, record the sampled and
successfully classified counts separately, and do not redraw the failed
chunk. Later cases still run. No extension or outcome-dependent stopping.

Each estimator subprocess receives only the original six-field order archive:
selection order, design, selected IDs, after/before bits and measurement N.
Use the unchanged full-covariance fifth-degree vacuum jet with ten free Weyl
and 24 free gradient components, fixed curvature-square corrections, four
starts, masks, cutoffs and tolerances. Repeated membership grouping is the
already audited exact integer identity. Do not use the rejected affine map.

Write the pure fit before importing any generator or truth-evaluation module;
record that import boundary. Only then evaluate it against withheld truth and
refit continuum clocks using the same mask and measured covariance. Keep the
95% nominal ball, 2% relative continuum-bias allowance, 1e-10 flat tolerance
and 32/64 quadrature check. A returned fit does not erase an evaluation failure.

## What these data can decide

For every case and level report the full-tidal error, all existing fit and
bias checks, nominal containment, radius, and `V=N*lambda_max`. Compare V and
errors with both the fresh lowest-count prefix and the original pilot, without
turning a favorable point comparison into an accuracy verdict. Retain the
flat cases' absolute errors; their relative errors are undefined.

One stream per fixed geometry cannot establish a population variance envelope,
95% coverage, or 90% power. Error need not decrease monotonically along a
nested sample. A stable V or smaller error is evidence about extrapolation on
these geometries only, not permission to treat the conditional 102-million
event price as certified. An adverse high-count bias, failure, or persistent
error can invalidate using that price. There is no threshold tuned to the
observed results and no pass/fail accuracy gate in this diagnostic. A later
accuracy campaign must still freeze its interval gate and power calculation
before independent final data.

## Execution and audit

Use one measurement generator and one sequential estimator subprocess, with
`OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` set before Python imports.
Other BLAS override variables must be unset. The old-data thread replay found
physical equivalence on all eight configurations, while its raw six-of-eight
comparison failure remains unchanged. Do not modify the active campaign's
worker count or environment. Report wall time, process CPU time, peak Windows
working set, packed archive bytes and landmark/event predicate counts. These
are descriptive costs under concurrent load, not a controlled speedup.

Before any fresh sample, replay old nested case 5 in a separate subprocess
with all BLAS overrides unset. Compare **every saved finite fit, inference,
continuum fit and evaluation field** with the old full wrapper at the existing
1e-8 absolute and relative numerical tolerance, with exact discrete flags.
This checks that splitting the old evaluator into an already-fitted entry
point changes no results. A failed replay stops before new data; do not relax
the rule or silently choose another case.

Bind the clean reachable implementation commit, source hashes, old input and
truth hashes, fixed counts, seed schedule and runtime before generation.
Preserve each prefix archive, pure fit, evaluated result, failures and costs.
Check the bound sources and old inputs again at completion. Independently
audit the new prefixes, assigned counts, errors and at least one pure fit
after completion. An interrupted run is incomplete and is never called a
completed diagnostic.
