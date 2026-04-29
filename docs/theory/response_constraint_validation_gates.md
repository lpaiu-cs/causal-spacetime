# Response Constraint Validation Gates

Milestone 30 defines validation gates for response-comparison constraint
pools. These gates are pre-embedding diagnostics: they decide whether a chosen
pool is stable enough to be considered as an input to future representability
tests, not whether it has spatial or metric meaning.

Pairwise response-profile comparison protocols from Milestone 29 produce
constraints of the form:

```text
(i, j) <_pi (k, l)
```

This means the selected protocol `pi` assigns a smaller response-profile
dissimilarity to target pair `(i, j)` than to target pair `(k, l)`. These are
response-comparison constraints. They are not spatial-distance claims.

Validated response-comparison constraints are not distance-order constraints.
They are pre-metric ordinal constraints that may become inputs to future
representability tests only after explicit additional assumptions.

## Gate Components

Held-out protocol agreement splits response-profile protocol columns into a
training subset and a held-out subset. A pool built from the training columns
is evaluated on held-out columns. This tests whether the response-comparison
ordering is specific to one subset of protocols or survives across shared
protocol columns.

Bootstrap stability resamples protocol columns with replacement and evaluates
a fixed pool across bootstrap profiles. Stable constraints are those that keep
agreeing across enough evaluable bootstrap samples.

Null-baseline separation compares a structured profile against profiles with
shuffled delays, shuffled reachability, permuted target profiles, or random
same-marginal structure. This tests whether agreement is above simple
marginal-preserving baselines. Passing a null baseline does not imply geometry.

Coverage diagnostics report whether enough targets and unordered target-pair
nodes appear in the pool. Low coverage warns that a future representation test
would be driven by a narrow subset of the target population.

## Missing Data

Missing-data policy affects constraint validity. A conservative policy may
discard many target pairs. A reachability-aware penalty can keep more pairs
but changes the protocol output. The policy must be declared before evaluating
held-out agreement, bootstrap stability, or null separation.

## Margins And Confidence

High-margin constraints compare target pairs with a larger gap in
response-profile dissimilarity. Increasing the margin threshold can improve
agreement but reduces pool size and coverage. Constraint confidence is
therefore a tradeoff among margin, coverage, held-out stability, bootstrap
stability, and null-baseline separation.

Milestone 30 does not fit an embedding. It does not reconstruct metric
structure, calibrated radar distance, finite-speed spatial geometry, or any
physical distance scale.

Milestone 31 adds the handoff step after validation. Passing a gate may allow
a manifest to be exported for future experiments, but the manifest remains a
no-fit input specification and does not assert representability.
