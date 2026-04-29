# Response Constraint Handoff Protocol

Milestone 31 defines predeclared handoff manifests. A handoff manifest is a
frozen input specification for future experiments, not an embedding result.
It records how a response-comparison constraint pool was generated, validated,
split, and either accepted or rejected.

Milestone 30 validation gates are necessary but not sufficient for
representation. A pool that passes held-out protocol agreement, bootstrap
stability, null-baseline separation, and coverage gates is only eligible for a
future representability test under explicit additional assumptions.

## Manifest Contents

A handoff manifest freezes:

- response-profile generation settings,
- pairwise response-comparison protocol,
- missing-data policy,
- margin threshold,
- constraint sampling seed,
- train and held-out protocol-column split,
- validation gate thresholds,
- null-baseline settings,
- accepted constraint pool,
- future frozen constraint split,
- failure reasons and stop rules,
- forbidden interpretations.

The manifest stores response-comparison constraints only. It does not contain
fitted embeddings, metric coordinates, calibrated radar distances, or spatial
geometry.

## Separation Rules

Train, held-out, and null-baseline information must remain separated. Future
representation procedures must not tune a fit using held-out validation
outcomes. Thresholds must be fixed before any future fit. Null baselines must
be declared before evaluation, not selected after seeing results.

## Stop Rules

A manifest is ineligible if it fails predeclared constraints on count,
held-out evaluability, held-out agreement, held-out inversion, bootstrap
confidence, null separation, target coverage, or pair-node coverage. Failed
handoff gates are reported as negative results, not silently removed.

Handoff eligibility does not imply metric representability. Passing handoff
does not imply distance, geometry, or calibrated scale.

Milestone 32 consumes these frozen manifests without changing their constraint
splits. The manifest train split is used for latent ordinal representation
fitting, while the manifest held-out split is evaluated only after fitting.
Held-out performance is a diagnostic of constraint generalization, not
physical geometry.
