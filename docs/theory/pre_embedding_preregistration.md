# Pre-Embedding Preregistration

Protocol preregistration freezes analysis choices before any future
representability experiment. It separates exploratory diagnostics from
predeclared no-fit validation.

## Frozen Settings

Before any future embedding attempt, the project must freeze:

- profile generation settings,
- target inclusion policy,
- pairwise response-comparison protocol,
- missing-data policy,
- margin threshold,
- constraint sampling seed,
- train and held-out constraint split,
- validation thresholds,
- null-baseline definitions,
- stop rules,
- forbidden interpretations.

These settings are recorded in the handoff manifest. The manifest is an input
specification for future experiments and does not contain fitted embeddings.

## Leakage Risks

Leakage occurs if an analysis choice is changed after inspecting validation or
future test outcomes. Examples include:

- selecting the comparison protocol after looking at held-out performance,
- changing the margin after null-baseline failure,
- tuning an embedding after seeing test error,
- excluding difficult targets after failure,
- adding easier null baselines after difficult nulls fail.

## Exploratory Versus Preregistered

Exploratory analyses may compare protocols, thresholds, and failure modes, but
they must be labeled exploratory. A preregistered handoff manifest freezes one
specific protocol and its stop rules. Future representability experiments must
report both eligible and failed manifests.

Handoff eligibility does not imply metric representability. It only states
that a response-comparison constraint pool has passed predeclared finite
diagnostics and may be used as a future input under explicit assumptions.
