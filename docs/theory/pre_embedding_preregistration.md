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

Milestone 32 is the first consumer of those preregistered inputs. It fits
latent ordinal representation models from frozen manifests, but it does not
change the manifest train/held-out split, tune thresholds using held-out
results, or exclude failed manifests after the fact.

Milestone 33 extends the preregistration discipline to family-level reports.
The family comparison uses fixed fit settings and fixed diagnostic thresholds;
it does not retune after seeing fit outcomes. Target-label permutation is
reported as a symmetry control, destructive nulls and marginal-preserving nulls
are separated, and failed or ineligible manifests remain visible.

Milestone 34 applies fixed carry-forward criteria to those reports. Threshold
sensitivity is recorded as sensitivity analysis only, not as a way to choose a
preferred threshold after the fact.

Milestone 35 adds a no-retuning stop-condition audit. If carry-forward gates
fail, the permitted continuation is failure decomposition and diagnostic
completeness reporting. It does not retune thresholds. It does not run stress
tests. Missing metrics are not success and not the same as measured failure.
Future remediation must be preregistered.

Milestone 36 records that future execution requires a new preregistered run.
The remediation plan freezes intended metrics, planned families, forbidden
interpretations, and the execution boundary before any future manifest
generation is attempted.

Milestone 37 executes that preregistered v2 generation run. It records fixed
settings, writes production v2 manifests, and preserves the no-retuning rule:
the resulting v2 diagnostics are not used inside M37 to alter carry-forward
criteria.

Milestone 38 preserves the same no-retuning boundary. It reads M37 outputs and
applies the fixed criteria without changing manifest generation, fit settings,
or thresholds after seeing decisions.

Milestone 39 preserves the boundary again. It audits v2 blocked decisions and
exports a planned-only v3 preregistration. V3 planned families are not current
results, and no thresholds are retuned.

Milestone 40 preserves the boundary before v3 execution. It audits
response-profile protocol invariance and patches the planned-only v3 design.
This is a pre-execution design correction, not threshold retuning. It does not
change M34, M38, or M39 decisions. It does not generate v3 manifests. It does
not run stress tests. It does not fit new representation models. A single
response profile must not mix measurement protocols.
## Milestone 41 Preregistration Boundary

M41 executes the protocol-invariant patched v3 manifest-generation run only
after the patched v3 preregistration is frozen. Top-down and hybrid handoff
provenance are admissible only when the design source, design digest, allowed
dependencies, forbidden dependencies, and evaluation boundary are recorded before
fit, carry-forward, or stress-test evaluation.

## Milestone 42 Boundary

Milestone 42 applies fixed carry-forward criteria to the M41 v3 protocol bundle.
Top-down/hybrid provenance is allowed only as preregistered manifest provenance,
not as evidence. M42 does not regenerate manifests, rerun fits, retune
thresholds, or run stress tests. Carry-forward is stress-test eligibility, not
geometry.

## Milestone 43 Boundary

Milestone 43 audits why v3 protocol families remain blocked and records a v4
upstream remediation design before execution. It does not change M42 decisions.
It does not retune thresholds. It does not run stress tests. It does not execute
v4. V4 design is planned-only and requires a later execution milestone.
