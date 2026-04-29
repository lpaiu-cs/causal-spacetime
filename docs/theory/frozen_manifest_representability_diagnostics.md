# Frozen Manifest Representability Diagnostics

Milestone 32 performs latent ordinal representation diagnostics from frozen manifests.
The manifest was produced before fitting by the Milestone 31 no-fit handoff
protocol, so it fixes the response-comparison constraint pool and the
manifest train/held-out split before any representation model is evaluated.

A latent ordinal representation is a finite mathematical representation whose
variables are adjusted to satisfy response-comparison constraints of the form
that one target pair has smaller response-profile dissimilarity than another
target pair. The fitted coordinates are not spatial coordinates. They are not
calibrated positions, not measured quantities, and not an interpretation of
physical geometry.

The manifest train/held-out split must not be changed. The train constraints
are the only constraints used during fitting. Held-out constraints are evaluated
after fitting as a diagnostic of constraint generalization, not physical
geometry. The held-out set is not a tuning set for protocol thresholds,
dimension choices, or post-hoc target removal.

Null representation baselines are also finite diagnostics. They compare the
frozen manifest fit to shuffled-side, random-constraint, and target-permutation
controls. Null baselines test artifacts, not geometry. Better-than-null behavior
can motivate further representability tests, but it does not infer spatial
distance or metric geometry.

Failed manifests must be reported, not silently omitted. Ineligible manifests
produce explicit no-fit rows unless an exploratory diagnostic deliberately opts
into fitting them. This protects the handoff procedure from survivorship bias.

The outcome of Milestone 32 is a manifest-level report: train and held-out
violation rates, hinge losses, restart stability, null comparisons, and
dimension-complexity diagnostics. These reports are finite response-comparison
representation diagnostics only.

Milestone 33 compares frozen manifest families under preregistered settings.
It keeps the Milestone 32 fitting procedure fixed while grouping structured,
ineligible, failed-control, destructive-null, symmetry-control, and
marginal-preserving-null cases. Target-label permutation is treated as a
symmetry control, not a destructive null. Family-level fit quality is a latent
ordinal diagnostic, not physical geometry.

Milestone 34 consumes those family-level reports to decide carry-forward
eligibility for future stress tests. It performs no new fitting and does not
retune thresholds after seeing fit outcomes.
