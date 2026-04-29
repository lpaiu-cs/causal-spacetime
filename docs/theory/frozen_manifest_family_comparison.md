# Frozen Manifest Family Comparison

Milestone 33 compares frozen manifest families under preregistered settings.
A manifest family is a declared group of handoff manifests with shared
metadata, such as comparison method, eligibility, margin policy, or failure
status. The goal is to compare families, not to tune individual latent ordinal
representation fits.

Families include eligible structured manifests, ineligible reported manifests,
failed controls, destructive nulls, symmetry controls, and marginal-preserving
nulls. Eligible manifests are fit by the same fixed diagnostic settings.
Ineligible and failed manifests are retained in accounting and produce no-fit
rows by default.

Target-label permutation is a symmetry control, not a destructive null. It
renames targets consistently and therefore can preserve the response-comparison
structure. Destructive nulls include shuffled sides and random constraints;
these are expected to disrupt response-comparison consistency. Marginal-
preserving profile nulls include delay shuffle, reachability shuffle, and
same-marginal random profiles; these preserve selected marginals and can be
harder baselines.

Family-level diagnostics include held-out violation, train/held-out
generalization gap, restart stability, latent-order stability, null separation,
no-fit rate, and failed-manifest coverage. Family-level fit quality is a latent
ordinal diagnostic, not physical geometry. The comparison does not retune
thresholds after seeing fit results.

Failed and ineligible manifests must be reported. A family comparison that
only reports successful fits would lose the negative evidence supplied by
handoff validation and failed-control manifests. No physical geometry is
inferred from family-level pass rates, null separation, or selected latent
representation dimension.
