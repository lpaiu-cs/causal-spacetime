# Diagnostic-Complete Manifest Requirements

A diagnostic-complete manifest family output bundle must include every metric
required by the cross-family carry-forward criteria before any future
carry-forward decision is attempted.

Required family-level metrics are:

- `manifest_count`,
- `fitted_fraction`,
- `no_fit_fraction`,
- `mean_heldout_violation`,
- `mean_generalization_gap`,
- `stricter_threshold_pass_fraction`,
- `destructive_null_gap`,
- `symmetry_control_gap`,
- `target_coverage_fraction`,
- `pair_node_coverage_fraction`,
- `restart_std`,
- `latent_order_disagreement`,
- `no_retuning_audit_pass`,
- `failed_accounting_present`.

Fit summary metrics should come from family-level representation diagnostic
summaries. Stricter-threshold metrics should come from fixed criteria
comparison outputs. Null-gap metrics should come from null-taxonomy outputs.
Coverage metrics should come from target and pair-node coverage diagnostics.
Restart and latent-order metrics should come from fit-stability outputs.
Accounting metrics should come from no-retuning and failed-family accounting
audits.

Missing metrics block or warn according to criterion type. Missing metrics are
not imputed. A future remediation run must predeclare metric production before
execution.

Milestone 37 executes the v2 manifest-generation run and is required to
produce every listed family-level metric in a v2 output bundle. Diagnostic
completeness is a reporting requirement, not a success criterion and not a
claim about metric geometry.

Milestone 38 uses those metrics as inputs to fixed carry-forward criteria. It
does not create replacement metrics after seeing the v2 decisions.
