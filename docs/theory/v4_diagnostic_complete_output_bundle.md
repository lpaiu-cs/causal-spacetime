# V4 Diagnostic-Complete Output Bundle

The v4 diagnostic-complete output bundle is the collection of production v4
handoff manifests and family-level diagnostic outputs produced by Milestone 44.

Required output files include:

- `outputs/manifests_v4/*.json`,
- `outputs/data/v4_protocol_manifest_generation.csv`,
- `outputs/data/v4_protocol_manifest_metadata_audit.csv`,
- `outputs/data/v4_protocol_manifest_family_fit_summary.csv`,
- `outputs/data/v4_protocol_manifest_family_null_taxonomy.csv`,
- `outputs/data/v4_protocol_manifest_family_stricter_criteria.csv`,
- `outputs/data/v4_protocol_manifest_family_failed_accounting.csv`,
- `outputs/data/v4_protocol_manifest_family_coverage_metrics.csv`,
- `outputs/data/v4_protocol_manifest_family_restart_stability.csv`,
- `outputs/data/v4_protocol_manifest_family_latent_order_stability.csv`,
- `outputs/data/v4_protocol_no_retuning_audit.csv`,
- `outputs/data/v4_protocol_cross_family_robustness_metrics.csv`,
- `outputs/data/v4_protocol_diagnostic_completeness_check.csv`,
- `outputs/data/v4_protocol_diagnostic_complete_bundle_report.csv`.

The required family-level metrics are:

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

Diagnostic-complete means all required metrics are present. It does not mean
carry-forward success. Carry-forward evaluation is deferred to M45.

