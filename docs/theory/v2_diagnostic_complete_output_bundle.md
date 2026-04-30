# V2 Diagnostic-Complete Output Bundle

The v2 output bundle is the set of Milestone 37 CSV and JSON artifacts needed
for a later carry-forward evaluation. It is a diagnostic bundle, not a
carry-forward decision and not a stress-test result.

Required output files include:

- `outputs/data/v2_manifest_generation.csv`
- `outputs/manifests_v2/*.json`
- `outputs/data/v2_manifest_family_fit_summary.csv`
- `outputs/data/v2_manifest_family_null_taxonomy.csv`
- `outputs/data/v2_manifest_family_stricter_criteria.csv`
- `outputs/data/v2_manifest_family_failed_accounting.csv`
- `outputs/data/v2_manifest_family_coverage_metrics.csv`
- `outputs/data/v2_manifest_family_restart_stability.csv`
- `outputs/data/v2_manifest_family_latent_order_stability.csv`
- `outputs/data/v2_no_retuning_audit.csv`
- `outputs/data/v2_cross_family_robustness_metrics.csv`
- `outputs/data/v2_diagnostic_completeness_check.csv`
- `outputs/data/v2_diagnostic_complete_bundle_report.csv`

Each v2 family row must expose the 14 required family-level metrics:

- `manifest_count`
- `fitted_fraction`
- `no_fit_fraction`
- `mean_heldout_violation`
- `mean_generalization_gap`
- `stricter_threshold_pass_fraction`
- `destructive_null_gap`
- `symmetry_control_gap`
- `target_coverage_fraction`
- `pair_node_coverage_fraction`
- `restart_std`
- `latent_order_disagreement`
- `no_retuning_audit_pass`
- `failed_accounting_present`

The fit-summary experiment produces the manifest, fitted, no-fit,
held-out-violation, and generalization-gap metrics. The null taxonomy
experiment produces destructive-null and symmetry-control gaps. The stricter
criteria experiment produces fixed-threshold pass fractions. The failed
accounting experiment keeps report-only and failed-control families visible.
The coverage experiment produces target and pair-node coverage. The stability
experiment produces restart and latent-order diagnostics. The no-retuning
audit supplies the fixed-settings accounting metric.

Missing metrics make the v2 bundle incomplete. Missing values are not
imputed as success. Failed-control families may carry conservative report-only
placeholder values where fitting is intentionally skipped, and those rows
remain separated by family kind.

Metric completeness is not carry-forward success. The bundle is a prerequisite
input for a future evaluation only.

Milestone 38 evaluates this bundle. Carry-forward means future stress-test
eligibility only, and carry-forward does not imply geometry.
It does not regenerate manifests. It does not rerun fits. It does not retune
thresholds. It does not run stress tests.

Milestone 39 uses the evaluated v2 bundle to separate structural and measured
blocking. Diagnostic completeness remains a property of the input bundle, not
a carry-forward outcome.

Milestone 40 does not alter the v2 bundle. It audits response-profile protocol
invariance before v3 execution and patches planned v3 family semantics. A
single response profile must not mix measurement protocols.
