# V5 Preregistered Remediation Design

Milestone 46 preregisters a v5 planned-only design. V5 design is planned-only
and requires a later execution milestone before any production manifests exist.

The v5 design goals are:

- retain fixed M34 criteria,
- retain protocol-invariant profile construction,
- retain parameter-complete measurement protocols,
- retain provenance-aware handoff manifests,
- preserve failed/report-only controls,
- target held-out stability without threshold changes,
- target latent-order stability without threshold changes,
- target destructive-null separation without changing null taxonomy post hoc,
- target pair-node coverage without removing hard cases post hoc,
- include explicit stopping criteria for future remediation cycles.

The planned-only v5 families are:

- `rank_gap_earliest_full_manifest_transfer_v5`,
- `rank_gap_high_coverage_strict_pair_v5`,
- `rank_gap_latent_stability_replicated_v5`,
- `combined_null_separation_v5`,
- `rank_gap_low_tie_reference_diverse_v5`,
- `failed_controls_v5`,
- `report_only_immediate_edge_v5`.

The v5 design is a no-execution remediation plan. It does not claim that the
next run will pass. It retains fixed M34 criteria and records remediation
iteration risk before any v5 production manifest generation.
