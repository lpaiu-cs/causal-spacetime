# V4 Preregistered Remediation Design

Milestone 43 preregisters a v4 upstream remediation design as planned-only
future work. It does not execute v4 and does not rescue v3 by changing fixed
criteria.

The v4 design goals are:

- retain fixed M34 criteria,
- retain protocol-invariant profile construction,
- retain parameter-complete measurement protocols,
- retain provenance-aware handoff manifests,
- increase manifest-level replication only if needed,
- improve held-out stability,
- improve latent-order stability,
- improve null taxonomy separation,
- improve response-profile strict-pair resolution,
- preserve failed and report-only controls,
- separate remediation by family/protocol, not by outcome cherry-picking.

The planned v4 families are:

- `rank_gap_earliest_full_stability_v4`,
- `rank_gap_earliest_retained_resolution_v4`,
- `rank_gap_gated_full_stability_v4`,
- `combined_earliest_full_stability_v4`,
- `rank_gap_high_resolution_reference_set_v4`,
- `failed_controls_v4`,
- `report_only_immediate_edge_v4`.

These are planned families, not current results. The v4 preregistration records
linked v3 root causes, required diagnostic metrics, forbidden actions, and
forbidden interpretations. It retains the fixed M34 criteria and rejects
threshold retuning.

Remediation design is not proof that a later run will pass. A future v4
execution would still need diagnostic production and a separate carry-forward
evaluation.

Milestone 44 executes the preregistered v4 protocol manifest-generation run from
this design. It does not evaluate carry-forward decisions. It does not run
stress tests. It does not retune thresholds. It does not infer metric geometry.
Diagnostic-complete v4 output is not carry-forward success.
