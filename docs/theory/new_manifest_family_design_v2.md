# New Manifest Family Design V2

Milestone 36 defines candidate v2 manifest families without running them.
Planned v2 families are not current results. Execution would be Milestone 37
or later and would require a new preregistered run.

The proposed families target observed or missing diagnostics:

- improve held-out agreement,
- reduce generalization gap,
- improve target and pair-node coverage,
- compute restart stability,
- compute latent-order stability,
- report null taxonomy completely.

Candidate planned families are:

- `rank_gap_more_protocol_columns_v2`,
- `rank_gap_coverage_enriched_v2`,
- `rank_gap_rank_resolution_enriched_v2`,
- `combined_diagnostic_complete_v2`,
- `report_only_failed_controls_v2`.

Each planned family must predeclare profile-column policy, target-inclusion
policy, pair-node coverage policy, null taxonomy policy, restart-stability
outputs, and latent-order stability outputs. These specifications are future
design proposals, not evidence that a family will pass carry-forward.

Milestone 37 is the first execution of this planned design. The execution
status is recorded in M37 output tables and v2 manifest files, while the M36
plan remains a planned preregistration artifact. The v2 families are not
carry-forward evaluated in Milestone 37.
