# V3 Preregistered Manifest-Family Design

Milestone 39 preregisters a v3 manifest-family design but does not execute it.
The design responds to the v2 blocked-decision audit by planning replicated
families and diagnostic-complete reporting before any future run.

The v3 design goals are:

- at least 3 manifests per structured family,
- all 14 diagnostic metrics,
- explicit family-level replication,
- fixed Milestone 34 criteria retained,
- no threshold retuning,
- failed and report-only controls retained.

The proposed planned families are:

- `rank_gap_multi_manifest_v3`,
- `rank_gap_more_protocol_columns_v3`,
- `rank_gap_rank_resolution_enriched_v3`,
- `rank_gap_coverage_enriched_v3`,
- `combined_diagnostic_complete_v3`,
- `failed_controls_v3`.

These are planned families, not current results. V3 planned families are not
current results. The v3 family design is preregistered but not executed.
Milestone 40 patches the planned design before execution. Execution would be
Milestone 41 or later and would require its own output boundary.

No thresholds are retuned. No stress tests are run. The v3 plan does not
reinterpret blocked v2 families as eligible.

Milestone 40 audits response-profile protocol invariance before v3 execution.
This is a pre-execution design correction, not threshold retuning. It does not
change M34, M38, or M39 decisions. It does not generate v3 manifests. It does
not run stress tests. It does not fit new representation models. A single
response profile must not mix measurement protocols.

A response profile used for pairwise response-profile dissimilarity may vary
reference chains inside one fixed measurement protocol. If emission, gate,
echo rule, spectrum type, subsampling, normalization, missing policy, tie
policy, or margin policy varies, those variants must form separate profile
families or be explicitly marked exploratory/report-only. The patched v3
families remain planned-only.
## Milestone 41 Protocol Metadata Patch

M41 executes the protocol-invariant patched v3 manifest-generation run from the
Milestone 40 patched preregistration. The original planned-only v3 design remains
part of the preregistration history; M41 adds parameter-complete protocol
metadata, response-profile metadata, and handoff provenance metadata before
production v3 manifests are written.
