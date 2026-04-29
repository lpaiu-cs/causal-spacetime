# Upstream Manifest Remediation Design

Milestone 35 audits why no family currently passes carry-forward. It does not
retune thresholds. It does not run stress tests. It only records possible
future upstream remediation directions.

## Possible Upstream Causes

Observed or missing diagnostics can point to several upstream causes:

- insufficient protocol columns,
- sparse reachability,
- high tie rate,
- weak null separation,
- low target coverage,
- low pair-node coverage,
- unstable restarts,
- missing output files,
- missing metric aggregation.

Missing metrics are not success and not the same as measured failure. A
missing coverage metric, for example, means the current output bundle is
incomplete for that criterion; it does not establish either coverage success
or measured coverage failure.

## Future Remediation Directions

Possible future designs include richer profile generation, more protocol
columns, stricter target inclusion before manifest creation, stronger coverage
diagnostics, better null taxonomy reporting, explicit restart-stability
outputs, and more complete pair-node coverage tracking.

Execution-changing proposals require new preregistration. They cannot be
applied retroactively to current Milestone 34 outputs. Future remediation must
be preregistered before a new manifest family is generated or selected.

Pure reporting proposals may improve diagnostic completeness in later runs,
but they still cannot convert a blocked family into a carried-forward family
inside Milestone 35.

Milestone 36 turns these report-only proposals into a preregistered remediation
plan. Planned v2 families are not current results. Missing-metric remediation
and measured-failure remediation are distinct. Future execution requires a new
preregistered run.
