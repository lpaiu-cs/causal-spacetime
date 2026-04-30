# Protocol-Invariant Response Profiles

Milestone 40 audits response-profile protocol invariance before v3 execution.
This is a pre-execution design correction, not threshold retuning. It does not
change M34, M38, or M39 decisions. It does not generate v3 manifests. It does
not run stress tests. It does not fit new representation models.

## Measurement Protocol

A measurement protocol pi is the fixed collection of rules used to turn a
target event and a reference chain into a response value. A measurement
protocol pi includes at least:

- echo_rule
- emission_rule
- gate_rule
- reference_subsampling_rule
- spectrum_type
- normalization_rule
- missing_policy
- tie_policy
- margin_policy

Reference-chain variation is allowed inside a single response profile only
when these measurement rules are fixed. A single response profile must not mix
measurement protocols.

## Protocol-Invariant Profile

A reference set R is the declared set of reference chains used by one fixed
measurement protocol. The admissible protocol-invariant profile is:

```text
P_i^pi = (
    D_{pi,R1}(e_i),
    D_{pi,R2}(e_i),
    ...,
    D_{pi,Rm}(e_i)
)
```

where pi is fixed and only R_a varies.

This profile may be used for pairwise response-profile dissimilarity because
each column has the same measurement semantics. It remains pre-metric and does
not infer metric geometry.

## Protocol-Mixed Profile

An inadmissible profile is:

```text
P_i = (
    D_{pi1,R1}(e_i),
    D_{pi2,R2}(e_i),
    D_{pi3,R3}(e_i)
)
```

where both pi and R vary.

A protocol-mixed profile is inadmissible for pairwise response-profile
dissimilarity unless explicitly marked exploratory/report-only. If emission,
gate, echo rule, spectrum type, subsampling, normalization, missing policy, tie
policy, or margin policy varies, those variants must form separate profile
families or be explicitly marked exploratory/report-only.

## Cross-Protocol Robustness

Cross-protocol robustness compares results across separately constructed
profile families. It is the correct place to compare measurement-protocol
variation. It must not be implemented by concatenating rule variants into one
production response profile.

## Exploratory Mixed-Context Profile

An exploratory mixed-context profile is a labeled report-only object used to
identify how sensitive results are to changing measurement rules. It is not a
production input for pairwise response-profile dissimilarity and cannot be used
to produce handoff manifests.
## Milestone 41 Extension

M41 executes the protocol-invariant patched v3 manifest-generation run. A
response profile used for pairwise response-profile dissimilarity may vary
reference chains inside one fixed measurement protocol. If emission, gate, echo
rule, spectrum type, subsampling, normalization, missing policy, tie policy, or
margin policy varies, those variants must form separate profile families or be
explicitly marked exploratory/report-only. Protocol metadata, profile metadata,
and provenance metadata are admissibility metadata, not physical interpretation.
