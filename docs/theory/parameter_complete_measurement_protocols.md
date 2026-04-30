# Parameter-Complete Measurement Protocols

Measurement protocol identity must include both rule names and relevant rule
parameters. The categorical rules identify what is being measured, and the
parameters identify the fixed choices inside those rules.

A measurement protocol includes at least the echo rule, emission rule, gate rule,
reference-subsampling rule, spectrum type, normalization rule, missing-data
policy, tie policy, and margin policy. The protocol hash also includes declared
parameters such as gate delay rank, emission position rule, reference stride,
retained-reference policy, normalization scope, and margin value.

Changing `reference_chain_ids` must not change `measurement_protocol_id`.
Changing `echo_rule`, `gate_rule`, `spectrum_type`, `missing_policy`,
`tie_policy`, or `margin_policy` must change `measurement_protocol_id`.
Changing declared gate delay, emission policy parameter, reference stride,
retained-reference policy, normalization parameter, or margin value must also
change `measurement_protocol_id`.

If required parameter values are omitted, the profile is underspecified and not
production-admissible for pairwise response-profile dissimilarity. Parameter
metadata is an admissibility check, not a physical interpretation.

