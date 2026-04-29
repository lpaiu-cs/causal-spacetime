"""Representability requirement ladder for response-order diagnostics."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class RepresentabilityRequirement:
    """One level in the pre-metric representability ladder."""

    level: str
    input_structure: str
    output_structure: str
    additional_assumptions: str
    what_it_does_not_imply: str


def default_response_representability_ladder() -> list[RepresentabilityRequirement]:
    """Return the default response-order representability ladder."""

    return [
        RepresentabilityRequirement(
            level="scalar_response_rank",
            input_structure="response-order signature",
            output_structure="topological scalar ranks",
            additional_assumptions="acyclicity",
            what_it_does_not_imply="target-target distance order or geometry",
        ),
        RepresentabilityRequirement(
            level="stable_scalar_response_core",
            input_structure="multiple protocol variants",
            output_structure="robust scalar order core",
            additional_assumptions="protocol stability threshold",
            what_it_does_not_imply="spatial distance",
        ),
        RepresentabilityRequirement(
            level="multi_reference_response_profile",
            input_structure="multiple reference chains or emission positions",
            output_structure="response-profile embedding candidate",
            additional_assumptions="aligned target identities and protocol labels",
            what_it_does_not_imply="metric space",
        ),
        RepresentabilityRequirement(
            level="pairwise_distance_order",
            input_structure="explicit target-target pair comparisons",
            output_structure="distance-order relation",
            additional_assumptions="pairwise comparison protocol",
            what_it_does_not_imply="ratios or metric tensor",
        ),
        RepresentabilityRequirement(
            level="ordinal_embedding_candidate",
            input_structure="pairwise distance-order constraints",
            output_structure="low-dimensional coordinates preserving many constraints",
            additional_assumptions=(
                "stable constraints, held-out validation, null baselines"
            ),
            what_it_does_not_imply="calibrated metric geometry",
        ),
        RepresentabilityRequirement(
            level="calibrated_metric_representation",
            input_structure=(
                "ordinal embedding plus calibration, measure, or atlas consistency"
            ),
            output_structure="effective metric representation",
            additional_assumptions=(
                "calibration, scale, consistency, dynamics-like constraints"
            ),
            what_it_does_not_imply="fundamental metric field or GR dynamics",
        ),
    ]


def representability_ladder_table() -> list[dict[str, str]]:
    """Return the default ladder as CSV-friendly rows."""

    return [
        {key: str(value) for key, value in asdict(requirement).items()}
        for requirement in default_response_representability_ladder()
    ]
