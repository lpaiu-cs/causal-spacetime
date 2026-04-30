"""Exact sanity checks for response-profile protocol invariance."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_measurement_protocol import (
    default_earliest_full_protocol,
    default_gated_full_protocol,
)
from causal_spacetime_lab.state_change_response_profile_metadata import (
    profile_metadata_from_protocols,
)


def run_experiment() -> list[dict[str, float | str]]:
    protocol = default_earliest_full_protocol()
    invariant = profile_metadata_from_protocols(
        "invariant_family",
        [protocol, protocol],
        ["r1", "r2"],
        "reference_set_a",
    )
    mixed = profile_metadata_from_protocols(
        "mixed_family",
        [protocol, default_gated_full_protocol()],
        ["r1", "r2"],
        "reference_set_a",
    )
    underspecified = profile_metadata_from_protocols(
        "underspecified_family",
        [],
        [],
        "",
    )
    exploratory = profile_metadata_from_protocols(
        "exploratory_mixed_family",
        [protocol, default_gated_full_protocol()],
        ["r1", "r2"],
        "reference_set_a",
        exploratory_mixed_context=True,
    )
    return [
        {
            "check": "protocol_invariant_profile_admissible",
            "passed": float(invariant.admissible_for_pairwise_dissimilarity),
        },
        {
            "check": "protocol_mixed_profile_blocked",
            "passed": float(
                mixed.profile_invariance_status == "protocol_mixed"
                and not mixed.admissible_for_pairwise_dissimilarity
            ),
        },
        {
            "check": "underspecified_profile_blocked",
            "passed": float(
                underspecified.profile_invariance_status == "underspecified"
                and not underspecified.admissible_for_pairwise_dissimilarity
            ),
        },
        {
            "check": "exploratory_mixed_context_report_only",
            "passed": float(
                exploratory.exploratory_mixed_context
                and exploratory.profile_invariance_status == "protocol_mixed"
                and not exploratory.admissible_for_pairwise_dissimilarity
            ),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/response_profile_invariance_exact_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote response profile invariance exact sanity: {path}")


if __name__ == "__main__":
    main()
