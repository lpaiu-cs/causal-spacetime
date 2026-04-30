"""Final exact sanity checks for M40 protocol-invariance patch."""

from __future__ import annotations

import json
from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    default_v3_protocol_invariant_family_patches,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preregistration import (
    build_v3_protocol_patched_preregistration,
    v3_protocol_patched_preregistration_to_jsonable,
)
from causal_spacetime_lab.state_change_measurement_protocol import (
    default_earliest_full_protocol,
    default_gated_full_protocol,
)
from causal_spacetime_lab.state_change_response_profile_metadata import (
    profile_metadata_from_protocols,
)


def run_experiment(output_dir: Path = Path("outputs")) -> list[dict[str, float | str]]:
    protocol = default_earliest_full_protocol()
    patches = default_v3_protocol_invariant_family_patches()
    structured_protocol_invariant = all(
        patch.execution_status == "planned_only"
        and patch.measurement_protocol_hash
        and patch.admissible_for_pairwise_dissimilarity
        for patch in patches
        if patch.family_kind == "structured"
    )
    mixed = profile_metadata_from_protocols(
        "mixed_report_only",
        [protocol, default_gated_full_protocol()],
        ["r1", "r2"],
        "refs",
        exploratory_mixed_context=True,
    )
    spec = build_v3_protocol_patched_preregistration(
        output_dir / "remediation" / "v3_preregistration_spec_m39.json",
        patches,
    )
    json.dumps(v3_protocol_patched_preregistration_to_jsonable(spec), sort_keys=True)
    no_carry_forward_decision = not (
        output_dir / "data" / "v3_cross_family_robustness_decisions.csv"
    ).exists() and not (
        output_dir / "carry_forward_v3" / "carry_forward_registry_v3.json"
    ).exists()
    return [
        {
            "check": "protocol_metadata_fields_exist",
            "passed": float(bool(protocol.echo_rule and protocol.missing_policy)),
        },
        {
            "check": "v3_patch_families_protocol_invariant",
            "passed": float(structured_protocol_invariant),
        },
        {
            "check": "protocol_mixed_context_report_only",
            "passed": float(
                mixed.profile_invariance_status == "protocol_mixed"
                and mixed.exploratory_mixed_context
                and not mixed.admissible_for_pairwise_dissimilarity
            ),
        },
        {
            "check": "patched_preregistration_json_serializes",
            "passed": 1.0,
        },
        {
            "check": "v3_execution_still_disallowed",
            "passed": float(not spec.execution_allowed_in_current_milestone),
        },
        {
            "check": "no_carry_forward_decision_produced",
            "passed": float(no_carry_forward_decision),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/protocol_invariance_final_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote protocol invariance final sanity: {path}")


if __name__ == "__main__":
    main()
