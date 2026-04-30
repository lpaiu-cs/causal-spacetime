from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    default_new_manifest_family_specs_v2,
)
from causal_spacetime_lab.state_change_manifest_v2_generation import (
    V2ManifestGenerationConfig,
    build_v2_handoff_manifest_for_family,
    comparison_protocol_for_v2_family,
    write_v2_handoff_manifests,
)
from causal_spacetime_lab.state_change_manifest_v2_profiles import (
    default_v2_profile_configs,
)
from causal_spacetime_lab.state_change_manifest_v2_spec import V2ManifestFamilySpec


def _specs() -> list[V2ManifestFamilySpec]:
    return [
        V2ManifestFamilySpec(**row)
        for row in default_new_manifest_family_specs_v2()
    ]


def test_comparison_protocol_for_v2_family() -> None:
    rank = comparison_protocol_for_v2_family("rank_gap_more_protocol_columns_v2")
    combined = comparison_protocol_for_v2_family("combined_diagnostic_complete_v2")

    assert rank.method == "rank_gap_mean"
    assert combined.method == "combined_gap_and_mismatch"


def test_build_v2_handoff_manifest_for_family_returns_manifest() -> None:
    spec = _specs()[0]
    profile_config = default_v2_profile_configs([spec], seed=0)[0]

    manifest = build_v2_handoff_manifest_for_family(
        spec,
        profile_config,
        V2ManifestGenerationConfig(
            max_constraints=100,
            bootstrap_count=1,
            null_repetitions=1,
        ),
    )

    assert manifest.created_by_milestone == "37"
    assert manifest.constraints.shape[1] == 4


def test_write_v2_handoff_manifests_writes_to_manifests_v2(tmp_path: Path) -> None:
    spec = _specs()[0]
    profile_config = default_v2_profile_configs([spec], seed=0)[0]
    manifest = build_v2_handoff_manifest_for_family(
        spec,
        profile_config,
        V2ManifestGenerationConfig(
            max_constraints=50,
            bootstrap_count=1,
            null_repetitions=1,
        ),
    )

    paths = write_v2_handoff_manifests([manifest], tmp_path / "manifests_v2")

    assert len(paths) == 1
    assert paths[0].parent.name == "manifests_v2"
    assert paths[0].exists()
