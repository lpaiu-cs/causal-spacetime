from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v4_execution_spec import (
    load_v4_protocol_execution_specs,
)
from causal_spacetime_lab.state_change_manifest_v4_generation import (
    V4ProtocolManifestGenerationConfig,
    build_v4_protocol_handoff_manifest_for_config,
    write_v4_protocol_handoff_manifests,
)
from causal_spacetime_lab.state_change_manifest_v4_profiles import (
    default_v4_protocol_profile_configs,
)
from causal_spacetime_lab.state_change_response_handoff import manifest_to_jsonable


def test_build_v4_protocol_handoff_manifest_includes_metadata(
    m43_v4_prereg_path,
) -> None:
    specs = load_v4_protocol_execution_specs(m43_v4_prereg_path)
    spec = [item for item in specs if item.family_kind == "structured"][0]
    config = default_v4_protocol_profile_configs([spec], seed=0)[0]

    manifest = build_v4_protocol_handoff_manifest_for_config(
        spec,
        config,
        V4ProtocolManifestGenerationConfig(
            max_constraints=40,
            bootstrap_count=1,
            null_repetitions=1,
        ),
    )
    jsonable = manifest_to_jsonable(manifest)

    assert isinstance(jsonable["measurement_protocol"], dict)
    assert isinstance(jsonable["profile_metadata"], dict)
    assert isinstance(jsonable["handoff_provenance"], dict)
    assert jsonable["profile_invariance_status"] == "protocol_invariant"
    assert jsonable["handoff_provenance_type"]


def test_v4_report_only_control_is_ineligible(m43_v4_prereg_path) -> None:
    specs = load_v4_protocol_execution_specs(m43_v4_prereg_path)
    spec = [item for item in specs if item.family_kind == "report_only"][0]
    config = default_v4_protocol_profile_configs([spec], seed=0)[0]

    manifest = build_v4_protocol_handoff_manifest_for_config(
        spec,
        config,
        V4ProtocolManifestGenerationConfig(max_constraints=20),
    )

    assert not manifest.handoff_decision.eligible
    assert manifest.handoff_provenance_type == "report_only_control"


def test_write_v4_protocol_handoff_manifests_only_allows_v4_directory(
    m43_v4_prereg_path,
    tmp_path: Path,
) -> None:
    specs = load_v4_protocol_execution_specs(m43_v4_prereg_path)
    spec = [item for item in specs if item.family_kind == "structured"][0]
    config = default_v4_protocol_profile_configs([spec], seed=0)[0]
    manifest = build_v4_protocol_handoff_manifest_for_config(
        spec,
        config,
        V4ProtocolManifestGenerationConfig(max_constraints=20),
    )

    paths = write_v4_protocol_handoff_manifests(
        [manifest],
        tmp_path / "manifests_v4",
        overwrite=True,
    )

    assert paths
    assert paths[0].parent.name == "manifests_v4"
