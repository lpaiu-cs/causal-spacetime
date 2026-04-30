from __future__ import annotations

import json

import pytest

from causal_spacetime_lab.state_change_manifest_v3_protocol_execution_spec import (
    load_v3_protocol_execution_specs,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_generation import (
    V3ProtocolManifestGenerationConfig,
    build_v3_protocol_handoff_manifest_for_config,
    build_v3_protocol_handoff_manifests,
    comparison_protocol_for_v3_protocol_family,
    write_v3_protocol_handoff_manifests,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_profiles import (
    default_v3_protocol_profile_configs,
)
from causal_spacetime_lab.state_change_response_handoff import manifest_to_jsonable
from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
    pairwise_response_dissimilarity_checked,
)
from causal_spacetime_lab.state_change_response_profile_metadata import (
    profile_metadata_from_protocols,
)
from causal_spacetime_lab.state_change_response_profiles import (
    EchoResponseProfileWithMetadata,
)


def test_comparison_protocol_for_v3_protocol_family(m40_prereg_path) -> None:
    specs = load_v3_protocol_execution_specs(m40_prereg_path)

    methods = {
        spec.family_name: comparison_protocol_for_v3_protocol_family(spec).method
        for spec in specs
    }

    assert methods["combined_earliest_full_reference_v3"] == (
        "combined_gap_and_mismatch"
    )
    assert methods["rank_gap_earliest_full_reference_v3"] == "rank_gap_mean"


def test_build_v3_protocol_handoff_manifest_includes_metadata(
    m40_prereg_path,
) -> None:
    specs = load_v3_protocol_execution_specs(m40_prereg_path)
    structured_spec = [spec for spec in specs if spec.family_kind == "structured"][0]
    config = [
        cfg
        for cfg in default_v3_protocol_profile_configs([structured_spec], seed=0)
        if cfg.family_name == structured_spec.family_name
    ][0]
    manifest = build_v3_protocol_handoff_manifest_for_config(
        structured_spec,
        config,
        V3ProtocolManifestGenerationConfig(
            max_constraints=50,
            bootstrap_count=1,
            null_repetitions=1,
        ),
    )
    jsonable = manifest_to_jsonable(manifest)

    assert jsonable["measurement_protocol_id"]
    assert jsonable["profile_metadata"]
    assert jsonable["handoff_provenance"]
    assert jsonable["handoff_provenance_type"] == (
        "hybrid_template_instantiated_from_profile"
    )


def test_pairwise_response_dissimilarity_checked_rejects_underspecified(
    m40_prereg_path,
) -> None:
    specs = load_v3_protocol_execution_specs(m40_prereg_path)
    config = default_v3_protocol_profile_configs(specs, seed=0)[0]
    wrapped = __import__(
        "causal_spacetime_lab.state_change_manifest_v3_protocol_profiles",
        fromlist=["build_v3_protocol_response_profile"],
    ).build_v3_protocol_response_profile(config)
    metadata = profile_metadata_from_protocols(
        "bad",
        [],
        [],
        "",
    )
    bad = EchoResponseProfileWithMetadata(wrapped.profile, metadata)

    with pytest.raises(ValueError):
        pairwise_response_dissimilarity_checked(
            bad,
            PairwiseResponseComparisonProtocol("p", "rank_gap_mean"),
        )


def test_write_v3_protocol_handoff_manifests_writes_to_manifest_dir(
    tmp_path,
    m40_prereg_path,
) -> None:
    specs = load_v3_protocol_execution_specs(m40_prereg_path)[:1]
    configs = default_v3_protocol_profile_configs(specs, seed=0)[:1]
    manifests = build_v3_protocol_handoff_manifests(
        specs,
        configs,
        V3ProtocolManifestGenerationConfig(
            max_constraints=30,
            bootstrap_count=1,
            null_repetitions=1,
        ),
    )
    paths = write_v3_protocol_handoff_manifests(
        manifests,
        tmp_path / "manifests_v3",
        overwrite=True,
    )

    assert paths
    assert all(path.parent.name == "manifests_v3" for path in paths)
    assert json.loads(paths[0].read_text(encoding="utf-8"))["handoff_provenance"]

