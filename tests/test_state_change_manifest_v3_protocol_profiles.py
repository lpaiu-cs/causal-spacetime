from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v3_protocol_execution_spec import (
    load_v3_protocol_execution_specs,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_profiles import (
    build_v3_protocol_response_profile,
    default_v3_protocol_profile_configs,
)


def test_default_v3_protocol_profile_configs_emit_planned_counts(
    m40_prereg_path,
) -> None:
    specs = load_v3_protocol_execution_specs(m40_prereg_path)
    configs = default_v3_protocol_profile_configs(specs, seed=0)

    for spec in specs:
        family_configs = [cfg for cfg in configs if cfg.family_name == spec.family_name]
        assert len(family_configs) == spec.planned_manifest_count


def test_build_v3_protocol_response_profile_is_protocol_invariant(
    m40_prereg_path,
) -> None:
    specs = load_v3_protocol_execution_specs(m40_prereg_path)
    config = default_v3_protocol_profile_configs(specs, seed=0)[0]
    wrapped = build_v3_protocol_response_profile(config)

    assert wrapped.metadata.profile_invariance_status == "protocol_invariant"
    assert wrapped.metadata.admissible_for_pairwise_dissimilarity
    assert wrapped.profile.delay_rank_matrix.shape[1] == len(config.reference_chain_ids)


def test_failed_control_profile_keeps_explicit_protocol_metadata(
    m40_prereg_path,
) -> None:
    specs = load_v3_protocol_execution_specs(m40_prereg_path)
    config = [
        cfg
        for cfg in default_v3_protocol_profile_configs(specs, seed=0)
        if cfg.family_kind == "failed_control"
    ][0]
    wrapped = build_v3_protocol_response_profile(config)

    assert wrapped.metadata.profile_invariance_status == "protocol_invariant"
    assert wrapped.metadata.measurement_protocol_id
