from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    default_new_manifest_family_specs_v2,
)
from causal_spacetime_lab.state_change_manifest_v2_profiles import (
    build_v2_response_profile,
    default_v2_profile_configs,
)
from causal_spacetime_lab.state_change_manifest_v2_spec import V2ManifestFamilySpec


def _specs() -> list[V2ManifestFamilySpec]:
    return [
        V2ManifestFamilySpec(**row)
        for row in default_new_manifest_family_specs_v2()
    ]


def test_default_v2_profile_configs_returns_all_planned_families() -> None:
    specs = _specs()

    configs = default_v2_profile_configs(specs, seed=12)

    assert {config.family_name for config in configs} == {
        spec.family_name for spec in specs
    }


def test_build_v2_response_profile_deterministic_shape() -> None:
    config = default_v2_profile_configs(_specs(), seed=7)[0]

    profile_a = build_v2_response_profile(config)
    profile_b = build_v2_response_profile(config)

    assert profile_a.delay_rank_matrix.shape[1] == config.protocol_column_count
    assert profile_a.target_event_ids.ndim == 1
    assert np.array_equal(profile_a.delay_rank_matrix, profile_b.delay_rank_matrix)
    assert np.array_equal(profile_a.reachable_matrix, profile_b.reachable_matrix)
