from __future__ import annotations

from pathlib import Path

import pytest

from causal_spacetime_lab.state_change_manifest_v3_protocol_execution_spec import (
    load_v3_protocol_execution_specs,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_generation import (
    V3ProtocolManifestGenerationConfig,
    build_v3_protocol_handoff_manifests,
    write_v3_protocol_handoff_manifests,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    default_v3_protocol_invariant_family_patches,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preregistration import (
    build_v3_protocol_patched_preregistration,
    write_v3_protocol_patched_preregistration,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_profiles import (
    default_v3_protocol_profile_configs,
)
from causal_spacetime_lab.state_change_manifest_v4_design import (
    default_v4_protocol_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v4_execution_spec import (
    load_v4_protocol_execution_specs,
)
from causal_spacetime_lab.state_change_manifest_v4_generation import (
    V4ProtocolManifestGenerationConfig,
    build_v4_protocol_handoff_manifests,
    write_v4_protocol_handoff_manifests,
)
from causal_spacetime_lab.state_change_manifest_v4_preregistration import (
    build_v4_protocol_preregistration_spec,
    write_v4_protocol_preregistration_spec,
)
from causal_spacetime_lab.state_change_manifest_v4_profiles import (
    default_v4_protocol_profile_configs,
)


@pytest.fixture()
def m40_prereg_path(tmp_path: Path) -> Path:
    """Write a deterministic M40 patched preregistration JSON for tests."""

    path = tmp_path / "remediation" / "v3_protocol_patched_preregistration_m40.json"
    spec = build_v3_protocol_patched_preregistration(
        tmp_path / "remediation" / "v3_preregistration_spec_m39.json",
        default_v3_protocol_invariant_family_patches(),
    )
    return write_v3_protocol_patched_preregistration(spec, path)


@pytest.fixture()
def m41_manifest_dir(tmp_path: Path, m40_prereg_path: Path) -> Path:
    """Write a small deterministic patched v3 manifest directory for tests."""

    specs = load_v3_protocol_execution_specs(m40_prereg_path)
    selected = [spec for spec in specs if spec.family_kind == "structured"][:1]
    profile_configs = default_v3_protocol_profile_configs(selected, seed=0)[:2]
    manifests = build_v3_protocol_handoff_manifests(
        selected,
        profile_configs,
        V3ProtocolManifestGenerationConfig(
            max_constraints=40,
            bootstrap_count=1,
            null_repetitions=1,
        ),
    )
    manifest_dir = tmp_path / "manifests_v3"
    write_v3_protocol_handoff_manifests(manifests, manifest_dir, overwrite=True)
    return manifest_dir


@pytest.fixture()
def m43_v4_prereg_path(tmp_path: Path) -> Path:
    """Write a deterministic M43 v4 preregistration JSON for tests."""

    path = tmp_path / "remediation" / "v4_protocol_preregistration_spec_m43.json"
    spec = build_v4_protocol_preregistration_spec(
        default_v4_protocol_family_designs()
    )
    return write_v4_protocol_preregistration_spec(spec, path)


@pytest.fixture()
def m44_manifest_dir(tmp_path: Path, m43_v4_prereg_path: Path) -> Path:
    """Write a small deterministic v4 manifest directory for tests."""

    specs = load_v4_protocol_execution_specs(m43_v4_prereg_path)
    selected = [spec for spec in specs if spec.family_kind == "structured"][:1]
    profile_configs = default_v4_protocol_profile_configs(selected, seed=0)[:2]
    manifests = build_v4_protocol_handoff_manifests(
        selected,
        profile_configs,
        V4ProtocolManifestGenerationConfig(
            max_constraints=40,
            bootstrap_count=1,
            null_repetitions=1,
        ),
    )
    manifest_dir = tmp_path / "manifests_v4"
    write_v4_protocol_handoff_manifests(manifests, manifest_dir, overwrite=True)
    return manifest_dir
