from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_remediation_plan import (
    default_new_manifest_family_specs_v2,
)
from causal_spacetime_lab.state_change_manifest_v2_diagnostics import (
    compute_v2_coverage_rows,
    compute_v2_failed_accounting_rows,
)
from causal_spacetime_lab.state_change_manifest_v2_generation import (
    V2ManifestGenerationConfig,
    build_v2_handoff_manifests,
    write_v2_handoff_manifests,
)
from causal_spacetime_lab.state_change_manifest_v2_profiles import (
    default_v2_profile_configs,
)
from causal_spacetime_lab.state_change_manifest_v2_spec import V2ManifestFamilySpec


def _write_test_manifests(tmp_path: Path) -> Path:
    specs = [
        V2ManifestFamilySpec(**row)
        for row in default_new_manifest_family_specs_v2()[:2]
    ]
    configs = default_v2_profile_configs(specs, seed=0)
    manifests = build_v2_handoff_manifests(
        specs,
        configs,
        V2ManifestGenerationConfig(
            max_constraints=80,
            bootstrap_count=1,
            null_repetitions=1,
        ),
    )
    manifest_dir = tmp_path / "manifests_v2"
    write_v2_handoff_manifests(manifests, manifest_dir)
    return manifest_dir


def test_compute_v2_coverage_rows_returns_target_and_pair_coverage(
    tmp_path: Path,
) -> None:
    rows = compute_v2_coverage_rows(_write_test_manifests(tmp_path))

    assert rows
    assert {"target_coverage_fraction", "pair_node_coverage_fraction"} <= set(rows[0])


def test_compute_v2_failed_accounting_rows(tmp_path: Path) -> None:
    rows = compute_v2_failed_accounting_rows(_write_test_manifests(tmp_path))

    assert rows
    assert any(row["row_type"] == "failed_accounting" for row in rows)
