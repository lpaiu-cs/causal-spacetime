from __future__ import annotations

import json
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v3_protocol_outputs import (
    load_v3_protocol_output_bundle,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preconditions import (
    evaluate_v3_protocol_preconditions,
)
from tests.v3_protocol_test_helpers import (
    family_from_manifest_dir,
    write_v3_protocol_bundle,
)


def test_evaluate_v3_protocol_preconditions_passes_valid_structured_manifest(
    tmp_path: Path,
    m41_manifest_dir: Path,
) -> None:
    family = family_from_manifest_dir(m41_manifest_dir)
    write_v3_protocol_bundle(tmp_path, family)
    bundle = load_v3_protocol_output_bundle(tmp_path)

    reports = evaluate_v3_protocol_preconditions(m41_manifest_dir, bundle)

    assert reports[0].preconditions_passed


def test_evaluate_v3_protocol_preconditions_blocks_missing_provenance(
    tmp_path: Path,
    m41_manifest_dir: Path,
) -> None:
    family = family_from_manifest_dir(m41_manifest_dir)
    broken_dir = tmp_path / "broken_manifests"
    broken_dir.mkdir()
    payload = json.loads(sorted(m41_manifest_dir.glob("*.json"))[0].read_text())
    payload.pop("handoff_provenance", None)
    (broken_dir / "broken.json").write_text(json.dumps(payload), encoding="utf-8")
    write_v3_protocol_bundle(tmp_path, family)

    reports = evaluate_v3_protocol_preconditions(
        broken_dir,
        load_v3_protocol_output_bundle(tmp_path),
    )

    assert "missing_handoff_provenance" in reports[0].failed_preconditions


def test_evaluate_v3_protocol_preconditions_blocks_protocol_mixed(
    tmp_path: Path,
    m41_manifest_dir: Path,
) -> None:
    family = family_from_manifest_dir(m41_manifest_dir)
    broken_dir = tmp_path / "mixed_manifests"
    broken_dir.mkdir()
    payload = json.loads(sorted(m41_manifest_dir.glob("*.json"))[0].read_text())
    payload["profile_invariance_status"] = "protocol_mixed"
    (broken_dir / "mixed.json").write_text(json.dumps(payload), encoding="utf-8")
    write_v3_protocol_bundle(tmp_path, family)

    reports = evaluate_v3_protocol_preconditions(
        broken_dir,
        load_v3_protocol_output_bundle(tmp_path),
    )

    assert "protocol_mixed" in reports[0].failed_preconditions

