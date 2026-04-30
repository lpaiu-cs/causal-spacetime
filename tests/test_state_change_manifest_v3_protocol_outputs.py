from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v3_protocol_outputs import (
    load_v3_protocol_output_bundle,
    missing_v3_protocol_bundle_inputs,
    v3_protocol_bundle_input_report,
)
from tests.v3_protocol_test_helpers import write_csv


def test_load_v3_protocol_output_bundle_handles_missing_files(tmp_path: Path) -> None:
    bundle = load_v3_protocol_output_bundle(tmp_path)

    assert bundle["metrics"] == []
    assert "metrics" in missing_v3_protocol_bundle_inputs(bundle)


def test_missing_v3_protocol_bundle_inputs_detects_absent_required_inputs(
    tmp_path: Path,
) -> None:
    write_csv(
        tmp_path / "data" / "v3_protocol_cross_family_robustness_metrics.csv",
        [{"family_name": "family"}],
    )

    bundle = load_v3_protocol_output_bundle(tmp_path)
    missing = missing_v3_protocol_bundle_inputs(bundle)

    assert "metrics" not in missing
    assert "completeness" in missing
    assert any(
        row["input_name"] == "metrics"
        for row in v3_protocol_bundle_input_report(bundle)
    )
