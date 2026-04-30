from __future__ import annotations

import csv
from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v2_outputs import (
    load_v2_output_bundle,
    missing_v2_bundle_inputs,
    v2_bundle_input_report,
)


def test_load_v2_output_bundle_handles_missing_files(tmp_path: Path) -> None:
    bundle = load_v2_output_bundle(tmp_path)

    assert bundle["metrics"] == []
    assert "metrics" in missing_v2_bundle_inputs(bundle)


def test_missing_v2_bundle_inputs_detects_absent_required_files(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    with (data_dir / "v2_cross_family_robustness_metrics.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["family_name"])
        writer.writeheader()
        writer.writerow({"family_name": "family"})

    bundle = load_v2_output_bundle(tmp_path)
    missing = missing_v2_bundle_inputs(bundle)

    assert "metrics" not in missing
    assert "completeness" in missing
    assert any(row["input_name"] == "metrics" for row in v2_bundle_input_report(bundle))
