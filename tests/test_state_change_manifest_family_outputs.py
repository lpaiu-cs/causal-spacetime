from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_outputs import (
    family_names_from_bundle,
    load_family_output_bundle,
    missing_family_output_inputs,
)


def test_load_family_output_bundle_handles_missing_files(tmp_path) -> None:  # type: ignore[no-untyped-def]
    bundle = load_family_output_bundle(tmp_path)

    assert bundle["fit_summary"] == []
    assert "fit_summary" in missing_family_output_inputs(bundle)


def test_family_names_from_bundle_collects_known_names() -> None:
    bundle = {
        "fit_summary": [{"family_name": "a"}],
        "stricter_criteria": [{"family_name": "b"}],
        "report_card": [{"family_name": "a"}],
    }

    assert family_names_from_bundle(bundle) == ["a", "b"]
