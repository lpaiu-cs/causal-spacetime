"""Exact sanity checks for handoff provenance metadata."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_handoff_provenance import (
    default_bottom_up_handoff_provenance,
    default_hybrid_handoff_provenance,
    default_report_only_handoff_provenance,
    default_top_down_handoff_provenance,
    validate_handoff_provenance,
)


def run_experiment() -> list[dict[str, float | str]]:
    bottom = default_bottom_up_handoff_provenance(
        design_source_label="bottom",
        design_digest="abc",
    )
    top = default_top_down_handoff_provenance(
        design_source_label="top",
        design_source_path="plan.json",
        design_digest="abc",
        template_id="t",
        template_hash="h",
    )
    hybrid = default_hybrid_handoff_provenance(
        design_source_label="hybrid",
        design_source_path="plan.json",
        design_digest="abc",
        constraint_selection_policy="profile_instantiated_rank_gap",
        template_id="t",
        template_hash="h",
    )
    fit_bad = replace(hybrid, fit_outputs_used=True)
    carry_bad = replace(hybrid, carry_forward_outputs_used=True)
    report = default_report_only_handoff_provenance(design_source_label="control")
    return [
        {
            "check": "bottom_up_validates",
            "passed": validate_handoff_provenance(bottom)["valid_provenance"],
        },
        {
            "check": "top_down_validates",
            "passed": validate_handoff_provenance(top)["valid_provenance"],
        },
        {
            "check": "hybrid_validates",
            "passed": validate_handoff_provenance(hybrid)["valid_provenance"],
        },
        {
            "check": "fit_outputs_used_blocked",
            "passed": float(
                validate_handoff_provenance(fit_bad)["valid_provenance"] == 0.0
            ),
        },
        {
            "check": "carry_forward_outputs_used_blocked",
            "passed": float(
                validate_handoff_provenance(carry_bad)["valid_provenance"] == 0.0
            ),
        },
        {
            "check": "report_only_control_requires_report_only",
            "passed": validate_handoff_provenance(report)["report_only"],
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/handoff_provenance_exact_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote handoff provenance sanity: {path}")


if __name__ == "__main__":
    main()
