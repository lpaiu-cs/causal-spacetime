"""Exact sanity checks for top-down handoff templates."""

from __future__ import annotations

from pathlib import Path

from manifest_family_experiment_helpers import write_csv

from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    default_v3_protocol_invariant_family_patches,
)
from causal_spacetime_lab.state_change_top_down_handoff_template import (
    template_from_v3_protocol_family_spec,
    top_down_template_digest,
    top_down_template_jsonable,
)


def run_experiment() -> list[dict[str, float | str]]:
    patch = [
        item
        for item in default_v3_protocol_invariant_family_patches()
        if item.family_kind == "structured"
    ][0]
    template = template_from_v3_protocol_family_spec(patch)
    jsonable = top_down_template_jsonable(template)
    digest = top_down_template_digest(template)
    forbidden_text = " ".join(template.forbidden_interpretations).lower()
    field_text = " ".join(jsonable).lower()
    return [
        {
            "check": "template_digest_stable",
            "passed": float(digest == top_down_template_digest(template)),
        },
        {
            "check": "template_includes_protocol_id",
            "passed": float(bool(template.measurement_protocol_id)),
        },
        {
            "check": "template_includes_comparison_method",
            "passed": float(bool(template.comparison_method)),
        },
        {
            "check": "template_includes_margin_policy_and_value",
            "passed": float(
                bool(template.margin_policy) and template.margin_value >= 0.0
            ),
        },
        {
            "check": "template_has_forbidden_interpretations",
            "passed": float(bool(template.forbidden_interpretations)),
        },
        {
            "check": "template_has_no_fit_or_metric_fields",
            "passed": float(
                "fit" not in field_text
                and "metric" not in field_text
                and "metric reconstruction" in forbidden_text
            ),
        },
    ]


def main() -> None:
    path = write_csv(
        run_experiment(),
        Path("outputs/data/top_down_handoff_template_exact_sanity.csv"),
        ["check", "passed"],
    )
    print(f"Wrote top-down handoff template sanity: {path}")


if __name__ == "__main__":
    main()
