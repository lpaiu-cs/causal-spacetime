from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v3_protocol_execution_spec import (
    load_v3_protocol_execution_specs,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_generation import (
    top_down_template_for_v3_protocol_family,
)
from causal_spacetime_lab.state_change_top_down_handoff_template import (
    top_down_template_digest,
    top_down_template_jsonable,
)


def test_top_down_template_digest_is_stable(m40_prereg_path) -> None:
    spec = load_v3_protocol_execution_specs(m40_prereg_path)[0]
    template = top_down_template_for_v3_protocol_family(spec)

    assert top_down_template_digest(template) == top_down_template_digest(template)


def test_top_down_template_includes_protocol_method_and_margin(m40_prereg_path) -> None:
    spec = load_v3_protocol_execution_specs(m40_prereg_path)[0]
    template = top_down_template_for_v3_protocol_family(spec)
    jsonable = top_down_template_jsonable(template)

    assert jsonable["measurement_protocol_id"]
    assert jsonable["measurement_protocol_hash"]
    assert jsonable["comparison_method"] in {
        "rank_gap_mean",
        "combined_gap_and_mismatch",
    }
    assert float(jsonable["margin_value"]) > 0.0
    assert jsonable["forbidden_interpretations"]

