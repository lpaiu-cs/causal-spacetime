from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_null_drilldown import (
    summarize_v3_protocol_null_taxonomy_failures,
)


def test_null_taxonomy_drilldown_separates_destructive_and_symmetry() -> None:
    rows = summarize_v3_protocol_null_taxonomy_failures(
        [
            {
                "family_name": "family",
                "taxonomy_class": "destructive_null",
                "mean_heldout_violation_rate": 0.4,
                "structured_heldout_violation_rate": 0.2,
            },
            {
                "family_name": "family",
                "taxonomy_class": "symmetry_control",
                "mean_heldout_violation_rate": 0.21,
                "structured_heldout_violation_rate": 0.2,
            },
        ],
        [{"family_name": "family", "mean_heldout_violation": 0.2}],
        default_cross_family_robustness_criteria(),
    )

    assert rows[0]["destructive_null_gap_pass"] == 1.0
    assert rows[0]["symmetry_control_gap_pass"] == 1.0
