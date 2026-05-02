from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v4_coverage_drilldown import (
    summarize_v4_coverage_failures,
)


def test_coverage_drilldown_reports_pair_node_coverage_failure() -> None:
    rows = summarize_v4_coverage_failures(
        [
            {
                "family_name": "family",
                "target_coverage_fraction": 1.0,
                "pair_node_coverage_fraction": 0.0,
            }
        ],
        default_cross_family_robustness_criteria(),
    )

    assert "pair_node_coverage_failure" in rows[0]["dominant_coverage_failure"]
