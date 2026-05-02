from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v4_stability_drilldown import (
    summarize_v4_stability_failures,
)


def test_stability_drilldown_reports_latent_order_failure() -> None:
    rows = summarize_v4_stability_failures(
        [{"family_name": "family", "restart_std": 0.01}],
        [{"family_name": "family", "latent_order_disagreement": 0.9}],
        default_cross_family_robustness_criteria(),
    )

    assert "latent_order_instability" in rows[0]["dominant_stability_failure"]
