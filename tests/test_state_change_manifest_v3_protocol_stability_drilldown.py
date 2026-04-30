from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_stability_drilldown import (
    summarize_v3_protocol_stability_failures,
)


def test_stability_drilldown_reports_latent_order_failure() -> None:
    rows = summarize_v3_protocol_stability_failures(
        [{"family_name": "family", "restart_std": 0.01}],
        [{"family_name": "family", "latent_order_disagreement": 0.6}],
        default_cross_family_robustness_criteria(),
    )

    assert rows[0]["latent_order_pass"] == 0.0
    assert "latent_order_instability" in str(rows[0]["dominant_stability_failure"])
