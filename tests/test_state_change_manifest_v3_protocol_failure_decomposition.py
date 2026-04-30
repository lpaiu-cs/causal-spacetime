from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_family_robustness import (
    default_cross_family_robustness_criteria,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_failure_decomposition import (  # noqa: E501
    decompose_v3_protocol_family_failures,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preconditions import (
    V3ProtocolPreconditionReport,
)
from tests.test_state_change_manifest_v3_protocol_carry_forward import _row


def test_decompose_v3_protocol_family_failures_separates_metric_and_provenance() -> (
    None
):
    records = decompose_v3_protocol_family_failures(
        [_row(heldout=0.5)],
        default_cross_family_robustness_criteria(),
        [
            V3ProtocolPreconditionReport(
                family_name="family",
                family_kind="structured",
                manifest_count=3,
                diagnostic_complete=True,
                all_manifests_have_measurement_protocol=True,
                all_manifests_have_profile_metadata=True,
                all_manifests_have_handoff_provenance=False,
                all_structured_protocol_invariant=True,
                all_structured_parameter_complete=True,
                all_structured_admissible_for_pairwise_dissimilarity=True,
                all_structured_valid_provenance=False,
                report_only_controls_ineligible=True,
                preconditions_passed=False,
                failed_preconditions=["missing_handoff_provenance"],
                warning_preconditions=[],
            )
        ],
    )
    roots = {record.root_cause_category for record in records}

    assert "heldout_failure" in roots
    assert "provenance_failure" in roots
