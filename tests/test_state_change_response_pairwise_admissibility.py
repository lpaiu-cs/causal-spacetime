from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_pairwise import (
    PairwiseResponseComparisonProtocol,
)
from causal_spacetime_lab.state_change_response_pairwise_admissibility import (
    compare_pairwise_protocols,
    pairwise_protocol_admissibility_report,
)
from causal_spacetime_lab.state_change_response_profiles import EchoResponseProfile


def _profile() -> EchoResponseProfile:
    delays = np.asarray(
        [[1, 2, 3], [1, 3, 4], [2, 3, 4], [1, -1, 3]],
        dtype=int,
    )
    return EchoResponseProfile(
        target_event_ids=np.arange(4, dtype=int),
        protocol_labels=["a", "b", "c"],
        delay_rank_matrix=delays,
        reachable_matrix=delays >= 0,
    )


def test_pairwise_protocol_admissibility_report_fields() -> None:
    report = pairwise_protocol_admissibility_report(
        _profile(),
        PairwiseResponseComparisonProtocol("gap", "rank_gap_mean"),
    )

    assert report["protocol_name"] == "gap"
    assert "valid_pair_fraction" in report
    assert report["symmetric"] == 1.0


def test_compare_pairwise_protocols_deterministic_case() -> None:
    rows = compare_pairwise_protocols(
        _profile(),
        [
            PairwiseResponseComparisonProtocol("gap", "rank_gap_mean"),
            PairwiseResponseComparisonProtocol("sep", "separation_fraction"),
        ],
    )

    assert len(rows) == 1
    assert rows[0]["protocol_a"] == "gap"
    assert "order_inversion_rate" in rows[0]
