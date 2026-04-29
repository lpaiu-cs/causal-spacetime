from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_constraint_coverage import (
    constraint_pair_node_coverage,
    constraint_pool_summary,
    constraint_target_coverage,
)
from causal_spacetime_lab.state_change_response_constraint_pool import (
    ResponseComparisonConstraintPool,
)


def _pool() -> ResponseComparisonConstraintPool:
    return ResponseComparisonConstraintPool(
        target_event_ids=np.asarray([10, 11, 12, 13], dtype=int),
        constraints=np.asarray(
            [
                [0, 1, 2, 3],
                [0, 2, 1, 3],
                [0, 1, 0, 3],
            ],
            dtype=int,
        ),
        margins=np.asarray([0.2, 0.4, 0.6], dtype=float),
        protocol_name="gap",
        method="rank_gap_mean",
        source_label="test",
    )


def test_constraint_target_coverage_deterministic_case() -> None:
    report = constraint_target_coverage(_pool(), target_count=4)

    assert report["touched_target_count"] == 4.0
    assert report["touched_target_fraction"] == 1.0


def test_constraint_pair_node_coverage_deterministic_case() -> None:
    report = constraint_pair_node_coverage(_pool(), target_count=4)

    assert report["possible_pair_node_count"] == 6.0
    assert report["touched_pair_node_count"] == 5.0


def test_constraint_pool_summary_fields() -> None:
    report = constraint_pool_summary(_pool(), target_count=4)

    assert report["constraint_count"] == 3.0
    assert report["mean_margin"] == 0.4000000000000001
    assert "touched_pair_node_fraction" in report
