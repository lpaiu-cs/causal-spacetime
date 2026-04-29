from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_response_underdetermination import (
    generate_rank_preserving_1d_layouts,
    pairwise_distance_order_disagreement,
    rank_preserving_scalar_values,
    response_order_preserved,
    single_reference_underdetermination_report,
)


def test_rank_preserving_scalar_values_preserve_order() -> None:
    delays = np.asarray([3, 1, 3, 2])
    reachable = np.ones(4, dtype=bool)

    values = rank_preserving_scalar_values(delays, reachable, spacing="exponential")

    assert response_order_preserved(delays, reachable, values)


def test_response_order_preserved_rejects_bad_values() -> None:
    delays = np.asarray([1, 2, 3])
    reachable = np.ones(3, dtype=bool)

    assert not response_order_preserved(delays, reachable, np.asarray([1.0, 3.0, 2.0]))


def test_generate_rank_preserving_1d_layouts_returns_multiple_layouts() -> None:
    delays = np.asarray([1, 2, 3, 4])
    reachable = np.ones(4, dtype=bool)

    layouts = generate_rank_preserving_1d_layouts(
        delays,
        reachable,
        layout_count=4,
        seed=0,
    )

    assert len(layouts) == 4
    assert all(
        response_order_preserved(delays, reachable, layout) for layout in layouts
    )


def test_pairwise_distance_order_disagreement_detects_difference() -> None:
    reachable = np.ones(4, dtype=bool)
    layout_a = np.asarray([1.0, 2.0, 3.0, 4.0])
    layout_b = np.asarray([1.0, 2.0, 100.0, 101.0])

    disagreement = pairwise_distance_order_disagreement(layout_a, layout_b, reachable)

    assert disagreement > 0.0


def test_single_reference_underdetermination_report_fields() -> None:
    report = single_reference_underdetermination_report(
        np.asarray([1, 2, 3, 4]),
        np.ones(4, dtype=bool),
        layout_count=4,
        seed=0,
    )

    assert report["reachable_count"] == 4.0
    assert report["response_order_preserved_fraction"] == 1.0
    assert report["max_pair_distance_order_disagreement"] > 0.0
