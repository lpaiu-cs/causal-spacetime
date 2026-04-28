from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_bracket_diagnostics import (
    bracket_coverage_summary,
    compare_bracket_rank_orders,
    rank_slice_summary,
)


def test_bracket_coverage_summary_fields_and_counts() -> None:
    predecessors = np.asarray([-1, 0, 0, -1])
    successors = np.asarray([-1, 2, 3, 1])
    accessible = np.asarray([False, True, True, False])

    summary = bracket_coverage_summary(
        predecessors,
        successors,
        accessible,
        reference_chain_length=4,
    )

    assert summary["target_count"] == 4.0
    assert summary["accessible_count"] == 2.0
    assert summary["accessible_fraction"] == 0.5
    assert summary["successor_only_count"] == 1.0
    assert summary["mean_bracket_width_rank"] == pytest.approx(2.5)


def test_rank_slice_summary_fields() -> None:
    summary = rank_slice_summary(np.asarray([-1, 0, 0, 1, 2]))

    assert summary["slice_count"] == 3.0
    assert summary["assigned_count"] == 4.0
    assert summary["max_slice_size"] == 2.0
    assert summary["singleton_slice_fraction"] == pytest.approx(2 / 3)


def test_compare_bracket_rank_orders_identical_has_zero_inversion() -> None:
    ranks = np.asarray([-1, 1, 2, 3])
    accessible = np.asarray([False, True, True, True])

    comparison = compare_bracket_rank_orders(ranks, ranks, accessible, accessible)

    assert comparison["common_accessible_count"] == 3.0
    assert comparison["order_inversion_rate"] == 0.0
    assert comparison["order_agreement_rate"] == 1.0
