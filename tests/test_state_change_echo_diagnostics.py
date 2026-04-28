from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.state_change_echo_diagnostics import (
    echo_delay_histogram,
    echo_order_resolution_summary,
    echo_reachability_summary,
)


def test_echo_reachability_summary_expected_fields() -> None:
    returns = np.asarray([-1, 2, 3, 3])
    delays = np.asarray([-1, 2, 3, 3])
    reachable = np.asarray([False, True, True, True])

    summary = echo_reachability_summary(
        returns,
        delays,
        reachable,
        emission_position=0,
        reference_chain_length=4,
    )

    assert summary["target_count"] == 4.0
    assert summary["reachable_count"] == 3.0
    assert summary["reachable_fraction"] == 0.75
    assert summary["mean_echo_delay_rank"] == pytest.approx(8 / 3)
    assert summary["distinct_delay_rank_count"] == 2.0


def test_echo_delay_histogram_deterministic_example() -> None:
    histogram = echo_delay_histogram(
        np.asarray([-1, 2, 3, 3]),
        np.asarray([False, True, True, True]),
    )

    assert histogram == {2: 1, 3: 2}


def test_echo_order_resolution_summary_deterministic_example() -> None:
    summary = echo_order_resolution_summary(
        np.asarray([-1, 2, 3, 3]),
        np.asarray([False, True, True, True]),
    )

    assert summary["reachable_count"] == 3.0
    assert summary["comparable_pair_count"] == 3.0
    assert summary["tied_pair_fraction"] == pytest.approx(1 / 3)
    assert summary["strict_order_pair_fraction"] == pytest.approx(2 / 3)
