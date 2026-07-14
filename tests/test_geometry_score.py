"""Tests for the P7 geometry order parameter."""

import pytest

from causal_spacetime_lab.positive_control.geometry_score import (
    clipped_gate_margin_score,
    geometry_order_parameter,
    minimum_gate_margin,
)


def test_score_boundary_and_saturation():
    assert clipped_gate_margin_score([0.0, 0.2]) == pytest.approx(0.5)
    assert clipped_gate_margin_score([0.5, 1.0]) == pytest.approx(1.0)
    assert clipped_gate_margin_score([-0.5, 1.0]) == pytest.approx(0.0)


def test_geometry_score_is_equivalent_to_frozen_gates_at_threshold():
    passing = geometry_order_parameter(
        status="ok", heldout=0.09, null_gap=0.11, truth_error=0.39
    )
    heldout_fail = geometry_order_parameter(
        status="ok", heldout=0.11, null_gap=0.20, truth_error=0.20
    )
    assert passing >= 0.5
    assert heldout_fail < 0.5


def test_structural_block_is_total_zero_without_metrics():
    assert geometry_order_parameter(status="structural_block: only 2 chains") == 0.0


def test_valid_score_requires_all_gate_quantities():
    with pytest.raises(ValueError, match="all three"):
        geometry_order_parameter(status="ok", heldout=0.02, null_gap=0.2)


def test_minimum_margin_supports_partial_source_specific_gates():
    assert minimum_gate_margin(
        heldout=0.025,
        heldout_max=0.05,
        truth_error=0.075,
        truth_error_max=0.15,
    ) == pytest.approx(0.5)


def test_geometry_score_is_raw_margin_with_fixed_offset_and_clipping():
    margin = minimum_gate_margin(
        heldout=0.09,
        heldout_max=0.10,
        null_gap=0.11,
        null_gap_min=0.10,
        truth_error=0.39,
        truth_error_max=0.40,
    )
    assert geometry_order_parameter(
        status="ok", heldout=0.09, null_gap=0.11, truth_error=0.39
    ) == pytest.approx(clipped_gate_margin_score((margin,)))
