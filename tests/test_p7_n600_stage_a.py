"""Focused tests for P7 N=600 summary statistics."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p7_n600_stage_a import _binder, _seed  # noqa: E402


def test_seed_formula_separates_beta_and_chain():
    assert _seed(900000, 14.0, 0) == 914000
    assert _seed(900000, 14.0, 1) == 914001
    assert _seed(900000, 18.0, 0) != _seed(900000, 14.0, 0)


def test_centered_binder_is_finite_for_nonconstant_series():
    value = _binder([1.0, 2.0, 3.0, 4.0])
    assert value == pytest.approx(0.45333333333333337)
