"""Unit tests for frozen P6b comparison statistics."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p6_diagnostics import _instrument_margin, roc_auc, spearman  # noqa: E402


def test_roc_auc_handles_perfect_reverse_and_ties():
    labels = [0, 0, 1, 1]
    assert roc_auc(labels, [0.0, 1.0, 2.0, 3.0]) == 1.0
    assert roc_auc(labels, [3.0, 2.0, 1.0, 0.0]) == 0.0
    assert roc_auc(labels, [1.0, 1.0, 1.0, 1.0]) == 0.5


def test_spearman_uses_average_tie_ranks():
    assert spearman([0.0, 1.0, 2.0], [2.0, 1.0, 0.0]) == pytest.approx(-1.0)
    assert spearman([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 1.0, 2.0]) > 0.9


def test_p1_instrument_margin_accepts_frozen_column_names():
    row = {
        "status": "ok",
        "heldout_violation": "0.025",
        "truth_order_error": "0.075",
    }
    assert _instrument_margin("P1", row) == pytest.approx(0.5)
