"""Statistics tests for the P6 near-critical characterization."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

pytest.importorskip("numba")

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p6_near_critical import integrated_autocorrelation  # noqa: E402


def test_iat_reports_full_ess_for_constant_or_alternating_series():
    tau, ess = integrated_autocorrelation([2.0] * 20)
    assert tau == 1.0
    assert ess == 20.0
    tau, ess = integrated_autocorrelation([0.0, 1.0] * 10)
    assert tau == 1.0
    assert ess == 20.0


def test_iat_detects_positive_serial_correlation():
    tau, ess = integrated_autocorrelation([float(index // 2) for index in range(20)])
    assert tau > 1.0
    assert ess < 20.0
