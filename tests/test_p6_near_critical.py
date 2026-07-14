"""Statistics tests for the P6 near-critical characterization."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from causal_spacetime_lab.positive_control.mcmc_diagnostics import (
    _initial_positive_pair_tau,
    integrated_autocorrelation,
)

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p6_near_critical import _validate_aggregate_inputs  # noqa: E402


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


def test_geyer_iat_stops_on_nonpositive_pairs_not_single_lags():
    tau = _initial_positive_pair_tau([1.0, 0.4, -0.1, 0.3, -0.4, -0.2])
    assert tau == pytest.approx(2.2)


def _rows(count: int, *, instrument: bool = False) -> list[dict[str, str]]:
    indices = [0, count // 2, count - 1] if instrument else range(count)
    return [
        {"beta": "12", "chain": "0", "sample": str(index)} for index in indices
    ]


def test_aggregate_rejects_partial_instrument_shard():
    with pytest.raises(SystemExit, match="instrument.*sample indices"):
        _validate_aggregate_inputs(
            {(12.0, 0)},
            _rows(45),
            _rows(45, instrument=True)[:-1],
            minimum_samples=45,
            nominal_samples=48,
        )


def test_aggregate_accepts_complete_chain_and_three_snapshots():
    _validate_aggregate_inputs(
        {(12.0, 0)},
        _rows(45),
        _rows(45, instrument=True),
        minimum_samples=45,
        nominal_samples=48,
    )
