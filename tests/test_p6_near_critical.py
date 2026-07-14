"""Statistics tests for the P6 near-critical characterization."""

from __future__ import annotations

from causal_spacetime_lab.positive_control.mcmc_diagnostics import (
    integrated_autocorrelation,
)


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
