from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.horizon_validation import (
    finite_coverage_mask_from_analytic_ticks,
    horizon_accessibility_confusion,
    radar_error_summary,
)


def test_horizon_accessibility_confusion_deterministic_masks() -> None:
    reconstructed = np.array([True, True, False, False])
    reference = np.array([True, False, True, False])

    summary = horizon_accessibility_confusion(reconstructed, reference)

    assert summary["true_positive"] == 1.0
    assert summary["false_positive"] == 1.0
    assert summary["true_negative"] == 1.0
    assert summary["false_negative"] == 1.0
    assert summary["precision"] == pytest.approx(0.5)
    assert summary["recall"] == pytest.approx(0.5)


def test_finite_coverage_mask_from_analytic_ticks_deterministic_arrays() -> None:
    tau_minus = np.array([-0.5, -2.0, np.nan, 0.1])
    tau_plus = np.array([0.5, 0.2, 0.3, 2.0])

    mask = finite_coverage_mask_from_analytic_ticks(
        tau_minus,
        tau_plus,
        tau_min=-1.0,
        tau_max=1.0,
    )

    assert mask.tolist() == [True, False, False, False]


def test_radar_error_summary_simple_arrays() -> None:
    summary = radar_error_summary(
        reconstructed_time=np.array([1.0, 2.0, np.nan]),
        reconstructed_distance=np.array([0.5, 1.5, np.nan]),
        analytic_time=np.array([1.1, 1.9, 0.0]),
        analytic_distance=np.array([0.4, 1.7, 0.0]),
        mask=np.array([True, True, False]),
    )

    assert summary["radar_time_mae"] == pytest.approx(0.1)
    assert summary["radar_distance_mae"] == pytest.approx(0.15)
    assert summary["radar_time_bias"] == pytest.approx(0.0)
    assert summary["radar_distance_bias"] == pytest.approx(-0.05)
