"""Validation helpers for horizon-limited reconstruction experiments."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from causal_spacetime_lab.estimators import bias, mae, rmse


def horizon_accessibility_confusion(
    reconstructed_accessible: ArrayLike,
    analytic_wedge_accessible: ArrayLike,
) -> dict[str, float]:
    """Return confusion metrics for reconstructed two-way accessibility."""

    reconstructed = np.asarray(reconstructed_accessible, dtype=bool)
    analytic = np.asarray(analytic_wedge_accessible, dtype=bool)
    if reconstructed.shape != analytic.shape:
        raise ValueError("accessibility masks must have the same shape")

    true_positive = float(np.count_nonzero(reconstructed & analytic))
    false_positive = float(np.count_nonzero(reconstructed & ~analytic))
    true_negative = float(np.count_nonzero(~reconstructed & ~analytic))
    false_negative = float(np.count_nonzero(~reconstructed & analytic))
    precision_denominator = true_positive + false_positive
    recall_denominator = true_positive + false_negative
    false_positive_denominator = false_positive + true_negative
    false_negative_denominator = false_negative + true_positive
    return {
        "true_positive": true_positive,
        "false_positive": false_positive,
        "true_negative": true_negative,
        "false_negative": false_negative,
        "precision": true_positive / precision_denominator
        if precision_denominator > 0.0
        else float("nan"),
        "recall": true_positive / recall_denominator
        if recall_denominator > 0.0
        else float("nan"),
        "false_positive_rate": false_positive / false_positive_denominator
        if false_positive_denominator > 0.0
        else float("nan"),
        "false_negative_rate": false_negative / false_negative_denominator
        if false_negative_denominator > 0.0
        else float("nan"),
    }


def finite_coverage_mask_from_analytic_ticks(
    tau_minus: ArrayLike,
    tau_plus: ArrayLike,
    tau_min: float,
    tau_max: float,
) -> np.ndarray:
    """Return mask for events covered by a finite observer-chain time range."""

    lower_ticks = np.asarray(tau_minus, dtype=float)
    upper_ticks = np.asarray(tau_plus, dtype=float)
    if lower_ticks.shape != upper_ticks.shape:
        raise ValueError("tau_minus and tau_plus must have the same shape")
    min_value = float(tau_min)
    max_value = float(tau_max)
    if max_value < min_value:
        raise ValueError("tau_max must be greater than or equal to tau_min")
    return (
        np.isfinite(lower_ticks)
        & np.isfinite(upper_ticks)
        & (lower_ticks >= min_value)
        & (upper_ticks <= max_value)
    )


def radar_error_summary(
    reconstructed_time: ArrayLike,
    reconstructed_distance: ArrayLike,
    analytic_time: ArrayLike,
    analytic_distance: ArrayLike,
    mask: ArrayLike,
) -> dict[str, float]:
    """Return radar-coordinate error metrics on a selected event mask."""

    time_hat = np.asarray(reconstructed_time, dtype=float)
    distance_hat = np.asarray(reconstructed_distance, dtype=float)
    time_true = np.asarray(analytic_time, dtype=float)
    distance_true = np.asarray(analytic_distance, dtype=float)
    selected = np.asarray(mask, dtype=bool)
    if not (
        time_hat.shape
        == distance_hat.shape
        == time_true.shape
        == distance_true.shape
        == selected.shape
    ):
        raise ValueError("all radar arrays and mask must have the same shape")

    selected &= (
        np.isfinite(time_hat)
        & np.isfinite(distance_hat)
        & np.isfinite(time_true)
        & np.isfinite(distance_true)
    )
    if not np.any(selected):
        return {
            "radar_time_mae": float("nan"),
            "radar_time_rmse": float("nan"),
            "radar_time_bias": float("nan"),
            "radar_distance_mae": float("nan"),
            "radar_distance_rmse": float("nan"),
            "radar_distance_bias": float("nan"),
        }

    return {
        "radar_time_mae": mae(time_true[selected], time_hat[selected]),
        "radar_time_rmse": rmse(time_true[selected], time_hat[selected]),
        "radar_time_bias": bias(time_true[selected], time_hat[selected]),
        "radar_distance_mae": mae(
            distance_true[selected],
            distance_hat[selected],
        ),
        "radar_distance_rmse": rmse(
            distance_true[selected],
            distance_hat[selected],
        ),
        "radar_distance_bias": bias(
            distance_true[selected],
            distance_hat[selected],
        ),
    }
