"""Utilities for validating Lorentz maps between reconstructed coordinates."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.lorentz import gamma


def _as_coord_array(coords: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(coords, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("coords must have shape (n, 2)")
    return array


def _validate_beta(beta: float) -> float:
    value = float(beta)
    if abs(value) >= 1.0:
        raise ValueError("beta must satisfy abs(beta) < 1 in units c = 1")
    return value


def true_lorentz_rest_coordinates_1p1(
    events_lab: ArrayLike,
    beta: float,
) -> NDArray[np.float64]:
    """Return moving-rest-frame coordinates ``(t', x')`` for lab events."""

    return lorentz_transform_coords_1p1(events_lab, beta)


def lorentz_transform_coords_1p1(
    coords: ArrayLike,
    beta: float,
) -> NDArray[np.float64]:
    """Transform lab coordinates ``(t, x)`` to moving rest-frame coordinates."""

    coord_array = _as_coord_array(coords)
    beta_value = _validate_beta(beta)
    gamma_value = float(gamma(beta_value))
    t = coord_array[:, 0]
    x = coord_array[:, 1]
    t_prime = gamma_value * (t - beta_value * x)
    x_prime = gamma_value * (x - beta_value * t)
    return np.column_stack((t_prime, x_prime)).astype(np.float64, copy=False)


def lorentz_residual_rmse(
    source_coords: ArrayLike,
    target_coords: ArrayLike,
    beta: float,
) -> float:
    """Return coordinate RMSE after applying the Lorentz map with ``beta``."""

    source = _as_coord_array(source_coords)
    target = _as_coord_array(target_coords)
    if source.shape != target.shape:
        raise ValueError("source_coords and target_coords must have the same shape")
    predicted = lorentz_transform_coords_1p1(source, beta)
    residual = predicted - target
    return float(np.sqrt(np.mean(residual * residual)))


def fit_lorentz_beta_grid(
    source_coords: ArrayLike,
    target_coords: ArrayLike,
    beta_min: float = -0.95,
    beta_max: float = 0.95,
    num_grid: int = 2001,
) -> tuple[float, float]:
    """Fit a Lorentz beta by grid search over residual RMSE."""

    source = _as_coord_array(source_coords)
    target = _as_coord_array(target_coords)
    if source.shape != target.shape:
        raise ValueError("source_coords and target_coords must have the same shape")
    if beta_max <= beta_min:
        raise ValueError("beta_max must be greater than beta_min")
    if abs(beta_min) >= 1.0 or abs(beta_max) >= 1.0:
        raise ValueError("beta bounds must lie inside (-1, 1)")
    grid_count = int(num_grid)
    if grid_count < 2:
        raise ValueError("num_grid must be at least 2")

    betas = np.linspace(beta_min, beta_max, grid_count)
    gammas = 1.0 / np.sqrt(1.0 - betas * betas)
    source_t = source[:, 0][None, :]
    source_x = source[:, 1][None, :]
    target_t = target[:, 0][None, :]
    target_x = target[:, 1][None, :]
    beta_column = betas[:, None]
    gamma_column = gammas[:, None]

    predicted_t = gamma_column * (source_t - beta_column * source_x)
    predicted_x = gamma_column * (source_x - beta_column * source_t)
    squared_residual = (predicted_t - target_t) ** 2 + (predicted_x - target_x) ** 2
    rmse_values = np.sqrt(np.mean(squared_residual, axis=1) / 2.0)
    best_index = int(np.argmin(rmse_values))
    return float(betas[best_index]), float(rmse_values[best_index])
