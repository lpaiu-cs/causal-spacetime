"""Affine Lorentz/Poincare map fitting utilities for 1+1D charts."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.lorentz_maps import lorentz_transform_coords_1p1


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


def _as_translation(translation: ArrayLike) -> NDArray[np.float64]:
    vector = np.asarray(translation, dtype=float)
    if vector.shape != (2,):
        raise ValueError("translation must have shape (2,)")
    return vector.astype(np.float64, copy=False)


def apply_affine_lorentz_map_1p1(
    coords: ArrayLike,
    beta: float,
    translation: tuple[float, float] | ArrayLike,
) -> NDArray[np.float64]:
    """Apply ``target = L(beta) source + translation``."""

    transformed = lorentz_transform_coords_1p1(coords, _validate_beta(beta))
    return transformed + _as_translation(translation)


def fit_affine_lorentz_beta_grid_1p1(
    source_coords: ArrayLike,
    target_coords: ArrayLike,
    beta_min: float = -0.95,
    beta_max: float = 0.95,
    num_grid: int = 2001,
) -> tuple[float, NDArray[np.float64], float]:
    """Fit an affine Lorentz map by grid search over beta."""

    source = _as_coord_array(source_coords)
    target = _as_coord_array(target_coords)
    if source.shape != target.shape:
        raise ValueError("source_coords and target_coords must have the same shape")
    if source.shape[0] == 0:
        raise ValueError("at least one coordinate pair is required")
    if beta_max <= beta_min:
        raise ValueError("beta_max must be greater than beta_min")
    if abs(beta_min) >= 1.0 or abs(beta_max) >= 1.0:
        raise ValueError("beta bounds must lie inside (-1, 1)")
    grid_count = int(num_grid)
    if grid_count < 2:
        raise ValueError("num_grid must be at least 2")

    betas = np.linspace(beta_min, beta_max, grid_count)
    gammas = 1.0 / np.sqrt(1.0 - betas * betas)
    beta_column = betas[:, None]
    gamma_column = gammas[:, None]
    source_t = source[:, 0][None, :]
    source_x = source[:, 1][None, :]
    target_t = target[:, 0][None, :]
    target_x = target[:, 1][None, :]

    transformed_t = gamma_column * (source_t - beta_column * source_x)
    transformed_x = gamma_column * (source_x - beta_column * source_t)
    translation_t = np.mean(target_t - transformed_t, axis=1)
    translation_x = np.mean(target_x - transformed_x, axis=1)
    residual_t = transformed_t + translation_t[:, None] - target_t
    residual_x = transformed_x + translation_x[:, None] - target_x
    rmse_values = np.sqrt(np.mean(residual_t**2 + residual_x**2, axis=1) / 2.0)
    best_index = int(np.argmin(rmse_values))
    translation = np.asarray(
        [translation_t[best_index], translation_x[best_index]],
        dtype=np.float64,
    )
    return float(betas[best_index]), translation, float(rmse_values[best_index])


def expected_relative_beta_1p1(source_beta: float, target_beta: float) -> float:
    """Return beta for ``target_coords = L(target) L(-source) source_coords``."""

    source = _validate_beta(source_beta)
    target = _validate_beta(target_beta)
    denominator = 1.0 - source * target
    if denominator == 0.0:
        raise ValueError("relative beta denominator is zero")
    return (target - source) / denominator


def compose_betas_1p1(beta_ab: float, beta_bc: float) -> float:
    """Compose boosts for maps ``A -> B -> C``."""

    first = _validate_beta(beta_ab)
    second = _validate_beta(beta_bc)
    denominator = 1.0 + first * second
    if denominator == 0.0:
        raise ValueError("composed beta denominator is zero")
    composed = (first + second) / denominator
    return _validate_beta(composed)


def compose_affine_lorentz_maps_1p1(
    beta_ab: float,
    translation_ab: ArrayLike,
    beta_bc: float,
    translation_bc: ArrayLike,
) -> tuple[float, NDArray[np.float64]]:
    """Compose affine maps ``A -> B -> C`` into a direct ``A -> C`` map."""

    composed_beta = compose_betas_1p1(beta_ab, beta_bc)
    translated_ab = apply_affine_lorentz_map_1p1(
        _as_translation(translation_ab)[None, :],
        beta_bc,
        (0.0, 0.0),
    )[0]
    composed_translation = translated_ab + _as_translation(translation_bc)
    return composed_beta, composed_translation.astype(np.float64, copy=False)
