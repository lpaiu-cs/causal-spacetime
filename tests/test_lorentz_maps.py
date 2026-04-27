from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.lorentz import gamma
from causal_spacetime_lab.lorentz_maps import (
    fit_lorentz_beta_grid,
    lorentz_residual_rmse,
    lorentz_transform_coords_1p1,
    true_lorentz_rest_coordinates_1p1,
)


def test_true_lorentz_rest_coordinates_beta_zero_is_identity() -> None:
    events = np.array([[0.0, 0.1], [0.4, -0.2], [-0.3, 0.0]])

    transformed = true_lorentz_rest_coordinates_1p1(events, beta=0.0)

    assert transformed == pytest.approx(events)


def test_true_lorentz_rest_coordinates_nonzero_beta_hand_coded() -> None:
    events = np.array([[1.0, 0.5]])
    beta = 0.6
    gamma_value = gamma(beta)

    transformed = true_lorentz_rest_coordinates_1p1(events, beta=beta)

    assert transformed[0, 0] == pytest.approx(gamma_value * (1.0 - beta * 0.5))
    assert transformed[0, 1] == pytest.approx(gamma_value * (0.5 - beta * 1.0))


def test_lorentz_residual_rmse_is_zero_for_exact_transformed_coordinates() -> None:
    source = np.array([[0.0, 0.2], [0.5, -0.1], [-0.4, 0.05]])
    target = lorentz_transform_coords_1p1(source, beta=0.3)

    residual = lorentz_residual_rmse(source, target, beta=0.3)

    assert residual == pytest.approx(0.0, abs=1e-14)


def test_fit_lorentz_beta_grid_recovers_known_beta_from_exact_coordinates() -> None:
    source = np.array(
        [
            [0.0, 0.2],
            [0.5, -0.1],
            [-0.4, 0.05],
            [0.2, 0.35],
        ]
    )
    target = lorentz_transform_coords_1p1(source, beta=0.4)

    best_beta, best_rmse = fit_lorentz_beta_grid(
        source,
        target,
        beta_min=-0.8,
        beta_max=0.8,
        num_grid=1601,
    )

    assert best_beta == pytest.approx(0.4, abs=1e-12)
    assert best_rmse == pytest.approx(0.0, abs=1e-14)
