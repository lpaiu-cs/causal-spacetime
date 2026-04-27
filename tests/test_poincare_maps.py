from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.observer_atlas import (
    affine_lab_to_rest_1p1,
    affine_rest_to_lab_1p1,
)
from causal_spacetime_lab.poincare_maps import (
    apply_affine_lorentz_map_1p1,
    compose_affine_lorentz_maps_1p1,
    compose_betas_1p1,
    expected_relative_beta_1p1,
    fit_affine_lorentz_beta_grid_1p1,
)


def test_apply_affine_lorentz_map_hand_coded_beta_zero() -> None:
    coords = np.array([[1.0, 0.5], [-0.5, 0.25]])

    transformed = apply_affine_lorentz_map_1p1(coords, 0.0, (0.2, -0.1))

    assert transformed == pytest.approx(coords + np.array([0.2, -0.1]))


def test_fit_affine_lorentz_beta_grid_recovers_exact_map() -> None:
    source = np.array([[-0.3, 0.2], [0.0, -0.1], [0.4, 0.25], [0.2, -0.3]])
    target = apply_affine_lorentz_map_1p1(source, 0.4, (0.2, -0.15))

    beta, translation, residual = fit_affine_lorentz_beta_grid_1p1(
        source,
        target,
        beta_min=-0.8,
        beta_max=0.8,
        num_grid=1601,
    )

    assert beta == pytest.approx(0.4)
    assert translation == pytest.approx(np.array([0.2, -0.15]))
    assert residual == pytest.approx(0.0, abs=1e-14)


def test_expected_relative_beta_agrees_with_exact_coordinate_transforms() -> None:
    source_coords = np.array([[-0.2, 0.1], [0.0, 0.0], [0.3, -0.15]])
    source_beta = 0.2
    target_beta = 0.5
    lab = affine_rest_to_lab_1p1(source_coords, source_beta, (0.0, 0.0))
    target_coords = affine_lab_to_rest_1p1(lab, target_beta, (0.0, 0.0))

    beta_rel = expected_relative_beta_1p1(source_beta, target_beta)
    predicted = apply_affine_lorentz_map_1p1(source_coords, beta_rel, (0.0, 0.0))

    assert predicted == pytest.approx(target_coords)


def test_compose_betas_agrees_with_successive_lorentz_transforms() -> None:
    coords = np.array([[-0.3, 0.1], [0.2, -0.05]])
    beta_ab = 0.2
    beta_bc = -0.4

    sequential = apply_affine_lorentz_map_1p1(
        apply_affine_lorentz_map_1p1(coords, beta_ab, (0.0, 0.0)),
        beta_bc,
        (0.0, 0.0),
    )
    composed = apply_affine_lorentz_map_1p1(
        coords,
        compose_betas_1p1(beta_ab, beta_bc),
        (0.0, 0.0),
    )

    assert composed == pytest.approx(sequential)


def test_compose_affine_lorentz_maps_agrees_with_sequential_application() -> None:
    coords = np.array([[-0.3, 0.1], [0.2, -0.05], [0.4, 0.2]])
    beta_ab = 0.2
    translation_ab = np.array([0.1, -0.2])
    beta_bc = -0.4
    translation_bc = np.array([-0.05, 0.08])

    sequential = apply_affine_lorentz_map_1p1(
        apply_affine_lorentz_map_1p1(coords, beta_ab, translation_ab),
        beta_bc,
        translation_bc,
    )
    beta_ac, translation_ac = compose_affine_lorentz_maps_1p1(
        beta_ab,
        translation_ab,
        beta_bc,
        translation_bc,
    )
    composed = apply_affine_lorentz_map_1p1(coords, beta_ac, translation_ac)

    assert composed == pytest.approx(sequential)
