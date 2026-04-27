from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.calibration import (
    affine_fit_1d,
    equal_step_deviation,
    ratio_error_under_transform,
)


def test_equal_step_deviation_arithmetic_sequence_is_zero() -> None:
    values = np.asarray([0.0, 2.0, 4.0, 6.0])

    assert equal_step_deviation(values) == pytest.approx(0.0)


def test_ratio_error_under_transform_zero_for_affine_scaling() -> None:
    values = np.asarray([0.0, 1.0, 3.0, 7.0])
    transformed = 2.5 * values + 4.0
    intervals = np.asarray([[0, 2, 0, 3], [1, 2, 0, 2]], dtype=int)

    error = ratio_error_under_transform(values, transformed, intervals)

    assert error == pytest.approx(0.0)


def test_ratio_error_under_transform_nonzero_for_square_transform() -> None:
    values = np.asarray([0.0, 1.0, 3.0, 7.0])
    transformed = values * values
    intervals = np.asarray([[0, 2, 0, 3], [1, 2, 0, 2]], dtype=int)

    error = ratio_error_under_transform(values, transformed, intervals)

    assert error > 0.0


def test_affine_fit_1d_recovers_scale_and_offset() -> None:
    source = np.asarray([0.0, 1.0, 2.0, 3.0])
    target = 3.0 * source - 2.0

    scale, offset, rmse = affine_fit_1d(source, target)

    assert scale == pytest.approx(3.0)
    assert offset == pytest.approx(-2.0)
    assert rmse == pytest.approx(0.0)
