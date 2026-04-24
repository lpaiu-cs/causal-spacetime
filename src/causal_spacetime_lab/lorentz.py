"""Lorentz transformations and length-contraction helpers in 1+1D."""

from __future__ import annotations

from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatOrArray: TypeAlias = float | NDArray[np.float64]


def _as_float_or_array(value: ArrayLike) -> FloatOrArray:
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        return float(array)
    return array.astype(np.float64, copy=False)


def _validate_beta(beta: ArrayLike) -> FloatOrArray:
    beta_value = _as_float_or_array(beta)
    if np.any(np.abs(beta_value) >= 1.0):
        raise ValueError("beta must satisfy abs(beta) < 1 in units c = 1")
    return beta_value


def _validate_rest_length(L0: float) -> float:
    rest_length = float(L0)
    if rest_length < 0.0:
        raise ValueError("L0 must be non-negative")
    return rest_length


def gamma(beta: ArrayLike) -> FloatOrArray:
    """Return the Lorentz factor for velocity ``beta = v/c`` with ``c = 1``."""

    beta_value = _validate_beta(beta)
    return _as_float_or_array(1.0 / np.sqrt(1.0 - np.asarray(beta_value) ** 2))


def lorentz_transform_1p1(
    t_prime: ArrayLike,
    x_prime: ArrayLike,
    beta: ArrayLike,
) -> tuple[FloatOrArray, FloatOrArray]:
    """Transform coordinates from rod rest frame ``S'`` to lab frame ``S``.

    The ``S'`` frame moves at velocity ``beta`` along the lab-frame x-axis:

    ``t = gamma * (t' + beta x')``
    ``x = gamma * (x' + beta t')``
    """

    beta_value = _validate_beta(beta)
    gamma_value = gamma(beta_value)
    t = np.asarray(gamma_value) * (
        np.asarray(t_prime, dtype=float) + np.asarray(beta_value) * x_prime
    )
    x = np.asarray(gamma_value) * (
        np.asarray(x_prime, dtype=float) + np.asarray(beta_value) * t_prime
    )
    return _as_float_or_array(t), _as_float_or_array(x)


def inverse_lorentz_transform_1p1(
    t: ArrayLike,
    x: ArrayLike,
    beta: ArrayLike,
) -> tuple[FloatOrArray, FloatOrArray]:
    """Transform coordinates from lab frame ``S`` to rod rest frame ``S'``."""

    beta_value = _validate_beta(beta)
    gamma_value = gamma(beta_value)
    t_prime = np.asarray(gamma_value) * (
        np.asarray(t, dtype=float) - np.asarray(beta_value) * x
    )
    x_prime = np.asarray(gamma_value) * (
        np.asarray(x, dtype=float) - np.asarray(beta_value) * t
    )
    return _as_float_or_array(t_prime), _as_float_or_array(x_prime)


def rod_endpoint_events_simultaneous_in_lab(
    L0: float,
    beta: float,
    lab_time: float = 0.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return rod endpoint events selected at the same lab-frame time.

    The rod is at rest in ``S'`` with endpoints ``x'=0`` and ``x'=L0``. The
    returned events are lab-frame ``(t, x)`` coordinates for the two endpoint
    events whose lab time is ``lab_time``.
    """

    rest_length = _validate_rest_length(L0)
    beta_value = float(_validate_beta(beta))
    gamma_value = float(gamma(beta_value))
    lab_time_value = float(lab_time)

    left_t_prime = lab_time_value / gamma_value
    right_t_prime = lab_time_value / gamma_value - beta_value * rest_length

    left_event = np.asarray(
        lorentz_transform_1p1(left_t_prime, 0.0, beta_value),
        dtype=float,
    )
    right_event = np.asarray(
        lorentz_transform_1p1(right_t_prime, rest_length, beta_value),
        dtype=float,
    )
    return left_event, right_event


def measured_rod_length_lab(L0: float, beta: float) -> float:
    """Measure a moving rod by selecting endpoint events simultaneous in the lab."""

    left_event, right_event = rod_endpoint_events_simultaneous_in_lab(
        L0,
        beta,
        lab_time=0.0,
    )
    return float(abs(right_event[1] - left_event[1]))

