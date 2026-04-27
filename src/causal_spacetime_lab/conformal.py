"""Conformal-factor utilities for controlled 1+1D ambiguity experiments."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class ConformalProfile:
    """Positive 1+1D conformal profile used as supplied measure structure."""

    name: str
    parameters: dict[str, float]

    def omega(self, events: ArrayLike) -> NDArray[np.float64]:
        """Return conformal factor values at events ``(t, x)``."""

        event_array = _as_event_array(events)
        t = event_array[:, 0]
        if self.name == "flat":
            return np.ones(event_array.shape[0], dtype=np.float64)
        if self.name == "constant":
            return np.full(
                event_array.shape[0],
                self.parameters["scale"],
                dtype=np.float64,
            )
        if self.name == "linear_time":
            amplitude = self.parameters["amplitude"]
            duration = self.parameters["T"]
            return (1.0 + amplitude * t / (duration / 2.0)).astype(
                np.float64,
                copy=False,
            )
        if self.name == "sinusoidal_time":
            amplitude = self.parameters["amplitude"]
            duration = self.parameters["T"]
            return (1.0 + amplitude * np.sin(np.pi * t / duration)).astype(
                np.float64,
                copy=False,
            )
        raise ValueError(f"unknown conformal profile: {self.name}")


def _as_event_array(events: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(events, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("events must have shape (n, 2)")
    return array


def _as_event(event: ArrayLike, name: str) -> NDArray[np.float64]:
    array = np.asarray(event, dtype=float)
    if array.shape != (2,):
        raise ValueError(f"{name} must have shape (2,)")
    return array


def flat_profile() -> ConformalProfile:
    """Return ``Omega = 1``."""

    return ConformalProfile(name="flat", parameters={})


def constant_profile(scale: float) -> ConformalProfile:
    """Return constant ``Omega = scale``."""

    value = float(scale)
    if value <= 0.0:
        raise ValueError("scale must be positive")
    return ConformalProfile(name="constant", parameters={"scale": value})


def linear_time_profile(amplitude: float, T: float) -> ConformalProfile:
    """Return ``Omega(t) = 1 + amplitude * t / (T / 2)``."""

    value = float(amplitude)
    duration = float(T)
    if duration <= 0.0:
        raise ValueError("T must be positive")
    if abs(value) >= 1.0:
        raise ValueError("abs(amplitude) must be less than 1")
    return ConformalProfile(
        name="linear_time",
        parameters={"amplitude": value, "T": duration},
    )


def sinusoidal_time_profile(amplitude: float, T: float) -> ConformalProfile:
    """Return ``Omega(t) = 1 + amplitude * sin(pi * t / T)``."""

    value = float(amplitude)
    duration = float(T)
    if duration <= 0.0:
        raise ValueError("T must be positive")
    if abs(value) >= 1.0:
        raise ValueError("abs(amplitude) must be less than 1")
    return ConformalProfile(
        name="sinusoidal_time",
        parameters={"amplitude": value, "T": duration},
    )


def omega_values(
    profile: ConformalProfile,
    events: ArrayLike,
) -> NDArray[np.float64]:
    """Return conformal factor values for ``events``."""

    values = profile.omega(events)
    if np.any(values <= 0.0):
        raise ValueError("conformal profile produced nonpositive Omega values")
    return values.astype(np.float64, copy=False)


def conformal_volume_density_1p1(
    profile: ConformalProfile,
    events: ArrayLike,
) -> NDArray[np.float64]:
    """Return 1+1D volume density weights ``Omega^2``."""

    omega = omega_values(profile, events)
    return (omega * omega).astype(np.float64, copy=False)


def conformal_clock_rate_1p1(
    profile: ConformalProfile,
    events: ArrayLike,
) -> NDArray[np.float64]:
    """Return clock rate ``Omega`` along a worldline."""

    return omega_values(profile, events)


def _validate_flat_timelike_pair(
    p: NDArray[np.float64],
    q: NDArray[np.float64],
) -> None:
    dt = q[0] - p[0]
    dx = q[1] - p[1]
    if dt <= 0.0 or dt * dt < dx * dx:
        raise ValueError("p must causally precede q in flat 1+1D Minkowski order")


def flat_interval_x_width_1p1(
    t_values: ArrayLike,
    p: ArrayLike,
    q: ArrayLike,
) -> NDArray[np.float64]:
    """Return flat Alexandrov interval x-width at each sampled time."""

    p_event = _as_event(p, "p")
    q_event = _as_event(q, "q")
    _validate_flat_timelike_pair(p_event, q_event)
    t = np.asarray(t_values, dtype=float)
    if t.ndim != 1:
        raise ValueError("t_values must be one-dimensional")

    lower = np.maximum(
        p_event[1] - (t - p_event[0]),
        q_event[1] - (q_event[0] - t),
    )
    upper = np.minimum(
        p_event[1] + (t - p_event[0]),
        q_event[1] + (q_event[0] - t),
    )
    return np.maximum(0.0, upper - lower).astype(np.float64, copy=False)


def integrate_conformal_interval_volume_1p1(
    p: ArrayLike,
    q: ArrayLike,
    profile: ConformalProfile,
    num_t: int = 4096,
) -> float:
    """Numerically integrate conformal physical volume over ``I(p, q)``."""

    p_event = _as_event(p, "p")
    q_event = _as_event(q, "q")
    _validate_flat_timelike_pair(p_event, q_event)
    sample_count = int(num_t)
    if sample_count < 2:
        raise ValueError("num_t must be at least 2")

    t_values = np.linspace(p_event[0], q_event[0], sample_count)
    widths = flat_interval_x_width_1p1(t_values, p_event, q_event)
    events = np.column_stack((t_values, np.zeros(sample_count, dtype=float)))
    weights = conformal_volume_density_1p1(profile, events)
    return float(np.trapezoid(widths * weights, t_values))


def central_observer_proper_time_1p1(
    t_min: float,
    t_max: float,
    profile: ConformalProfile,
    x: float = 0.0,
    num_t: int = 4096,
) -> float:
    """Integrate ``tau = integral Omega(t, x) dt`` for a fixed-x worldline."""

    lower = float(t_min)
    upper = float(t_max)
    if upper <= lower:
        raise ValueError("t_max must be greater than t_min")
    sample_count = int(num_t)
    if sample_count < 2:
        raise ValueError("num_t must be at least 2")
    t_values = np.linspace(lower, upper, sample_count)
    events = np.column_stack((t_values, np.full(sample_count, float(x))))
    rates = conformal_clock_rate_1p1(profile, events)
    return float(np.trapezoid(rates, t_values))


def _as_causal_matrix(causal_matrix: ArrayLike) -> NDArray[np.bool_]:
    matrix = np.asarray(causal_matrix, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("causal_matrix must be a square matrix")
    return matrix


def _validate_interval_indices(matrix: NDArray[np.bool_], i: int, j: int) -> None:
    n_events = matrix.shape[0]
    if not 0 <= i < n_events or not 0 <= j < n_events:
        raise IndexError("interval endpoint index out of range")
    if not matrix[i, j]:
        raise ValueError("i must causally precede j")


def weighted_interval_volume_estimate_1p1(
    events: ArrayLike,
    causal_matrix: ArrayLike,
    i: int,
    j: int,
    profile: ConformalProfile,
    coordinate_density: float,
) -> float:
    """Estimate conformal volume from supplied ``Omega^2`` measure weights."""

    event_array = _as_event_array(events)
    matrix = _as_causal_matrix(causal_matrix)
    if event_array.shape[0] != matrix.shape[0]:
        raise ValueError("events and causal_matrix must contain the same count")
    if coordinate_density <= 0.0:
        raise ValueError("coordinate_density must be positive")
    start = int(i)
    end = int(j)
    _validate_interval_indices(matrix, start, end)

    inside = matrix[start] & matrix[:, end]
    weights = conformal_volume_density_1p1(profile, event_array[inside])
    return float(np.sum(weights) / coordinate_density)


def unweighted_interval_coordinate_volume_estimate_1p1(
    causal_matrix: ArrayLike,
    i: int,
    j: int,
    coordinate_density: float,
) -> float:
    """Estimate flat coordinate volume as interval count divided by density."""

    matrix = _as_causal_matrix(causal_matrix)
    if coordinate_density <= 0.0:
        raise ValueError("coordinate_density must be positive")
    start = int(i)
    end = int(j)
    _validate_interval_indices(matrix, start, end)
    return float(np.count_nonzero(matrix[start] & matrix[:, end]) / coordinate_density)
