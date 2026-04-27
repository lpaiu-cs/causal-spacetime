"""Observer-atlas utilities for oriented affine inertial protocols."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from causal_spacetime_lab.discrete_radar import discrete_radar_coordinates_from_order
from causal_spacetime_lab.lorentz import gamma
from causal_spacetime_lab.observer import (
    inertial_chain_inside_diamond_mask,
    observer_chain_indices,
)
from causal_spacetime_lab.oriented_radar import (
    oriented_radar_coordinates_from_two_chains,
)


@dataclass(frozen=True)
class ObserverProtocolSpec:
    """Supplied oriented affine inertial observer protocol."""

    name: str
    beta: float
    origin_lab_time: float
    origin_lab_position: float
    beacon_separation: float

    @property
    def origin_lab(self) -> tuple[float, float]:
        """Return the lab-frame origin event for this protocol."""

        return (self.origin_lab_time, self.origin_lab_position)


@dataclass(frozen=True)
class ReconstructedChart:
    """Coordinates reconstructed for one observer protocol."""

    protocol_name: str
    target_indices: NDArray[np.int_]
    accessible: NDArray[np.bool_]
    reconstructed_coords: NDArray[np.float64]


@dataclass(frozen=True)
class ProtocolChainIndices:
    """Primary and beacon chain indices for one protocol."""

    primary: NDArray[np.int_]
    beacon: NDArray[np.int_]
    clock_times: NDArray[np.float64]


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


def _validate_origin(origin_lab: tuple[float, float]) -> tuple[float, float]:
    if len(origin_lab) != 2:
        raise ValueError("origin_lab must contain two values")
    return float(origin_lab[0]), float(origin_lab[1])


def affine_rest_to_lab_1p1(
    coords_rest: ArrayLike,
    beta: float,
    origin_lab: tuple[float, float],
) -> NDArray[np.float64]:
    """Map protocol rest-frame coordinates ``(T, X)`` to lab coordinates."""

    rest = _as_coord_array(coords_rest)
    beta_value = _validate_beta(beta)
    origin_t, origin_x = _validate_origin(origin_lab)
    gamma_value = float(gamma(beta_value))
    rest_t = rest[:, 0]
    rest_x = rest[:, 1]
    lab_t = origin_t + gamma_value * (rest_t + beta_value * rest_x)
    lab_x = origin_x + gamma_value * (rest_x + beta_value * rest_t)
    return np.column_stack((lab_t, lab_x)).astype(np.float64, copy=False)


def affine_lab_to_rest_1p1(
    coords_lab: ArrayLike,
    beta: float,
    origin_lab: tuple[float, float],
) -> NDArray[np.float64]:
    """Map lab coordinates to affine protocol rest-frame coordinates."""

    lab = _as_coord_array(coords_lab)
    beta_value = _validate_beta(beta)
    origin_t, origin_x = _validate_origin(origin_lab)
    gamma_value = float(gamma(beta_value))
    shifted_t = lab[:, 0] - origin_t
    shifted_x = lab[:, 1] - origin_x
    rest_t = gamma_value * (shifted_t - beta_value * shifted_x)
    rest_x = gamma_value * (shifted_x - beta_value * shifted_t)
    return np.column_stack((rest_t, rest_x)).astype(np.float64, copy=False)


def make_affine_inertial_chain_1p1(
    beta: float,
    num_ticks: int,
    tau_min: float,
    tau_max: float,
    origin_lab: tuple[float, float] = (0.0, 0.0),
    x_prime: float = 0.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return an affine inertial chain in lab coordinates with rest clock labels."""

    _validate_beta(beta)
    ticks = int(num_ticks)
    if ticks < 2:
        raise ValueError("num_ticks must be at least 2")
    tau_min_value = float(tau_min)
    tau_max_value = float(tau_max)
    if tau_max_value <= tau_min_value:
        raise ValueError("tau_max must be greater than tau_min")

    clock_times = np.linspace(tau_min_value, tau_max_value, ticks, dtype=np.float64)
    rest_coords = np.column_stack(
        (
            clock_times,
            np.full(ticks, float(x_prime), dtype=np.float64),
        )
    )
    return affine_rest_to_lab_1p1(rest_coords, beta, origin_lab), clock_times


def safe_tau_range_for_affine_inertial_chain_1p1(
    T_global: float,
    beta: float,
    origin_lab: tuple[float, float],
    x_prime: float = 0.0,
    margin: float = 0.98,
    samples: int = 20_001,
) -> tuple[float, float]:
    """Return a conservative rest-clock range inside the global causal diamond."""

    if T_global <= 0:
        raise ValueError("T_global must be positive")
    _validate_beta(beta)
    origin_t, origin_x = _validate_origin(origin_lab)
    if not 0.0 < margin <= 1.0:
        raise ValueError("margin must satisfy 0 < margin <= 1")
    sample_count = int(samples)
    if sample_count < 3:
        raise ValueError("samples must be at least 3")

    gamma_value = float(gamma(beta))
    broad_radius = (
        T_global
        + abs(origin_t)
        + abs(origin_x)
        + gamma_value * (1.0 + abs(beta)) * (abs(float(x_prime)) + T_global)
    )
    broad_radius = max(broad_radius, T_global)
    candidates = np.linspace(-broad_radius, broad_radius, sample_count)
    rest_coords = np.column_stack(
        (
            candidates,
            np.full(sample_count, float(x_prime), dtype=np.float64),
        )
    )
    lab_events = affine_rest_to_lab_1p1(rest_coords, beta, (origin_t, origin_x))
    allowed = inertial_chain_inside_diamond_mask(lab_events, T_global)
    if not np.any(allowed):
        raise ValueError("no safe tau range exists for this chain in the diamond")

    allowed_tau = candidates[allowed]
    tau_min = float(np.min(allowed_tau))
    tau_max = float(np.max(allowed_tau))
    midpoint = 0.5 * (tau_min + tau_max)
    half_width = 0.5 * (tau_max - tau_min) * float(margin)
    if half_width <= 0.0:
        raise ValueError("safe tau range has zero width")
    return midpoint - half_width, midpoint + half_width


def common_safe_tau_range_for_oriented_protocol_1p1(
    T_global: float,
    spec: ObserverProtocolSpec,
    margin: float = 0.98,
) -> tuple[float, float]:
    """Return a common safe clock range for a protocol's primary and beacon."""

    if spec.beacon_separation <= 0.0:
        raise ValueError("beacon_separation must be positive")
    primary = safe_tau_range_for_affine_inertial_chain_1p1(
        T_global,
        spec.beta,
        spec.origin_lab,
        x_prime=0.0,
        margin=margin,
    )
    beacon = safe_tau_range_for_affine_inertial_chain_1p1(
        T_global,
        spec.beta,
        spec.origin_lab,
        x_prime=spec.beacon_separation,
        margin=margin,
    )
    tau_min = max(primary[0], beacon[0])
    tau_max = min(primary[1], beacon[1])
    if tau_max <= tau_min:
        raise ValueError("primary and beacon chains have no common safe tau range")
    return tau_min, tau_max


def reconstruct_oriented_chart_from_order(
    causal_matrix: ArrayLike,
    target_indices: ArrayLike,
    primary_indices: ArrayLike,
    beacon_indices: ArrayLike,
    clock_times: ArrayLike,
    beacon_separation: float,
    protocol_name: str,
) -> ReconstructedChart:
    """Reconstruct one oriented chart from order data and supplied clocks."""

    targets = np.asarray(target_indices, dtype=int)
    if targets.ndim != 1:
        raise ValueError("target_indices must be one-dimensional")

    primary = discrete_radar_coordinates_from_order(
        causal_matrix,
        primary_indices,
        targets,
        clock_times,
    )
    beacon = discrete_radar_coordinates_from_order(
        causal_matrix,
        beacon_indices,
        targets,
        clock_times,
    )
    oriented = oriented_radar_coordinates_from_two_chains(
        primary,
        beacon,
        beacon_separation,
    )

    coords = np.full((targets.size, 2), np.nan, dtype=np.float64)
    accessible = np.zeros(targets.size, dtype=bool)
    for index, result in enumerate(oriented):
        if result.accessible:
            accessible[index] = True
            coords[index, 0] = float(result.oriented_time)
            coords[index, 1] = float(result.signed_position)

    return ReconstructedChart(
        protocol_name=str(protocol_name),
        target_indices=targets.astype(int, copy=False),
        accessible=accessible,
        reconstructed_coords=coords,
    )


def append_oriented_protocol_chains(
    chain_events: list[NDArray[np.float64]],
    spec: ObserverProtocolSpec,
    tick_count: int,
    tau_range: tuple[float, float],
    start_index: int,
) -> tuple[ProtocolChainIndices, int]:
    """Append primary and beacon chain events and return their combined indices."""

    primary_events, clocks = make_affine_inertial_chain_1p1(
        spec.beta,
        tick_count,
        tau_range[0],
        tau_range[1],
        origin_lab=spec.origin_lab,
        x_prime=0.0,
    )
    beacon_events, beacon_clocks = make_affine_inertial_chain_1p1(
        spec.beta,
        tick_count,
        tau_range[0],
        tau_range[1],
        origin_lab=spec.origin_lab,
        x_prime=spec.beacon_separation,
    )
    if not np.allclose(clocks, beacon_clocks):
        raise RuntimeError("primary and beacon clocks are not synchronized")

    primary_indices = observer_chain_indices(start_index, tick_count)
    beacon_indices = observer_chain_indices(start_index + tick_count, tick_count)
    chain_events.extend([primary_events, beacon_events])
    return (
        ProtocolChainIndices(
            primary=primary_indices,
            beacon=beacon_indices,
            clock_times=clocks,
        ),
        start_index + 2 * tick_count,
    )
