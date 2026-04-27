"""Observer-slice reconstruction pipeline helpers."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.observer_atlas import (
    ObserverProtocolSpec,
    append_oriented_protocol_chains,
    common_safe_tau_range_for_oriented_protocol_1p1,
    reconstruct_oriented_chart_from_order,
)
from causal_spacetime_lab.radar_order import radar_tick_brackets_from_order
from causal_spacetime_lab.spatial_slices import (
    assign_slices_from_radar_time_rank,
    radar_distance_rank_from_tick_brackets,
    radar_time_rank_from_tick_brackets,
)


def reconstruct_stationary_oriented_slices_1p1(
    support_events: NDArray[np.float64],
    T: float,
    tick_count: int,
    beacon_separation: float,
    bin_width: int = 2,
) -> dict[str, NDArray[np.float64] | NDArray[np.int_] | NDArray[np.bool_]]:
    """Reconstruct oriented coordinates and radar-time slices from order data."""

    events = np.asarray(support_events, dtype=float)
    if events.ndim != 2 or events.shape[1] != 2:
        raise ValueError("support_events must have shape (n, 2)")
    if T <= 0:
        raise ValueError("T must be positive")
    ticks = int(tick_count)
    if ticks < 2:
        raise ValueError("tick_count must be at least 2")
    if beacon_separation <= 0.0:
        raise ValueError("beacon_separation must be positive")

    spec = ObserverProtocolSpec(
        name="stationary_oriented",
        beta=0.0,
        origin_lab_time=0.0,
        origin_lab_position=0.0,
        beacon_separation=float(beacon_separation),
    )
    tau_range = common_safe_tau_range_for_oriented_protocol_1p1(T, spec)
    chain_events: list[NDArray[np.float64]] = []
    indices, _ = append_oriented_protocol_chains(
        chain_events,
        spec,
        ticks,
        tau_range,
        events.shape[0],
    )
    combined = np.vstack((events, *chain_events))
    causal_order = causal_matrix_1p1(combined)
    targets = np.arange(events.shape[0], dtype=int)
    chart = reconstruct_oriented_chart_from_order(
        causal_order,
        targets,
        indices.primary,
        indices.beacon,
        indices.clock_times,
        float(beacon_separation),
        spec.name,
    )
    predecessor, successor, bracket_accessible = radar_tick_brackets_from_order(
        causal_order,
        indices.primary,
        targets,
    )
    accessible = chart.accessible & bracket_accessible
    radar_time_rank = radar_time_rank_from_tick_brackets(
        predecessor,
        successor,
        accessible,
    )
    radar_distance_rank = radar_distance_rank_from_tick_brackets(
        predecessor,
        successor,
        accessible,
    )
    slice_labels = assign_slices_from_radar_time_rank(
        radar_time_rank,
        accessible,
        bin_width=bin_width,
    )
    return {
        "accessible": accessible,
        "reconstructed_T": chart.reconstructed_coords[:, 0],
        "reconstructed_X": chart.reconstructed_coords[:, 1],
        "predecessor_tick_positions": predecessor,
        "successor_tick_positions": successor,
        "radar_time_rank": radar_time_rank,
        "radar_distance_rank": radar_distance_rank,
        "slice_labels": slice_labels,
    }
