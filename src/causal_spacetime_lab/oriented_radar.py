"""Oriented two-chain radar reconstruction utilities."""

from __future__ import annotations

from dataclasses import dataclass

from causal_spacetime_lab.discrete_radar import DiscreteRadarCoordinates


@dataclass(frozen=True)
class OrientedRadarCoordinates:
    """Signed coordinate reconstructed relative to a two-chain protocol."""

    target_index: int
    accessible: bool
    primary_radar_time: float | None
    primary_radar_distance: float | None
    beacon_radar_time: float | None
    beacon_radar_distance: float | None
    oriented_time: float | None
    signed_position: float | None


def signed_position_from_two_distances(
    primary_distance: float,
    beacon_distance: float,
    beacon_separation: float,
) -> float:
    """Return signed coordinate from distances to two synchronized chains."""

    primary = float(primary_distance)
    beacon = float(beacon_distance)
    separation = float(beacon_separation)
    if separation <= 0.0:
        raise ValueError("beacon_separation must be positive")
    if primary < 0.0 or beacon < 0.0:
        raise ValueError("distances must be nonnegative")
    return (primary * primary - beacon * beacon + separation * separation) / (
        2.0 * separation
    )


def oriented_radar_coordinates_from_two_chains(
    primary_results: list[DiscreteRadarCoordinates],
    beacon_results: list[DiscreteRadarCoordinates],
    beacon_separation: float,
) -> list[OrientedRadarCoordinates]:
    """Reconstruct signed coordinates relative to a supplied two-chain protocol."""

    if len(primary_results) != len(beacon_results):
        raise ValueError("primary_results and beacon_results must have equal length")

    oriented: list[OrientedRadarCoordinates] = []
    for primary, beacon in zip(primary_results, beacon_results, strict=True):
        if primary.target_index != beacon.target_index:
            raise ValueError("primary and beacon target indices must match")
        if (
            not primary.accessible
            or not beacon.accessible
            or primary.radar_time is None
            or primary.radar_distance is None
            or beacon.radar_time is None
            or beacon.radar_distance is None
        ):
            oriented.append(
                OrientedRadarCoordinates(
                    target_index=primary.target_index,
                    accessible=False,
                    primary_radar_time=primary.radar_time,
                    primary_radar_distance=primary.radar_distance,
                    beacon_radar_time=beacon.radar_time,
                    beacon_radar_distance=beacon.radar_distance,
                    oriented_time=None,
                    signed_position=None,
                )
            )
            continue

        signed_position = signed_position_from_two_distances(
            primary.radar_distance,
            beacon.radar_distance,
            beacon_separation,
        )
        oriented_time = 0.5 * (primary.radar_time + beacon.radar_time)
        oriented.append(
            OrientedRadarCoordinates(
                target_index=primary.target_index,
                accessible=True,
                primary_radar_time=primary.radar_time,
                primary_radar_distance=primary.radar_distance,
                beacon_radar_time=beacon.radar_time,
                beacon_radar_distance=beacon.radar_distance,
                oriented_time=oriented_time,
                signed_position=signed_position,
            )
        )
    return oriented

