from __future__ import annotations

import pytest

from causal_spacetime_lab.discrete_radar import DiscreteRadarCoordinates
from causal_spacetime_lab.oriented_radar import (
    oriented_radar_coordinates_from_two_chains,
    signed_position_from_two_distances,
)


@pytest.mark.parametrize("x", [-0.4, 0.0, 0.2, 0.7])
def test_signed_position_from_two_distances_recovers_hand_coded_x(x: float) -> None:
    beacon_separation = 0.5
    primary_distance = abs(x)
    beacon_distance = abs(x - beacon_separation)

    signed_x = signed_position_from_two_distances(
        primary_distance,
        beacon_distance,
        beacon_separation,
    )

    assert signed_x == pytest.approx(x)


def test_signed_position_from_two_distances_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="positive"):
        signed_position_from_two_distances(0.1, 0.2, 0.0)

    with pytest.raises(ValueError, match="nonnegative"):
        signed_position_from_two_distances(-0.1, 0.2, 0.5)


def test_oriented_radar_coordinates_from_two_chains_accessible_case() -> None:
    primary = [
        DiscreteRadarCoordinates(
            target_index=7,
            accessible=True,
            tau_minus=-0.3,
            tau_plus=0.3,
            radar_time=0.0,
            radar_distance=0.3,
        )
    ]
    beacon = [
        DiscreteRadarCoordinates(
            target_index=7,
            accessible=True,
            tau_minus=-0.2,
            tau_plus=0.4,
            radar_time=0.1,
            radar_distance=0.3,
        )
    ]

    result = oriented_radar_coordinates_from_two_chains(primary, beacon, 0.6)[0]

    assert result.accessible
    assert result.oriented_time == pytest.approx(0.05)
    assert result.signed_position == pytest.approx(0.3)


def test_oriented_radar_coordinates_from_two_chains_inaccessible_case() -> None:
    primary = [
        DiscreteRadarCoordinates(
            target_index=7,
            accessible=True,
            tau_minus=-0.3,
            tau_plus=0.3,
            radar_time=0.0,
            radar_distance=0.3,
        )
    ]
    beacon = [
        DiscreteRadarCoordinates(
            target_index=7,
            accessible=False,
            tau_minus=None,
            tau_plus=None,
            radar_time=None,
            radar_distance=None,
        )
    ]

    result = oriented_radar_coordinates_from_two_chains(primary, beacon, 0.6)[0]

    assert not result.accessible
    assert result.oriented_time is None
    assert result.signed_position is None
