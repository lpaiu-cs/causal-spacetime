from __future__ import annotations

import numpy as np
import pytest

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.observer_atlas import (
    ObserverProtocolSpec,
    affine_lab_to_rest_1p1,
    affine_rest_to_lab_1p1,
    common_safe_tau_range_for_oriented_protocol_1p1,
    make_affine_inertial_chain_1p1,
    reconstruct_oriented_chart_from_order,
    safe_tau_range_for_affine_inertial_chain_1p1,
)


def test_affine_rest_to_lab_and_lab_to_rest_are_inverse_maps() -> None:
    rest = np.array([[-0.2, 0.1], [0.0, 0.0], [0.4, -0.2]])
    beta = 0.35
    origin = (0.05, -0.04)

    lab = affine_rest_to_lab_1p1(rest, beta, origin)
    round_trip = affine_lab_to_rest_1p1(lab, beta, origin)

    assert round_trip == pytest.approx(rest)


def test_make_affine_inertial_chain_validates_inputs() -> None:
    with pytest.raises(ValueError, match=r"abs\(beta\) < 1"):
        make_affine_inertial_chain_1p1(1.0, 4, -0.5, 0.5)

    with pytest.raises(ValueError, match="at least 2"):
        make_affine_inertial_chain_1p1(0.0, 1, -0.5, 0.5)


def test_safe_tau_range_for_affine_chain_lies_inside_diamond() -> None:
    tau_min, tau_max = safe_tau_range_for_affine_inertial_chain_1p1(
        T_global=2.0,
        beta=0.3,
        origin_lab=(0.05, -0.05),
        x_prime=0.1,
    )

    events, _ = make_affine_inertial_chain_1p1(
        beta=0.3,
        num_ticks=32,
        tau_min=tau_min,
        tau_max=tau_max,
        origin_lab=(0.05, -0.05),
        x_prime=0.1,
    )

    assert np.all(np.abs(events[:, 1]) <= 1.0 - np.abs(events[:, 0]) + 1e-12)


def test_common_safe_tau_range_covers_primary_and_beacon() -> None:
    spec = ObserverProtocolSpec(
        name="B",
        beta=0.3,
        origin_lab_time=0.05,
        origin_lab_position=-0.05,
        beacon_separation=0.15,
    )

    tau_min, tau_max = common_safe_tau_range_for_oriented_protocol_1p1(2.0, spec)
    primary, _ = make_affine_inertial_chain_1p1(
        spec.beta,
        32,
        tau_min,
        tau_max,
        origin_lab=spec.origin_lab,
        x_prime=0.0,
    )
    beacon, _ = make_affine_inertial_chain_1p1(
        spec.beta,
        32,
        tau_min,
        tau_max,
        origin_lab=spec.origin_lab,
        x_prime=spec.beacon_separation,
    )

    assert np.all(np.abs(primary[:, 1]) <= 1.0 - np.abs(primary[:, 0]) + 1e-12)
    assert np.all(np.abs(beacon[:, 1]) <= 1.0 - np.abs(beacon[:, 0]) + 1e-12)


def test_reconstruct_oriented_chart_from_order_deterministic_example() -> None:
    clock_times = np.array([-1.0, -0.4, -0.2, 0.0, 0.2, 0.4, 1.0])
    primary = np.column_stack((clock_times, np.zeros(clock_times.size)))
    beacon = np.column_stack((clock_times, np.full(clock_times.size, 0.6)))
    target = np.array([[0.0, 0.2]])
    events = np.vstack((primary, beacon, target))
    causal_matrix = causal_matrix_1p1(events)

    chart = reconstruct_oriented_chart_from_order(
        causal_matrix,
        target_indices=np.array([14]),
        primary_indices=np.arange(0, 7),
        beacon_indices=np.arange(7, 14),
        clock_times=clock_times,
        beacon_separation=0.6,
        protocol_name="A",
    )

    assert chart.protocol_name == "A"
    assert chart.accessible.tolist() == [True]
    assert chart.reconstructed_coords[0, 0] == pytest.approx(0.0)
    assert chart.reconstructed_coords[0, 1] == pytest.approx(0.2)
