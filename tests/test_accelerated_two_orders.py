"""Trajectory equivalence for the optional Numba replay accelerator."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("numba")

from causal_spacetime_lab.positive_control.accelerated_two_orders import (
    mcmc_2d_order_replay_accelerated,
)
from causal_spacetime_lab.positive_control.two_orders import mcmc_2d_order_fast


def test_accelerated_replay_matches_validated_sampler_trajectory():
    pi0 = np.random.default_rng(7).permutation(30)
    expected, expected_acceptance, expected_perms = mcmc_2d_order_fast(
        pi0,
        beta=0.8,
        eps=0.2,
        steps=8_000,
        seed=8,
        sample_every=200,
        collect_perms=True,
    )
    actual, actual_acceptance, actual_perms = mcmc_2d_order_replay_accelerated(
        pi0,
        beta=0.8,
        eps=0.2,
        steps=8_000,
        seed=8,
        sample_every=200,
        collect_perms=True,
    )
    assert actual_acceptance == expected_acceptance
    assert len(actual) == len(expected)
    for expected_perm, actual_perm in zip(expected_perms, actual_perms, strict=True):
        assert np.array_equal(actual_perm, expected_perm)
    for expected_row, actual_row in zip(expected, actual, strict=True):
        assert actual_row["S"] == pytest.approx(expected_row["S"], abs=1e-10)
        assert actual_row["n0"] == expected_row["n0"]
        assert actual_row["height"] == expected_row["height"]
