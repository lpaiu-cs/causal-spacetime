"""Trajectory equivalence for the optional Numba replay accelerator."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("numba")

from causal_spacetime_lab.positive_control.accelerated_two_orders import (
    mcmc_2d_order_replay_accelerated,
)
from causal_spacetime_lab.positive_control.two_orders import mcmc_2d_order_fast


def _assert_replays_match(expected, actual, expected_perms, actual_perms):
    assert len(actual) == len(expected)
    for expected_perm, actual_perm in zip(expected_perms, actual_perms, strict=True):
        assert np.array_equal(actual_perm, expected_perm)
    for expected_row, actual_row in zip(expected, actual, strict=True):
        assert actual_row["S"] == pytest.approx(expected_row["S"], abs=1e-10)
        assert actual_row["n0"] == expected_row["n0"]
        assert actual_row["height"] == expected_row["height"]


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
    _assert_replays_match(expected, actual, expected_perms, actual_perms)


@pytest.mark.slow
def test_accelerated_replay_matches_reference_at_production_n_across_resync():
    pi0 = np.random.default_rng(17).permutation(600)
    kwargs = {
        "beta": 16.0,
        "eps": 0.02,
        "steps": 1_001,
        "seed": 18,
        "sample_every": 1_000,
        "burn_frac": 0.0,
        "collect_perms": True,
    }
    expected, expected_acceptance, expected_perms = mcmc_2d_order_fast(pi0, **kwargs)
    actual, actual_acceptance, actual_perms = mcmc_2d_order_replay_accelerated(
        pi0, resync_every=1_000, **kwargs
    )
    assert actual_acceptance == expected_acceptance
    _assert_replays_match(expected, actual, expected_perms, actual_perms)


@pytest.mark.slow
def test_production_n_trajectory_is_stable_over_50k_steps_between_resyncs():
    pi0 = np.random.default_rng(27).permutation(600)
    kwargs = {
        "beta": 16.0,
        "eps": 0.02,
        "steps": 50_000,
        "seed": 28,
        "sample_every": 10_000,
        "burn_frac": 0.0,
        "collect_perms": True,
    }
    frequent, frequent_acceptance, frequent_perms = (
        mcmc_2d_order_replay_accelerated(pi0, resync_every=1_000, **kwargs)
    )
    deferred, deferred_acceptance, deferred_perms = (
        mcmc_2d_order_replay_accelerated(pi0, resync_every=50_001, **kwargs)
    )
    assert frequent_acceptance == deferred_acceptance
    _assert_replays_match(frequent, deferred, frequent_perms, deferred_perms)
