"""Regression tests for G4b: what D determines in 2+1D.

The claims these pin are structural, not statistical: a validated
Jacobian, an exact flex count, and an explicit same-D counterexample.
A failure means a statement in the theory document is false, not that a
sample was unlucky.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "theory"
sys.path.insert(0, str(EXPERIMENT_DIR))

from t1_g4b_unlabeled_2plus1d import (  # noqa: E402
    check_frozen_scene_is_in_the_rigid_regime,
    check_r3_explicit_counterexample,
    dissimilarity,
    flatten,
    gauge_dimension,
    jacobian,
    nullity,
    rigid_motion_gauge,
    scene,
    scene_edm,
)


def test_rigid_motions_are_null_directions_and_scale_is_not():
    """Validates the instrument before any verdict rests on it. D is
    invariant under rigid motion of the whole scene and homogeneous of
    degree 1 under scaling, so J @ (scale direction) reproduces D
    exactly -- which is why the absolute scale is recoverable."""

    X, P = scene(10, 5, seed=5)
    theta = flatten(X, P)
    J = jacobian(theta, 10, 5, 2)
    base = dissimilarity(theta, 10, 5, 2)

    gauge = rigid_motion_gauge(theta, 2)
    assert gauge.shape[1] == gauge_dimension(2) == 3
    for k in range(gauge.shape[1]):
        assert np.linalg.norm(J @ gauge[:, k]) < 1e-10 * np.linalg.norm(base)

    # Euler's identity for a degree-1 homogeneous map: J @ theta == D
    assert np.allclose(J @ theta, base, rtol=1e-9, atol=1e-12)


def test_1p1d_is_never_rigid():
    """The control. Lemma 4f says D carries order and nothing metric in
    1+1D, so extra flexes must persist however much data is added."""

    for n, R in ((6, 3), (16, 6), (24, 8)):
        X, P = scene(n, R, seed=100 + n + 7 * R, d=1)
        theta = flatten(X, P)
        null, _ = nullity(theta, n, R, 1)
        assert null > gauge_dimension(1), (n, R, null)


def test_three_observers_keep_exactly_six_extra_flexes():
    """With R = 3 the centered profiles span the 2-dimensional mean-zero
    subspace of R^3, so the profile surface cannot curve inside its
    ambient space. The deficiency is structural: it is the same 6 at
    every target count, so no amount of data removes it."""

    for n in (6, 14, 34):
        X, P = scene(n, 3, seed=200 + n)
        theta = flatten(X, P)
        null, _ = nullity(theta, n, 3, 2)
        assert null - gauge_dimension(2) == 6, (n, null)


def test_four_observers_reach_rigidity_where_three_never_do():
    """The threshold that matters: at a target count where R = 3 is
    still 6-fold flexible, R >= 4 has nullity exactly equal to the
    rigid-motion gauge."""

    n = 14
    flexible, _ = nullity(flatten(*scene(n, 3, seed=200 + n)), n, 3, 2)
    assert flexible - gauge_dimension(2) == 6

    for R in (4, 5, 6, 8):
        X, P = scene(n, R, seed=500 + n)
        null, _ = nullity(flatten(X, P), n, R, 2)
        assert null == gauge_dimension(2), (R, null)


def test_explicit_same_dissimilarity_counterexample_at_r3():
    """The decisive form of 'R = 3 does not suffice': a concrete scene
    with the same D to machine precision and a different shape. This is
    the 2+1D analogue of Lemma 4f's same-D counterexample."""

    result = check_r3_explicit_counterexample(n=12, steps=60)

    assert result["max_D_drift"] < 1e-12
    # the two scenes must differ by an amount comparable to the
    # configuration itself, not by a numerical whisker
    assert result["scene_shape_gap"] > 0.05 * result["target_spread"]
    before = np.array(result["observer_triangle_before"])
    after = np.array(result["observer_triangle_after"])
    assert np.max(np.abs(before - after)) > 1e-3
    assert result["passed"]


def test_counterexample_really_preserves_the_observable():
    """Independent of the harness's own bookkeeping: take the two scenes
    it returns, recompute D for each from scratch, and confirm the
    observable agrees while the scenes themselves do not."""

    result = check_r3_explicit_counterexample(n=12, steps=60)
    assert result["arc_length"] > 0.0

    before = np.array(result["scene_before"])
    after = np.array(result["scene_after"])
    d_before = dissimilarity(before, 12, 3, 2)
    d_after = dissimilarity(after, 12, 3, 2)

    # the observable is identical to machine precision ...
    assert np.max(np.abs(d_before - d_after)) < 1e-12
    assert np.all(d_before > 0.0)
    # ... and the scenes are not congruent, by a wide margin
    assert np.max(np.abs(scene_edm(before) - scene_edm(after))) > 1e-3


def test_frozen_2plus1d_instrument_sits_in_the_rigid_regime():
    """The configuration the instrument actually uses -- 8 chains on the
    frozen non-collinear ring, its own selected targets -- must be
    rigid, with a clear margin rather than a marginal one."""

    result = check_frozen_scene_is_in_the_rigid_regime()

    assert result["n_observers"] == 8
    assert result["n_targets"] >= 30
    assert result["extra_flexes"] == 0
    assert result["nullity"] == gauge_dimension(2)
    # the gap between the true zeros and the smallest real singular
    # value is what makes the rank decision safe
    assert result["smallest_nonzero_singular_value"] > 1e-4
    assert result["passed"]


def test_scene_edm_separates_congruence_from_genuine_change():
    """The fingerprint used to judge the counterexample: a rigid motion
    must leave it untouched, a real deformation must not."""

    X, P = scene(8, 4, seed=1)
    theta = flatten(X, P)
    base = scene_edm(theta)

    points = theta.reshape(-1, 2)
    angle = 0.7
    rot = np.array([[np.cos(angle), -np.sin(angle)],
                    [np.sin(angle), np.cos(angle)]])
    moved = (points @ rot.T) + np.array([0.3, -0.2])
    assert np.allclose(scene_edm(moved.ravel()), base, atol=1e-12)

    reflected = points * np.array([1.0, -1.0])
    assert np.allclose(scene_edm(reflected.ravel()), base, atol=1e-12)

    stretched = points * np.array([1.1, 1.0])
    assert not np.allclose(scene_edm(stretched.ravel()), base, atol=1e-6)
