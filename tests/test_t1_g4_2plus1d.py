"""Regression tests for G4a: T1 statements that survive in 2+1D.

The checks are deterministic (seed-pinned) or exact, so a failure means
a claim in the theory document is false in 2+1D, not that a sample was
unlucky. Reduced sizes keep CI fast; the full grid lives in
docs/theory/t1_g4_2plus1d_results.json.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "theory"
sys.path.insert(0, str(EXPERIMENT_DIR))

from t1_g4_2plus1d import (  # noqa: E402
    _ring,
    _sample_targets,
    check_builder_density_invariance,
    check_fold_structure,
    check_model_p_same_slice,
    check_multilateration,
    check_pipeline_band,
    measure_widths,
    multilaterate,
    radial_distances,
    residuals_in_band,
)


def test_quantization_band_holds_in_2plus1d():
    """Lemma 2's band never used the dimension: with |dx| read as a
    Euclidean norm, every measured width must sit in [-1, 1) around
    2|dx|/delta + 1."""

    ticks, span = 96, 1.4
    delta = span / (ticks - 1)
    rng = np.random.default_rng(0)
    positions = _ring(4, 0.25)
    targets = _sample_targets(rng, 60, 0.10, 0.22)

    widths = measure_widths(targets, positions, ticks)
    distances = radial_distances(targets, positions)
    reachable = ~np.isnan(widths)

    assert reachable.all()
    residuals = widths[reachable] - (2.0 * distances[reachable] / delta + 1.0)
    assert residuals_in_band(residuals).all()


def test_single_observer_cannot_separate_a_circle_of_equal_radius():
    """The 2+1D fold is a circle, not a mirror pair: equal radius at
    equal time gives EXACTLY equal integer widths."""

    result = check_fold_structure()
    assert result["single_observer_folds_the_circle"]
    assert result["two_observers_fold_the_mirror_pair"]
    # a third observer off the mirror line separates them decisively --
    # far beyond the +/-1 quantization spread of each width
    assert result["third_observer_width_gap"] > 2.0
    assert result["passed"]


def test_multilateration_recovers_position_within_the_proved_bound():
    """Theorem 1's labeled clause in 2+1D: three non-collinear observers
    determine the target's 2D position, with error under the O(delta)
    bound from the squared-distance linear system."""

    result = check_multilateration(seed=7, n_targets=40)
    good = result["non_collinear"]
    assert good["n_solved"] == 40
    assert good["all_within_proved_bound"]
    assert good["max_position_error"] < 0.05


def test_near_collinear_observers_degrade_exactly_as_the_bound_says():
    """'Non-collinear' is quantitative in 2+1D: the error scales with
    ||A^+||, so a nearly degenerate layout is measurably worse while
    still respecting its own (much larger) bound. This has no 1+1D
    analogue."""

    result = check_multilateration(seed=7, n_targets=40)
    good, bad = result["non_collinear"], result["near_collinear"]

    assert bad["pinv_norm"] > 10.0 * good["pinv_norm"]
    assert bad["median_position_error"] > 10.0 * good["median_position_error"]
    assert bad["all_within_proved_bound"]


def test_multilateration_is_exact_without_quantization():
    """Sanity on the solver itself: fed exact distances (delta -> 0 in
    the estimator), the linear system returns the true position to
    machine precision, so the measured error is quantization, not a
    reconstruction defect."""

    positions = _ring(3, 0.25)
    truth = np.array([0.07, -0.11])
    distances = np.linalg.norm(positions - truth, axis=1)
    # invert d_hat = (W - 1) * delta / 2 exactly
    delta = 0.01
    widths = 2.0 * distances / delta + 1.0

    estimate, pinv_norm = multilaterate(widths, positions, delta)
    assert np.allclose(estimate, truth, atol=1e-12)
    assert pinv_norm > 0.0


def test_same_slice_pairs_stay_pathwise_monotone_in_2plus1d():
    """Theorem 2 in 2+1D: concentric brackets, so ties happen and
    strict inversions do not. The tie constant counts observers --
    exp(-2 lambda g) for one, exp(-4 lambda g) for two flanking."""

    result = check_model_p_same_slice(trials=1500)
    for block in ("single_observer", "two_observer_flanking"):
        assert result[block]["strict_inversions"] == 0
        assert abs(
            result[block]["tie_rate"] - result[block]["tie_rate_predicted"]
        ) < 0.04


def test_frozen_2plus1d_builder_keeps_the_clock_independent_of_density():
    """The Model-D falsifier through the frozen builder: chain
    worldlines must not move when the sprinkling density changes."""

    assert check_builder_density_invariance()["passed"]


def test_band_holds_end_to_end_through_the_frozen_builder():
    result = check_pipeline_band(seed=0)
    assert result["violations"] == 0
    assert result["n_measurements"] > 0
    assert result["passed"]
