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
    measure_widths_full_order,
    multilaterate,
    multilateration_matrix,
    multilateration_pinv_norm,
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


def test_near_collinear_observers_stay_inside_their_own_larger_bound():
    """'Non-collinear' is quantitative in 2+1D: the bound scales with
    ||A^+||, so a nearly degenerate layout gets a much larger bound and
    must still respect it. The degradation RATIO is recorded, not
    asserted -- the bound is an upper bound, so theory does not predict
    how much worse the layout actually performs."""

    result = check_multilateration(seed=7, n_targets=40)
    good, bad = result["non_collinear"], result["near_collinear"]

    assert bad["pinv_norm"] > 10.0 * good["pinv_norm"]
    assert bad["all_within_proved_bound"]
    assert bad["max_proved_bound"] > good["max_proved_bound"]

    recorded = result["conditioning_recorded"]
    assert recorded["gating"] is False
    assert recorded["median_error_ratio"] > 1.0


def test_multilateration_is_exact_without_quantization():
    """Sanity on the solver itself: fed exact distances (delta -> 0 in
    the estimator), the linear system returns the true position to
    machine precision, so the measured error is quantization, not a
    reconstruction defect. Also: passing the prebuilt design matrix must
    not change the answer."""

    positions = _ring(3, 0.25)
    truth = np.array([0.07, -0.11])
    distances = np.linalg.norm(positions - truth, axis=1)
    # invert d_hat = (W - 1) * delta / 2 exactly
    delta = 0.01
    widths = 2.0 * distances / delta + 1.0

    estimate = multilaterate(widths, positions, delta)
    assert np.allclose(estimate, truth, atol=1e-12)

    matrix = multilateration_matrix(positions)
    assert np.allclose(
        estimate, multilaterate(widths, positions, delta, matrix), atol=0.0
    )
    assert multilateration_pinv_norm(matrix) > 0.0


def test_fast_and_library_width_paths_agree_exactly():
    """measure_widths evaluates only the target-to-tick block; the
    library path builds the full causal matrix and calls
    find_radar_ticks_from_order. They must agree entry for entry,
    including which entries are unreachable -- otherwise the cheap path
    is measuring something else."""

    rng = np.random.default_rng(3)
    positions = _ring(4, 0.25)
    targets = _sample_targets(rng, 25, 0.10, 0.22)
    bulk = np.column_stack((
        rng.uniform(-0.5, 0.5, size=200),
        rng.uniform(-0.3, 0.3, size=200),
        rng.uniform(-0.3, 0.3, size=200),
    ))

    for extra in (None, bulk):
        fast = measure_widths(targets, positions, 96, extra_events=extra)
        library = measure_widths_full_order(
            targets, positions, 96, extra_events=extra
        )
        assert np.array_equal(np.isnan(fast), np.isnan(library))
        assert np.array_equal(
            np.nan_to_num(fast, nan=-1.0), np.nan_to_num(library, nan=-1.0)
        )


def test_same_slice_pairs_stay_pathwise_monotone_in_2plus1d():
    """Theorem 2 in 2+1D: concentric brackets, so ties happen and
    strict inversions do not. The tie constant counts observers --
    exp(-2 lambda g) for one, exp(-4 lambda g) for two flanking."""

    result = check_model_p_same_slice(trials=1500)
    for block in ("single_observer", "two_observer_flanking"):
        assert result[block]["strict_inversions"] == 0
        # no draw may enter the statistics without a real bracket
        assert result[block]["unreachable_trials"] == 0
        assert abs(
            result[block]["tie_rate"] - result[block]["tie_rate_predicted"]
        ) < 0.04


def test_unbracketed_draws_are_counted_not_silently_used():
    """The reachability guard must actually bite: at a rate low enough
    that draws routinely have no tick outside the radar interval, the
    check must report unreachable trials and fail rather than feeding
    fabricated widths into the tie statistics."""

    result = check_model_p_same_slice(lam=0.4, gap=0.02, trials=300)
    assert result["single_observer"]["unreachable_trials"] > 0
    assert not result["passed"]


def test_frozen_2plus1d_builder_keeps_the_clock_independent_of_density():
    """The Model-D falsifier through the frozen builder: chain
    worldlines must not move when the sprinkling density changes."""

    assert check_builder_density_invariance()["passed"]


def test_band_holds_end_to_end_through_the_frozen_builder():
    result = check_pipeline_band(seed=0)
    assert result["violations"] == 0
    assert result["n_measurements"] > 0
    assert result["passed"]
