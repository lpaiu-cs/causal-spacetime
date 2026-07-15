"""Regression tests pinning the T1 [PROVED] statements to the instrument.

Every assertion here mirrors a claim tagged [PROVED] in
docs/theory/t1_parallax_identifiability.md (v0.2.1). These are not
statistical checks with tolerances tuned to pass: the quantization band,
the fold, and the density invariance are exact consequences of Lemma 1,
Lemma 2 (Model D) and the null-inclusive order convention. A failure here
means the theory document and the code have diverged -- one of them is
wrong, and the test does not say which.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "theory"
sys.path.insert(0, str(EXPERIMENT_DIR))

import t1_verification  # noqa: E402
from t1_verification import (  # noqa: E402
    BAND_TOL,
    bracket_widths_ranks,
    check_builder_density_invariance,
    check_centered_residue,
    check_coincident_tick_orphan,
    check_density_invariance,
    check_fold,
    check_null_aligned_tick,
    check_pipeline_band,
    check_quantization_band,
    check_resolution_scaling,
    pc_v1_centered_profile,
    predicted_center,
    residuals_in_band,
)


def test_quantization_band_holds_on_a_controlled_chain():
    """Lemma 2 (Model D): W - (2|dx|/delta + 1) in [-1, 1), no exceptions."""

    result = check_quantization_band(seed=0, n_targets=150)
    assert result["n_reachable"] > 0
    assert result["violations"] == 0, result
    assert result["residual_min"] >= -1.0 - BAND_TOL
    assert result["residual_max"] < 1.0


def test_band_predicate_rejects_a_true_residual_of_plus_one():
    """The open upper edge must stay open through float fuzz: a target
    spacetime-coincident with a tick (dt = 0, unrelated under dt > 0)
    genuinely produces W = 2 and residual exactly +1 -- outside [-1, 1).
    A symmetric tolerance would have admitted it; the predicate must not."""

    result = check_coincident_tick_orphan()
    assert result["width"] == 2.0
    assert result["residual"] == 1.0
    assert result["band_rejects_it"], result
    # And directly: +1 is out, the closed lower edge -1 is in.
    assert not residuals_in_band(np.array([1.0]))[0]
    assert residuals_in_band(np.array([-1.0]))[0]


def test_null_aligned_example_from_lemma_2_is_pinned():
    """The document's own load-bearing example: target (0, 0.5) against
    ticks -0.75..0.75. Null-inclusive order: W = 4 = N + 1, residual
    exactly -1 (theta = -1, the closed band edge, in band). Strict order:
    W = 6, residual +1, out of band. The one deterministic configuration
    where the two conventions measure different widths."""

    result = check_null_aligned_tick()
    assert result["width_inclusive"] == 4.0
    assert result["residual"] == -1.0
    assert result["in_band"]
    assert result["width_strict_counterfactual"] == 6.0
    assert result["passed"], result


def test_swapping_in_the_strict_order_is_detected(monkeypatch):
    """The review's swap experiment, kept as a regression: replace
    causal_matrix_1p1 with the strict relation dt > |dx| at runtime.
    Before check_null_aligned_tick existed, every harness check stayed
    green under this swap -- all sampled targets sit in general position,
    where the conventions agree -- so the suite did not pin the convention
    it declares load-bearing. Now at least one check must fail."""

    def strict_causal_matrix(events, *, atol=1e-12):
        events = np.asarray(events, dtype=float)
        dt = events[None, :, 0] - events[:, None, 0]
        dx = events[None, :, 1] - events[:, None, 1]
        return (dt > 0.0) & (dt * dt - dx * dx > atol)

    monkeypatch.setattr(t1_verification, "causal_matrix_1p1", strict_causal_matrix)
    result = check_null_aligned_tick()
    assert result["width_inclusive"] == 6.0  # the strict width leaks through
    assert not result["passed"], result


def test_builder_density_invariance_exercises_the_scene_builder():
    """The falsifier at the level the theory doc talks about: scenes built
    by build_positive_control_scene() at different n_events must place
    bit-identical chains and measure bit-identical widths for one fixed
    appended target set. The direct-order variant cannot catch a builder
    that derives its clock from n_events; this can."""

    result = check_builder_density_invariance(
        n_targets=12, n_events_grid=(300, 900)
    )
    assert result["chains_identical"], result
    assert result["widths_identical"], result
    assert result["n_reachable"] > 0
    assert result["passed"]


def test_quantization_band_holds_through_the_scene_pipeline():
    """Same band, measured end-to-end through build_positive_control_scene()
    and measure_bracket_echo_profiles() -- divergence cannot hide in glue."""

    result = check_pipeline_band(seed=0)
    assert result["n_measurements"] > 0
    assert result["violations"] == 0, result


def test_one_observer_folds_mirrored_targets_and_a_second_separates_them():
    """Theorem 1 sharpness: the fold is exact, and parallax breaks it."""

    result = check_fold()
    assert result["fold_exact"], result
    assert result["separation_measured"] > 2.0  # far beyond quantization
    assert (
        abs(result["separation_measured"] - result["separation_predicted"])
        <= 2.0
    )


def test_bulk_sprinkling_density_never_moves_a_single_width():
    """Model D's clock does not know rho: the direct falsifier from the
    review that exposed the v0.1 chain-model error."""

    result = check_density_invariance(seed=1, n_targets=40, bulk_sizes=(0, 500))
    assert result["identical_across_density"], result


def test_position_error_is_bounded_by_half_a_tick_spacing():
    """Section 5 (Model D): |dx|_hat errs by at most delta/2, pointwise,
    at every rung of the tick ladder; RMSE falls like 1/K."""

    result = check_resolution_scaling(
        seed=2, n_targets=120, tick_ladder=(12, 48, 192)
    )
    for row in result["ladder"]:
        assert row["within_bound"], row
    assert -1.3 < result["rmse_loglog_slope"] < -0.7


def test_review_counterexample_is_reproduced_exactly():
    """Lemma 3b: the centered residue exists (a pure time translation moves
    the profile) and stays within the proved 2-rank bound. The two profiles
    are pinned to the reviewer's own numbers from the PR #4 thread."""

    result = check_centered_residue(x=0.10, t_shift=0.003)
    assert result["time_shift_moved_profile"], result
    assert result["max_abs_residue"] < 2.0 + BAND_TOL

    reviewer_t0 = [77 / 3, 35 / 3, -7 / 3, -49 / 3, -49 / 3, -7 / 3]
    reviewer_shifted = [76 / 3, 34 / 3, -8 / 3, -47 / 3, -47 / 3, -8 / 3]
    assert result["profile_t0"] == pytest.approx(reviewer_t0)
    assert result["profile_shifted"] == pytest.approx(reviewer_shifted)


def test_widths_are_integers_in_rank_units():
    """W = s - p over integer ranks: any non-integer would mean the harness
    is not measuring what the theory defines."""

    rng = np.random.default_rng(7)
    targets = np.column_stack([
        rng.uniform(-0.3, 0.3, size=25),
        rng.uniform(-0.3, 0.3, size=25),
    ])
    widths = bracket_widths_ranks(targets, 0.0, 1.4, 96)
    finite = widths[~np.isnan(widths)]
    assert finite.size > 0
    assert np.array_equal(finite, np.round(finite))


def test_predicted_center_is_affine_in_distance():
    """The Lemma 2 center: slope 2/delta, intercept exactly +1."""

    delta = 1.4 / 95
    assert predicted_center(np.array([0.0]), delta)[0] == pytest.approx(1.0)
    lo, hi = predicted_center(np.array([0.1, 0.2]), delta)
    assert hi - lo == pytest.approx(0.2 / delta)


def test_centered_profile_helper_matches_direct_measurement():
    """pc_v1_centered_profile is just W minus its row mean on the PC-V1 grid."""

    profile = pc_v1_centered_profile(0.0, 0.10)
    assert profile.shape == (6,)
    assert profile.sum() == pytest.approx(0.0)
