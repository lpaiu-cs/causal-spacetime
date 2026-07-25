"""Regressions for G4c: the general-dimension reading of the G4b threshold.

These pin what was MEASURED, not what was predicted. The predictions
live in `docs/theory/t1_g4c_predictions.json` and were committed before
the harness existed; wiring CI to them instead would mean a wrong
prediction turns the build red, which is quiet pressure to edit the
prediction. The scorecard is allowed to say MISS.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "theory"
sys.path.insert(0, str(EXPERIMENT_DIR))

from t1_g4b_unlabeled_2plus1d import (  # noqa: E402
    flatten,
    gauge_dimension,
    rigid_motion_gauge,
    scene,
    nullity,
)
from t1_g4c_general_dimension import (  # noqa: E402
    affine_rank,
    check_machinery_at_d3,
    consecutive_differences,
    counting_bound,
    nullity_by_n,
)


def test_generalized_gauge_still_emits_the_old_plane_columns():
    """The pair loop replaced a hand-written ``(-y, x)`` special case.

    If it ever stops reproducing that column exactly, every 2+1D number
    in G4b shifts underneath us without any test in that file noticing,
    because they all measure nullity against this basis.
    """

    theta = np.array([0.3, -0.7, 1.1, 0.2, -0.5, 0.9])
    gauge = rigid_motion_gauge(theta, 2)
    points = theta.reshape(-1, 2)

    assert gauge.shape[1] == gauge_dimension(2) == 3
    np.testing.assert_allclose(gauge[:, 0], np.array([1.0, 0.0] * 3))
    np.testing.assert_allclose(gauge[:, 1], np.array([0.0, 1.0] * 3))
    np.testing.assert_allclose(
        gauge[:, 2], np.column_stack((-points[:, 1], points[:, 0])).ravel()
    )


def test_gauge_dimension_grows_as_d_times_d_plus_one_over_two():
    theta = np.zeros(4 * 7)
    for d, expected in ((1, 1), (2, 3), (3, 6), (4, 10)):
        columns = rigid_motion_gauge(theta[: d * 7], d).shape[1]
        assert columns == expected == gauge_dimension(d)


def test_machinery_validates_in_three_dimensions():
    """Six gauge directions must be exact nulls, scale must not be, and
    the observers must actually span the space they are placed in."""

    outcome = check_machinery_at_d3()
    assert outcome["passed"]
    assert outcome["gauge_columns"] == 6
    assert max(outcome["gauge_residuals"]) < 1e-10
    assert outcome["scale_response"] > 0.5
    assert outcome["observer_affine_rank"] == 3


def test_observers_placed_off_the_plane_really_span_three_dimensions():
    _, P = scene(10, 5, seed=5, d=3)
    assert affine_rank(P) == 3


def test_three_observers_in_three_dimensions_leave_a_fibre_per_target():
    """R - 1 = 2 < 3 = d, so each target slides along a 1-dimensional
    fibre and the flex count grows one per target."""

    rows = nullity_by_n(range(6, 13), R=3, d=3, seed_base=800)
    assert [row["nullity"] for row in rows] == [18, 19, 20, 21, 22, 23, 24]
    assert set(consecutive_differences(rows)) == {1}


def test_two_observers_in_the_plane_leave_a_fibre_per_target():
    """The same regime one dimension down -- the case G4b never ran.

    This is the retro-check: the law claims an under-observed regime
    exists at d = 2 as well, and if it did not, the law would be wrong
    in the dimension already understood.
    """

    rows = nullity_by_n(range(6, 13), R=2, d=2, seed_base=700)
    assert [row["nullity"] for row in rows] == [11, 12, 13, 14, 15, 16, 17]
    assert set(consecutive_differences(rows)) == {1}


def test_four_observers_in_three_dimensions_saturate_at_eighteen():
    """R - 1 = 3 = d: the profile surface fills its ambient, so targets
    stop buying information. Flat from 6 to 34 targets."""

    rows = nullity_by_n((6, 10, 14, 20, 34), R=4, d=3, seed_base=900)
    assert {row["nullity"] for row in rows} == {18}
    assert {row["extra_flexes"] for row in rows} == {12}


def test_five_observers_in_three_dimensions_reach_rigidity():
    """R >= d + 2 = 5: the surface curves and the nullity falls to the
    six-dimensional rigid-motion gauge. Measured threshold n = 19."""

    for n, expect_rigid in ((18, False), (19, True), (25, True)):
        X, P = scene(n, 5, seed=1000 + 13 * 0 + n, d=3)
        null, _ = nullity(flatten(X, P), n, 5, 3)
        assert (null == gauge_dimension(3)) is expect_rigid, (n, null)


def _is_rigid(n: int, R: int, d: int, seed: int, shell_variant: int = 0) -> bool:
    X, P = scene(n, R, seed=seed, d=d, shell_variant=shell_variant)
    null, _ = nullity(flatten(X, P), n, R, d, )
    return null == gauge_dimension(d)


def test_fibre_slope_tracks_d_minus_r_plus_one_not_just_one():
    """The out-of-sample cell that separates the formula from a pattern.

    Every cell measured before round 2 -- d=2 R=2, d=3 R=3 -- had fibre
    dimension exactly 1, so "slope 1" fitted all of it. At d=4, R=3 the
    formula says 2 and mere pattern-matching says 1.
    """

    steep = nullity_by_n(range(6, 13), R=3, d=4, seed_base=1100)
    assert [row["nullity"] for row in steep] == [27, 29, 31, 33, 35, 37, 39]
    assert set(consecutive_differences(steep)) == {2}

    shallow = nullity_by_n(range(6, 13), R=4, d=4, seed_base=1200)
    assert [row["nullity"] for row in shallow] == [28, 29, 30, 31, 32, 33, 34]
    assert set(consecutive_differences(shallow)) == {1}


def test_counting_law_holds_in_four_and_five_dimensions():
    """nullity = d*R + d(d+1)/2 in the saturated regime, now at its
    third and fourth confirmation (9 at d=2, 18 at d=3)."""

    d4 = nullity_by_n((6, 10, 14, 20, 30), R=5, d=4, seed_base=1300)
    assert {row["nullity"] for row in d4} == {4 * 5 + gauge_dimension(4)} == {30}

    d5 = nullity_by_n((8, 12, 16, 20), R=6, d=5, seed_base=1700)
    assert {row["nullity"] for row in d5} == {5 * 6 + gauge_dimension(5)} == {45}


def test_rigidity_thresholds_at_r_equals_d_plus_two():
    """Measured 11, 19, 29, 41 for d = 2, 3, 4, 5.

    The d = 5 point was predicted at 41 by n = d^2 + 3d + 1 fitted to the
    first three, and committed before the scan ran; it is the only thing
    that makes the formula more than interpolation.
    """

    assert not _is_rigid(28, 6, 4, seed=1400 + 8)
    assert _is_rigid(29, 6, 4, seed=1400 + 29)
    assert not _is_rigid(40, 7, 5, seed=1500 + 40)
    assert _is_rigid(41, 7, 5, seed=1500 + 41)

    for d, measured in ((2, 11), (3, 19), (4, 29), (5, 41)):
        assert d * d + 3 * d + 1 == measured


def test_threshold_survives_a_redrawn_observer_shell():
    """A threshold is a claim about where a rank transition sits, so it
    could have been an artifact of one arbitrary observer placement.
    Redrawing the shell leaves d = 3, R = 5 at 19."""

    assert not _is_rigid(18, 5, 3, seed=1600 + 18, shell_variant=1)
    assert _is_rigid(19, 5, 3, seed=1600 + 19, shell_variant=1)


def test_shell_variant_zero_leaves_the_plane_construction_alone():
    """The variant knob must not perturb any existing G4b number."""

    plain = scene(9, 4, seed=42, d=2)
    same = scene(9, 4, seed=42, d=2, shell_variant=0)
    np.testing.assert_array_equal(plain[0], same[0])
    np.testing.assert_array_equal(plain[1], same[1])
    assert not np.allclose(scene(9, 4, seed=42, d=2, shell_variant=1)[1], plain[1])


def test_counting_bound_is_necessary_but_was_not_tight_anywhere():
    """It gave 7 where the plane measured 11, and 9 where three
    dimensions measured 19. Recorded so the gap is visible rather than
    quietly assumed away."""

    assert counting_bound(2, 4) == 7
    assert counting_bound(3, 5) == 9
    assert counting_bound(3, 6) == 10
    assert counting_bound(3, 8) == 11
