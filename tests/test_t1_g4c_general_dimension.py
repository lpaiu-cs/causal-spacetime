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


def test_counting_bound_is_necessary_but_was_not_tight_anywhere():
    """It gave 7 where the plane measured 11, and 9 where three
    dimensions measured 19. Recorded so the gap is visible rather than
    quietly assumed away."""

    assert counting_bound(2, 4) == 7
    assert counting_bound(3, 5) == 9
    assert counting_bound(3, 6) == 10
    assert counting_bound(3, 8) == 11
