"""Regressions for the G4c proof program (`docs/theory/t1_g4c_proof.md`).

These pin the LEMMAS, not only the conclusion. A closed form can come out
right while the argument behind it is wrong, and that is the failure mode
worth catching early: the conclusion is measured in half a dozen other
places, the derivation is not.

As elsewhere in the G4c series, CI pins measurements rather than the
frozen predictions; a missed prediction belongs on the scorecard.
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
    jacobian,
    nullity,
    scene,
)
from t1_g4c_general_dimension import nullity_by_n  # noqa: E402
from t1_g4c_proof_check import (  # noqa: E402
    MEASURED_THRESHOLDS,
    closed_form_nullity,
    cloud_flex_space,
    intersection_dim,
    null_basis,
    operator_L,
    rank_of,
    threshold_count,
)


def _cell(d: int, R: int, n: int) -> np.ndarray:
    X, P = scene(n, R, seed=4242 + n, d=d)
    return flatten(X, P)


def test_lemma_a_every_flex_moves_the_cloud_isometrically():
    """The keystone. If a flex of D induced a cloud motion that was not
    an infinitesimal isometry, the whole derivation would be about a
    different object than the one the instrument measures.
    """

    for d, R, n in ((2, 3, 6), (3, 4, 6), (2, 4, 11)):
        theta = _cell(d, R, n)
        flexes = null_basis(jacobian(theta, n, R, d))
        L = operator_L(theta, n, R, d)
        flex_cloud, _, _ = cloud_flex_space(theta, n, R, d)
        for k in range(flexes.shape[1]):
            image = L @ flexes[:, k]
            norm = np.linalg.norm(image)
            if norm < 1e-9:
                continue
            residual = image - flex_cloud @ (flex_cloud.T @ image)
            assert np.linalg.norm(residual) / norm < 1e-8, (d, R, n, k)


def test_lemma_b_flex_dimension_formula():
    """dim F = n(m - q) + q(q+1)/2 against the measured affine span."""

    for d, R, n in ((2, 2, 6), (2, 3, 6), (3, 4, 6), (3, 5, 12), (4, 5, 6)):
        theta = _cell(d, R, n)
        flex, q, m = cloud_flex_space(theta, n, R, d)
        assert flex.shape[1] == n * (m - q) + q * (q + 1) // 2, (d, R, n)


def test_lemma_c_decomposition_holds_on_rigid_cells_too():
    """nullity = dim ker L + dim(Im L cap F). The rigid cell is the one
    that matters: there the intersection is 0 and the identity is doing
    real work rather than restating a coincidence."""

    for d, R, n, expect in ((2, 3, 6, 9), (2, 4, 11, 3), (3, 5, 19, 6)):
        theta = _cell(d, R, n)
        measured, _ = nullity(theta, n, R, d)
        L = operator_L(theta, n, R, d)
        flex, _, _ = cloud_flex_space(theta, n, R, d)
        ker_L = theta.size - rank_of(L)
        assert measured == expect
        assert measured == ker_L + intersection_dim(L, flex), (d, R, n)


def test_theorem_1_closed_form_matches_measurement():
    """nullity = dR + n(d - R + 1) + R(R-1)/2 for R <= d + 1."""

    for d, R, n in ((2, 2, 6), (2, 3, 10), (3, 3, 9), (3, 4, 12),
                    (4, 3, 6), (4, 4, 6), (4, 5, 6)):
        theta = _cell(d, R, n)
        measured, _ = nullity(theta, n, R, d)
        assert measured == closed_form_nullity(d, R, n), (d, R, n)


def test_theorem_2a_bound_is_respected_and_attained_everywhere_measured():
    """Necessity is proved; attainment is the open half. Both are
    recorded so a future cell that respects the bound without attaining
    it shows up as a change here rather than as a silent weakening."""

    for (d, R), measured in MEASURED_THRESHOLDS.items():
        bound = threshold_count(d, R)
        assert measured >= bound, (d, R, measured, bound)
        assert measured == bound, (d, R, measured, bound)


def test_threshold_count_collapses_to_the_quadratic_at_r_equals_d_plus_two():
    """Round 3 fitted d^2 + 3d + 1 to three points and called it
    speculative. It is the m - d = 1 case of the derived count."""

    for d in range(2, 8):
        assert threshold_count(d, d + 2) == d * d + 3 * d + 1


def test_counting_alone_does_not_imply_rigidity_in_one_dimension():
    """The standing counterexample to the sufficiency step.

    At d = 1 the necessary count is satisfied with room to spare and the
    configuration is still never rigid, because the profile surface is a
    polyline of zero curvature. This is why Theorem 2b is [CONJECTURED]
    and not [PROVABLE] by counting.
    """

    for n, R in ((6, 3), (10, 4), (16, 6)):
        X, P = scene(n, R, seed=100 + n + 7 * R, d=1)
        theta = flatten(X, P)
        null, _ = nullity(theta, n, R, 1)
        flex, _, m = cloud_flex_space(theta, n, R, 1)
        counting_bound = 1 * (n + R) + flex.shape[1] - m * n
        assert counting_bound <= gauge_dimension(1)   # count permits rigidity
        assert null > gauge_dimension(1)              # reality refuses it


def test_round4_slope_three_and_its_intercept():
    """First slope-3 cell ever run, and the first under-observed cell
    whose intercept was predicted rather than only its slope."""

    rows = nullity_by_n(range(6, 13), R=4, d=6, seed_base=1800)
    assert [row["nullity"] for row in rows] == [48, 51, 54, 57, 60, 63, 66]
    for row in rows:
        assert row["nullity"] == 3 * row["n"] + 30


def test_round4_threshold_at_six_dimensions():
    """d = 6, R = 8: the derived count says 55, in a scene of 55 targets
    and 8 observers -- the largest in the series."""

    assert threshold_count(6, 8) == 55
    for n, expect_rigid in ((54, False), (55, True)):
        X, P = scene(n, 8, seed=2000 + 100 * 6 + 8 + n, d=6)
        null, _ = nullity(flatten(X, P), n, 8, 6)
        assert (null == gauge_dimension(6)) is expect_rigid, (n, null)


def test_round4_non_exact_divisions_round_up_and_are_attained():
    """P18 and P19 both had non-integer counts (39/2, 53/2). Attainment
    at the first integer above is a claim in its own right."""

    assert threshold_count(4, 7) == 20
    assert threshold_count(5, 8) == 27
    for n, R, d, expect_rigid in ((19, 7, 4, False), (20, 7, 4, True)):
        X, P = scene(n, R, seed=2000 + 100 * d + R + n, d=d)
        null, _ = nullity(flatten(X, P), n, R, d)
        assert (null == gauge_dimension(d)) is expect_rigid, (n, null)
