"""Independent identities for the new response instrument, not fitted results."""

import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(
    0, str(Path(__file__).resolve().parents[1] / "experiments" / "positive_control")
)
import window_response as wr  # noqa: E402
from p14_plane_wave import PlaneWaveGeometry, Slab, causal_relation  # noqa: E402
from scipy.integrate import quad  # noqa: E402


@pytest.mark.parametrize("region", wr.BOXES)
def test_coefficients_against_unreduced_integral(region):
    for order in (24, 40):
        np.testing.assert_allclose(
            wr.coefficient_quadrature(region, order),
            wr.box_coefficients(region),
            atol=2e-13,
            rtol=0,
        )


@pytest.mark.parametrize("amplitude", [0, 1 / 64, 1 / 4, 1, 2])
def test_order_matrix_against_certified_scalar_predicate(amplitude):
    rng = np.random.default_rng([20260921, 6101, 1])
    points = wr.sample_points("R+", amplitude, 24, rng)
    points = points[np.argsort(points[:, 0])]
    matrix, _ = wr.order_matrix(points, amplitude)
    geometry = PlaneWaveGeometry(Slab(1, 8, 16, 16), amplitude**0.5)
    for i, j in itertools.combinations(range(len(points)), 2):
        reference = causal_relation(geometry, points[i], points[j])
        assert bool(matrix[i, j]) == reference.related
    brute = sum(
        matrix[i, j] and matrix[j, k]
        for i, j, k in itertools.combinations(range(len(points)), 3)
    )
    assert brute == np.dot(matrix.sum(axis=0), matrix.sum(axis=1))
    swapped = wr.order_matrix(points[:, [0, 1, 3, 2]], -amplitude)[0]
    np.testing.assert_array_equal(matrix, swapped)
    permuted = points[rng.permutation(len(points))]
    assert wr.statistics(permuted, amplitude) == wr.statistics(points, amplitude)


def test_clipped_lines_against_independent_adaptive_integral():
    lo, hi = np.array([-0.5] * 4), np.array([0.5] * 4)
    a, b = np.array([1, -2, 0.5, 10.0]), np.array([0.1, 0.4, -2, 0.2])
    c, d = np.array([2, 1, -3, -5.0]), np.array([0.5, -0.2, 3, 0.7])
    p, both = wr.integrate_v(lo, hi, a, b, c, d)
    for i in range(4):
        roots = [-b[i] / a[i], (1 - b[i]) / a[i], -d[i] / c[i], (1 - d[i]) / c[i]]
        points = [x for x in roots if lo[i] < x < hi[i]]

        def f(x, i=i):
            return np.clip(a[i] * x + b[i], 0, 1)

        def g(x, i=i):
            return np.clip(c[i] * x + d[i], 0, 1)

        expected_p = quad(lambda x: (f(x) + g(x)) / 2, lo[i], hi[i], points=points)[0]
        expected_both = quad(lambda x: f(x) * g(x), lo[i], hi[i], points=points)[0]
        assert abs(p[i] - expected_p) < 1e-12
        assert abs(both[i] - expected_both) < 1e-12


def test_flat_diamond_volume_and_membership():
    import math

    assert abs(wr.diamond_cdf(0)[1] - math.pi / 6) < 1e-13
    for amplitude in (0, 1, 2):
        rng = np.random.default_rng([20260921, 6101, 2])
        points = wr.sample_points("D", amplitude, 128, rng)
        geometry = PlaneWaveGeometry(Slab(1, 8, 16, 16), amplitude**0.5)
        for point in points:
            assert causal_relation(geometry, np.array([0, 0.5, 0, 0]), point).related
            assert causal_relation(geometry, point, np.array([1, -0.5, 0, 0])).related


def test_invalid_inputs_and_degenerate_small_poset():
    with pytest.raises(ValueError):
        wr.action(0.5, 0, 1, 0, 1, 3)
    with pytest.raises(ValueError):
        wr.sample_points("D", 1, 1, np.random.default_rng(0))
    # Same u, separated x is spacelike; same-u null generator is handled by legacy.
    points = np.array([[0.2, 0.2, 0, 0], [0.2, 0.1, 1, 0], [0.2, 0, 0, 0]])
    matrix, count = wr.order_matrix(points, 1)
    assert count == 3 and matrix[0, 2] and not matrix[0, 1] and not matrix[1, 2]


def test_diamond_conditional_degree_integrals():
    import math

    flat = wr.diamond_moments(0)
    np.testing.assert_allclose(
        [flat["p"], flat["zeta"], flat["volume"]],
        [0.1, 0.01, math.pi / 6],
        rtol=0,
        atol=1e-14,
    )
    for amplitude in (0.015625, 0.25, 1, 2):
        coarse = wr.diamond_moments(amplitude, 32)
        fine = wr.diamond_moments(amplitude, 64)
        np.testing.assert_allclose(
            list(coarse.values()), list(fine.values()), rtol=0, atol=1e-13
        )
        assert abs(fine["volume"] - wr.diamond_cdf(amplitude)[1]) < 1e-11
    small = wr.diamond_moments(1 / 64)
    assert abs((small["p"] - 0.1) / (1 / 64) ** 2 + 23 / 75600) < 1e-8
