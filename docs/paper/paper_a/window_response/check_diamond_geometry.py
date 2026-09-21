"""Post-study equation audit, 22 September 2026. No sprinkling data are read.

Requires SymPy (and its mpmath dependency). Recomputes the interval-volume
coefficient in Brinkmann and Riemann normal coordinates, contracts curvature
from the metric, and compares Wang arXiv:1904.01034v2, equation (79).
This is a deterministic mathematical check, not a new confirmatory campaign.
"""

import hashlib
import itertools
import json
import platform
from pathlib import Path

import mpmath as mp
import sympy as sp


def curvature_invariants():
    u, v, x, y, a = sp.symbols("u v x y A", real=True)
    coordinates = (u, v, x, y)
    g = sp.Matrix(
        [[a * (x * x - y * y), 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    )
    inverse = g.inv()
    gamma = sp.MutableDenseNDimArray.zeros(4, 4, 4)
    for i, j, k in itertools.product(range(4), repeat=3):
        gamma[i, j, k] = sp.simplify(
            sum(
                inverse[i, m]
                * (
                    sp.diff(g[m, k], coordinates[j])
                    + sp.diff(g[m, j], coordinates[k])
                    - sp.diff(g[j, k], coordinates[m])
                )
                / 2
                for m in range(4)
            )
        )
    mixed = sp.MutableDenseNDimArray.zeros(4, 4, 4, 4)
    lower = sp.MutableDenseNDimArray.zeros(4, 4, 4, 4)
    for i, j, k, m in itertools.product(range(4), repeat=4):
        mixed[i, j, k, m] = sp.simplify(
            sp.diff(gamma[i, j, m], coordinates[k])
            - sp.diff(gamma[i, j, k], coordinates[m])
            + sum(
                gamma[i, k, b] * gamma[b, j, m] - gamma[i, m, b] * gamma[b, j, k]
                for b in range(4)
            )
        )
    for i, j, k, m in itertools.product(range(4), repeat=4):
        lower[i, j, k, m] = sp.simplify(
            sum(g[i, b] * mixed[b, j, k, m] for b in range(4))
        )
    ricci = sp.Matrix(4, 4, lambda j, m: sum(mixed[i, j, i, m] for i in range(4)))
    assert ricci == sp.zeros(4)
    assert lower[0, 2, 0, 2] == -a and lower[0, 3, 0, 3] == a
    # This constant-profile plane wave is locally symmetric: no derivative
    # terms are being dropped when comparing the curvature-squared coefficient.
    for b, i, j, k, m in itertools.product(range(4), repeat=5):
        indices = [i, j, k, m]
        covariant = sp.diff(lower[i, j, k, m], coordinates[b])
        for position in range(4):
            for c in range(4):
                shifted = indices.copy()
                shifted[position] = c
                covariant -= gamma[c, b, indices[position]] * lower[tuple(shifted)]
        assert sp.simplify(covariant) == 0
    root = sp.sqrt(2)
    frame = (
        (1 / root, -1 / root, 0, 0),
        (0, 0, 1, 0),
        (0, 0, 0, 1),
        (1 / root, 1 / root, 0, 0),
    )
    assert sp.Matrix(frame) * g.subs({x: 0, y: 0}) * sp.Matrix(frame).T == sp.diag(
        -1, 1, 1, 1
    )
    nonzero = [
        (indices, lower[indices])
        for indices in itertools.product(range(4), repeat=4)
        if lower[indices] != 0
    ]

    def component(i, j, k, m):
        return sp.simplify(
            sum(
                value * frame[i][b] * frame[j][c] * frame[k][d] * frame[m][e]
                for (b, c, d, e), value in nonzero
            )
        )

    electric = sum(
        component(0, i, 0, j) ** 2 for i, j in itertools.product(range(1, 4), repeat=2)
    )
    magnetic = sum(
        component(0, i, j, k) ** 2
        for i, j, k in itertools.product(range(1, 4), repeat=3)
    )
    spatial = sum(
        component(i, j, k, m) ** 2
        for i, j, k, m in itertools.product(range(1, 4), repeat=4)
    )
    assert electric == a**2 / 2 and magnetic == a**2 and spatial == 2 * a**2
    d = sp.Integer(4)
    # Eq. (79), divided by its flat volume; tip half-length l=1/sqrt(2).
    wang = sp.simplify(
        d
        * (d - 1)
        / 8
        * (
            (7 * d**3 + 58 * d**2 + 146 * d + 108) * electric
            + 6 * (d + 2) * (d + 6) * magnetic
            - sp.Rational(3, 2) * (d + 2) * spatial
        )
        / (90 * (d - 1) * (d + 1) * (d + 2) * (d + 3) * (d + 4))
        / a**2
    )
    assert wang == sp.Rational(43, 9450)
    return {
        "E_squared": str(electric),
        "H_squared": str(magnetic),
        "D_squared": str(spatial),
        "Wang_relative_A2": str(wang),
        "Wang_74_76_reintegrated_relative_A2": str(reintegrate_wang(component, a)),
    }


def reintegrate_wang(component, amplitude):
    """Directly specialize Eqs. (74),(76), separately from final Eq. (79)."""
    t, radius = sp.symbols("t r", real=True)
    directions = sp.symbols("nx ny nz", real=True)

    def angular_average(expression):
        value = 0
        for powers, coefficient in sp.Poly(sp.expand(expression), directions).terms():
            if any(power % 2 for power in powers):
                continue
            moment = sp.prod(sp.factorial2(power - 1) for power in powers)
            value += coefficient * moment / sp.factorial2(sum(powers) + 1)
        return sp.simplify(value)

    electric = sp.Matrix(3, 3, lambda i, j: component(0, i + 1, 0, j + 1))
    n = sp.Matrix(directions)
    enn = (n.T * electric * n)[0]
    e2nn = (n.T * electric * electric * n)[0]
    hh = sum(
        sum(
            directions[i] * directions[j] * component(0, i + 1, m + 1, j + 1)
            for i, j in itertools.product(range(3), repeat=2)
        )
        ** 2
        for m in range(3)
    )
    b1 = enn / 6
    b2 = enn**2 / 24 - e2nn * t / 45 + (1 - t) ** 2 * (e2nn + hh - enn**2) / 90
    boundary = (1 - t) ** 3 * angular_average(b2 + b1**2)
    point = (t, *(radius * direction for direction in directions))
    signs = (-1, 1, 1, 1)
    contraction = 0
    for i, j in itertools.product(range(4), repeat=2):
        tidal = sum(
            component(i, k, j, m) * point[k] * point[m]
            for k, m in itertools.product(range(4), repeat=2)
        )
        contraction += signs[i] * signs[j] * tidal**2
    assert (
        sp.simplify(contraction - amplitude**2 * (t + radius * directions[2]) ** 4 / 2)
        == 0
    )
    measure = -sp.integrate(
        radius**2 * angular_average(contraction) / 180, (radius, 0, 1 - t)
    )
    result = sp.simplify(3 * sp.integrate(boundary + measure, (t, 0, 1)) / amplitude**2)
    assert result == sp.Rational(1, 252)
    return result


def independent_volume_expansions():
    u, w = sp.symbols("u w", real=True)
    # Direct expansion of the exact Brinkmann slice, without importing J.
    product = sp.series(sp.sin(w) * sp.sinh(w) / w**2, w, 0, 8).removeO()
    assert sp.expand(product) == 1 - w**4 / 90
    coefficient = sp.integrate(
        2 * u * (1 - u) * (1 - u**4 - (1 - u) ** 4) / 180, (u, 0, 1)
    )
    assert 3 * coefficient == sp.Rational(1, 252)

    # Normal coordinates at the midpoint. Set the tip half-length l=1
    # while integrating; dimensional scaling restores l^4 afterwards.
    t, z, x, y, a = sp.symbols("t z x y A", real=True)
    s = (t + z) / sp.sqrt(2)
    remaining = 1 / sp.sqrt(2) - s
    bx = x * (1 + a * s**2 / 6 + a**2 * s**4 / 120)
    by = y * (1 - a * s**2 / 6 + a**2 * s**4 / 120)
    bv = (
        (-t + z) / sp.sqrt(2)
        - a * s * (x * x - y * y) / 3
        - a**2 * s**3 * (x * x + y * y) / 15
    )
    null = sp.series(
        2 * remaining * (-1 / sp.sqrt(2) - bv)
        + bx**2
        + by**2
        + a * remaining**2 * (bx**2 - by**2) / 3
        - a**2 * remaining**4 * (bx**2 + by**2) / 45,
        a,
        0,
        3,
    ).removeO()
    expected = -((1 - t) ** 2) + x * x + y * y + z * z + a * (x * x - y * y) / 6
    expected += a**2 * (x * x + y * y) * (-((t + z) ** 2) + 4 * (t + z) - 1) / 180
    assert sp.simplify(null - expected) == 0
    jacobian = sp.Matrix([s, bv, bx, by]).jacobian([t, z, x, y]).det()
    jacobian = sp.series(jacobian, a, 0, 3).removeO()
    assert sp.simplify(jacobian - (1 - a**2 * (t + z) ** 4 / 360)) == 0

    # Solve the null equation for radius r=r0*(1+A*b1+A^2*b2).
    # n_z=k, n_x^2-n_y^2=difference, n_x^2+n_y^2=1-k^2.
    k, difference, r = sp.symbols("k difference r", real=True)
    r0 = 1 - t
    b1 = -difference / 12
    b2 = (
        difference**2 / 96
        + (1 - k * k) * ((t + r0 * k) ** 2 - 4 * (t + r0 * k) + 1) / 360
    )
    radial_null = -((1 - t) ** 2) + r * r + a * r * r * difference / 6
    radial_null += (
        a**2 * r * r * (1 - k * k) * (-((t + r * k) ** 2) + 4 * (t + r * k) - 1) / 180
    )
    residual = sp.series(
        radial_null.subs(r, r0 * (1 + a * b1 + a * a * b2)), a, 0, 3
    ).removeO()
    assert sp.simplify(residual) == 0
    # Angular averages: <difference^2>=4/15, <1-k^2>=2/3,
    # <(1-k^2)k^2>=2/15. Odd terms vanish in the angular integral.
    average_b2 = (
        sp.Rational(1, 360)
        + (sp.Rational(2, 3) * (t * t - 4 * t + 1) + sp.Rational(2, 15) * r0 * r0) / 360
    )
    boundary = r0**3 * (average_b2 + sp.Rational(1, 540))
    measure = -(t**4 * r0**3 / 3 + 2 * t * t * r0**5 / 5 + r0**7 / 35) / 360
    # Upper and lower cones give 8*pi; V_flat=2*pi*l^4/3.
    normal_coefficient = sp.integrate(boundary + measure, (t, 0, 1)) * 12
    assert normal_coefficient == sp.Rational(1, 63)
    fixed_span_coefficient = normal_coefficient / 4  # l=1/sqrt(2)
    assert fixed_span_coefficient == 3 * coefficient
    return {
        "Brinkmann_relative_A2": str(3 * coefficient),
        "RNC_relative_A2": str(fixed_span_coefficient),
        "Wang_minus_exact": str(sp.Rational(43, 9450) - fixed_span_coefficient),
    }


def high_precision_volume():
    mp.mp.dps = 60
    rows = []
    for amplitude in ("0.1", "0.05", "0.025", "0.0125"):
        a = mp.mpf(amplitude)
        w = mp.sqrt(a)

        def slice_volume(u, w=w):
            return (
                2
                / w
                * mp.sqrt(
                    mp.sin(w * u)
                    * mp.sinh(w * u)
                    * mp.sin(w * (1 - u))
                    * mp.sinh(w * (1 - u))
                    / (mp.sin(w) * mp.sinh(w))
                )
            )

        scaled = (3 * mp.quad(slice_volume, [0, mp.mpf("0.5"), 1]) - 1) / (a * a)
        assert abs(scaled - mp.mpf(1) / 252) < 1e-6
        rows.append({"A": amplitude, "(V/V0-1)/A2": mp.nstr(scaled, 24)})
    assert abs(scaled - mp.mpf(1) / 252) < 1e-8
    return rows


def mean_inverted_dimension():
    d = sp.symbols("d", positive=True)
    flat = sp.gamma(d + 1) * sp.gamma(d / 2) / (2 * sp.gamma(3 * d / 2))
    derivative = sp.simplify(sp.diff(flat, d).subs(d, 4))
    shift = sp.simplify(-sp.Rational(23, 75600) / derivative)
    assert derivative == -sp.Rational(101, 1200)
    assert shift == sp.Rational(23, 6363)
    return {
        "flat_p_prime_at_4": str(derivative),
        "dimension_A2_coefficient": str(shift),
        "scope": "inverse of population mean; not mean of the finite-N inverse",
    }


if __name__ == "__main__":
    print(
        json.dumps(
            {
                "kind": "Post-study deterministic audit; no new confirmation data",
                "date": "2026-09-22",
                "python": platform.python_version(),
                "sympy": sp.__version__,
                "mpmath": mp.__version__,
                "source_sha256_lf": hashlib.sha256(
                    Path(__file__).read_bytes().replace(b"\r\n", b"\n")
                ).hexdigest(),
                "curvature": curvature_invariants(),
                "volume": independent_volume_expansions(),
                "high_precision": high_precision_volume(),
                "dimension": mean_inverted_dimension(),
                "literature_status": (
                    "Disagreement with arXiv:1904.01034v2 Eq. (79); "
                    "no general correction derived"
                ),
            },
            indent=2,
        )
    )
