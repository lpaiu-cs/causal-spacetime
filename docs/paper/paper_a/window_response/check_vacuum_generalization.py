"""Post-study local vacuum-volume derivation checks; no sampled data are read.

The proof and its scope are in vacuum_generalization_derivation.md. This file
checks Jacobi coefficients, arbitrary-dimensional algebraic curvature channels,
and separate null-arrival calculations in actual, time-dependent Kasner metrics.
It is not a numerical survey or a machine-checked proof of differential geometry.
Requires the already installed SymPy; outputs a source-bound JSON record.
"""

import hashlib
import itertools
import json
import platform
import time
from pathlib import Path

import sympy as sp


def sphere_moment(powers):
    if any(power % 2 for power in powers):
        return sp.S.Zero
    degree = sum(powers) // 2
    return sp.prod(sp.factorial2(power - 1) for power in powers) / sp.prod(
        len(powers) + 2 * j for j in range(degree)
    )


def sphere_average(expression, directions):
    return sp.expand(
        sum(
            coefficient * sphere_moment(powers)
            for powers, coefficient in sp.Poly(
                sp.expand(expression), directions
            ).terms()
        )
    )


def relative_coefficient(dimension, electric, magnetic, spatial):
    d = sp.Integer(dimension)
    numerator = (
        (7 * d**3 + 58 * d**2 + 146 * d + 116) * electric
        + 3 * (d + 2) * (d + 6) * magnetic
        - sp.Rational(3, 2) * (d + 2) * spatial
    )
    return sp.factor(d * numerator / (180 * (d + 1) * (d + 2) * (d + 3) * (d + 4)))


def jacobi_and_action_coefficients():
    s, eps, k0, k1, k2 = sp.symbols("s eps K0 K1 K2")
    # At this order the only curvature product is K0*K0. No commutativity
    # assumption about different curvature matrices is used.
    j = s - eps**2 * s**3 * k0 / 6 - eps**3 * s**4 * k1 / 12
    j += eps**4 * s**5 * (k0**2 / 120 - k2 / 40)
    residual = sp.expand(
        sp.diff(j, s, 2) + (eps**2 * k0 + eps**3 * s * k1 + eps**4 * s**2 * k2 / 2) * j
    )
    assert all(residual.coeff(eps, order) == 0 for order in range(5))
    metric = sp.series(j.subs(s, 1) ** 2, eps, 0, 5).removeO()
    assert metric.coeff(eps, 4) == 2 * k0**2 / 45 - k2 / 20
    assert sp.Rational(1, 120) - sp.Rational(1, 72) == -sp.Rational(1, 180)
    aa, av, vv = sp.symbols("AA AV VV")
    # Metric along the straight segment, minus the stationary path correction.
    action = sp.integrate(
        sp.Rational(2, 45) * (aa + 2 * s * av + s**2 * vv)
        - sp.Rational(1, 9) * (2 * s - 1) ** 2 * vv,
        (s, 0, 1),
    )
    assert sp.expand(action - (2 * aa + 2 * av - vv) / 45) == 0
    return {
        "metric_series": str(metric),
        "curvature_squared_action": str(action),
        "vacuum_measure_K_squared": "-1/180",
    }


def algebraic_weyl(dimension, e_weight, h_weight, w_weight):
    """Independent E/H/spatial-Weyl channels, represented on two-forms."""
    m = dimension - 1
    electric = sp.diag(1, 2, -3, *([0] * (m - 3))) * e_weight
    b = sp.Matrix([[0, 1, 0], [1, 0, 2], [0, 2, 0]]) * h_weight
    pairs = list(itertools.combinations(range(dimension), 2))
    lookup = {pair: index for index, pair in enumerate(pairs)}
    matrix = sp.zeros(len(pairs))
    spatial_sections = {
        (1, 2): 1,
        (1, 3): 2,
        (1, 4): -3,
        (2, 3): -3,
        (2, 4): 2,
        (3, 4): 1,
    }
    for p, (i, j) in enumerate(pairs):
        for q, (k, a) in enumerate(pairs):
            if i == 0 and k == 0:
                value = electric[j - 1, a - 1]
            elif i == 0:
                value = (
                    sum(sp.LeviCivita(k - 1, a - 1, c) * b[j - 1, c] for c in range(3))
                    if max(j, k, a) <= 3
                    else 0
                )
            elif k == 0:
                value = matrix[q, p]
            else:
                value = (
                    (i == k) * electric[j - 1, a - 1]
                    - (i == a) * electric[j - 1, k - 1]
                    - (j == k) * electric[i - 1, a - 1]
                    + (j == a) * electric[i - 1, k - 1]
                ) / (m - 2)
                if (i, j) == (k, a) and m >= 4:
                    value += w_weight * spatial_sections.get((i, j), 0)
            matrix[p, q] = value

    def component(i, j, k, a):
        if i == j or k == a:
            return sp.S.Zero
        sign = (1 if i < j else -1) * (1 if k < a else -1)
        return (
            sign * matrix[lookup[tuple(sorted((i, j)))], lookup[tuple(sorted((k, a)))]]
        )

    for i, j, k, a in itertools.product(range(dimension), repeat=4):
        assert (
            component(i, j, k, a) + component(i, k, a, j) + component(i, a, j, k) == 0
        )
    for i, j in itertools.product(range(dimension), repeat=2):
        assert (
            sum((-1 if k == 0 else 1) * component(k, i, k, j) for k in range(dimension))
            == 0
        )
    return component


def check_boundary_and_volume(dimension, weights, kappa=sp.S.Zero):
    weyl = algebraic_weyl(dimension, *weights)
    m = dimension - 1
    t, r = sp.symbols("t r")
    directions = sp.symbols(f"n1:{dimension}")
    n = sp.Matrix(directions)
    delta = (t - 1, *(r * entry for entry in n))
    signs = (-1, *([1] * m))

    def c(i, j, k, a):
        constant = signs[i] * signs[j] * ((i == k and j == a) - (i == a and j == k))
        return weyl(i, j, k, a) + kappa * constant

    # A_a=R_{a q Delta q}, V_a=R_{a Delta Delta q}, with q=(1,0,...).
    avec = [
        sum(c(i, 0, j, 0) * delta[j] for j in range(dimension))
        for i in range(dimension)
    ]
    vvec = [
        sum(
            c(i, j, k, 0) * delta[j] * delta[k]
            for j, k in itertools.product(range(dimension), repeat=2)
        )
        for i in range(dimension)
    ]
    action2 = sp.expand(
        sum(
            signs[i] * (2 * avec[i] ** 2 + 2 * avec[i] * vvec[i] - vvec[i] ** 2)
            for i in range(dimension)
        )
        / 45
    )
    electric = sp.Matrix(m, m, lambda i, j: c(0, i + 1, 0, j + 1))
    e = (n.T * electric * n)[0]
    e2n = (n.T * electric * electric * n)[0]
    h = sp.Matrix(
        [
            sum(
                c(0, j + 1, i + 1, k + 1) * n[j] * n[k]
                for j, k in itertools.product(range(m), repeat=2)
            )
            for i in range(m)
        ]
    )
    hh, eh = (h.T * h)[0], (n.T * electric * h)[0]
    expected_action = (
        r**2
        * ((4 * t - t**2 - 1) * e2n + 2 * r * (t - 2) * eh + r**2 * (e**2 - hh))
        / 45
    )
    assert sp.expand(action2 - expected_action) == 0
    r0 = 1 - t
    b1 = e / 6
    b2 = sp.cancel(e**2 / 24 - action2.subs(r, r0) / (2 * r0**2))
    boundary = r0**m * sphere_average(b2 + sp.Rational(m - 1, 2) * b1**2, directions)
    point = (t, *(r * entry for entry in n))
    q = sum(
        signs[i]
        * signs[j]
        * sum(
            c(i, k, j, a) * point[k] * point[a]
            for k, a in itertools.product(range(dimension), repeat=2)
        )
        ** 2
        for i, j in itertools.product(range(dimension), repeat=2)
    )
    average_q = sphere_average(q, directions)
    e2 = sum(
        weyl(0, i, 0, j) ** 2
        for i, j in itertools.product(range(1, dimension), repeat=2)
    )
    h2 = sum(
        weyl(0, i, j, k) ** 2
        for i, j, k in itertools.product(range(1, dimension), repeat=3)
    )
    d2 = sum(
        weyl(i, j, k, a) ** 2
        for i, j, k, a in itertools.product(range(1, dimension), repeat=4)
    )
    assert (
        sp.expand(
            average_q
            - e2 * t**4
            - 3 * h2 * t**2 * r**2 / m
            - (3 * e2 - 3 * h2 + sp.Rational(3, 2) * d2) * r**4 / (m * (m + 2))
            - m * kappa**2 * (r * r - t * t) ** 2
        )
        == 0
    )
    measure2 = -m * kappa * (r * r - t * t) / 6
    average_b1 = sphere_average(b1, directions)
    second = (
        dimension
        * m
        * sp.integrate(
            r0**m * average_b1 + sp.integrate(r ** (m - 1) * measure2, (r, 0, r0)),
            (t, 0, 1),
        )
    )
    assert sp.cancel(second + dimension * m * kappa / (3 * (dimension + 2))) == 0
    # The linear measure also changes when its upper radial boundary moves.
    boundary += r0**m * average_b1 * measure2.subs(r, r0)
    measure4 = m * m * kappa**2 * (r * r - t * t) ** 2 / 72 - average_q / 180
    measure = sp.integrate(r ** (m - 1) * measure4, (r, 0, r0))
    computed = sp.factor(dimension * m * sp.integrate(boundary + measure, (t, 0, 1)))
    expected = relative_coefficient(dimension, e2, h2, d2)
    expected += (
        sp.Rational(dimension * m * (5 * dimension + 2), 90 * (dimension + 4))
        * kappa**2
    )
    assert sp.cancel(computed - expected) == 0
    return {
        "dimension": dimension,
        "channel_weights": weights,
        "constant_sectional_curvature": str(kappa),
        "E2_H2_D2": list(map(str, (e2, h2, d2))),
        "relative_l4": str(computed),
        "relative_l2": str(sp.factor(second)),
    }


def kasner_arrival_volume(exponents):
    """Solve the null eikonal in original Kasner coordinates, not RNC.

    tau=1 at the midpoint, tips tau=1+-l, sqrt(-g)=tau.
    No Weyl-volume coefficient is used until the final comparison.
    """
    assert sum(exponents) == sum(p * p for p in exponents) == 1
    r, z1, z2, ell, rho = sp.symbols("r z1 z2 ell rho")
    z = (z1, z2, 1 - z1 - z2)
    f = {1: sp.S.One}
    for order in range(2, 6):
        front = 1 + sum(r**k * value for k, value in f.items())
        residual = 0
        for i, p in enumerate(exponents):
            gradient = 0
            for k, value in f.items():
                radial = z1 * sp.diff(value, z1) + z2 * sp.diff(value, z2)
                derivative = sp.diff(value, (z1, z2)[i]) if i < 2 else 0
                gradient += r ** (k - 1) * (k * value + 2 * derivative - 2 * radial)
            factor = sp.series(front ** (-2 * p), r, 0, order).removeO()
            gradient_square = sp.Poly(sp.expand(gradient**2), r)
            residual += z[i] * sum(
                factor.coeff(r, j)
                * gradient_square.coeff_monomial(r ** (order - 1 - j))
                for j in range(order)
            )
        f[order] = sp.factor(-residual / (2 * order))
    # Along each principal axis, the null equation has an independent closed
    # solution: T=[1+(1-p)*r]**(1/(1-p)), or exp(r) if p=1.
    for axis, p in enumerate(exponents):
        closed = sp.exp(r) if p == 1 else (1 + (1 - p) * r) ** (1 / (1 - p))
        closed = sp.series(closed, r, 0, 6).removeO()
        point = {z1: int(axis == 0), z2: int(axis == 1)}
        assert all(
            sp.simplify(f[k].subs(point) - closed.coeff(r, k)) == 0 for k in range(1, 6)
        )
    if exponents == (1, 0, 0):
        # Exact flat Milne light cone: cosh(r*n1) + sqrt(sinh(r*n1)^2+r^2*(1-z1)).
        exact_flat = {
            1: 1,
            2: z1 / 2,
            3: z1**2 / 6,
            4: z1**2 / 24,
            5: z1**3 / 45 - z1**4 / 72,
        }
        assert all(sp.expand(f[k] - exact_flat[k]) == 0 for k in f)
    tip = 1 - ell
    factors = [sp.series(tip ** (2 * (p - 1)), ell, 0, 5).removeO() for p in exponents]
    radial_factor = sp.expand(sum(z[i] * factors[i] for i in range(3)))
    lower = tip
    for k in range(1, 6):
        for (j1, j2), coefficient in sp.Poly(f[k], z1, z2).terms():
            exponent = 1 + 2 * j1 * (exponents[0] - 1) + 2 * j2 * (exponents[1] - 1)
            term = sp.series(tip**exponent, ell, 0, 6 - k).removeO()
            term *= sp.series(
                radial_factor ** (sp.Rational(k, 2) - j1 - j2), ell, 0, 6 - k
            ).removeO()
            lower += coefficient * z1**j1 * z2**j2 * rho**k * ell**k * sp.expand(term)
    coefficients = sp.Poly(sp.expand(lower), ell)
    a = {k: coefficients.coeff_monomial(ell**k) for k in range(1, 6)}
    assert a[1] == rho - 1
    height2 = -2 * (a[3] + a[1] * a[2])
    height4 = -2 * (a[5] + a[2] * a[3] + a[1] * a[4])
    edge2 = -a[3].subs(rho, 1)

    def average(expression):
        return sp.factor(
            sum(
                coefficient * sphere_moment((2 * i, 2 * j, 0))
                for (i, j), coefficient in sp.Poly(
                    sp.expand(expression), z1, z2
                ).terms()
            )
        )

    coefficient2 = 6 * average(sp.integrate(rho**2 * height2, (rho, 0, 1)))
    coefficient4 = 6 * average(sp.integrate(rho**2 * height4, (rho, 0, 1)) + edge2**2)
    e2 = sum((p * (1 - p)) ** 2 for p in exponents)
    expected = relative_coefficient(4, e2, 0, 4 * e2)
    assert coefficient2 == 0 and coefficient4 == expected
    return {
        "Kasner_exponents": list(map(str, exponents)),
        "E2": str(e2),
        "direct_relative_l2": str(coefficient2),
        "direct_relative_l4": str(coefficient4),
        "predicted_relative_l4": str(expected),
        "method": "null eikonal and moving edge in original Kasner coordinates",
    }


def varying_plane_wave():
    """Exact Jacobi propagator expansion for A(u)=a0+a1*u+a2*u^2.

    u=L*x with tips x=+-1 and proper half-length l=sqrt(2)*L.
    The two transverse equations have opposite signs, so Ricci is zero.
    """
    x, y, s, length, a0, a1, a2 = sp.symbols("x y s L a0 a1 a2")
    first = [sp.integrate((x - s) * s**k * (s - y), (s, y, x)) for k in range(3)]
    repeated = sp.integrate((x - s) * first[0].subs(x, s), (s, y, x))
    propagators = [
        (x - y)
        + sign
        * (
            a0 * length**2 * first[0]
            + a1 * length**3 * first[1]
            + a2 * length**4 * first[2]
        )
        + a0 * a0 * length**4 * repeated
        for sign in (-1, 1)
    ]
    product = sp.expand(propagators[0] * propagators[1])
    assert product.coeff(length, 2) == product.coeff(length, 3) == 0
    quartic = sp.factor(product.coeff(length, 4))
    assert quartic == -a0 * a0 * (x - y) ** 6 / 90
    # Wronskian identity gives the diamond slice from the three propagators.
    correction = (1 - x * x) / 2 * a0 * a0 * (16 - (x + 1) ** 4 - (1 - x) ** 4) / 180
    coefficient_L4 = sp.Rational(3, 2) * sp.integrate(correction, (x, -1, 1))
    coefficient_l4 = sp.factor(coefficient_L4 / 4)
    assert coefficient_l4 == a0 * a0 / 63
    return {
        "profile": "A(u)=a0+a1*u+a2*u^2",
        "relative_l4": str(coefficient_l4),
        "linear_curvature_derivative_terms": "cancel in the propagator product",
    }


def einstein_constant_curvature():
    # Exact static-patch radius tanh(sqrt(kappa)*s)/sqrt(kappa).
    # Analytic continuation kappa<0 gives the AdS result as well.
    z, m = sp.symbols("z m", positive=True)
    radius = sp.series(sp.tanh(z) / z, z, 0, 6).removeO()
    series = sp.series(radius**m, z, 0, 6).removeO().expand()
    second = sp.factor(series.coeff(z, 2) * (m + 1) / (m + 3))
    fourth = sp.factor(series.coeff(z, 4) * (m + 1) / (m + 5))
    assert sp.cancel(second + m * (m + 1) / (3 * (m + 3))) == 0
    assert sp.cancel(fourth - m * (m + 1) * (5 * m + 7) / (90 * (m + 5))) == 0
    return {
        "spatial_dimension": "m=d-1",
        "kappa_l2": str(second),
        "kappa2_l4": str(fourth),
        "four_dimensional": {
            "kappa_l2": str(second.subs(m, 3)),
            "kappa2_l4": str(fourth.subs(m, 3)),
        },
    }


if __name__ == "__main__":
    started = time.perf_counter()
    tensor_rows = []
    for dimension in (4, 5, 6):
        channels = [(1, 0, 0), (0, 1, 0), (1, 2, 3 if dimension >= 5 else 0)]
        if dimension >= 5:
            channels.append((0, 0, 1))
        for weights in channels:
            tensor_rows.append(check_boundary_and_volume(dimension, weights))
    for dimension in (4, 5):
        tensor_rows.append(
            check_boundary_and_volume(dimension, (1, 2, 0), sp.Symbol("kappa"))
        )
    kasner = [
        kasner_arrival_volume(tuple(map(sp.Rational, exponents)))
        for exponents in [
            ("1", "0", "0"),
            ("-1/3", "2/3", "2/3"),
            ("-2/7", "3/7", "6/7"),
        ]
    ]
    record = {
        "kind": "Post-study deterministic derivation audit; not peer review",
        "date": "2026-09-22",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "source_sha256_lf": hashlib.sha256(
            Path(__file__).read_bytes().replace(b"\r\n", b"\n")
        ).hexdigest(),
        "scope": (
            "Local small diamonds, fixed central timelike unit vector, "
            "smooth Ricci-flat neighborhood; Einstein extension stated separately"
        ),
        "jacobi_and_action": jacobi_and_action_coefficients(),
        "algebraic_tensor_checks": tensor_rows,
        "independent_Kasner": kasner,
        "constant_curvature": einstein_constant_curvature(),
        "varying_profile_plane_wave": varying_plane_wave(),
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }
    print(json.dumps(record, indent=2))
