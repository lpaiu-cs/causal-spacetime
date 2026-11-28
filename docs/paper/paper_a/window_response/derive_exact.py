"""Reconstruct the manuscript's exact coefficients without simulation or fitting.

Optional symbolic audit: Python + SymPy. The experimental runtime does not
depend on SymPy. This reconstructs already frozen formulas after confirmation.
"""

import sympy as sp


def box(a, b):
    x, y, s = sp.symbols("x y s", positive=True)
    mx, my = (a - x) ** 2 / 12, (b - y) ** 2 / 12
    c0 = (x * x + y * y) / (2 * s)
    c1 = s * ((mx - my) / 2 + (x * x - y * y) / 24)
    c2 = -(s**3) * ((mx + my) / 24 + (x * x + y * y) / 1440)
    c11 = s**2 * (
        ((a - x) ** 4 + (b - y) ** 4) / 320
        - mx * my / 2
        + (x * x - y * y) * (mx - my) / 24
        + (x * x - y * y) ** 2 / 576
    )
    weight = (a - x) * (b - y) / (a * a * b * b)
    coefficients = []
    for integrand in ((1 - c0) ** 2 / 2, -(1 - c0) * c1, -(1 - c0) * c2 + c11 / 2):
        polynomial = sp.Poly(sp.expand(s**2 * weight * integrand), x, y, s)
        value = 0
        for (i, j, shifted_k), coefficient in polynomial.terms():
            h = sp.Rational(i + j + 2, 2)
            t = shifted_k - 2 + h
            beta = (
                sp.gamma(sp.Rational(i + 1, 2))
                * sp.gamma(sp.Rational(j + 1, 2))
                / sp.gamma(h)
            )
            value += coefficient * 4 * beta * 2**h / ((i + j + 2) * (t + 1) * (t + 2))
        coefficients.append(sp.simplify(value))
    return coefficients


def diamond():
    u, e = sp.symbols("u e", real=True)
    # e=A^2; discard e^2 only after forming the normalized moments.
    j = lambda s: s**2 / 3 * (1 + e * s**4 / 252)  # noqa: E731
    slice_weight = 2 * u * (1 - u) * (1 + e * (1 - u**4 - (1 - u) ** 4) / 180)
    j1 = j(sp.Integer(1))
    p = sp.integrate(sp.expand(slice_weight * (j(u) + j(1 - u))), (u, 0, 1)) / (
        6 * j1**2
    )
    second = (
        sp.integrate(
            sp.expand(
                slice_weight * ((j(u) ** 2 + j(1 - u) ** 2) / 15 + j(u) * j(1 - u) / 45)
            ),
            (u, 0, 1),
        )
        / j1**3
    )
    return (
        sp.series(p, e, 0, 2).removeO(),
        sp.series(second - p**2, e, 0, 2).removeO(),
        e,
    )


if __name__ == "__main__":
    r = sp.Rational
    plus = [
        sp.pi / 108 + r(1, 5184) - 128 * sp.sqrt(2) / 33075,
        sp.pi / 54 + r(1, 540) - 352 * sp.sqrt(2) / 25515,
        1613 * sp.pi / 48600 + r(16733, 2177280) - 1316008 * sp.sqrt(2) / 30405375,
    ]
    square = [
        sp.pi / 108 + r(1, 5184) - 64 * sp.sqrt(6) / 33075,
        0,
        257 * sp.pi / 48600 + r(73, 435456) - 137204 * sp.sqrt(6) / 30405375,
    ]
    for name, a, b, expected in (
        ("R+", sp.Integer(2), sp.Integer(6), plus),
        ("R-", sp.Integer(6), sp.Integer(2), [plus[0], -plus[1], plus[2]]),
        ("R0", sp.sqrt(12), sp.sqrt(12), square),
    ):
        actual = box(a, b)
        assert all(
            sp.simplify(v - w) == 0 for v, w in zip(actual, expected, strict=True)
        )
        print(name, [str(v) for v in actual])
    p, zeta, e = diamond()
    assert sp.simplify(p - (r(1, 10) - r(23, 75600) * e)) == 0
    assert sp.simplify(zeta - (r(1, 100) - r(211, 4158000) * e)) == 0
    print("D", p, zeta)
    print("All exact coefficient identities verified; no poset data read.")
