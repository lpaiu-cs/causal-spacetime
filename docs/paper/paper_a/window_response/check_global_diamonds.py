"""Post-study exact global examples; no frozen inputs or sampled posets are read.

See global_diamond_exact_cases.md for proofs and scope. Floating-point quadrature
and polygon refinement below are diagnostics, not certified interval enclosures.
The Kasner bounds are separate analytic bounds, not errors on the quadrature.
"""

import hashlib
import json
import platform
import time
from fractions import Fraction
from math import isqrt
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp


def cylinder_checks():
    ell, a, radius, y = sp.symbols("ell a r y", positive=True)
    small = 2 * sp.pi * ell**4 / 3
    ball_large = 2 * sp.pi * a * radius**2 - 2 * sp.pi * a**3 / 3
    by_slices = small.subs(ell, a) + 2 * sp.integrate(ball_large, (radius, a, ell))
    # Independent time-integrated optical formula, after the transverse integral.
    by_optical_distance = (
        8 * sp.pi * sp.integrate(ell**3 / 6 - ell * y**2 / 2 + y**3 / 3, (y, 0, a))
    )
    assert sp.expand(by_slices - by_optical_distance) == 0
    assert sp.simplify((by_slices - small).subs(ell, a)) == 0
    assert sp.simplify(sp.diff(by_slices - small, ell).subs(ell, a)) == 0
    circumference = sp.Symbol("C", positive=True)
    large = sp.expand(by_slices.subs(a, circumference / 2))
    ratios = []
    for scale in (sp.Rational(1, 4), sp.Rational(1, 2), 1, 2, 5):
        volume = small if scale <= sp.Rational(1, 2) else large
        ratio = sp.factor((volume / small).subs({ell: scale, circumference: 1}))
        ratios.append({"ell_over_C": str(scale), "V_over_Minkowski": str(ratio)})
    assert ratios[2]["V_over_Minkowski"] == "13/16"
    return {"large_branch": str(large), "ratios": ratios}


def plane_ratio(z):
    if z <= 0 or z >= mp.pi:
        raise ValueError("Finite-volume quadrature requires 0 < omega*L < pi.")

    def integrand(x):
        return mp.sqrt(
            mp.sinh(x)
            * mp.sinh(z - x)
            * mp.sin(x)
            * mp.sin(z - x)
            / (mp.sinh(z) * mp.sin(z))
        )

    return 6 * mp.quad(integrand, [0, z / 2, z]) / z**2


def plane_wave_checks():
    x, z, length, omega = sp.symbols("x z L omega", positive=True)
    assert (
        sp.trigsimp(sp.sin(z - x) * sp.cos(x) + sp.cos(z - x) * sp.sin(x) - sp.sin(z))
        == 0
    )
    assert (
        sp.simplify(
            sp.sinh(z - x) * sp.cosh(x) + sp.cosh(z - x) * sp.sinh(x) - sp.sinh(z)
        )
        == 0
    )
    # Dirichlet oscillator test mode gives the first loss of positivity.
    mode = sp.sin(sp.pi * x / length)
    energy = sp.integrate(sp.diff(mode, x) ** 2 - omega**2 * mode**2, (x, 0, length))
    assert sp.simplify(energy - length * ((sp.pi / length) ** 2 - omega**2) / 2) == 0
    assert sp.simplify(energy.subs(length, sp.pi / omega)) == 0
    assert sp.trigsimp(sp.cot(x) + sp.cot(sp.pi - x)) == 0
    parameters = ("0.5", "1", "2", "3", "pi-0.01", "pi-0.001")
    values = []
    for precision in (60, 90):
        with mp.workdps(precision):
            rows = []
            for parameter in parameters:
                zvalue = (
                    mp.pi - mp.mpf(parameter[3:])
                    if parameter.startswith("pi-")
                    else mp.mpf(parameter)
                )
                exact = plane_ratio(zvalue)
                local = 1 + zvalue**4 / 252
                rows.append((exact, 100 * (local / exact - 1)))
            values.append(rows)
    with mp.workdps(90):
        agreement = max(abs(a[0] / b[0] - 1) for a, b in zip(*values, strict=True))
        assert agreement < mp.mpf("1e-55")
        constant = (
            6
            / (mp.pi**2 * mp.sqrt(mp.sinh(mp.pi)))
            * mp.quad(
                lambda x: mp.sin(x) * mp.sqrt(mp.sinh(x) * mp.sinh(mp.pi - x)),
                [0, mp.pi / 2, mp.pi],
            )
        )
        near_critical = []
        for gap in ("0.001", "0.00001", "0.0000001"):
            ratio = plane_ratio(mp.pi - mp.mpf(gap)) * mp.sqrt(mp.mpf(gap)) / constant
            near_critical.append(
                {"pi_minus_z": gap, "V_ratio_sqrt_gap_over_C": mp.nstr(ratio, 18)}
            )
        assert abs(ratio - 1) < mp.mpf("1e-6")
        return {
            "finite_domain": "0 < z=omega*L < pi",
            "at_and_above_first_conjugate_span": "infinite volume; proof in note",
            "dirichlet_mode_energy": str(energy),
            "precision_agreement_relative": mp.nstr(agreement, 8),
            "rows": [
                {
                    "z": parameter,
                    "V_over_V0": mp.nstr(row[0], 20),
                    "quartic_relative_error_percent": mp.nstr(row[1], 12),
                }
                for parameter, row in zip(parameters, values[1], strict=True)
            ],
            "divergence_C": mp.nstr(constant, 22),
            "divergence_checks": near_critical,
        }


def constant_curvature_checks():
    z = sp.Symbol("z", positive=True)
    ds = 4 * sp.pi / 3 * (1 / sp.cosh(z) ** 2 + 2 * sp.log(sp.cosh(z)) - 1)
    ads = 4 * sp.pi / 3 * (sp.tan(z) ** 2 + 2 * sp.log(sp.cos(z)))
    assert sp.simplify(sp.diff(ds, z) - 8 * sp.pi * sp.tanh(z) ** 3 / 3) == 0
    assert sp.trigsimp(sp.diff(ads, z) - 8 * sp.pi * sp.tan(z) ** 3 / 3) == 0
    flat = 2 * sp.pi * z**4 / 3
    ds_series = sp.series(ds / flat, z, 0, 6)
    ads_series = sp.series(ads / flat, z, 0, 6)
    assert sp.expand(ds_series.removeO()) == 1 - 2 * z**2 / 3 + 11 * z**4 / 30
    assert sp.expand(ads_series.removeO()) == 1 + 2 * z**2 / 3 + 11 * z**4 / 30
    rows = []
    with mp.workdps(60):
        for kind, points in (
            ("dS", ("0.5", "1", "3", "5")),
            ("AdS", ("0.5", "1", "1.5")),
        ):
            for point in points:
                t = mp.mpf(point)
                if kind == "dS":
                    closed = (
                        4 * mp.pi / 3 * (mp.sech(t) ** 2 + 2 * mp.log(mp.cosh(t)) - 1)
                    )
                    numeric = 8 * mp.pi / 3 * mp.quad(lambda x: mp.tanh(x) ** 3, [0, t])
                else:
                    closed = 4 * mp.pi / 3 * (mp.tan(t) ** 2 + 2 * mp.log(mp.cos(t)))
                    numeric = 8 * mp.pi / 3 * mp.quad(lambda x: mp.tan(x) ** 3, [0, t])
                assert abs(closed / numeric - 1) < mp.mpf("1e-55")
                rows.append(
                    {
                        "geometry": kind,
                        "sqrt_abs_kappa_ell": point,
                        "V_over_V0": mp.nstr(closed / (2 * mp.pi * t**4 / 3), 18),
                    }
                )
    return {
        "dS_relative_series": str(ds_series),
        "AdS_relative_series": str(ads_series),
        "dS_finite_domain": "all finite ell > 0",
        "AdS_finite_domain": (
            "sqrt(-kappa)*ell < pi/2, universal cover without boundary points"
        ),
        "AdS_at_and_above_critical": "infinite volume; proof in note",
        "rows": rows,
    }


def intersection_of_revolution(front, back):
    """Exact integral of two piecewise-linear meridian profiles, mirrored in x."""
    x1, r1 = front
    x2, r2 = back
    end = min(x1[-1], x2[-1])
    knots = np.unique(
        np.concatenate(
            ([0.0, end], x1[(x1 > 0) & (x1 < end)], x2[(x2 > 0) & (x2 < end)])
        )
    )
    delta = np.interp(knots, x1, r1) - np.interp(knots, x2, r2)
    switch = delta[:-1] * delta[1:] < 0
    crossings = (
        knots[:-1][switch]
        - delta[:-1][switch] * np.diff(knots)[switch] / np.diff(delta)[switch]
    )
    knots = np.sort(np.concatenate((knots, crossings)))
    radius = np.minimum(np.interp(knots, x1, r1), np.interp(knots, x2, r2))
    return (
        2
        * np.pi
        / 3
        * np.sum(
            np.diff(knots)
            * (radius[:-1] ** 2 + radius[:-1] * radius[1:] + radius[1:] ** 2)
        )
    )


def kasner_limit_numeric(angles, ray_nodes, time_nodes):
    """Independent full null-front integration for p=(-1/3,2/3,2/3), b=1.

    tau=w^3 removes the integrable singularity of the ray integrals at tau=0.
    This is a refinement diagnostic, not interval arithmetic or a proof of error.
    """
    angle = np.pi / 2 * np.linspace(0, 1, angles + 1) ** 2
    kx, kr = np.cos(angle)[:, None], np.sin(angle)[:, None]
    nodes, weights = np.polynomial.legendre.leggauss(ray_nodes)

    def front(left, right):
        w = left + (right - left) * (nodes + 1) / 2
        weight = weights * (right - left) / 2
        denominator = np.sqrt(kr * kr + kx * kx * w[None, :] ** 6)
        axis = 3 * np.sum(weight * kx * w**6 / denominator, axis=1)
        radius = 3 * np.sum(weight * kr / denominator, axis=1)
        axis[0], radius[0] = 3 * (right**4 - left**4) / 4, 0
        axis[-1], radius[-1] = 0, 3 * (right - left)
        assert np.all(np.diff(axis) <= 1e-14)
        return axis[::-1], radius[::-1]

    # Flat-ball control exercises the independent intersection integrator.
    control = (np.sin(angle), np.cos(angle))
    flat = intersection_of_revolution(control, control)
    assert abs(flat / (4 * np.pi / 3) - 1) < 10 / angles**2
    tnodes, tweights = np.polynomial.legendre.leggauss(time_nodes)
    total = 0.0
    breaks = (0, 1 / 8, 2 ** (-3 / 4), 1)
    for left, right in zip(breaks[:-1], breaks[1:], strict=True):
        times = left + (right - left) * (tnodes + 1) / 2
        for tau, weight in zip(times, tweights * (right - left) / 2, strict=True):
            split = tau ** (1 / 3)
            total += (
                weight
                * tau
                * intersection_of_revolution(front(0, split), front(split, 1))
            )
    return {
        "angle_intervals": angles,
        "ray_nodes": ray_nodes,
        "time_nodes_per_piece": time_nodes,
        "V_0_to_1_approx": float(total),
        "flat_ball_relative_error": float(flat / (4 * np.pi / 3) - 1),
    }


def kasner_lower_certificate(start_w, intervals=256):
    """Rational lower enclosure from inscribed ellipsoid intersections.

    No floating-point number participates in the inequality. Machin's identity
    and alternating arctangent bounds enclose pi; integer square roots bound
    the lens subtraction outward. start_w is the cube root of the lower tip.
    """

    def atan_bounds(inverse):
        x = Fraction(1, inverse)
        partial = sum(
            ((-1) ** k * x ** (2 * k + 1) / (2 * k + 1) for k in range(12)),
            Fraction(0),
        )
        return partial, partial + x**25 / 25

    a_bounds, b_bounds = atan_bounds(5), atan_bounds(239)
    pi_lower = 16 * a_bounds[0] - 4 * b_bounds[1]
    pi_upper = 16 * a_bounds[1] - 4 * b_bounds[0]

    def sqrt_upper(value):
        scale = 10**18
        upper = Fraction(
            isqrt(value.numerator * scale**2 // value.denominator) + 1, scale
        )
        assert upper * upper >= value
        return upper

    def lens_lower(axis1, radius1, axis2, radius2):
        if min(axis1, radius1, axis2, radius2) == 0:
            return Fraction(0)
        if axis1 > axis2:
            axis1, radius1, axis2, radius2 = axis2, radius2, axis1, radius1
        bracket = axis1 * radius1**2
        if radius1 > radius2:
            gap = radius1**2 - radius2**2
            crossing_squared = gap / (radius1**2 / axis1**2 - radius2**2 / axis2**2)
            bracket -= gap * sqrt_upper(crossing_squared)
        return max(Fraction(0), 4 * pi_lower * bracket / 3)

    lower = Fraction(0)
    for j in range(intervals):
        w0 = start_w + (1 - start_w) * Fraction(j, intervals)
        w1 = start_w + (1 - start_w) * Fraction(j + 1, intervals)
        volume = lens_lower(
            Fraction(3, 4) * (w0**4 - start_w**4),
            3 * (w0 - start_w),
            Fraction(3, 4) * (1 - w1**4),
            3 * (1 - w1),
        )
        lower += (w1**6 - w0**6) * volume / 2
    start = start_w**3
    half, midpoint = (1 - start) / 2, (1 + start) / 2
    quartic_upper = (
        2 * pi_upper * half**4 / 3 * (1 + Fraction(68, 8505) * (half / midpoint) ** 4)
    )
    assert lower > quartic_upper
    scale = 10**12
    lower_floor = lower.numerator * scale // lower.denominator
    upper_ceiling = -((-quartic_upper.numerator * scale) // quartic_upper.denominator)
    assert lower_floor > upper_ceiling
    return {
        "lower_tip_tau": str(start),
        "upper_tip_tau": "1",
        "w_partition_intervals": intervals,
        "volume_lower_floor_1e12": lower_floor,
        "quartic_upper_ceiling_1e12": upper_ceiling,
        "strict_disagreement": True,
        "arithmetic": (
            "Exact Fraction; Machin identity with 12-term alternating bounds; "
            "integer sqrt upper enclosure at denominator 10^18"
        ),
    }


def kasner_checks():
    w = sp.Symbol("w", positive=True)
    q = 2 ** (-sp.Rational(1, 4))
    integral = sp.Rational(81, 4) * (
        sp.integrate(w**11, (w, 0, sp.Rational(1, 2)))
        + sp.integrate(w**9 * (1 - w) ** 2, (w, sp.Rational(1, 2), q))
        + sp.integrate(w**5 * (1 - w**4) * (1 - w) ** 2, (w, q, 1))
    )
    integral = sp.simplify(integral)
    lower, upper = 4 * sp.pi * integral / 3, 8 * integral
    # Check the closed hypergeometric support against a nonsingular ray integral.
    with mp.workdps(60):
        for kx, kr, lo, hi in ((2, 1, 0, 1), (3, 2, mp.mpf("0.2"), mp.mpf("0.7"))):

            def primitive(t, kx=kx, kr=kr):
                return (
                    3
                    * kr
                    * mp.root(t, 3)
                    * mp.hyp2f1(
                        -mp.mpf("0.5"),
                        mp.mpf(1) / 6,
                        mp.mpf(7) / 6,
                        -((mp.mpf(kx) / kr * t) ** 2),
                    )
                )

            closed = primitive(hi) - primitive(lo)
            direct = 3 * mp.quad(
                lambda x, kr=kr, kx=kx: mp.sqrt(kr * kr + kx * kx * x**6),
                [mp.root(lo, 3), mp.root(hi, 3)],
            )
            assert abs(closed / direct - 1) < mp.mpf("1e-55")
    rows = [
        kasner_limit_numeric(*counts)
        for counts in ((128, 64, 48), (256, 96, 72), (512, 128, 96), (1024, 192, 160))
    ]
    assert all(float(lower) < row["V_0_to_1_approx"] < float(upper) for row in rows)
    last_change = abs(rows[-1]["V_0_to_1_approx"] - rows[-2]["V_0_to_1_approx"])
    return {
        "exponents": ["-1/3", "2/3", "2/3"],
        "singularity": "Kretschmann=64/(27*tau^4); tau=0 is not a spacetime point",
        "support_primitive": "3*k_perp*tau^(1/3)*2F1(-1/2,1/6;7/6;-(kx*tau/k_perp)^2)",
        "bound_integral_exact": str(integral),
        "analytic_lower_exact": str(lower),
        "analytic_upper_exact": str(upper),
        "analytic_bounds_decimal": [str(sp.N(value, 22)) for value in (lower, upper)],
        "numeric_refinements_not_certified": rows,
        "last_refinement_change": last_change,
        "rational_lower_certificates": [
            kasner_lower_certificate(Fraction(0)),
            kasner_lower_certificate(Fraction(1, 256)),
        ],
        "quartic_extrapolation_at_boundary": str(
            sp.N(sp.pi / 24 * (1 + sp.Rational(68, 8505)), 20)
        ),
        "scaling": (
            "V(a,b)=b^4*F(a/b); the finite a->0 limit is not a diamond "
            "with a singular endpoint"
        ),
    }


if __name__ == "__main__":
    started = time.perf_counter()
    result = {
        "kind": (
            "Post-study exact global examples and deterministic diagnostics; "
            "not peer review"
        ),
        "date": "2026-09-22",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "mpmath": mp.__version__,
        "numpy": np.__version__,
        "source_sha256_lf": hashlib.sha256(
            Path(__file__).read_bytes().replace(b"\r\n", b"\n")
        ).hexdigest(),
        "flat_spatial_cylinder": cylinder_checks(),
        "plane_wave": plane_wave_checks(),
        "constant_curvature": constant_curvature_checks(),
        "Kasner_singular_limit": kasner_checks(),
        "numerical_scope": (
            "Precision/refinement agreement is diagnostic, not a certified "
            "error bound. Analytic counterexamples and Kasner bounds "
            "are proved separately in the note."
        ),
    }
    result["elapsed_seconds"] = round(time.perf_counter() - started, 3)
    print(json.dumps(result, indent=2))
