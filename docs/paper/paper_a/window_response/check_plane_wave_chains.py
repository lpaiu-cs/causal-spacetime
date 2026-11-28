"""Deterministic checks of plane_wave_chain_hierarchy.md; no sampled posets.

Run with the existing NumPy/SciPy experiment environment. Grid refinements are
diagnostics, not certified error enclosures or a confirmatory statistical gate.
"""

import hashlib
import json
import math
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, simpson, solve_ivp
from scipy.optimize import brentq
from scipy.special import gammaln, spherical_jn

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/positive_control"))
from window_response import diamond_moments  # noqa: E402


def transfer(profile, start, end, evaluation):
    """Fundamental matrix of x''=H(u)x, integrated independently per call."""

    def rhs(u, flat):
        f = flat.reshape(4, 4)
        return np.vstack((f[2:], profile(u) @ f[:2])).ravel()

    sol = solve_ivp(
        rhs,
        (start, end),
        np.eye(4).ravel(),
        t_eval=evaluation,
        method="DOP853",
        rtol=2e-12,
        atol=2e-14,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T.reshape(-1, 4, 4)


def rotating(u):
    theta = 2 * math.pi * u + 0.3 * math.sin(2 * math.pi * u)
    amp = 4 * (1 + 0.2 * math.sin(4 * math.pi * u))
    c, s = math.cos(2 * theta), math.sin(2 * theta)
    return amp * np.array([[c, s], [s, -c]])


def changing(u):
    return np.diag([5 * math.cos(2 * math.pi * u), -5 * math.cos(2 * math.pi * u)])


def matrix_kernel(profile, panels):
    u = np.linspace(0, 1, panels + 1)
    f = transfer(profile, 0, 1, u)
    inv = np.linalg.inv(f)
    kernel = np.zeros((len(u), len(u)))
    for j in range(1, len(u)):
        segments = f[j] @ inv[:j]
        det = np.linalg.det(segments[:, :2, 2:])
        if np.any(det <= 0):
            raise ValueError(
                "Nonpositive Jacobi determinant: finite-diamond scope failed"
            )
        kernel[j, :j] = np.sqrt(det)
    symplectic = np.block(
        [[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]]
    )
    residual = np.max(np.abs(f.transpose(0, 2, 1) @ symplectic @ f - symplectic))
    return kernel, float(residual)


def constant_kernel(z, panels):
    if not 0 <= z < math.pi:
        raise ValueError("Requires 0 <= omega*L < pi; the limiting kernel is separate")
    u = np.linspace(0, 1, panels + 1)
    s = np.maximum(u[:, None] - u[None, :], 0)
    if z == 0:
        return s
    return np.sqrt(np.sin(z * s) * np.sinh(z * s)) / z


def convolution(kernel, f):
    # ponytail: O(M^2) quadrature suffices for this audit; no general solver framework.
    step = 1 / (len(f) - 1)
    return np.array(
        [0.0]
        + [simpson(kernel[j, : j + 1] * f[: j + 1], dx=step) for j in range(1, len(f))]
    )


def chain_statistics(kernel, max_k=6):
    """d=4, [a,b]=[0,1], effective h=1. C_r are Volterra powers."""
    f = kernel[:, 0]
    powers = [float(f[-1])]
    c2 = None
    for r in range(2, max_k + 2):
        f = convolution(kernel, f)
        powers.append(float(f[-1]))
        if r == 2:
            c2 = f
    kfull, cfull = powers[:2]
    gamma = cfull / kfull
    volume = math.pi * gamma
    ordered = [
        (2 * math.pi) ** k / math.factorial(2 * k) * powers[k] / kfull
        for k in range(1, max_k + 1)
    ]
    probabilities = [
        math.factorial(k) * ordered[k - 1] / volume**k for k in range(1, max_k + 1)
    ]
    c2right = convolution(kernel[::-1, ::-1].T, kernel[-1, ::-1])[::-1]
    alpha = np.divide(c2, kernel[:, 0], out=np.zeros_like(c2), where=kernel[:, 0] > 0)
    beta = np.divide(c2right, kernel[-1], out=np.zeros_like(c2), where=kernel[-1] > 0)
    w = kernel[:, 0] * kernel[-1] / kfull
    step = 1 / (len(f) - 1)
    qmean = simpson(w * (alpha + beta), dx=step) / (6 * gamma**2)
    qsecond = (
        simpson(w * ((alpha**2 + beta**2) / 15 + alpha * beta / 45), dx=step) / gamma**3
    )
    return {
        "volume": volume,
        "P_k": probabilities,
        "I_k": ordered,
        "C_r": powers,
        "q_mean": float(qmean),
        "q_second": float(qsecond),
        "zeta": float(qsecond - probabilities[1] ** 2),
        "P3_over_P2_squared": probabilities[2] / probabilities[1] ** 2,
    }


def log_flat_probability(d, k):
    return (
        gammaln(k)
        + (k - 1) * (gammaln(d + 1) - math.log(2))
        + gammaln(d / 2)
        + gammaln(d)
        - gammaln(k * d / 2)
        - gammaln((k + 1) * d / 2)
    )


def flat_checks():
    stats = chain_statistics(constant_kernel(0, 1024))
    expected = [
        Fraction(
            math.factorial(k) * 12**k, math.factorial(2 * k) * math.factorial(2 * k + 1)
        )
        for k in range(1, 7)
    ]
    error = max(
        abs(a / float(b) - 1) for a, b in zip(stats["P_k"], expected, strict=True)
    )
    assert error < 2e-8, error
    assert abs(stats["zeta"] - 1 / 100) < 1e-10
    for d in (3, 4, 5, 6, 8):
        nu = d / 2
        for k in range(1, 7):
            # Independent flat Dirichlet-integral simplification of the theorem.
            log_p = (
                gammaln(k + 1)
                + k * (gammaln(nu + 1) + gammaln(2 * nu))
                - gammaln(k * nu + 1)
                - gammaln((k + 1) * nu)
                + (1 - k) * gammaln(nu)
            )
            assert abs(log_p - log_flat_probability(d, k)) < 1e-12
    return {
        "P_1_through_P_6_exact": [str(x) for x in expected],
        "max_relative_grid_error": error,
        "zeta": stats["zeta"],
        "P3_over_P2_squared_exact": "2/7",
    }


def weak_profile_coefficients():
    """Exact simplex/Beta moments, independent of numerical quadrature."""
    rows = []
    for k in range(1, 7):
        r = k + 1
        # Under Dirichlet(2,...,2), E sum_i gap_i^4 = r*(2)_4/(2r)_4.
        beta_coefficient = Fraction(
            r * math.factorial(5),
            180 * math.prod(range(2 * r, 2 * r + 4)),
        )
        closed = Fraction(1, 3 * (2 * r + 1) * (2 * r + 2) * (2 * r + 3))
        assert beta_coefficient == closed
        b = Fraction(1, 180) - closed - Fraction(k, 252)
        rows.append({"k": k, "relative_P_k_z4_coefficient": str(b)})
    assert Fraction(rows[1]["relative_P_k_z4_coefficient"]) / 10 == -Fraction(23, 75600)
    ratio = Fraction(rows[2]["relative_P_k_z4_coefficient"]) - 2 * Fraction(
        rows[1]["relative_P_k_z4_coefficient"]
    )
    assert ratio == -Fraction(5, 8316)
    return {"P_coefficients": rows, "relative_chain_ratio_z4_coefficient": str(ratio)}


def determinant_checks():
    rows = []
    for name, profile in (
        ("rotating_polarization", rotating),
        ("changing_sign", changing),
    ):
        for k in range(1, 6):
            times = np.r_[0, np.linspace(0.09, 0.91, k), 1]
            fs = [
                transfer(profile, a, b, [b])[0]
                for a, b in zip(times[:-1], times[1:], strict=True)
            ]
            full = transfer(profile, 0, 1, [1])[0]
            m = np.zeros((2 * k, 2 * k))
            for j in range(k):
                left, right = fs[j], fs[j + 1]
                l_inv, r_inv = np.linalg.inv(left[:2, 2:]), np.linalg.inv(right[:2, 2:])
                m[2 * j : 2 * j + 2, 2 * j : 2 * j + 2] = (
                    left[2:, 2:] @ l_inv + r_inv @ right[:2, :2]
                )
                if j < k - 1:
                    m[2 * j : 2 * j + 2, 2 * j + 2 : 2 * j + 4] = -r_inv
                    m[2 * j + 2 : 2 * j + 4, 2 * j : 2 * j + 2] = -r_inv.T
            predicted = np.linalg.det(full[:2, 2:]) / math.prod(
                np.linalg.det(f[:2, 2:]) for f in fs
            )
            rel = abs(np.linalg.det(m) / predicted - 1)
            eigen_min = float(np.linalg.eigvalsh((m + m.T) / 2).min())
            assert eigen_min > 0 and rel < 2e-10
            assert np.max(np.abs(m - m.T)) < 1e-10
            rows.append(
                {
                    "profile": name,
                    "k": k,
                    "relative_error": float(rel),
                    "minimum_Hessian_eigenvalue": eigen_min,
                }
            )
    commutator = rotating(0.1) @ rotating(0.37) - rotating(0.37) @ rotating(0.1)
    assert np.linalg.norm(commutator) > 1
    return {
        "partition_checks": rows,
        "noncommuting_profile_commutator_norm": float(np.linalg.norm(commutator)),
        "analytic_positivity_bound": "sup ||H|| <= 5 < pi^2 on length 1 (Poincare)",
    }


def profile_checks():
    rows = []
    for name, profile in (
        ("rotating_polarization", rotating),
        ("changing_sign", changing),
    ):
        refinements = []
        for panels in (256, 512, 1024, 2048):
            kernel, residual = matrix_kernel(profile, panels)
            stats = chain_statistics(kernel)
            refinements.append(
                {"panels": panels, "symplectic_residual": residual, **stats}
            )
        p_old, p_new = (np.array(row["P_k"]) for row in refinements[-2:])
        change = float(np.max(np.abs(p_new / p_old - 1)))
        assert change < 1e-7
        assert abs(refinements[-1]["q_mean"] - p_new[1]) < 2e-10
        reversed_kernel, _ = matrix_kernel(
            lambda u, profile=profile: profile(1 - u), 2048
        )
        reverse = chain_statistics(reversed_kernel)
        reversal_error = float(np.max(np.abs(np.array(reverse["P_k"]) / p_new - 1)))
        assert reversal_error < 2e-8
        rows.append(
            {
                "profile": name,
                "refinements": refinements,
                "last_relative_P_change": change,
                "time_reversal_relative_error": reversal_error,
            }
        )
    return rows


def old_formula_checks():
    rows = []
    for amplitude in (0, 1 / 64, 1, 2):
        new = chain_statistics(constant_kernel(math.sqrt(amplitude), 1024))
        old = diamond_moments(amplitude, order=96)
        errors = {
            "volume": new["volume"] - old["volume"],
            "p": new["P_k"][1] - old["p"],
            "zeta": new["zeta"] - old["zeta"],
        }
        assert max(abs(x) for x in errors.values()) < 2e-10
        rows.append(
            {"amplitude": amplitude, "differences_from_old_Gauss_formula": errors}
        )
    return rows


def multiplicity_two_checks():
    """Independent closed Laplace/Bessel control, H=-omega^2 I (nonvacuum)."""
    u = np.linspace(0, 1, 2049)
    gaps = np.maximum(u[:, None] - u[None, :], 0)
    rows = []
    for z in (1.0, math.pi - 0.01, math.pi - 0.000001):
        stats = chain_statistics(np.sin(z * gaps) / z)
        closed = [
            float(spherical_jn(r - 1, z) / ((2 * z) ** (r - 1) * math.factorial(r - 1)))
            for r in range(1, 8)
        ]
        rel = max(abs(x / y - 1) for x, y in zip(stats["C_r"], closed, strict=True))
        assert rel < 2e-8
        delta = math.pi - z
        rows.append(
            {
                "pi_minus_z": delta,
                "max_relative_error_to_closed_C_r": rel,
                "P_k_divided_by_delta_power_k_minus_1": [
                    p / delta ** (k - 1) for k, p in enumerate(stats["P_k"], 1)
                ],
                "chain_ratio": stats["P3_over_P2_squared"],
            }
        )
    limit = 2 * (15 - math.pi**2) / 45
    assert abs(rows[-1]["chain_ratio"] - limit) < 2e-6
    return {
        "scope": "Nonvacuum isotropic plane wave, determinant zero of order two",
        "rows": rows,
        "exact_limiting_ratio": "2*(15-pi^2)/45",
        "limit": limit,
    }


def caustic_checks():
    rows = []
    for delta in (0.1, 0.01, 0.001, 0.0001, 0.000001):
        z = math.pi - delta
        coarse = chain_statistics(constant_kernel(z, 1024), max_k=4)
        fine = chain_statistics(constant_kernel(z, 2048), max_k=4)
        p = fine["P_k"]
        d_eff = brentq(
            lambda d, p2=p[1]: log_flat_probability(d, 2) - math.log(p2), 2, 100
        )
        flat_r = math.exp(
            log_flat_probability(d_eff, 3) - 2 * log_flat_probability(d_eff, 2)
        )
        change = max(
            abs(x / y - 1) for x, y in zip(fine["P_k"], coarse["P_k"], strict=True)
        )
        assert change < 2e-7
        rows.append(
            {
                "pi_minus_z": delta,
                "volume": fine["volume"],
                "P_k": p,
                "I_k_times_sqrt_delta": [x * math.sqrt(delta) for x in fine["I_k"]],
                "P_k_divided_by_predicted_power": [
                    x / delta ** ((k - 1) / 2) for k, x in enumerate(p, 1)
                ],
                "flat_calibrated_dimension": d_eff,
                "wave_P3_over_P2_squared": fine["P3_over_P2_squared"],
                "flat_ratio_at_same_P2": flat_r,
                "expected_3_chains_at_N1200": math.comb(1200, 3) * p[2],
                "expected_flat_3_chains_at_N1200_same_P2": (
                    math.comb(1200, 3) * flat_r * p[1] ** 2
                ),
                "last_relative_P_change": change,
            }
        )

    # Independent nested adaptive quadrature at the limiting kernel, no grid.
    def kc(s):
        return (
            math.sqrt(max(0, math.sin(math.pi * s) * math.sinh(math.pi * s))) / math.pi
        )

    def c2(t):
        return quad(lambda s: kc(s) * kc(t - s), 0, t, epsabs=2e-12, epsrel=2e-11)[0]

    def c3(t):
        return quad(lambda s: c2(s) * kc(t - s), 0, t, epsabs=2e-12, epsrel=2e-10)[0]

    # C4= C2*C2, reducing the nested adaptive integral's cost.
    limits = [c2(1), c3(1), quad(lambda s: c2(s) * c2(1 - s), 0, 1, epsabs=2e-12)[0]]
    ratio = 3 / 5 * limits[0] * limits[2] / limits[1] ** 2
    assert abs(rows[-1]["wave_P3_over_P2_squared"] - ratio) < 2e-6
    return {
        "rows": rows,
        "limiting_C2_C3_C4_adaptive": limits,
        "limiting_chain_ratio": ratio,
        "flat_large_d_log_rate_for_ratio": -1.5 * math.log(4 / 3),
    }


if __name__ == "__main__":
    started = time.perf_counter()
    result = {
        "kind": (
            "Post-study deterministic theorem checks; "
            "no peer review or confirmatory samples"
        ),
        "date": "2026-09-22",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_sha256_lf": hashlib.sha256(
            Path(__file__).read_bytes().replace(b"\r\n", b"\n")
        ).hexdigest(),
        "old_crosscheck_source_sha256_lf": hashlib.sha256(
            (ROOT / "experiments/positive_control/window_response.py")
            .read_bytes()
            .replace(b"\r\n", b"\n")
        ).hexdigest(),
        "flat": flat_checks(),
        "weak_profile_coefficients": weak_profile_coefficients(),
        "noncommuting_determinants": determinant_checks(),
        "arbitrary_profiles": profile_checks(),
        "old_formula_crosscheck": old_formula_checks(),
        "multiplicity_two_closed_control": multiplicity_two_checks(),
        "caustic": caustic_checks(),
        "numerical_scope": (
            "Refinement and cross-method agreement are not certified error enclosures."
        ),
    }
    result["elapsed_seconds"] = round(time.perf_counter() - started, 3)
    destination = Path(__file__).with_name("plane_wave_chain_audit.json")
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "written": str(destination),
                "all_checks_passed": True,
                "elapsed_seconds": result["elapsed_seconds"],
                "caustic_limiting_chain_ratio": result["caustic"][
                    "limiting_chain_ratio"
                ],
            },
            indent=2,
        )
    )
