"""Challenge a nonlocal wave-correlation inverse from multiscale pair fractions.

--self-check verifies exact algebra and analytic inverse controls. --output DIR
runs the frozen continuum diagnostic. No finite-event accuracy claim is made.
"""

import argparse
import hashlib
import json
import subprocess
import time
from fractions import Fraction
from pathlib import Path
from statistics import NormalDist

import numpy as np
from numpy.polynomial import Chebyshev
from numpy.polynomial.chebyshev import chebvander
from order_wave_energy import (
    relation_probability,
    rosen_kernel_from_coefficients,
    rule,
    transfer,
)

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "docs/paper/paper_a/order_wave_energy/correlation_challenge.md"
SIZES = (1 - np.cos(np.arange(65) * np.pi / 64)) / 2
LAGS = np.linspace(0.15, 0.85, 29)
DEGREES = (12, 18, 24)


def inverse_from_difference(difference, degree):
    """Input is only 0.1-P2 at SIZES; no spectrum/energy is accepted."""
    vandermonde = SIZES[:, None] ** 4 * chebvander(2 * SIZES**2 - 1, degree)
    coefficients = np.linalg.lstsq(vandermonde, difference, rcond=1e-14)[0]
    x = Chebyshev([0.5, 0.5], domain=[0, 1])
    e = 262.5 * x * Chebyshev(coefficients, domain=[0, 1])
    current = e
    local = np.zeros_like(LAGS)
    for coefficient in (384, 74, 59, 10, 1):
        local += coefficient * current(LAGS**2)
        current = 2 * x * current.deriv()
    nodes, weights = rule(128)
    a, b = 9.5, np.sqrt(23) / 2
    angles = b * np.log(1 / nodes)
    green = nodes ** (a - 1) * (864 * np.cos(angles) + 18576 / b * np.sin(angles))
    integral = np.sum(weights * green * e((LAGS[:, None] * nodes) ** 2), axis=1)
    return (local - integral) / 105


def inputs(name):
    if name == "separated":
        n = np.array([1, 4, 9])
        weights = np.array([0.25, 0.35, 0.40])
        ellipticity, rotations = np.zeros(3), np.zeros(3)
    elif name == "rotating":
        n = np.array([3, 7, 11])
        weights = np.array([0.20, 0.35, 0.45])
        ellipticity, rotations = np.array([0.7, -0.4, 1]), np.array([0.2, 0.8, 1.3])
    elif name == "broad":
        n = np.arange(1, 9)
        weights = 1 / n
        weights /= weights.sum()
        ellipticity, rotations = 0.6 * (-1.0) ** (n + 1), 0.17 * n
    else:
        raise ValueError(name)
    return n, weights, ellipticity, rotations


def coefficients_at(u, name, q, phase):
    n, weights, eta, rotations = inputs(name)
    omega = 2 * n
    amplitudes = np.sqrt(2 * q * weights / (1 + eta**2)) / omega
    z = omega * u + n * phase + 0.31 * n
    sine, cosine = np.sin(z), np.cos(z)
    c, s = np.cos(rotations), np.sin(rotations)
    p = np.sum(amplitudes * (c * sine - s * eta * cosine))
    cross = np.sum(amplitudes * (s * sine + c * eta * cosine))
    pp = np.sum(amplitudes * omega * (c * cosine + s * eta * sine))
    cp = np.sum(amplitudes * omega * (s * cosine - c * eta * sine))
    radius2 = p * p + cross * cross
    radius = np.sqrt(radius2)
    sinhc = np.sinh(2 * radius) / (2 * radius) if radius else 1.0
    sinv = np.array(
        [
            np.cosh(2 * radius) - 2 * sinhc * p,
            -2 * sinhc * cross,
            np.cosh(2 * radius) + 2 * sinhc * p,
        ]
    )
    if radius2 < 1e-24:
        forcing = pp * pp + cp * cp
    else:
        forcing = (
            (p * pp + cross * cp) ** 2 + sinhc**2 * (p * cp - cross * pp) ** 2
        ) / radius2
    return sinv, forcing


def exact_curve(name, q, phases=24, quadrature=64):
    n = inputs(name)[0]
    curves = []
    for phase in np.arange(phases) * 2 * np.pi / phases:

        def coefficients(u, phase=phase):
            return coefficients_at(u, name, q, phase)

        kernel, _ = rosen_kernel_from_coefficients(
            q, 1 / (2 * max(n)), coefficients, tighter=True
        )
        values = [0.1]
        for length in SIZES[1:]:

            def scaled(s, t, length=length, kernel=kernel):
                return kernel(length * s, length * t) / length

            values.append(relation_probability(scaled, quadrature))
        curves.append(values)
    return np.array(curves)


def target(name, q):
    n, weights, _, _ = inputs(name)
    return q * np.sum(weights * (1 - np.cos(2 * n * LAGS[:, None])), axis=1)


def analytic_difference(name, q):
    n, weights, _, _ = inputs(name)
    e = np.array(
        [
            q
            * sum(w * transfer(2 * k * length) for k, w in zip(n, weights, strict=True))
            for length in SIZES
        ]
    )
    return e * SIZES**2 / 262.5


def scores(difference, name, q):
    truth = target(name, q)
    return [
        {
            "degree": degree,
            "recovered": (
                prediction := inverse_from_difference(difference, degree)
            ).tolist(),
            "normalized_rmse": float(np.sqrt(np.mean(((prediction - truth) / q) ** 2))),
            "normalized_max_error": float(np.max(np.abs((prediction - truth) / q))),
        }
        for degree in DEGREES
    ]


def self_check():
    for z in range(2, 25):
        lhs = (
            1
            + Fraction(360, (z + 2) * (z + 3) * (z + 4) * (z + 5))
            - Fraction(24, (z + 2) * (z + 3))
        )
        rhs = Fraction(
            z * (z - 1) * (z * z + 15 * z + 62), (z + 2) * (z + 3) * (z + 4) * (z + 5)
        )
        assert lhs == rhs
        numerator = (z + 3) * (z + 4) ** 2 * (z + 5) * (z + 6) * (z + 7)
        quotient = z**4 + 10 * z**3 + 59 * z * z + 74 * z + 384
        assert numerator == quotient * (z * z + 19 * z + 96) - (864 * z + 26784)
    # Monomials check theta derivatives, scale normalization, and Green integral.
    for power in (2, 4, 6, 8):
        difference = SIZES ** (power + 2) / 262.5
        multiplier = (
            (power + 3)
            * (power + 4) ** 2
            * (power + 5)
            * (power + 6)
            * (power + 7)
            / (105 * (power * power + 19 * power + 96))
        )
        error = inverse_from_difference(difference, 12) - multiplier * LAGS**power
        assert np.max(np.abs(error)) < 2e-7
    # Check exact Q using a distinct matrix derivative, including both polarizations.
    u, step = 0.19, 1e-6

    def matrix(at):
        values, _ = coefficients_at(at, "rotating", 0.07, 0.21)
        return np.array([[values[0], values[1]], [values[1], values[2]]])

    sinv = matrix(u)
    derivative = (matrix(u + step) - matrix(u - step)) / (2 * step)
    product = np.linalg.solve(sinv, derivative)
    _, q = coefficients_at(u, "rotating", 0.07, 0.21)
    assert abs(np.trace(product @ product) / 8 - q) < 2e-10
    print(
        "Controls passed: rational identity, inverse monomials, exact polarization Q."
    )


def run(output):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Freeze source/design and run from a clean checkout")
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    stamp = (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    )
    paths = [Path(__file__), Path(__file__).with_name("order_wave_energy.py"), DESIGN]
    manifest = {
        "code_version": stamp,
        "status": "in_progress",
        "files": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths
        },
        "sizes": SIZES.tolist(),
        "lags": LAGS.tolist(),
        "scope": "quadratic-order inverse; deterministic continuum diagnostic",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    started = time.time()
    with (output / "rows.jsonl").open("w", encoding="utf-8") as stream:

        def record(row):
            stream.write(json.dumps(row) + "\n")
            stream.flush()

        for name in ("separated", "rotating", "broad"):
            difference = analytic_difference(name, 1)
            record(
                {
                    "stage": "analytic_inverse",
                    "case": name,
                    "difference": difference.tolist(),
                    "target": target(name, 1).tolist(),
                    "scores": scores(difference, name, 1),
                }
            )
        for name in ("separated", "rotating", "broad"):
            for q in (0.0025, 0.00125):
                curves = exact_curve(name, q)
                difference = 0.1 - np.mean(curves, axis=0)
                record(
                    {
                        "stage": "exact_inverse",
                        "case": name,
                        "q": q,
                        "phases": 24,
                        "quadrature": 64,
                        "curves": curves.tolist(),
                        "difference": difference.tolist(),
                        "target": target(name, q).tolist(),
                        "scores": scores(difference, name, q),
                    }
                )
                print(f"Exact case complete: {name}, q={q}", flush=True)
        curves = exact_curve("rotating", 0.00125, phases=48, quadrature=128)
        difference = 0.1 - np.mean(curves, axis=0)
        record(
            {
                "stage": "refinement",
                "case": "rotating",
                "q": 0.00125,
                "phases": 48,
                "quadrature": 128,
                "curves": curves.tolist(),
                "difference": difference.tolist(),
                "target": target("rotating", 0.00125).tolist(),
                "scores": scores(difference, "rotating", 0.00125),
            }
        )
        power = NormalDist().inv_cdf(0.975) + NormalDist().inv_cdf(0.95)
        for degree in DEGREES:
            matrix = np.column_stack(
                [inverse_from_difference(row, degree) for row in np.eye(len(SIZES))]
            )
            # L=0 has exactly P2=.1, so there is no measured random value there.
            noise_norm = np.linalg.norm(matrix[:, 1:], axis=1)
            n = (power * 0.2 * noise_norm / (0.25 * 0.00125)) ** 2
            record(
                {
                    "stage": "conditioning",
                    "degree": degree,
                    "row_norms": noise_norm.tolist(),
                    "normal_planning_N_per_scale_for_each_lag": n.tolist(),
                    "assumptions": (
                        "Independent scale diamonds; flat variance; no bias; "
                        "not a lower bound"
                    ),
                }
            )
    manifest.update(
        status="completed",
        elapsed_seconds=time.time() - started,
        rows_sha256=hashlib.sha256((output / "rows.jsonl").read_bytes()).hexdigest(),
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(
        json.dumps(
            {
                k: manifest[k]
                for k in ("code_version", "status", "elapsed_seconds", "rows_sha256")
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.self_check:
        self_check()
    if args.output:
        run(args.output.resolve())
    if not (args.self_check or args.output):
        parser.error("choose --self-check or --output")
