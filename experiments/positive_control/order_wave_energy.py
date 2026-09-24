"""Deterministic challenge of an order-only high-frequency wave-energy readout.

Run --self-check for analytic controls; --output DIR runs the frozen diagnostic.
This is a continuum calculation, not a finite-sprinkling accuracy campaign.
Requires the already installed NumPy and SciPy; no GPU is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import time
from functools import cache
from pathlib import Path

import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.special import hyp0f1, sici

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "docs/paper/paper_a/order_wave_energy/design.md"


@cache
def rule(n):
    x, w = leggauss(n)
    return (x + 1) / 2, w / 2


def shape(u, q, lam, profile, phase):
    """Return beta, theta, and Q with a''=-Q*a (all cases are exact vacuum)."""
    z = u / lam + phase
    if profile == "circular":
        return 0.5 * np.arcsinh(2 * lam * np.sqrt(q)), z / 2, q
    if profile == "sine":
        f, fp = np.sin(z), np.cos(z)
    elif profile == "mixed":
        f = (np.sin(z) + 0.25 * np.sin(2 * z)) / np.sqrt(1.25)
        fp = (np.cos(z) + 0.5 * np.cos(2 * z)) / np.sqrt(1.25)
    else:
        raise ValueError(profile)
    return lam * np.sqrt(2 * q) * f, 0.0, 2 * q * fp**2


def inverse_shape(beta, theta):
    c, s = np.cosh(2 * beta), np.sinh(2 * beta)
    return np.array(
        [c - s * np.cos(2 * theta), -s * np.sin(2 * theta), c + s * np.cos(2 * theta)]
    )


def rosen_kernel(q, lam, profile, phase=0.37, *, tighter=False):
    """K(s,t)=a(s)a(t)sqrt(det integral_s^t gamma^-1), on [-1/2,1/2]."""

    def coefficients(u):
        beta, theta, forcing = shape(u, q, lam, profile, phase)
        return inverse_shape(beta, theta), forcing

    return rosen_kernel_from_coefficients(q, lam, coefficients, tighter=tighter)


def rosen_kernel_from_coefficients(q, lam, coefficients, *, tighter=False):
    """Solve the vacuum Rosen ODE from (inverse unimodular metric, Q)."""
    if not (0 <= q < math.pi**2 and lam > 0):
        raise ValueError("q must be in [0,pi^2); lambda must be positive")
    k = np.sqrt(q)

    def rhs(u, state):
        a, ap = state[:2]
        sinv, forcing = coefficients(u)
        return np.r_[ap, -forcing * a, sinv / a**2]

    initial = [np.cos(k / 2), k * np.sin(k / 2), 0, 0, 0]
    sol = solve_ivp(
        rhs,
        (-0.5, 0.5),
        initial,
        method="DOP853",
        dense_output=True,
        rtol=2e-13 if tighter else 2e-11,
        atol=2e-15 if tighter else 2e-13,
        max_step=min(0.02, lam / (24 if tighter else 12)),
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    minimum = float(np.min(sol.sol(np.linspace(-0.5, 0.5, 4097))[0]))
    if minimum <= 0:
        raise RuntimeError("Rosen chart has reached a caustic")

    def kernel(s, t):
        s, t = np.broadcast_arrays(s, t)
        if np.any(t < s):
            raise ValueError("kernel requires s <= t")
        left = sol.sol(s.ravel()).reshape((5,) + s.shape)
        right = sol.sol(t.ravel()).reshape((5,) + t.shape)
        delta = right[2:] - left[2:]
        det = delta[0] * delta[2] - delta[1] ** 2
        if np.any((t > s) & (det <= 0)):
            raise RuntimeError("nonpositive segment determinant")
        return left[0] * right[0] * np.sqrt(np.maximum(det, 0))

    return kernel, {"minimum_a": minimum, "ode_evaluations": sol.nfev}


def relation_probability(kernel, n):
    """Exact chain hierarchy reduced to 1D and 2D Gaussian quadratures."""
    x, w = rule(n)
    v = x - 0.5
    c2 = np.sum(w * kernel(-0.5, v) * kernel(v, 0.5))
    u = -0.5 + x[:, None] * x[None, :]
    c3 = np.sum(
        w[None, :]
        * x[None, :]
        * kernel(v, 0.5)[None, :]
        * np.sum(w[:, None] * kernel(-0.5, u) * kernel(u, v[None, :]), axis=0)
    )
    return float(kernel(-0.5, 0.5) * c3 / (3 * c2**2))


def checked_probability(kernel):
    values = {n: relation_probability(kernel, n) for n in (96, 192)}
    if abs(values[192] - values[96]) > 2e-11:
        values[384] = relation_probability(kernel, 384)
    last = list(values.values())
    return last[-1], {
        "quadrature_values": values,
        "last_quadrature_difference": abs(last[-1] - last[-2]),
    }


def background_probability(q):
    # Laplace transform C_m=(p^2+q)^(-m); hyp0f1 avoids small-q cancellation.
    c = [hyp0f1(m + 0.5, -q / 4) / math.factorial(2 * m - 1) for m in (1, 2, 3)]
    return float(c[0] * c[2] / (3 * c[1] ** 2))


def readout(p2):
    # L=1 here. For arbitrary L this is L^2 times the energy-like readout.
    return 52.5 * (1 - 10 * p2)


def kernel_response(z):
    """D(ell)/(q*ell^3); stable at z=0 (second order in wave amplitude)."""
    z = np.asarray(z, dtype=float)
    result = np.zeros_like(z)
    small = np.abs(z) < 1
    for n in range(3, 13):
        result[small] += 4 * (-1) ** n * z[small] ** (2 * n - 4) / math.factorial(2 * n)
    large = z[~small]
    result[~small] = -1 / 6 + 2 / large**2 - 8 * np.sin(large / 2) ** 2 / large**4
    return result


def transfer(r):
    """Phase-averaged weak-wave energy readout divided by 2q."""
    x, w = rule(max(96, int(r) + 32))
    d = x**3 * kernel_response(r * x)
    log_response = kernel_response(r) + np.sum(
        w * (60 * (1 - x) ** 3 - 24 * (1 - x)) * d
    )
    return float(-105 / 4 * log_response)


def spectrum_prediction(r, profile):
    if profile == "sine":
        return transfer(r)
    if profile == "mixed":
        return 0.8 * transfer(r) + 0.2 * transfer(2 * r)
    raise ValueError(profile)


def self_check():
    flat = lambda s, t: np.asarray(t) - np.asarray(s)  # noqa: E731
    assert abs(relation_probability(flat, 24) - 0.1) < 2e-15
    for q in (0.0001, 0.25, 1.0, 4.0):
        k = np.sqrt(q)
        kernel = lambda s, t, k=k: np.sin(k * (t - s)) / k  # noqa: E731
        assert abs(relation_probability(kernel, 48) - background_probability(q)) < 2e-14
    derivative = (background_probability(-1e-4) - background_probability(1e-4)) / 2e-4
    assert abs(derivative - 2 / 525) < 1e-10
    assert abs(transfer(0.001) / 0.001**2 - 23 / 288) < 1e-8
    # Verify the rotating-polarization Q identity through a separate derivative.
    q, lam, u, step = 0.7, 0.13, 0.19, 1e-6

    def matrix(at):
        beta, theta, _ = shape(at, q, lam, "circular", 0.21)
        entries = inverse_shape(-beta, theta)
        return np.array([[entries[0], entries[1]], [entries[1], entries[2]]])

    s = matrix(u)
    sp = (matrix(u + step) - matrix(u - step)) / (2 * step)
    product = np.linalg.solve(s, sp)
    assert abs(np.trace(product @ product) / 8 - q) < 1e-9
    kernel, _ = rosen_kernel(0, 0.13, "sine")
    assert abs(kernel(-0.4, 0.3) - 0.7) < 2e-13
    # Independent Brinkmann Jacobi integration, using trace-free H directly.
    q, lam, phase = 0.7, 0.13, 0.37
    k = np.sqrt(q)

    def jacobi(u, y):
        a, ap, bx, bxp, by, byp = y
        z = u / lam + phase
        bp = np.sqrt(2 * q) * np.cos(z)
        bpp = -np.sqrt(2 * q) * np.sin(z) / lam
        h = bpp + 2 * ap * bp / a
        return [ap, -(bp**2) * a, bxp, h * bx, byp, -h * by]

    independent = solve_ivp(
        jacobi,
        (-0.5, 0.5),
        [np.cos(k / 2), k * np.sin(k / 2), 0, 1, 0, 1],
        method="DOP853",
        rtol=2e-13,
        atol=2e-15,
        max_step=lam / 24,
    )
    assert independent.success
    reference = np.sqrt(independent.y[2, -1] * independent.y[4, -1])
    kernel, _ = rosen_kernel(q, lam, "sine", phase, tighter=True)
    assert abs(kernel(-0.5, 0.5) - reference) < 2e-11
    for r in (0.3, 3.0, 8.0, 64.0):
        a = (
            36 * (np.euler_gamma + np.log(r) - sici(r)[1])
            - 85
            - np.cos(r)
            + 36 * np.sin(r) / r
            + (180 - 60 * np.cos(r)) / r**2
            - 120 * np.sin(r) / r**3
        )
        assert abs(transfer(r) - (1 + 105 * a / r**4)) < 2e-8
    print("Controls passed: flat, Ricci, circular Q, Brinkmann, transfer identity.")


def run(output):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit the design and source; run from a clean worktree")
    stamp = (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    )
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    started = time.time()
    manifest = {
        "status": "diagnostic_in_progress",
        "code_version": stamp,
        "design_sha256": hashlib.sha256(DESIGN.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "utc_started": time.strftime("%FT%TZ", time.gmtime()),
        "scope": "deterministic continuum diagnostic, no sampling confidence claim",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    with (output / "rows.jsonl").open("w", encoding="utf-8") as stream:

        def record(row):
            stream.write(json.dumps(row) + "\n")
            stream.flush()

        for q in (0.25, 1.0, 4.0):
            target = background_probability(q)
            for profile in ("sine", "mixed", "circular"):
                for denominator in (8, 16, 32, 64, 128):
                    lam = 1 / denominator
                    kernel, details = rosen_kernel(q, lam, profile)
                    p2, quadrature = checked_probability(kernel)
                    record(
                        {
                            "stage": "limit",
                            "q": q,
                            "lambda": lam,
                            "profile": profile,
                            "phase": 0.37,
                            "p2": p2,
                            "background_p2": target,
                            "residual": p2 - target,
                            **details,
                            **quadrature,
                        }
                    )
                print(f"Limit family complete: q={q}, {profile}", flush=True)

        for profile in ("sine", "mixed"):
            for r in (0.25, 0.5, 1.0, 2.0, 4.0, 2 * math.pi, 8.0, 16.0, 64.0):
                for q in (0.0025, 0.00125):
                    values = []
                    errors = []
                    for phase in np.arange(12) * 2 * math.pi / 12:
                        kernel, _ = rosen_kernel(q, 1 / r, profile, float(phase))
                        p2, details = checked_probability(kernel)
                        values.append(readout(p2) / (2 * q))
                        errors.append(details["last_quadrature_difference"])
                    mean = float(np.mean(values))
                    record(
                        {
                            "stage": "transfer",
                            "profile": profile,
                            "r": r,
                            "q": q,
                            "phases": 12,
                            "phase_values": values,
                            "mean": mean,
                            "prediction": spectrum_prediction(r, profile),
                            "residual": mean - spectrum_prediction(r, profile),
                            "max_quadrature_difference": max(errors),
                        }
                    )
                print(f"Transfer family complete: {profile}, r={r:.6g}", flush=True)

        # Predeclared solver/quadrature/phase refinement; retain both runs.
        for profile, q, r in (
            ("mixed", 4.0, 128.0),
            ("sine", 0.00125, 0.25),
            ("mixed", 0.00125, 8.0),
        ):
            phases = [0.37] if q == 4 else np.arange(24) * 2 * math.pi / 24
            p2_values = []
            for phase in phases:
                kernel, _ = rosen_kernel(q, 1 / r, profile, float(phase), tighter=True)
                p2_values.append(relation_probability(kernel, 384))
            record(
                {
                    "stage": "refinement",
                    "profile": profile,
                    "q": q,
                    "r": r,
                    "phases": len(phases),
                    "p2_values": p2_values,
                    "mean_readout_ratio": float(
                        np.mean([readout(p) / (2 * q) for p in p2_values])
                    ),
                }
            )

        for r in np.geomspace(0.01, 128, 160):
            record(
                {
                    "stage": "analytic_transfer",
                    "r": float(r),
                    "sine": transfer(r),
                    "mixed": spectrum_prediction(r, "mixed"),
                }
            )
    manifest.update(
        status="diagnostic_completed",
        elapsed_seconds=time.time() - started,
        rows_sha256=hashlib.sha256((output / "rows.jsonl").read_bytes()).hexdigest(),
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2), flush=True)


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
        parser.error("select --self-check or --output")
