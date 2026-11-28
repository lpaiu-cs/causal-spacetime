"""Frozen independent coverage validation; evaluator truth never enters the fitter."""

import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
import scipy
from intrinsic_weyl_joint import (
    TIDAL_GRAM,
    WEYL,
    axis_signs,
    confidence_ball,
    fit_joint,
    fit_order_archive,
    gauge,
    tidal_distance,
)
from intrinsic_weyl_joint_audit import (
    aligned_tensor_error,
    normal_coordinates,
    volume_table,
)
from intrinsic_weyl_joint_pilot import ROOT, generate_archive
from scipy.stats import beta, binom

N_CASES = 256
POPULATION_SEED = 2026092250
DATA_SEED = 2026092300
MEASUREMENT_N = 262144
BIAS = 0.02


def selfcheck():
    c = np.linspace(-0.3, 0.4, 10)
    tensor = np.einsum("k,kabcd->abcd", c, WEYL)
    for signs in ((-1, 1, 1), (1, -1, 1), (1, 1, -1), (-1, -1, -1)):
        s = np.r_[1, signs]
        transformed = tensor * np.einsum("a,b,c,d->abcd", s, s, s, s)
        expected = np.einsum("k,kabcd->abcd", c * axis_signs(signs), WEYL)
        assert np.max(abs(transformed - expected)) < 1e-14
        d = np.diag(axis_signs(signs))
        assert np.allclose(d @ TIDAL_GRAM @ d, TIDAL_GRAM)
    assert tidal_distance(c, c * axis_signs([1, 1, -1])) < 1e-14
    assert tidal_distance(c, 0.5 * c) > 0
    fit = {
        "optimizer_success": True,
        "boundary_excess": 0,
        "coordinate_rank": 5,
        "coordinate_parameters": 5,
        "projected_curvature_rank": 10,
        "coordinates_inferred_from_counts": [
            [-0.5, 0, 0, 0],
            [0, 0.15, 0, 0],
            [0, 0.1, 0.12, 0],
            [0.5, 0, 0, 0],
        ],
        "coefficients": c,
        "nominal_component_covariance": np.eye(10),
    }
    before = confidence_ball(fit)
    signs = [-1, 1, -1]
    transform = axis_signs(signs)
    reflected = fit | {
        "coordinates_inferred_from_counts": np.array(
            fit["coordinates_inferred_from_counts"]
        )
        * np.r_[1, signs],
        "coefficients": transform * c,
        "nominal_component_covariance": transform[:, None] * np.eye(10) * transform,
    }
    after = confidence_ball(reflected)
    assert (
        tidal_distance(
            before["coefficients_xy_positive"], after["coefficients_xy_positive"]
        )
        < 1e-14
    )
    assert abs(before["radius"] - after["radius"]) < 1e-14
    minimum = next(
        k for k in range(1, N_CASES + 1) if beta.ppf(0.025, k, N_CASES - k + 1) >= 0.95
    )
    assert minimum == 251
    assert binom.sf(minimum - 1, N_CASES, 0.99) > 0.90
    return {
        "tensor_parity": "passed",
        "minimum_valid_covered": minimum,
        "power_at_099": float(binom.sf(minimum - 1, N_CASES, 0.99)),
    }


def evaluate(path, truth_path):
    # The estimator entrypoint reads exactly the order archive. The separate
    # truth archive is opened only after it has returned.
    fit, clocks = fit_order_archive(path, with_clocks=True)
    inference = confidence_ball(fit, relative_bias=BIAS)
    with np.load(truth_path, allow_pickle=False) as truth:
        points, amplitude = truth["points"], float(truth["amplitude"])
    normal = normal_coordinates(points, amplitude)
    canonical = gauge(normal)
    expected = aligned_tensor_error(
        {"coefficients": np.zeros(10), "coordinates_inferred_from_counts": canonical},
        normal,
        amplitude,
    )["aligned_truth_coefficients"]
    counts, mask, covariance, _, _ = clocks
    volume = volume_table(points, amplitude)
    coarse = volume_table(points, amplitude, 32)
    refinement = float(np.max(abs(coarse[mask] / volume[mask] - 1)))
    if refinement >= 1e-10:
        raise ArithmeticError("Continuum quadrature did not stabilize")
    mean = fit_joint(
        counts[0, -1] * volume / volume[0, -1],
        clock_mask=mask,
        clock_covariance=covariance,
    )
    mean_inference = confidence_ball(mean, relative_bias=BIAS)
    model_bias = tidal_distance(mean_inference["coefficients_xy_positive"], expected)
    error = tidal_distance(inference["coefficients_xy_positive"], expected)
    norm = 2 * amplitude
    bias_ok = model_bias <= BIAS * norm + 1e-10
    covered = error <= inference["radius"]
    return {
        "fit": fit,
        "inference": inference,
        "evaluation": {
            "truth_coefficients_xy_positive": expected,
            "truth_norm": norm,
            "full_tidal_error_mod_reflection": error,
            "continuum_bias_mod_reflection": model_bias,
            "quadrature_relative_change": refinement,
            "bias_condition_met": bool(bias_ok),
            "covered": bool(covered),
            "valid_and_covered": bool(bias_ok and covered),
        },
    }


def run(destination):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit protocol and code before validation data")
    if destination.exists():
        raise FileExistsError("Preserve all prior validation data")
    checked = selfcheck()
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    sources = [
        Path(__file__),
        *[
            Path(__file__).with_name(p)
            for p in (
                "intrinsic_weyl_joint.py",
                "intrinsic_weyl_joint_pilot.py",
                "intrinsic_weyl_joint_audit.py",
                "intrinsic_weyl.py",
                "window_response.py",
                "p14_plane_wave.py",
            )
        ],
        ROOT / "docs/paper/paper_a/intrinsic_weyl/coverage_protocol.md",
    ]
    manifest = {
        "code_version": stamp,
        "selfcheck": checked,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "source_sha256_lf": {
            str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
    }
    destination.mkdir(parents=True)
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    populations = np.random.default_rng(POPULATION_SEED).choice(
        [0.0, 0.25, 0.75], N_CASES
    )
    rows = []
    for case, amplitude in enumerate(populations):
        path, truth = (
            destination / f"case_{case}.npz",
            destination / f"truth_{case}.npz",
        )
        start = time.perf_counter()
        row = {"case": case, "amplitude": float(amplitude), "seed": DATA_SEED + case}
        try:
            row["generation"] = generate_archive(
                path,
                float(amplitude),
                DATA_SEED + case,
                "connected",
                measurement_n=MEASUREMENT_N,
                truth_path=truth,
            )
            result = evaluate(path, truth)
            (destination / f"case_{case}_full.json").write_text(
                json.dumps(result, indent=2) + "\n"
            )
            row.update(result["evaluation"])
            row["confidence_radius"] = result["inference"]["radius"]
            row["largest_tidal_variance"] = result["inference"][
                "largest_tidal_variance"
            ]
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            row["failure"] = f"{type(error).__name__}: {error}"
            row["valid_and_covered"] = False
        row["seconds"] = time.perf_counter() - start
        row["artifact_sha256"] = {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (
                path,
                truth,
                destination / f"case_{case}_full.json",
                destination / f"case_{case}_selection_failed.npz",
            )
            if p.exists()
        }
        rows.append(row)
        (destination / f"case_{case}_result.json").write_text(
            json.dumps(row, indent=2) + "\n"
        )
        print(
            json.dumps(
                {
                    k: v
                    for k, v in row.items()
                    if k not in ("artifact_sha256", "truth_coefficients_xy_positive")
                }
            ),
            flush=True,
        )
    count = sum(r["valid_and_covered"] for r in rows)
    interval = [
        float(beta.ppf(0.025, count, N_CASES - count + 1)) if count else 0.0,
        float(beta.ppf(0.975, count + 1, N_CASES - count)) if count < N_CASES else 1.0,
    ]
    summary = manifest | {
        "cases": rows,
        "valid_covered": count,
        "n": N_CASES,
        "exact_95_two_sided_interval": interval,
        "coverage_gate": "SUPPORTED" if interval[0] >= 0.95 else "NOT_SUPPORTED",
        "scope": (
            "End-to-end validity/coverage in the frozen iid constant-wave mixture; "
            "not accuracy or varying-profile validation"
        ),
    }
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k not in ("cases", "source_sha256_lf")}
        ),
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/intrinsic_weyl_coverage_v1"
    )
    args = parser.parse_args()
    if args.check:
        print(json.dumps(selfcheck()))
    else:
        run(args.output)
