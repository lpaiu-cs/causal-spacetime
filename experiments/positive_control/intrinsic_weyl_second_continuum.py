"""Continuum diagnostic: does the sixth-degree jet remove the half-cycle bias?

Development diagnostic, not a test with a pre-specified criterion. Everything is
deterministic: exact plane-wave geometry from profile_plane_wave, the landmark
configurations and clock covariances of the Ibrahim-Bombelli benchmark run, and
no sampling noise. Two parts:

1. World-function check. On the exact normal coordinates of the benchmark
   landmarks, the exact squared interval of the generator (Jacobi propagator)
   is compared with the curvature terms of the fourth-, fifth- and sixth-degree
   models by linear least squares. For a profile quadratic in u every
   curvature derivative beyond the second vanishes, so at linear order the
   sixth-degree model must be exact and the fitted midpoint tensor must equal
   the one of the constant profile with the same midpoint amplitude.
2. Continuum volume fit. The estimator itself (interval volumes as clocks, the
   run's clock mask and covariance, curvature-square terms) is applied with the
   fifth- and sixth-degree jets to profiles with the same midpoint amplitude.
"""

import argparse
import json
import platform
import subprocess
import sys
import time
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PROFILES = ("constant", "linear", "quadratic", "half_cycle")


def build(profile, amplitude):
    """Profiles on u in [0, 1] with a(1/2) = amplitude."""
    from profile_plane_wave import PRECISION, ProfileWave, decimal_pi, half_cycle

    with localcontext() as context:
        context.prec = PRECISION
        a = Decimal(amplitude)
        if profile == "constant":
            return ProfileWave([a], abs(float(a)), profile)
        if profile == "linear":
            # a(u) = 2 a u: slope 2a, same midpoint value.
            return ProfileWave([Decimal(0), 2 * a], 2 * abs(float(a)), profile)
        if profile == "quadratic":
            # Same value and second derivative at u = 1/2 as a sin(pi u).
            half = decimal_pi() ** 2 * a / 2
            taylor = [a - half / 4, half, -half]
            return ProfileWave(taylor, abs(float(a)), profile)
        if profile == "half_cycle":
            return half_cycle(amplitude)
    raise ValueError("Unknown profile")


def oriented_pairs(points):
    """All landmark pairs with distinct u, ordered so that u increases."""
    i, j = np.triu_indices(len(points), 1)
    keep = abs(points[i, 0] - points[j, 0]) > 1e-6
    i, j = i[keep], j[keep]
    swap = points[i, 0] > points[j, 0]
    return np.column_stack((np.where(swap, j, i), np.where(swap, i, j)))


def exact_squared_interval(wave, points, pairs):
    """Omega = 2 du (dv + S) in the Brinkmann chart, S the transverse action."""
    left, right = points[pairs[:, 0]], points[pairs[:, 1]]
    cost, _, unusable = wave._pair(left, right)
    if np.any(unusable):
        raise ArithmeticError("Unusable Jacobi block among landmark pairs")
    return 2 * (right[:, 0] - left[:, 0]) * (right[:, 1] - left[:, 1] + cost)


def world_function_check(points, amplitudes):
    from intrinsic_weyl import ETA
    from intrinsic_weyl_joint import TIDAL_ROOT, model

    pairs = oriented_pairs(points)
    flat_wave = build("constant", "0")
    normal = flat_wave.normal_coordinates(points)
    delta = normal[pairs[:, 1]] - normal[pairs[:, 0]]
    flat_value = np.einsum("pa,ab,pb->p", delta, ETA, delta)
    scale = float(
        np.median(exact_squared_interval(flat_wave, points, pairs) / flat_value)
    )
    if abs(scale - 2) > 1e-9:
        raise ArithmeticError("Normal coordinates are not twice-Minkowski scaled")
    rows = []
    for amplitude in amplitudes:
        reference = None
        for profile in PROFILES:
            wave = build(profile, amplitude)
            normal = wave.normal_coordinates(points)
            delta = normal[pairs[:, 1]] - normal[pairs[:, 0]]
            target = exact_squared_interval(wave, points, pairs) / scale - np.einsum(
                "pa,ab,pb->p", delta, ETA, delta
            )
            row = {"amplitude": amplitude, "profile": profile}
            for order, size in ((4, 10), (5, 34), (6, 76)):
                design = model(normal, np.zeros(size), pairs)[3]
                solution, *_ = np.linalg.lstsq(design, target, rcond=None)
                residual = target - design @ solution
                row[f"order_{order}"] = {
                    "weyl": solution[:10].tolist(),
                    "relative_residual": float(
                        np.linalg.norm(residual) / np.linalg.norm(target)
                    ),
                    "rank": int(np.linalg.matrix_rank(design)),
                    "columns": size,
                }
            if profile == "constant":
                reference = np.array(row["order_4"]["weyl"])
            norm = float(np.linalg.norm(TIDAL_ROOT @ reference))
            for order in (4, 5, 6):
                fitted = np.array(row[f"order_{order}"]["weyl"])
                row[f"order_{order}"]["relative_weyl_error"] = float(
                    np.linalg.norm(TIDAL_ROOT @ (fitted - reference)) / norm
                )
            rows.append(row)
    return {"pairs": len(pairs), "normal_metric_scale": scale, "rows": rows}


def continuum_fits(archive, truth, amplitude):
    from intrinsic_weyl_joint import (
        confidence_ball,
        fit_joint,
        gauge,
        read_order_clocks,
        tidal_distance,
    )
    from intrinsic_weyl_joint_audit import aligned_tensor_error

    started = time.time()
    counts, mask, covariance, _, n = read_order_clocks(archive)
    clock_seconds = time.time() - started
    with np.load(truth, allow_pickle=False) as saved:
        points = saved["points"]
    pairs = np.column_stack(np.nonzero(np.triu(mask, 1)))
    rows = []
    for profile in PROFILES:
        wave = build(profile, amplitude)
        volume = wave.volumes(points)
        usable = volume[pairs[:, 0], pairs[:, 1]] > 0
        local = np.zeros_like(mask)
        kept = pairs[usable]
        local[kept[:, 0], kept[:, 1]] = local[kept[:, 1], kept[:, 0]] = True
        weights = covariance[np.ix_(usable, usable)]
        normal = wave.normal_coordinates(points)
        expected = np.asarray(
            aligned_tensor_error(
                {
                    "coefficients": np.zeros(10),
                    "coordinates_inferred_from_counts": gauge(normal),
                },
                normal,
                float(amplitude),
            )["aligned_truth_coefficients"]
        )
        row = {"profile": profile, "dropped_clocks": int((~usable).sum())}
        for order in (5, 6):
            begin = time.time()
            fit = fit_joint(
                counts[0, -1] * volume / volume[0, -1],
                clock_mask=local,
                clock_covariance=weights,
                model_order=order,
                square_correction=True,
            )
            ball = confidence_ball(fit)
            truth_norm = 2 * float(amplitude)
            row[f"order_{order}"] = {
                "relative_continuum_bias": tidal_distance(
                    ball["coefficients_xy_positive"], expected
                )
                / truth_norm,
                "relative_statistical_radius_95": ball["statistical_radius"]
                / truth_norm,
                "relative_radius_95_with_allowance": ball["radius"] / truth_norm,
                "continuum_tidal_norm_over_truth": ball["tidal_norm"] / truth_norm,
                "projected_curvature_rank": fit["projected_curvature_rank"],
                "nuisance_rank": fit["nuisance_rank"],
                "nuisance_parameters": fit["nuisance_parameters"],
                "optimizer_success": fit["optimizer_success"],
                "seconds": round(time.time() - begin, 1),
            }
        rows.append(row)
    return {
        "archive": str(archive),
        "measurement_n": n,
        "landmarks": len(points) - 2,
        "clocks": len(pairs),
        "clock_read_seconds": round(clock_seconds, 1),
        "rows": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", type=Path)
    parser.add_argument("--level", type=int, default=16777216)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude
    from intrinsic_weyl_joint import selfcheck_second

    amplitude = str(benchmark_amplitude())
    commit = subprocess.run(
        ["git", "rev-parse", "--short=7", "HEAD"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    ).stdout.strip()
    result = {
        "diagnostic": "sixth-degree jet, continuum; no sampling, no criterion",
        "code_version": commit,
        "worktree_clean": not dirty,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "selfcheck_second": selfcheck_second(),
        "benchmark_amplitude": amplitude,
        "world_function": world_function_check(
            np.load(args.benchmark / "truth_0.npz")["points"], ["0.001", "0.01"]
        ),
        "continuum": [
            continuum_fits(
                args.benchmark / f"case_{case}_n{args.level}.npz",
                args.benchmark / f"truth_{case}.npz",
                amplitude,
            )
            for case in (0, 1)
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["world_function"]["rows"][-4:], indent=1)[:4000])
    for block in result["continuum"]:
        for row in block["rows"]:
            print(
                row["profile"],
                {
                    k: {
                        m: round(v, 4) if isinstance(v, float) else v
                        for m, v in row[k].items()
                        if m.startswith("relative") or m.endswith("rank")
                    }
                    for k in ("order_5", "order_6")
                },
            )


if __name__ == "__main__":
    sys.exit(main())
