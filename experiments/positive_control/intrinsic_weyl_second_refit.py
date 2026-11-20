"""Post hoc refit of the Ibrahim-Bombelli benchmark archives with the sixth-degree jet.

Descriptive only. The archives were generated and judged under the frozen
fifth-degree protocol; refitting them with a model chosen afterwards is a re-
analysis, not a test, and no outcome here changes the recorded verdicts. The
order-only boundary is kept: clocks are read from the saved order bits, both
jets are fitted, and only then is the generator truth loaded for evaluation.
"""

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CASES = (("wave", 0), ("wave", 1), ("flat", 2), ("flat", 3))
LEVELS = (1048576, 4194304, 16777216, 67108864)


def refit(archive, truth_path):
    import intrinsic_weyl_joint as joint
    from intrinsic_weyl_gpu_counts import clock_second_moment
    from intrinsic_weyl_joint import (
        confidence_ball,
        fit_joint,
        gauge,
        read_order_clocks,
        tidal_distance,
    )

    # The exact integer-weighted second moment on GPU, as in the benchmark worker.
    joint._clock_second_moment = clock_second_moment
    begin = time.time()
    counts, mask, covariance, _, n = read_order_clocks(archive)
    clock_seconds = time.time() - begin
    fits = {}
    for order in (5, 6):
        start = time.time()
        fit = fit_joint(
            counts,
            clock_mask=mask,
            clock_covariance=covariance,
            model_order=order,
            square_correction=True,
        )
        fits[order] = (fit, time.time() - start)
    # Truth is read only after both order-only fits exist.
    from intrinsic_weyl_bombelli_benchmark import build_wave
    from intrinsic_weyl_joint_audit import aligned_tensor_error

    with np.load(truth_path, allow_pickle=False) as saved:
        points, kind = saved["points"], str(saved["kind"].item())
        amplitude = str(saved["amplitude"].item())
    wave, value = build_wave(kind, amplitude), float(amplitude)
    normal = wave.normal_coordinates(points)
    expected = np.asarray(
        aligned_tensor_error(
            {
                "coefficients": np.zeros(10),
                "coordinates_inferred_from_counts": gauge(normal),
            },
            normal,
            value,
        )["aligned_truth_coefficients"]
    )
    rows = {}
    for order, (fit, seconds) in fits.items():
        c = np.asarray(fit["coefficients"])
        s = np.asarray(fit["nominal_component_covariance"])
        ball = confidence_ball(fit)
        row = {
            "q": float(c @ np.linalg.solve(s, c)),
            "tidal_norm_estimate": ball["tidal_norm"],
            "statistical_radius_95": ball["statistical_radius"],
            "radius_95_with_allowance": ball["radius"],
            "electric_eigenvalues": fit["electric_eigenvalues"],
            "projected_curvature_rank": fit["projected_curvature_rank"],
            "nuisance_rank": fit["nuisance_rank"],
            "nuisance_parameters": fit["nuisance_parameters"],
            "optimizer_success": fit["optimizer_success"],
            "curved_cost": fit["curved_cost"],
            "measured_intervals": fit["measured_intervals"],
            "fit_seconds": round(seconds, 1),
        }
        error = tidal_distance(ball["coefficients_xy_positive"], expected)
        row["full_tidal_error"] = error
        if value:
            norm = 2 * value
            row["relative_error"] = error / norm
            row["relative_statistical_radius_95"] = ball["statistical_radius"] / norm
            row["relative_radius_95_with_allowance"] = ball["radius"] / norm
        rows[f"order_{order}"] = row
    from scipy.stats import chi2

    # Whitened costs are half squared norms, so twice the drop is the nominal
    # likelihood-ratio statistic for the 42 second-derivative parameters.
    ratio = 2 * (rows["order_5"]["curved_cost"] - rows["order_6"]["curved_cost"])
    return {
        "measurement_n": n,
        "kind": kind,
        "truth_norm": 2 * value,
        "clock_read_seconds": round(clock_seconds, 1),
        "second_derivative_likelihood_ratio": ratio,
        "second_derivative_nominal_p": float(chi2.sf(ratio, 42)),
        "fits": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--levels", type=int, nargs="*", default=list(LEVELS))
    args = parser.parse_args()
    commit = subprocess.run(
        ["git", "rev-parse", "--short=7", "HEAD"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    ).stdout.strip()
    dirty = bool(
        subprocess.run(
            ["git", "status", "--porcelain"],
            capture_output=True,
            text=True,
            cwd=ROOT,
            check=True,
        ).stdout.strip()
    )
    result = {
        "analysis": "post hoc sixth-degree refit of saved archives; descriptive",
        "code_version": commit,
        "worktree_clean": not dirty,
        "rows": [],
    }
    for level in args.levels:
        for kind, case in CASES:
            row = refit(
                args.benchmark / f"case_{case}_n{level}.npz",
                args.benchmark / f"truth_{case}.npz",
            )
            row["case"] = case
            result["rows"].append(row)
            args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
            o5, o6 = row["fits"]["order_5"], row["fits"]["order_6"]
            extra = (
                f" err5={o5['relative_error']:.3f} err6={o6['relative_error']:.3f}"
                f" rad5={o5['relative_statistical_radius_95']:.3f}"
                f" rad6={o6['relative_statistical_radius_95']:.3f}"
                if kind == "wave"
                else ""
            )
            print(
                f"n={level} case={case} {kind}: q5={o5['q']:.2f} q6={o6['q']:.2f}"
                f" lr={row['second_derivative_likelihood_ratio']:.1f}{extra}",
                flush=True,
            )


if __name__ == "__main__":
    main()
