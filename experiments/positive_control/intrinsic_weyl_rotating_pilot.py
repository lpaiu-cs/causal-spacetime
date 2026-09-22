"""Prospective varying-profile finite-order pilot; no accuracy certification."""

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
    confidence_ball,
    fit_joint,
    fit_order_archive,
    gauge,
    selection,
    selfcheck_jet,
    tidal_distance,
)
from intrinsic_weyl_joint_audit import aligned_tensor_error
from intrinsic_weyl_joint_pilot import ROOT
from rotating_plane_wave import RotatingWave, selfcheck

CASES = [
    (0, 0),
    (0, 0.5),
    (0.25, 0),
    (0.75, 0),
    (0.25, 0.5),
    (0.75, 0.5),
    (0.25, -0.5),
    (0.75, -0.5),
]
SELECTION_N, MEASUREMENT_N, CHUNK = 4096, 262144, 4096


def generate(path, truth_path, amplitude, rotation_rate, seed):
    wave = RotatingWave(amplitude, rotation_rate)
    rng = np.random.default_rng([seed, 0])
    points = np.vstack(([0, 0.5, 0, 0], wave.sample(SELECTION_N, rng), [1, -0.5, 0, 0]))
    # Labels are randomized independently of the coordinates; the order routine
    # preserves those labels. Only supplied outer tips have distinguished IDs.
    permutation = np.r_[
        0, rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
    ]
    points = points[permutation]
    order = wave.order(points)
    try:
        ids, _ = selection(order, design="connected")
    except ValueError:
        np.savez_compressed(
            path.with_name(path.stem + "_selection_failed.npz"),
            selection_order=order,
            selection_design="connected",
        )
        raise
    landmarks = points[ids]
    np.savez_compressed(
        truth_path, points=landmarks, amplitude=amplitude, rotation_rate=rotation_rate
    )
    after = np.empty((len(ids), MEASUREMENT_N // 8), dtype=np.uint8)
    before = np.empty_like(after)
    rng = np.random.default_rng([seed, 1])
    for start in range(0, MEASUREMENT_N, CHUNK):
        events = wave.sample(CHUNK, rng)
        future, past = wave.relations(landmarks[:, None], events[None])
        after[:, start // 8 : (start + CHUNK) // 8] = np.packbits(future, axis=1)
        before[:, start // 8 : (start + CHUNK) // 8] = np.packbits(past, axis=1)
    np.savez_compressed(
        path,
        selection_order=order,
        selection_design="connected",
        selected_ids=ids,
        after=after,
        before=before,
        measurement_n=MEASUREMENT_N,
    )
    return {
        "landmarks": len(ids) - 2,
        "selection_n": SELECTION_N,
        "measurement_n": MEASUREMENT_N,
    }


def evaluate(path, truth_path, *, model_order=4, square_correction=False):
    fit, clocks = fit_order_archive(
        path,
        with_clocks=True,
        model_order=model_order,
        square_correction=square_correction,
    )
    inference = confidence_ball(fit)
    # The estimator has returned before oracle inputs are read. The fit remains
    # a general ten-component Weyl fit, with no profile or type-N prior.
    with np.load(truth_path, allow_pickle=False) as saved:
        points = saved["points"]
        amplitude, rotation = float(saved["amplitude"]), float(saved["rotation_rate"])
    wave = RotatingWave(amplitude, rotation)
    normal = wave.normal_coordinates(points)
    canonical = gauge(normal)
    expected = aligned_tensor_error(
        {"coefficients": np.zeros(10), "coordinates_inferred_from_counts": canonical},
        normal,
        amplitude,
    )["aligned_truth_coefficients"]
    counts, mask, covariance, _, _ = clocks
    volume, coarse = wave.volumes(points), wave.volumes(points, 32)
    refinement = float(np.max(abs(coarse[mask] / volume[mask] - 1)))
    if refinement >= 1e-10:
        raise ArithmeticError("Rotating interval quadrature did not stabilize")
    continuum = fit_joint(
        counts[0, -1] * volume / volume[0, -1],
        clock_mask=mask,
        clock_covariance=covariance,
        model_order=model_order,
        square_correction=square_correction,
    )
    mean = confidence_ball(continuum)
    bias = tidal_distance(mean["coefficients_xy_positive"], expected)
    error = tidal_distance(inference["coefficients_xy_positive"], expected)
    return {
        "fit": fit,
        "inference": inference,
        "continuum_fit": continuum,
        "evaluation": {
            "truth_coefficients_xy_positive": expected,
            "truth_norm": 2 * amplitude,
            "full_tidal_error_mod_reflection": error,
            "relative_full_tidal_error": error / (2 * amplitude) if amplitude else None,
            "continuum_bias_mod_reflection": bias,
            "relative_continuum_bias": bias / (2 * amplitude) if amplitude else None,
            "bias_condition_met": bool(bias <= 0.02 * (2 * amplitude) + 1e-10),
            "nominal_ball_contains_truth": bool(error <= inference["radius"]),
            "confidence_radius": inference["radius"],
            "largest_tidal_variance": inference["largest_tidal_variance"],
            "quadrature_relative_change": refinement,
        },
    }


def run(
    destination,
    *,
    model_order=4,
    data_seed=2026092600,
    protocol="rotating_pilot_protocol.md",
    square_correction=False,
):
    if square_correction and model_order != 5:
        raise ValueError("Curvature-square pilot requires the fifth-degree model")
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit protocol and code before new rotating orders")
    if destination.exists():
        raise FileExistsError("Preserve all previous pilot outputs")
    checked = selfcheck()
    if model_order == 5:
        checked["jet_inverse"] = selfcheck_jet()
    if square_correction:
        from intrinsic_weyl_square_diagnostic import physics_check

        checked["curvature_square"] = physics_check()
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    sources = [
        Path(__file__),
        *[
            Path(__file__).with_name(name)
            for name in (
                "rotating_plane_wave.py",
                "intrinsic_weyl_joint.py",
                "intrinsic_weyl.py",
                "intrinsic_weyl_joint_audit.py",
                "intrinsic_weyl_joint_pilot.py",
                "window_response.py",
                "p14_plane_wave.py",
            )
        ],
        ROOT / "docs/paper/paper_a/intrinsic_weyl" / protocol,
        ROOT / "docs/paper/paper_a/intrinsic_weyl/rotating_generator.md",
    ]
    if square_correction:
        sources.extend(
            [
                Path(__file__).with_name("intrinsic_weyl_square.py"),
                Path(__file__).with_name("intrinsic_weyl_square_diagnostic.py"),
                Path(__file__).with_name("intrinsic_weyl_gradient_diagnostic.py"),
                ROOT / "docs/paper/paper_a/intrinsic_weyl/square_correction.md",
                ROOT
                / "docs/paper/paper_a/window_response"
                / "vacuum_generalization_derivation.md",
            ]
        )
    manifest = {
        "code_version": stamp,
        "selfcheck": checked,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "cases_fixed_before_data": CASES,
        "selection_n": SELECTION_N,
        "measurement_n": MEASUREMENT_N,
        "model_order": model_order,
        "data_seed": data_seed,
        "source_sha256_lf": {
            str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
    }
    if square_correction:
        manifest["curvature_square_correction"] = True
    destination.mkdir(parents=True)
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    rows = []
    for i, (a, w) in enumerate(CASES):
        started = time.perf_counter()
        row = {"case": i, "amplitude": a, "rotation_rate": w, "seed": data_seed + i}
        path, truth = destination / f"case_{i}.npz", destination / f"truth_{i}.npz"
        try:
            row["generation"] = generate(path, truth, a, w, row["seed"])
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            row["failure"] = f"{type(error).__name__}: {error}"
        if "failure" not in row:
            variants = [(model_order, square_correction, "")]
            if square_correction:
                variants.append((5, False, "fifth"))
            elif model_order == 5:
                variants.append((4, False, "quartic"))
            for order, corrected, reference in variants:
                label = f"{reference}_full" if reference else "full"
                try:
                    result = evaluate(
                        path, truth, model_order=order, square_correction=corrected
                    )
                    (destination / f"case_{i}_{label}.json").write_text(
                        json.dumps(result, indent=2) + "\n"
                    )
                    if reference:
                        row[f"{reference}_reference"] = result["evaluation"]
                    else:
                        row.update(result["evaluation"])
                except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
                    key = f"{reference}_reference_failure" if reference else "failure"
                    row[key] = f"{type(error).__name__}: {error}"
        row["seconds"] = time.perf_counter() - started
        row["artifact_sha256"] = {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in destination.glob(f"*_{i}*.npz")
        }
        full = destination / f"case_{i}_full.json"
        if full.exists():
            row["artifact_sha256"][full.name] = hashlib.sha256(
                full.read_bytes()
            ).hexdigest()
        reference_label = "fifth" if square_correction else "quartic"
        baseline_path = destination / f"case_{i}_{reference_label}_full.json"
        if baseline_path.exists():
            row["artifact_sha256"][baseline_path.name] = hashlib.sha256(
                baseline_path.read_bytes()
            ).hexdigest()
        rows.append(row)
        (destination / f"case_{i}_result.json").write_text(
            json.dumps(row, indent=2) + "\n"
        )
        print(
            json.dumps(
                {
                    k: v
                    for k, v in row.items()
                    if k not in ("truth_coefficients_xy_positive", "artifact_sha256")
                }
            ),
            flush=True,
        )
    summary = manifest | {
        "cases": rows,
        "scope": (
            "Prospective varying-profile development pilot; "
            "no accuracy or coverage verdict"
        ),
    }
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Written {destination / 'summary.json'}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/intrinsic_weyl_rotating_pilot_v1"
    )
    parser.add_argument("--model-order", type=int, choices=(4, 5), default=4)
    parser.add_argument("--data-seed", type=int, default=2026092600)
    parser.add_argument("--protocol", default="rotating_pilot_protocol.md")
    parser.add_argument("--square-correction", action="store_true")
    args = parser.parse_args()
    run(
        args.output,
        model_order=args.model_order,
        data_seed=args.data_seed,
        protocol=args.protocol,
        square_correction=args.square_correction,
    )
