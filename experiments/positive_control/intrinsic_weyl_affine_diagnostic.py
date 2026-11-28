"""Fixed eight-case diagnostic of a pilot-frozen affine count inverse.

No events are generated and no confidence or power verdict is produced.
All maps are saved before any truth module is imported or truth file opened.
"""

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from intrinsic_weyl_gradient_information import affine_map, information, linearization
from intrinsic_weyl_joint import (
    TIDAL_ROOT,
    axis_signs,
    gauge,
    model,
    read_order_clocks,
    tidal_distance,
)

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as output:
        output.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run(folder, destination):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit code and the fixed diagnostic plan before running")
    destination.mkdir(parents=True, exist_ok=False)
    sources = list((ROOT / "experiments/positive_control").glob("intrinsic_weyl*.py"))
    sources += [
        ROOT / "experiments/positive_control" / name
        for name in (
            "rotating_plane_wave.py",
            "window_response.py",
            "p14_plane_wave.py",
        )
    ]
    sources.append(
        ROOT / "docs/paper/paper_a/intrinsic_weyl/local_consistency_scope.md"
    )
    source_hashes = {
        p.relative_to(ROOT).as_posix(): hashlib.sha256(
            p.read_bytes().replace(b"\r\n", b"\n")
        ).hexdigest()
        for p in sources
    }
    rows = []
    inputs = {}
    for case in range(8):
        archive = folder / f"case_{case}_nested.npz"
        saved = folder / f"case_{case}_nested_full.json"
        inputs[archive.name], inputs[saved.name] = digest(archive), digest(saved)
        # Only this subobject is supplied to the inverse; the enclosing file's
        # continuum fit, truth coefficients and evaluation fields are not inputs.
        fit = json.loads(saved.read_text())["fit"]
        row = {
            "case": case,
            "order_archive_sha256": inputs[archive.name],
            "saved_wrapper_sha256": inputs[saved.name],
            "fit_subobject_sha256": hashlib.sha256(
                json.dumps(fit, sort_keys=True, allow_nan=False).encode()
            ).hexdigest(),
        }
        start = time.perf_counter()
        try:
            counts, mask, covariance, ids, n = read_order_clocks(archive)
            matrix, offset, pairs, details = affine_map(
                fit, counts, mask, covariance, n
            )
            map_path = destination / f"case_{case}_map.npz"
            np.savez(
                map_path, matrix=matrix, offset=offset, pairs=pairs, selected_ids=ids
            )
            row.update(details, map_sha256=digest(map_path))
        except (
            ValueError,
            ArithmeticError,
            AssertionError,
            np.linalg.LinAlgError,
        ) as error:
            row["pure_failure"] = f"{type(error).__name__}: {error}"
        row["pure_seconds"] = time.perf_counter() - start
        rows.append(row)
        print(json.dumps({"phase": "pure", **row}, allow_nan=False), flush=True)
    forbidden = ("rotating_plane_wave", "window_response", "p14_plane_wave")
    if any(name in sys.modules for name in forbidden) or any(
        name.startswith(("intrinsic_weyl_joint_audit", "intrinsic_weyl_joint_pilot"))
        for name in sys.modules
    ):
        raise RuntimeError(
            "An oracle or generator was imported before maps were frozen"
        )
    manifest = {
        "scope": "Eight existing pilot cases; fixed affine-map diagnostic only",
        "code_version": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256_lf": source_hashes,
        "pilot_manifest_sha256": digest(folder / "manifest.json"),
        "input_sha256": inputs,
        "all_maps_frozen_before_truth": True,
        "oracle_or_generator_imported_in_pure_stage": False,
        "fresh_events": 0,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "variance_scope": "Same-pilot empirical scores; not independent calibration",
        "cases": rows,
    }
    pure_path = destination / "pure.json"
    write_json(pure_path, manifest)
    pure_hash = digest(pure_path)
    # Truth is confined to evaluation after all maps (including failures) freeze.
    from intrinsic_weyl_joint_audit import aligned_tensor_error
    from rotating_plane_wave import RotatingWave

    for row in rows:
        if "pure_failure" in row:
            continue
        case = row["case"]
        map_path = destination / f"case_{case}_map.npz"
        if digest(map_path) != row["map_sha256"]:
            raise RuntimeError("A frozen map changed before evaluation")
        with np.load(map_path, allow_pickle=False) as saved:
            matrix, offset, pairs = saved["matrix"], saved["offset"], saved["pairs"]
        truth_path = folder / f"truth_{case}_nested.npz"
        row["truth_archive_sha256"] = digest(truth_path)
        start = time.perf_counter()
        try:
            with np.load(truth_path, allow_pickle=False) as truth:
                points = truth["points"]
                amplitude = float(truth["amplitude"])
                rotation = float(truth["rotation_rate"])
            wave = RotatingWave(amplitude, rotation)
            normal = wave.normal_coordinates(points)
            expected = aligned_tensor_error(
                {
                    "coefficients": np.zeros(10),
                    "coordinates_inferred_from_counts": gauge(normal),
                },
                normal,
                amplitude,
            )["aligned_truth_coefficients"]
            volume, coarse = wave.volumes(points, 64), wave.volumes(points, 32)
            a, b = pairs.T
            refinement = float(
                max(
                    np.max(abs(coarse[a, b] / volume[a, b] - 1)),
                    abs(coarse[0, -1] / volume[0, -1] - 1),
                )
            )
            if not np.isfinite(refinement) or refinement >= 1e-10:
                raise ArithmeticError("Interval-volume quadrature did not stabilize")
            probabilities = volume[a, b] / volume[0, -1]
            if np.any(probabilities <= 0) or np.any(probabilities > 1):
                raise ArithmeticError("Invalid exact interval-volume fractions")
            mean = offset + matrix @ probabilities
            coefficients = np.linalg.solve(TIDAL_ROOT, mean)
            bias = tidal_distance(coefficients, expected)
            truth_norm = float(np.linalg.norm(TIDAL_ROOT @ expected))
            row["evaluation"] = {
                "amplitude": amplitude,
                "rotation_rate": rotation,
                "truth_norm": truth_norm,
                "truth_coefficients_xy_positive": expected,
                "conditional_mean_tidal": mean.tolist(),
                "conditional_mean_coefficients_xy_positive": coefficients.tolist(),
                "mean_bias_mod_reflection": bias,
                "relative_mean_bias": bias / truth_norm if truth_norm else None,
                "old_bias_condition_met": bool(bias <= 0.02 * truth_norm + 1e-10),
                "quadrature_relative_change": refinement,
            }
        except (
            ValueError,
            ArithmeticError,
            AssertionError,
            np.linalg.LinAlgError,
        ) as error:
            row["evaluation_failure"] = f"{type(error).__name__}: {error}"
        row["evaluation_seconds"] = time.perf_counter() - start
        print(
            json.dumps(
                {
                    "phase": "evaluation",
                    "case": case,
                    "evaluation": row.get("evaluation"),
                    "failure": row.get("evaluation_failure"),
                },
                allow_nan=False,
            ),
            flush=True,
        )
    if digest(pure_path) != pure_hash or any(
        digest(folder / name) != value for name, value in inputs.items()
    ):
        raise RuntimeError("A pure artifact or pilot input changed during evaluation")
    if any(
        hashlib.sha256((ROOT / name).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        != value
        for name, value in source_hashes.items()
    ):
        raise RuntimeError("Diagnostic sources changed during the run")
    for row in rows:
        if (
            "map_sha256" in row
            and digest(destination / f"case_{row['case']}_map.npz") != row["map_sha256"]
        ):
            raise RuntimeError("A frozen map changed during evaluation")
        if (
            "truth_archive_sha256" in row
            and digest(folder / f"truth_{row['case']}_nested.npz")
            != row["truth_archive_sha256"]
        ):
            raise RuntimeError("A truth archive changed during evaluation")
    manifest["pure_manifest_sha256"] = pure_hash
    manifest["bias_condition_met_cases"] = sum(
        row.get("evaluation", {}).get("old_bias_condition_met", False) for row in rows
    )
    manifest["cases_attempted"] = len(rows)
    manifest["pure_maps_returned"] = sum("map_sha256" in row for row in rows)
    manifest["evaluations_returned"] = sum("evaluation" in row for row in rows)
    manifest["result_scope"] = (
        "Conditional-mean and empirical-range diagnostic; no finite-sample coverage, "
        "power, recovery or universal bias certification"
    )
    write_json(destination / "summary.json", manifest)


def selfcheck():
    """One deterministic synthetic check of derivatives, gauge and count mapping."""
    k = np.arange(48)
    z = 1 - 2 * (k + 0.5) / len(k)
    phi = k * np.pi * (3 - np.sqrt(5))
    points = np.zeros((50, 4))
    points[0, 0], points[-1, 0] = -0.5, 0.5
    points[1:-1, 0] = 0.27 * ((k % 3) - 1) + 0.015 * np.sin(k)
    points[1:-1, 1:] = (0.12 + 0.08 * np.sin(k * 2.1) ** 2)[:, None] * np.column_stack(
        (np.sqrt(1 - z * z) * np.cos(phi), np.sqrt(1 - z * z) * np.sin(phi), z)
    )
    points = gauge(points)
    difference = points[:, None] - points[None]
    tau2 = difference[..., 0] ** 2 - np.sum(difference[..., 1:] ** 2, axis=-1)
    mask = tau2 > 0.001
    mask[0, -1] = mask[-1, 0] = False
    coefficients = np.linspace(-0.01, 0.01, 34)
    fit = {
        "coordinates_inferred_from_counts": points.tolist(),
        "coefficients": coefficients[:10].tolist(),
        "gradient_coefficients": coefficients[10:].tolist(),
        "model_order": 5,
        "curvature_square_correction": True,
        "optimizer_success": True,
        "boundary_excess": 0,
        "projected_curvature_rank": 10,
    }
    value, jacobian, pairs, dimension = linearization(fit, mask)
    assert np.linalg.matrix_rank(jacobian) == jacobian.shape[1]
    inverse = np.linalg.pinv(jacobian)[dimension : dimension + 10]
    covariance = np.eye(len(pairs)) * 1e-8
    fit.update(
        coordinate_rank=dimension,
        coordinate_parameters=dimension,
        nuisance_rank=dimension + 24,
        nuisance_parameters=dimension + 24,
        nominal_component_covariance=(1e-8 * inverse @ inverse.T).tolist(),
    )
    counts = np.zeros_like(tau2)
    counts[pairs[:, 0], pairs[:, 1]] = 262144 * value**2
    counts += counts.T
    counts[0, -1] = counts[-1, 0] = 262144
    matrix, offset, _, details = affine_map(fit, counts, mask, covariance, 262144)
    signs = axis_signs([np.sign(points[1, 1]), np.sign(points[2, 2]), 1])
    transform = TIDAL_ROOT * signs
    assert (
        np.max(abs(offset + matrix @ value**2 - transform @ coefficients[:10])) < 1e-10
    )
    maximum = 0.0
    for component in (0, 7, 10, 33):
        step = np.eye(34)[component] * 1e-6
        plus = model(points, coefficients + step, pairs, square_correction=True)[0]
        minus = model(points, coefficients - step, pairs, square_correction=True)[0]
        mapped = matrix @ (plus**2 - minus**2) / 2e-6
        expected = transform[:, component] if component < 10 else np.zeros(10)
        maximum = max(maximum, float(np.max(abs(mapped - expected))))
    for landmark, component in ((7, 2), (4, 0)):
        step = np.zeros_like(points)
        step[landmark, component] = 1e-6
        plus = model(points + step, coefficients, pairs, square_correction=True)[0]
        minus = model(points - step, coefficients, pairs, square_correction=True)[0]
        maximum = max(maximum, float(np.max(abs(matrix @ (plus**2 - minus**2) / 2e-6))))
    assert maximum < 1e-5
    quartic = dict(fit)
    del quartic["gradient_coefficients"]
    quartic.update(model_order=4, curvature_square_correction=False)
    _, quartic_jacobian, _, _ = linearization(quartic, mask)
    quartic_inverse = np.linalg.pinv(quartic_jacobian)[dimension : dimension + 10]
    quartic["nominal_component_covariance"] = (
        1e-8 * quartic_inverse @ quartic_inverse.T
    ).tolist()
    old_diagnostic = information(quartic, mask, covariance)
    return {
        "finite_difference_map_max_error": maximum,
        "left_inverse_max_error": details["local_left_inverse_max_error"],
        "covariance_replay_error": details["saved_covariance_relative_replay_error"],
        "quartic_information_replay_error": old_diagnostic[
            "old_covariance_relative_replay_difference"
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--folder", type=Path, default=ROOT / "outputs/intrinsic_weyl_nested_pilot_v1"
    )
    parser.add_argument(
        "--destination",
        type=Path,
        default=ROOT / "outputs/intrinsic_weyl_affine_diagnostic_v1",
    )
    parser.add_argument("--selfcheck", action="store_true")
    args = parser.parse_args()
    if args.selfcheck:
        print(json.dumps(selfcheck()))
    else:
        run(args.folder, args.destination)
