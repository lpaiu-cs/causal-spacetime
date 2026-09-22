"""Curvature-square physical selfcheck and labelled post-result diagnostics."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from intrinsic_weyl import WEYL
from intrinsic_weyl_gradient_diagnostic import corrections, local_tensor
from intrinsic_weyl_joint import (
    confidence_ball,
    fit_joint,
    model,
    read_order_clocks,
    tidal_distance,
)
from intrinsic_weyl_joint_pilot import ROOT
from intrinsic_weyl_square import correction, selfcheck
from rotating_plane_wave import RotatingWave


def proper_squares(wave, points):
    result = np.zeros((len(points), len(points)))
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            p, q = points[i], points[j]
            if p[0] > q[0]:
                p, q = q, p
            du = q[0] - p[0]
            if du == 0:
                raise ValueError("The fixed diagnostic has equal-u pairs")
            value = du * (q[1] - p[1] + wave.action(du, p[2:], q[2:])[0])
            result[i, j] = result[j, i] = value
    return result


def true_coefficients(amplitude):
    tensor = local_tensor(np.diag([amplitude, -amplitude]))
    c = np.linalg.lstsq(WEYL.reshape(10, -1).T, tensor.ravel(), rcond=None)[0]
    assert np.max(abs(np.einsum("k,kabcd->abcd", c, WEYL) - tensor)) < 1e-12
    return np.r_[c, np.zeros(24)]


def physics_check():
    points = np.array(
        [
            [0, 0.5, 0, 0],
            [0.22, 0.20, 0.08, -0.06],
            [0.46, 0.02, -0.06, 0.07],
            [0.71, -0.18, 0.09, 0.04],
            [1, -0.5, 0, 0],
        ]
    )
    pairs = np.column_stack(np.triu_indices(len(points), 1))
    rows = []
    for scale in (1, 0.5, 0.25, 0.125):
        wave = RotatingWave(0.75 * scale**2, 0)
        normal, c = wave.normal_coordinates(points), true_coefficients(wave.amplitude)
        base = model(normal, c, pairs)[0]
        world = correction(normal, c, pairs, volume_clock=False)[0]
        total = correction(normal, c, pairs)[0]
        exact = proper_squares(wave, points)[pairs[:, 0], pairs[:, 1]]
        volumes = wave.volumes(points)
        clocks = -np.sqrt(volumes[pairs[:, 0], pairs[:, 1]] / volumes[0, -1])
        row = {
            "scale": scale,
            "world_before": float(np.max(abs(exact - base))),
            "world_after": float(np.max(abs(exact - base - world))),
            "volume_before": float(np.max(abs(clocks - base))),
            "volume_after": float(np.max(abs(clocks - base - total))),
        }
        rows.append(row)
        assert (
            row["world_after"] < row["world_before"]
            and row["volume_after"] < row["volume_before"]
        )
    for name in ("world_after", "volume_after"):
        assert all(rows[i][name] / rows[i + 1][name] > 30 for i in (0, 1))
    return {"derivatives": selfcheck(), "constant_wave_physical_scaling": rows}


def run(folder, output):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit this diagnostic before execution")
    if output.exists():
        raise FileExistsError("Preserve previous diagnostic outputs")
    checked = physics_check()
    rows, inputs = [], []
    for case in (27, 29, 49):
        archive, truth, full = (
            folder / f"case_{case}.npz",
            folder / f"truth_{case}.npz",
            folder / f"case_{case}_full.json",
        )
        inputs.extend((archive, truth, full))
        original = json.loads(full.read_text())
        counts, mask, covariance, _, _ = read_order_clocks(archive)
        with np.load(truth, allow_pickle=False) as saved:
            points = saved["points"]
            a, w = float(saved["amplitude"]), float(saved["rotation_rate"])
        assert (a, w) == (0.75, -0.5)
        expected = original["evaluation"]["truth_coefficients_xy_positive"]
        wave = RotatingWave(a, w)
        normal, c = wave.normal_coordinates(points), true_coefficients(a)
        pairs = np.column_stack(np.triu_indices(len(points), 1))
        world, total = np.zeros_like(counts), np.zeros_like(counts)
        for matrix, volume in ((world, False), (total, True)):
            matrix[pairs[:, 0], pairs[:, 1]] = correction(
                normal, c, pairs, volume_clock=volume
            )[0]
            matrix += matrix.T.copy()
        volumes = wave.volumes(points)
        clocks = -np.sqrt(volumes / volumes[0, -1])
        delta = normal[:, None] - normal[None]
        flat = -(delta[..., 0] ** 2) + np.sum(delta[..., 1:] ** 2, axis=-1)
        fourth, fifth = corrections(normal, a, w)
        variants = {
            "oracle_exact_proper_time": proper_squares(wave, points),
            "oracle_exact_fifth_polynomial": flat + fourth + fifth,
            "oracle_subtract_square_world": clocks - world,
            "oracle_subtract_square_volume": clocks - (total - world),
            "oracle_subtract_both": clocks - total,
            "blind_square_continuum": clocks,
            "blind_square_finite": None,
        }
        row = {"case": case, "original": original["evaluation"], "variants": {}}
        for label, target in variants.items():
            try:
                if target is not None and np.any(target[mask] >= 0):
                    raise ValueError("A diagnostic changed a selected clock sign")
                fit = fit_joint(
                    counts
                    if target is None
                    else counts[0, -1] * np.minimum(target, 0) ** 2,
                    clock_mask=mask,
                    clock_covariance=covariance,
                    model_order=5,
                    square_correction=label.startswith("blind_square"),
                )
                inference = confidence_ball(fit)
                error = tidal_distance(
                    inference["coefficients_xy_positive"], expected
                ) / (2 * a)
                row["variants"][label] = {
                    "relative_tensor_error": error,
                    "largest_tidal_variance": inference["largest_tidal_variance"],
                    "fit": fit,
                }
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
                row["variants"][label] = {"failure": f"{type(error).__name__}: {error}"}
            print(
                json.dumps(
                    {
                        "case": case,
                        "variant": label,
                        **{
                            k: v
                            for k, v in row["variants"][label].items()
                            if k != "fit"
                        },
                    }
                ),
                flush=True,
            )
        rows.append(row)
    sources = (
        [Path(__file__)]
        + [
            Path(__file__).with_name(name)
            for name in (
                "intrinsic_weyl_square.py",
                "intrinsic_weyl_joint.py",
                "intrinsic_weyl.py",
                "intrinsic_weyl_gradient_diagnostic.py",
                "rotating_plane_wave.py",
            )
        ]
        + [
            ROOT / "docs/paper/paper_a/intrinsic_weyl/square_correction.md",
            ROOT
            / "docs/paper/paper_a/window_response/vacuum_generalization_derivation.md",
        ]
    )
    result = {
        "code_version": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "scope": (
            "Post-result diagnosis on selected failures; oracle paths labelled "
            "separately; no independent validation or changed coverage verdict"
        ),
        "selfcheck": checked,
        "cases": rows,
        "source_sha256_lf": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
        "input_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--folder", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "outputs/intrinsic_weyl_square_diagnostic_v1.json",
    )
    args = parser.parse_args()
    if args.check:
        print(json.dumps(physics_check(), indent=2))
    elif args.folder is None:
        parser.error("--folder is required for the post-result diagnostic")
    else:
        run(args.folder, args.output)
