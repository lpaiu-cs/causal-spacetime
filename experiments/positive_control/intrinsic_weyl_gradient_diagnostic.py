"""Post-pilot ORACLE gradient subtraction, never an order-only recovery claim.

The fifth-order world-function term follows Brewin 0903.2087, section 9.
This diagnoses existing varying-profile bias using known generator truth.
"""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from intrinsic_weyl_joint import (
    confidence_ball,
    fit_joint,
    read_order_clocks,
    tidal_distance,
)
from intrinsic_weyl_joint_pilot import ROOT
from rotating_plane_wave import RotatingWave


def local_tensor(transverse):
    brinkmann = np.zeros((4, 4, 4, 4))
    for i in range(2):
        for j in range(2):
            value = -transverse[i, j]
            brinkmann[0, i + 2, 0, j + 2] = value
            brinkmann[i + 2, 0, j + 2, 0] = value
            brinkmann[i + 2, 0, 0, j + 2] = -value
            brinkmann[0, i + 2, j + 2, 0] = -value
    # Normal coordinates are divided by the tip time sqrt(2): M=sqrt(2)*tetrad.
    # Metric and covariant curvature in these coordinates are divided by T^2=2.
    m = np.array(
        [[1, 0, 0, 1], [-1, 0, 0, 1], [0, np.sqrt(2), 0, 0], [0, 0, np.sqrt(2), 0]]
    )
    return 0.5 * np.einsum("ia,jb,kc,ld,ijkl->abcd", m, m, m, m, brinkmann)


def corrections(normal, amplitude, rotation):
    c = local_tensor(np.diag([amplitude, -amplitude]))
    dc_du = local_tensor(
        np.array([[0, 2 * amplitude * rotation], [2 * amplitude * rotation, 0]])
    )
    quartic = -np.einsum("abcd,ia,jb,ic,jd->ij", c, normal, normal, normal, normal) / 3
    u = normal[:, 0] + normal[:, 3]
    fifth = (
        -np.einsum("abcd,ia,jb,ic,jd->ij", dc_du, normal, normal, normal, normal)
        * (u[:, None] + u[None])
        / 12
    )
    np.fill_diagonal(quartic, 0)
    np.fill_diagonal(fifth, 0)
    return quartic, fifth


def selfcheck():
    points = np.array(
        [
            [0, 0.5, 0, 0],
            [0.22, 0.20, 0.08, -0.06],
            [0.46, 0.02, -0.06, 0.07],
            [0.71, -0.18, 0.09, 0.04],
            [1, -0.5, 0, 0],
        ]
    )
    rows = []
    for scale in (1, 0.5, 0.25, 0.125):
        a, w = 0.75 * scale**2, 0.5 * scale
        wave = RotatingWave(a, w)
        normal = wave.normal_coordinates(points)
        delta = normal[:, None] - normal[None]
        flat = -(delta[..., 0] ** 2) + np.sum(delta[..., 1:] ** 2, axis=-1)
        quartic, fifth = corrections(normal, a, w)
        exact = np.zeros_like(flat)
        for i in range(len(points)):
            for j in range(i + 1, len(points)):
                p, q = points[i], points[j]
                length = q[0] - p[0]
                cost = wave.action(length, p[2:], q[2:])[0]
                exact[i, j] = exact[j, i] = length * (q[1] - p[1] + cost)
        e4 = float(np.max(abs(exact - flat - quartic)))
        e5 = float(np.max(abs(exact - flat - quartic - fifth)))
        assert e5 < e4
        rows.append({"scale": scale, "quartic_error": e4, "fifth_order_error": e5})
    assert all(
        rows[i]["fifth_order_error"] / rows[i + 1]["fifth_order_error"] > 12
        for i in range(3)
    )
    return {
        "exact_proper_time_expansion": rows,
        "scope": "Deterministic scaling check; oracle term only",
    }


def run(folder, output):
    if output.exists():
        raise FileExistsError("Preserve earlier diagnostic outputs")
    checked = selfcheck()
    summary = json.loads((folder / "summary.json").read_text())
    rows = []
    for old in summary["cases"]:
        if not old["amplitude"] or not old["rotation_rate"]:
            continue
        case = old["case"]
        row = {
            "case": case,
            "amplitude": old["amplitude"],
            "rotation_rate": old["rotation_rate"],
            "original_relative_continuum_bias": old["relative_continuum_bias"],
        }
        try:
            counts, mask, covariance, _, _ = read_order_clocks(
                folder / f"case_{case}.npz"
            )
            with np.load(folder / f"truth_{case}.npz") as truth:
                points = truth["points"]
            wave = RotatingWave(old["amplitude"], old["rotation_rate"])
            normal = wave.normal_coordinates(points)
            fifth = corrections(normal, old["amplitude"], old["rotation_rate"])[1]
            volumes = wave.volumes(points)
            clocks = -np.sqrt(volumes / volumes[0, -1])
            corrected = np.where(volumes > 0, clocks - fifth, 0)
            if np.any(corrected[mask] >= 0):
                raise ValueError("Oracle correction changed a selected clock sign")
            fitted = fit_joint(
                counts[0, -1] * corrected**2,
                clock_mask=mask,
                clock_covariance=covariance,
            )
            inference = confidence_ball(fitted)
            error = tidal_distance(
                inference["coefficients_xy_positive"],
                old["truth_coefficients_xy_positive"],
            )
            row["oracle_subtracted_relative_bias"] = error / (2 * old["amplitude"])
            row["projected_curvature_rank"] = fitted["projected_curvature_rank"]
            row["oracle_fit"] = fitted
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            row["failure"] = f"{type(error).__name__}: {error}"
        rows.append(row)
        print(
            json.dumps({k: v for k, v in row.items() if k != "oracle_fit"}), flush=True
        )
    result = {
        "scope": (
            "Post-pilot oracle gradient subtraction: uses TRUE coordinates/profile; "
            "NOT operational reconstruction"
        ),
        "selfcheck": checked,
        "cases": rows,
        "original_summary_sha256": hashlib.sha256(
            (folder / "summary.json").read_bytes()
        ).hexdigest(),
        "code_version": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256_lf": hashlib.sha256(
            Path(__file__).read_bytes().replace(b"\r\n", b"\n")
        ).hexdigest(),
    }
    output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--folder", type=Path, default=ROOT / "outputs/intrinsic_weyl_rotating_pilot_v1"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "outputs/intrinsic_weyl_gradient_diagnostic_v1.json",
    )
    args = parser.parse_args()
    if args.check:
        print(json.dumps(selfcheck()))
    else:
        run(args.folder, args.output)
