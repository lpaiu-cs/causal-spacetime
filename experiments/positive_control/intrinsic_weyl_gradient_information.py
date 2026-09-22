"""Order-only local information price when all vacuum Weyl gradients are free.

This is a Jacobian diagnostic at the returned quartic fit, not a gradient fit
or a finite-sample accuracy claim. No generator, profile or truth is imported.
"""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from intrinsic_weyl import gradient_basis
from intrinsic_weyl_joint import TIDAL_GRAM, TIDAL_ROOT, fit_order_archive, model


def information(fit, mask, covariance):
    points = np.array(fit["coordinates_inferred_from_counts"])
    pairs = np.column_stack(np.nonzero(np.triu(mask, 1)))
    _, dl, dr, curvature = model(points, np.array(fit["coefficients"]), pairs)
    free = np.ones_like(points, dtype=bool)
    free[[0, -1]] = False
    free[1, 2:] = False
    free[2, 3] = False
    indices = np.full(points.shape, -1)
    indices[free] = np.arange(free.sum())
    geometry = np.zeros((len(pairs), free.sum()))
    for side, derivative in enumerate((dl, dr)):
        for component in range(4):
            index = indices[pairs[:, side], component]
            keep = index >= 0
            geometry[np.flatnonzero(keep), index[keep]] = derivative[keep, component]
    # -C(x,y,x,y)/3 becomes -(nabla_e C)(x,y,x,y)(x+y)^e/12.
    gradient = (
        np.einsum(
            "pe,pk->pek", (points[pairs[:, 0]] + points[pairs[:, 1]]) / 4, curvature
        ).reshape(len(pairs), 40)
        @ gradient_basis()
    )
    scales = np.sqrt(np.diag(covariance))
    values, vectors = np.linalg.eigh(covariance / scales[:, None] / scales)
    keep = values > values[-1] * 1e-8
    whiten = (vectors[:, keep] / np.sqrt(values[keep])).T / scales
    geometry, curvature, gradient = (
        whiten @ x for x in (geometry, curvature, gradient)
    )
    result = {}
    for label, nuisance in (
        ("quartic", geometry),
        ("free_vacuum_gradient", np.column_stack((geometry, gradient))),
    ):
        projected = (
            curvature - nuisance @ np.linalg.lstsq(nuisance, curvature, rcond=1e-11)[0]
        )
        _, singular, vt = np.linalg.svd(projected, full_matrices=False)
        threshold = 1e-9 * np.linalg.norm(curvature, 2)
        rank = int(np.sum(singular > threshold))
        row = {
            "nuisance_columns": nuisance.shape[1],
            "nuisance_rank": int(np.linalg.matrix_rank(nuisance)),
            "projected_curvature_rank": rank,
            "rank_absolute_threshold": float(threshold),
            "projected_curvature_singular_values": singular.tolist(),
        }
        if rank == 10:
            cov = (vt.T / singular**2) @ vt
            row["largest_tidal_variance"] = float(
                np.linalg.eigvalsh(TIDAL_ROOT @ cov @ TIDAL_ROOT.T)[-1]
            )
            row["electric_tidal_variance"] = float(
                np.trace(TIDAL_GRAM[:5, :5] @ cov[:5, :5])
            )
            row["magnetic_tidal_variance"] = float(
                np.trace(TIDAL_GRAM[5:, 5:] @ cov[5:, 5:])
            )
        result[label] = row
    original = float(
        np.linalg.eigvalsh(
            TIDAL_ROOT @ np.array(fit["nominal_component_covariance"]) @ TIDAL_ROOT.T
        )[-1]
    )
    discrepancy = abs(result["quartic"]["largest_tidal_variance"] / original - 1)
    assert discrepancy < 1e-8
    result["old_covariance_relative_replay_difference"] = discrepancy
    if result["free_vacuum_gradient"]["projected_curvature_rank"] == 10:
        result["gradient_to_quartic_ratios"] = {
            key: result["free_vacuum_gradient"][key] / result["quartic"][key]
            for key in (
                "largest_tidal_variance",
                "electric_tidal_variance",
                "magnetic_tidal_variance",
            )
        }
    return result


def run(folder, output):
    root = Path(__file__).resolve().parents[2]
    if output.exists():
        raise FileExistsError("Preserve earlier information diagnostics")
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=root).strip():
        raise RuntimeError("Commit code before recording the information diagnostic")
    rows = []
    # Fixed existing cases: two constant and two rotating, weak and strong.
    for case in (2, 3, 4, 5):
        path = folder / f"case_{case}.npz"
        fit, clocks = fit_order_archive(path, with_clocks=True)
        row = {
            "case": case,
            "order_archive_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "information": information(fit, clocks[1], clocks[2]),
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
    sources = [
        Path(__file__),
        Path(__file__).with_name("intrinsic_weyl_joint.py"),
        Path(__file__).with_name("intrinsic_weyl.py"),
    ]
    result = {
        "scope": (
            "Local Jacobian information at order-only quartic fits; no truth inputs, "
            "no gradient recovery or accuracy claim"
        ),
        "code_version": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        "source_sha256_lf": {
            str(p.relative_to(root)).replace("\\", "/"): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
        "gradient_dimension": 24,
        "cases": rows,
    }
    output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--folder", type=Path, default=Path("outputs/intrinsic_weyl_rotating_pilot_v1")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/intrinsic_weyl_gradient_information_v1.json"),
    )
    args = parser.parse_args()
    run(args.folder, args.output)
