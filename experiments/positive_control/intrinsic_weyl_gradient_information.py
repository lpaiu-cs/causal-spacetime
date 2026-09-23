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
from intrinsic_weyl_joint import (
    TIDAL_GRAM,
    TIDAL_ROOT,
    axis_signs,
    confidence_ball,
    fit_order_archive,
    model,
)


def linearization(fit, mask):
    """Reuse the fitted model and its gauge, without any oracle inputs."""
    points = np.array(fit["coordinates_inferred_from_counts"])
    pairs = np.column_stack(np.nonzero(np.triu(mask, 1)))
    coefficients = np.r_[fit["coefficients"], fit.get("gradient_coefficients", [])]
    value, dl, dr, curvature = model(
        points,
        coefficients,
        pairs,
        square_correction=fit.get("curvature_square_correction", False),
    )
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
    return value, np.column_stack((geometry, curvature)), pairs, int(free.sum())


def whitening(covariance):
    scales = np.sqrt(np.diag(covariance))
    if not np.isfinite(covariance).all() or np.any(scales <= 0):
        raise ValueError("Finite covariance with positive variances required")
    values, vectors = np.linalg.eigh(covariance / scales[:, None] / scales)
    keep = values > values[-1] * 1e-8
    return (vectors[:, keep] / np.sqrt(values[keep])).T / scales


def affine_map(fit, counts, mask, covariance, n):
    """Freeze a one-step count map; its conditional mean bias is still unknown."""
    confidence_ball(fit)  # Reuse rank, optimizer, frame and covariance guards.
    if fit.get("model_order") != 5 or not fit.get("curvature_square_correction"):
        raise ValueError("The fixed diagnostic requires fifth-degree plus square terms")
    value, jacobian, pairs, dimension = linearization(fit, mask)
    if not np.isfinite(value).all() or np.any(value >= 0) or np.any(value**2 > 1):
        raise ValueError(
            "Pilot model clocks outside the negative probability domain: "
            f"nonnegative={int(np.sum(value >= 0))}, "
            f"squared_above_one={int(np.sum(value**2 > 1))}"
        )
    whiten = whitening(covariance)
    influence = np.linalg.pinv(whiten @ jacobian)[dimension : dimension + 10]
    clock_map = influence @ whiten
    count_map = clock_map / (2 * value)
    selector = np.eye(jacobian.shape[1])[dimension : dimension + 10]
    left_error = float(
        np.max(abs(count_map @ (2 * value[:, None] * jacobian) - selector))
    )
    if not np.isfinite(left_error) or left_error > 1e-7:
        raise ValueError(f"Local curvature left inverse failed: {left_error}")
    expected_covariance = np.asarray(fit["nominal_component_covariance"])
    covariance_error = float(
        np.linalg.norm(clock_map @ covariance @ clock_map.T - expected_covariance)
        / np.linalg.norm(expected_covariance)
    )
    if not np.isfinite(covariance_error) or covariance_error > 1e-8:
        raise ValueError(f"Saved covariance replay failed: {covariance_error}")
    points = np.asarray(fit["coordinates_inferred_from_counts"])
    signs = axis_signs([np.sign(points[1, 1]), np.sign(points[2, 2]), 1])
    transform = TIDAL_ROOT * signs
    matrix = transform @ count_map
    offset = transform @ (np.asarray(fit["coefficients"]) - count_map @ value**2)
    probabilities = counts[pairs[:, 0], pairs[:, 1]] / n
    # Undo the empirical square-root-clock delta method. The same sample
    # chose this map, so these score variances are descriptive, not calibration.
    score_clock_map = matrix * (-2 * np.sqrt(probabilities))
    score_covariance = n * score_clock_map @ covariance @ score_clock_map.T
    if not np.isfinite(score_covariance).all() or np.any(
        np.diag(score_covariance) <= 0
    ):
        raise ValueError("Nonpositive empirical score variance")
    details = {
        "local_left_inverse_max_error": left_error,
        "saved_covariance_relative_replay_error": covariance_error,
        "model_probability_minimum": float(np.min(value**2)),
        "model_probability_maximum": float(np.max(value**2)),
        "component_range_bounds": np.sum(abs(matrix), axis=1).tolist(),
        "descriptive_component_score_sd": np.sqrt(np.diag(score_covariance)).tolist(),
        "descriptive_score_covariance": score_covariance.tolist(),
        "same_sample_affine_estimate_tidal": (offset + matrix @ probabilities).tolist(),
        "measurement_n": int(n),
        "measured_intervals": len(pairs),
        "retained_covariance_rank": len(whiten),
    }
    return matrix, offset, pairs, details


def information(fit, mask, covariance):
    if fit.get("model_order", 4) != 4:
        raise ValueError(
            "The original gradient-price diagnostic requires a quartic fit"
        )
    points = np.array(fit["coordinates_inferred_from_counts"])
    _, jacobian, pairs, dimension = linearization(fit, mask)
    geometry, curvature = jacobian[:, :dimension], jacobian[:, dimension:]
    # -C(x,y,x,y)/3 becomes -(nabla_e C)(x,y,x,y)(x+y)^e/12.
    gradient = (
        np.einsum(
            "pe,pk->pek", (points[pairs[:, 0]] + points[pairs[:, 1]]) / 4, curvature
        ).reshape(len(pairs), 40)
        @ gradient_basis()
    )
    whiten = whitening(covariance)
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
