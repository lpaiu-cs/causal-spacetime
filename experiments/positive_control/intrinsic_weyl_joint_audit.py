"""Post-pilot continuum-bias audit on the SAME selected landmark populations.

This evaluator uses generator truth. It is not a second blind recovery result
or independent confirmation; the sealed finite-order pilot is kept unchanged.
"""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from intrinsic_weyl import ETA, WEYL
from intrinsic_weyl_joint import fit_joint, read_order_clocks
from intrinsic_weyl_joint_pilot import ROOT
from scipy.integrate import solve_ivp
from scipy.stats import chi2
from window_response import (
    action,
    diamond_moments,
    diamond_slice,
    order_matrix,
    sample_points,
)


def normal_coordinates(points, amplitude):
    """Evaluator-only inverse exponential map at the central tip geodesic.

    x(s)=X*s*sinhc(sqrt(A)*U*s), y(s)=Y*s*sinc(sqrt(A)*U*s),
    and v(1)=V-U*A*integral(x(s)^2-y(s)^2,0,1). Divide by tip time sqrt(2).
    """
    u = points[:, 0] - 0.5
    z = np.sqrt(amplitude) * u / np.pi
    x = points[:, 2] / np.sinc(1j * z).real
    y = points[:, 3] / np.sinc(z)
    s, w = np.polynomial.legendre.leggauss(32)
    s, w = (s + 1) / 2, w / 2
    integral = x**2 * (np.sinc(1j * z[:, None] * s).real ** 2 @ (w * s**2)) - y**2 * (
        np.sinc(z[:, None] * s) ** 2 @ (w * s**2)
    )
    v = points[:, 1] + u * amplitude * integral
    return np.column_stack(((u - v) / 2, x / np.sqrt(2), y / np.sqrt(2), (u + v) / 2))


def aligned_tensor_error(fit, normal, amplitude):
    """Truth-only comparison of all ten components; no spectra-only shortcut.

    Spatial orientation is arbitrary in the order. Align coordinate clouds by
    O(3), then transform the full covariant tensor, including reflection parity.
    """
    inferred = np.array(fit["coordinates_inferred_from_counts"])
    left, _, right = np.linalg.svd(normal[1:-1, 1:].T @ inferred[1:-1, 1:])
    rotation = left @ right
    frame = np.eye(4)
    frame[1:, 1:] = rotation
    brinkmann = np.zeros((4, 4, 4, 4))
    for axis, value in ((2, -amplitude), (3, amplitude)):
        brinkmann[0, axis, 0, axis] = brinkmann[axis, 0, axis, 0] = value
        brinkmann[axis, 0, 0, axis] = brinkmann[0, axis, axis, 0] = -value
    tetrad = np.array(
        [[1, 0, 0, 1], [-1, 0, 0, 1], [0, np.sqrt(2), 0, 0], [0, 0, np.sqrt(2), 0]]
    ) / np.sqrt(2)
    transform = tetrad @ frame
    tensor = 2 * np.einsum(
        "ia,jb,kc,ld,ijkl->abcd", transform, transform, transform, transform, brinkmann
    )
    expected = np.linalg.lstsq(WEYL.reshape(10, -1).T, tensor.ravel(), rcond=None)[0]
    assert np.max(abs(np.einsum("k,kabcd->abcd", expected, WEYL) - tensor)) < 1e-12
    difference = np.array(fit["coefficients"]) - expected
    electric = np.einsum("k,kij->ij", difference[:5], WEYL[:5, 0, 1:, 0, 1:])
    magnetic = np.einsum("k,kij->ij", difference[5:], WEYL[:5, 0, 1:, 0, 1:])
    error = float(np.sqrt(np.sum(electric**2) + np.sum(magnetic**2)))
    return {
        "electric_magnetic_frobenius_error": error,
        "relative_full_tidal_error": error / (2 * amplitude) if amplitude else None,
        "aligned_truth_coefficients": expected.tolist(),
        "alignment_parity": int(round(np.linalg.det(rotation))),
        "scope": (
            "Oracle coordinate alignment for evaluation only; not a confidence bound"
        ),
    }


def selfcheck():
    """Independent geodesic ODE and action checks for the truth-only map."""
    points = np.array(
        [
            [0, 0.5, 0, 0],
            [0.25, 0.12, 0.08, -0.04],
            [0.5, 0.03, 0.06, 0.01],
            [0.73, -0.11, -0.06, 0.03],
            [1, -0.5, 0, 0],
        ]
    )
    for amplitude in (0, 0.25, 0.75):
        normal = normal_coordinates(points, amplitude)
        for point, n in zip(points, normal, strict=True):
            u, v = n[0] + n[3], n[3] - n[0]
            x, y = np.sqrt(2) * n[1:3]

            def rhs(_, a, amplitude=amplitude):
                _, _, xx, yy, du, dv, dx, dy = a
                return [
                    du,
                    dv,
                    dx,
                    dy,
                    0,
                    -2 * amplitude * du * (xx * dx - yy * dy),
                    amplitude * du**2 * xx,
                    -amplitude * du**2 * yy,
                ]

            end = solve_ivp(
                rhs, [0, 1], [0.5, 0, 0, 0, u, v, x, y], rtol=1e-12, atol=1e-14
            ).y[:4, -1]
            assert np.max(abs(end - point)) < 1e-11
            if u:
                cost, _ = action(abs(u), 0, point[2], 0, point[3], amplitude)
                assert abs(n @ ETA @ n - (u * point[1] + abs(u) * cost)) < 1e-12
        fit = {"coefficients": np.zeros(10), "coordinates_inferred_from_counts": normal}
        comparison = aligned_tensor_error(fit, normal, amplitude)
        assert (
            abs(comparison["electric_magnetic_frobenius_error"] - 2 * amplitude) < 1e-12
        )
    return "inverse exponential map and tidal normalization independently checked"


def magnetic_spectrum(fit):
    magnetic = np.einsum("k,kij->ij", fit["coefficients"][5:], WEYL[:5, 0, 1:, 0, 1:])
    return np.linalg.eigvalsh(magnetic).tolist()


def volume_table(points, amplitude, quadrature_order=64):
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    nodes, weights = (nodes + 1) / 2, weights / 2
    volumes = np.zeros((len(points), len(points)))
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            p, q = points[i], points[j]
            if p[0] > q[0]:
                p, q = q, p
            du = q[0] - p[0]
            if du == 0:
                continue
            cost, _ = action(du, p[2], q[2], p[3], q[3], amplitude)
            h = p[1] - q[1] - cost
            if h <= 0:
                continue
            standard = weights @ diamond_slice(nodes, amplitude * du**2)
            volumes[i, j] = volumes[j, i] = du**2 * h**2 * standard
    return volumes


def nominal_eigen_intervals(fit, magnetic=False):
    c = np.array(fit["coefficients"])
    covariance = np.array(fit["nominal_component_covariance"])
    basis = WEYL[:, 0, 1:, 0, 1:]
    if magnetic:
        basis = np.concatenate((np.zeros_like(basis[:5]), basis[:5]))
    electric = np.einsum("k,kij->ij", c, basis)
    values, vectors = np.linalg.eigh(electric)
    se = []
    for vector in vectors.T:
        gradient = np.einsum("i,kij,j->k", vector, basis, vector)
        se.append(float(np.sqrt(max(0, gradient @ covariance @ gradient))))
    return {
        "eigenvalues": values.tolist(),
        "nominal_eigen_se": se,
        "nominal_marginal_95_intervals": [
            [float(v - 1.959963984540054 * s), float(v + 1.959963984540054 * s)]
            for v, s in zip(values, se, strict=True)
        ],
        "scope": "Delta-method approximation; not valid at degenerate flat eigenvalues",
        "full_tensor_flat_wald": float(c @ np.linalg.pinv(covariance) @ c),
        "nominal_flat_chi2_95_df10": float(chi2.ppf(0.95, 10)),
    }


def run(folder, target):
    if target.exists():
        raise FileExistsError("Preserve prior audits; choose a new output path")
    checked = selfcheck()
    original = json.loads((folder / "summary.json").read_text())
    rows = []
    for row in original["cases"]:
        if "fit" not in row:
            rows.append({"case": row["case"], "original_failure": row.get("failure")})
            continue
        path = folder / f"case_{row['case']}.npz"
        counts, mask, covariance, ids, _ = read_order_clocks(path)
        with np.load(path, allow_pickle=False) as saved:
            selection_n = len(saved["selection_order"]) - 2
        amplitude, seed = row["generator_amplitude"], row["seed"]
        rng = np.random.default_rng([seed, 0])
        points = np.vstack(
            (
                [0, 0.5, 0, 0],
                sample_points("D", amplitude, selection_n, rng),
                [1, -0.5, 0, 0],
            )
        )
        ordered = points[np.lexsort((-points[:, 1], points[:, 0]))]
        order, _ = order_matrix(points, amplitude)
        permutation = np.r_[
            0, rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
        ]
        ordered, order = ordered[permutation], order[permutation][:, permutation]
        with np.load(path, allow_pickle=False) as saved:
            assert np.array_equal(order, saved["selection_order"])
        points = ordered[ids]
        normal = normal_coordinates(points, amplitude)
        volume = volume_table(points, amplitude)
        coarse = volume_table(points, amplitude, 32)
        relative_quadrature_change = float(np.max(abs(coarse[mask] / volume[mask] - 1)))
        assert relative_quadrature_change < 1e-10
        assert abs(volume[0, -1] / diamond_moments(amplitude)["volume"] - 1) < 1e-10
        mean_counts = counts[0, -1] * volume / volume[0, -1]
        fit = fit_joint(mean_counts, clock_mask=mask, clock_covariance=covariance)
        expected = np.array(row["expected_electric_eigenvalues"])
        bias = np.array(fit["electric_eigenvalues"]) - expected
        result = {
            "case": row["case"],
            "amplitude": amplitude,
            "quadrature_relative_change": relative_quadrature_change,
            "continuum_eigenvalues": fit["electric_eigenvalues"],
            "continuum_eigen_bias": bias.tolist(),
            "continuum_magnetic_eigenvalues": magnetic_spectrum(fit),
            "finite_magnetic_eigenvalues": magnetic_spectrum(row["fit"]),
            "continuum_optimizer_success": fit["optimizer_success"],
            "continuum_projected_curvature_rank": fit["projected_curvature_rank"],
            "finite_order_uncertainty": nominal_eigen_intervals(row["fit"]),
            "finite_magnetic_uncertainty": nominal_eigen_intervals(
                row["fit"], magnetic=True
            ),
            "finite_full_tensor_error": aligned_tensor_error(
                row["fit"], normal, amplitude
            ),
            "continuum_full_tensor_error": aligned_tensor_error(fit, normal, amplitude),
        }
        rows.append(result)
        print(json.dumps(result), flush=True)
    sources = [
        Path(__file__),
        Path(__file__).with_name("intrinsic_weyl_joint.py"),
        Path(__file__).with_name("intrinsic_weyl_joint_pilot.py"),
        Path(__file__).with_name("intrinsic_weyl.py"),
        Path(__file__).with_name("window_response.py"),
    ]
    output = {
        "kind": "Post-pilot continuum bias diagnostic on previously selected landmarks",
        "scope": "No frozen accuracy verdict or independent coverage certification",
        "selfcheck": checked,
        "audit_code_version": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "original_pilot_code_version": original["code_version"],
        "original_summary_sha256": hashlib.sha256(
            (folder / "summary.json").read_bytes()
        ).hexdigest(),
        "cases": rows,
        "source_sha256_lf": {
            str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(
                p.read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
            for p in sources
        },
    }
    target.write_text(json.dumps(output, indent=2) + "\n")
    print(f"Written {target}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.folder, args.output)
