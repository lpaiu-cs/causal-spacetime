"""Post-pilot continuum-bias audit on the SAME selected landmark populations.

This evaluator uses generator truth. It is not a second blind recovery result
or independent confirmation; the sealed finite-order pilot is kept unchanged.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
from intrinsic_weyl import WEYL
from intrinsic_weyl_joint import fit_joint, read_order_clocks
from intrinsic_weyl_joint_pilot import ROOT, SELECTION_N
from scipy.stats import chi2
from window_response import (
    action,
    diamond_moments,
    diamond_slice,
    order_matrix,
    sample_points,
)


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


def nominal_eigen_intervals(fit):
    c = np.array(fit["coefficients"])
    covariance = np.array(fit["nominal_component_covariance"])
    basis = WEYL[:, 0, 1:, 0, 1:]
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


def run():
    folder = ROOT / "outputs/intrinsic_weyl_joint_pilot_v1"
    original = json.loads((folder / "summary.json").read_text())
    rows = []
    for row in original["cases"]:
        if "fit" not in row:
            rows.append({"case": row["case"], "original_failure": row.get("failure")})
            continue
        path = folder / f"case_{row['case']}.npz"
        counts, mask, covariance, ids, _ = read_order_clocks(path)
        amplitude, seed = row["generator_amplitude"], row["seed"]
        rng = np.random.default_rng([seed, 0])
        points = np.vstack(
            (
                [0, 0.5, 0, 0],
                sample_points("D", amplitude, SELECTION_N, rng),
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
    target = ROOT / "docs/paper/paper_a/intrinsic_weyl/joint_pilot_bias_audit.json"
    target.write_text(json.dumps(output, indent=2) + "\n")
    print(f"Written {target}")


if __name__ == "__main__":
    run()
