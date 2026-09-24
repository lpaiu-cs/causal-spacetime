"""Generator-coordinate diagnostic of joint cone constraints, not an observer."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import time
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy
from order_wave_null_coordinate import ROOT, WaveFixture, series
from scipy.linalg import lstsq
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsmr

DESIGN = ROOT / "docs/paper/paper_a/order_wave_energy/cone_constraint_challenge.md"


def background_points(n, seed):
    """Rejection sample the q=1 volume law; columns are (u,m,v,y1,y2)."""
    rng = np.random.default_rng(seed)
    total = 2 * math.tan(0.5)
    kept = []
    size = 0
    while size < n:
        values = rng.random((4 * n, 5))
        m, v = values[:, 0] * total, values[:, 1]
        y = (2 * values[:, 2:4] - 1) * math.sqrt(total / 2)
        u = np.arctan(m - math.tan(0.5))
        inside = np.sum(y**2, axis=1) < 2 * np.minimum(m * v, (total - m) * (1 - v))
        accept = inside & (values[:, 4] < np.cos(u) ** 4)
        batch = np.column_stack((u[accept], m[accept], v[accept], y[accept]))
        kept.append(batch)
        size += len(batch)
    return np.concatenate(kept)[:n]


def cone_system(points):
    n = len(points)
    total = 2 * math.tan(0.5)
    x = np.column_stack(
        (
            (points[:, 1] + points[:, 2]) / math.sqrt(2),
            (points[:, 1] - points[:, 2]) / math.sqrt(2),
            points[:, 3:],
        )
    )
    i, j = np.triu_indices(n, 1)
    swap = x[i, 0] > x[j, 0]
    i, j = np.where(swap, j, i), np.where(swap, i, j)
    delta = x[j] - x[i]
    radius = np.linalg.norm(delta[:, 1:], axis=1)
    width = math.sqrt(2 * total) * 4 * math.log(n) / n
    select = (
        (abs(delta[:, 0] - radius) <= width)
        & (points[j, 1] - points[i, 1] >= total / 10)
        & (points[j, 2] - points[i, 2] >= 0.1)
    )
    i, j, delta, radius = i[select], j[select], delta[select], radius[select]
    first = np.column_stack((-np.ones(len(i)), delta[:, 1:] / radius[:, None]))
    columns = np.column_stack(
        (4 * i[:, None] + np.arange(4), 4 * j[:, None] + np.arange(4))
    )
    matrix = csr_matrix(
        (
            np.column_stack((first, -first)).ravel(),
            (np.repeat(np.arange(len(i)), 8), columns.ravel()),
        ),
        shape=(len(i), 4 * n),
    )
    return x, i, j, delta[:, 0] + radius, matrix, width


def project(matrix, forcing, *, tighter=False):
    tolerance = 1e-13 if tighter else 1e-11
    fit = lsmr(
        matrix,
        forcing,
        atol=tolerance,
        btol=tolerance,
        conlim=1e12,
        maxiter=40000 if tighter else 20000,
    )
    residual = forcing - matrix @ fit[0]
    norm = np.linalg.norm(forcing)
    return residual, {
        "stop": fit[1],
        "iterations": fit[2],
        "relative_residual": float(np.linalg.norm(residual) / norm) if norm else 0.0,
        "relative_stationarity": float(np.linalg.norm(matrix.T @ residual) / norm)
        if norm
        else 0.0,
        "adjustment_norm": float(np.linalg.norm(fit[0])),
        "condition_estimate": float(fit[6]),
    }


def forcing(points, i, j, denominator, wave):
    _, dm, integral = wave.parts(points[i, 0], points[j, 0])
    z = wave.factor * integral
    e = np.stack((z.real, z.imag, z.imag, -z.real), axis=-1).reshape(-1, 2, 2)
    dy = points[j, 3:] - points[i, 3:]
    return action_changes(dm, dy, e, denominator)


def action_changes(dm, dy, e, denominator):
    ey = np.einsum("nij,nj->ni", e, dy)
    linear = -np.sum(dy * ey, axis=1) / (dm * denominator)
    exact = (
        -np.sum(
            dy
            * np.linalg.solve(dm[:, None, None] * np.eye(2) - e, ey[..., None])[..., 0],
            axis=1,
        )
        / denominator
    )
    return linear, exact


def self_check():
    points = background_points(64, 71)
    x, i, j, _, matrix, _ = cone_system(points)
    shift = np.random.default_rng(82).normal(size=x.shape)

    def margins(values):
        change = values[j] - values[i]
        return change[:, 0] - np.linalg.norm(change[:, 1:], axis=1)

    step = 1e-6
    finite = (margins(x + step * shift) - margins(x - step * shift)) / (2 * step)
    assert np.max(abs(finite - matrix @ shift.ravel())) < 2e-8
    for axis in range(4):
        translation = np.zeros_like(x)
        translation[:, axis] = 1
        assert np.max(abs(matrix @ translation.ravel())) == 0
    residual, control = project(matrix, matrix @ shift.ravel())
    assert control["stop"] in (1, 2, 4, 5) and control["relative_residual"] < 1e-7
    zero, _ = project(matrix, np.zeros(len(i)))
    assert not np.any(zero)
    dm = points[j, 1] - points[i, 1]
    dy = points[j, 3:] - points[i, 3:]
    delta = x[j] - x[i]
    denominator = delta[:, 0] + np.linalg.norm(delta[:, 1:], axis=1)
    e = dm[:, None, None] * np.array([[0.01, 0.02], [0.02, -0.01]])
    linear, _ = action_changes(dm, dy, e, denominator)
    identity = dm[:, None, None] * np.eye(2)
    inverse_change = (
        np.linalg.solve(identity - step * e, dy[..., None])[..., 0]
        - np.linalg.solve(identity + step * e, dy[..., None])[..., 0]
    ) / (2 * step)
    independent = -dm * np.sum(dy * inverse_change, axis=1) / denominator
    assert np.max(abs(independent - linear)) < 2e-8
    assert not any(np.any(item) for item in action_changes(dm, dy, 0 * e, denominator))
    rotated = points.copy()
    rotated[:, 3], rotated[:, 4] = -points[:, 4], points[:, 3]
    _, ri, rj, rden, rmatrix, _ = cone_system(rotated)
    assert np.array_equal(i, ri) and np.array_equal(j, rj)
    transformed = matrix.toarray().reshape(len(i), len(points), 4).copy()
    original = transformed.copy()
    transformed[:, :, 2], transformed[:, :, 3] = -original[:, :, 3], original[:, :, 2]
    assert (
        np.max(abs(rmatrix.toarray().reshape(transformed.shape) - transformed)) < 1e-14
    )
    opposite, _ = action_changes(dm, rotated[j, 3:] - rotated[i, 3:], e, rden)
    assert np.max(abs(opposite + linear)) < 1e-14
    print(
        json.dumps(
            {
                "jacobian_max_error": float(
                    np.max(abs(finite - matrix @ shift.ravel()))
                ),
                "coordinate_control": control,
                "zero_control": True,
                "quarter_turn_odd_forcing_and_equal_coordinate_range": True,
            },
            indent=2,
        )
    )


def run(output):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit design and source and run from a clean tree")
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    stamp = (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    )
    files = [Path(__file__), DESIGN] + [
        Path(__file__).with_name(name)
        for name in (
            "order_wave_null_coordinate.py",
            "order_wave_resolution.py",
            "order_wave_energy.py",
        )
    ]
    manifest = {
        "status": "in_progress",
        "code_version": stamp,
        "utc_started": time.strftime("%FT%TZ", time.gmtime()),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "files": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files
        },
        "scope": (
            "known-coordinate linear obstruction; not an order-only test "
            "or a new resolution bound"
        ),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    started = time.time()
    waves = {r: WaveFixture(1, r) for r in (32, 128, 512)}
    with (output / "rows.jsonl").open("w", encoding="utf-8") as stream:
        for n in (128, 256, 512):
            for seed in (118, 819):
                points = background_points(n, seed)
                _, i, j, denominator, matrix, width = cone_system(points)
                control_shift = np.random.default_rng(seed + 1000).normal(size=4 * n)
                _, control = project(matrix, matrix @ control_shift)
                valid = (
                    control["stop"] in (1, 2, 4, 5)
                    and control["relative_residual"] < 1e-7
                )
                degree = np.bincount(np.r_[i, j], minlength=n)
                for r, wave in waves.items():
                    linear, exact = forcing(points, i, j, denominator, wave)
                    residual, fit = project(matrix, linear)
                    row = {
                        "n": n,
                        "seed": seed,
                        "r": r,
                        "edges": len(i),
                        "width": width,
                        "minimum_degree": int(degree.min()),
                        "mean_degree": float(degree.mean()),
                        "coordinate_control": control,
                        "control_valid": bool(valid),
                        "fit": fit,
                        "scaled_forcing_rms": float(np.sqrt(np.mean(linear**2)) * r**2),
                        "scaled_residual_rms": float(
                            np.sqrt(np.mean(residual**2)) * r**2
                        ),
                        "exact_linear_relative_change": float(
                            np.linalg.norm(exact - linear) / np.linalg.norm(linear)
                        ),
                    }
                    if (n, seed, r) == (128, 118, 512):
                        dense = matrix.toarray()
                        fitted, _, rank, _ = lstsq(
                            dense, linear, cond=1e-12, lapack_driver="gelsd"
                        )
                        dense_residual = linear - dense @ fitted
                        row["dense_svd"] = {
                            "rank": int(rank),
                            "relative_residual": float(
                                np.linalg.norm(dense_residual) / np.linalg.norm(linear)
                            ),
                            "residual_vector_relative_change": float(
                                np.linalg.norm(dense_residual - residual)
                                / np.linalg.norm(linear)
                            ),
                        }
                    if (n, seed, r) == (512, 819, 512):
                        other, refined = project(matrix, linear, tighter=True)
                        row["refinement"] = {
                            "fit": refined,
                            "residual_vector_relative_change": float(
                                np.linalg.norm(other - residual)
                                / np.linalg.norm(linear)
                            ),
                        }
                    stream.write(json.dumps(row) + "\n")
                    stream.flush()
                    print(json.dumps(row), flush=True)
    manifest.update(
        status="completed",
        elapsed_seconds=time.time() - started,
        rows_sha256=hashlib.sha256((output / "rows.jsonl").read_bytes()).hexdigest(),
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))


def signed_audit(folder):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit signed-audit plan before running")
    target = folder / "signed_information_audit.json"
    if target.exists():
        raise FileExistsError(target)
    raw = (folder / "rows.jsonl").read_bytes()
    manifest = json.loads((folder / "manifest.json").read_text())
    assert hashlib.sha256(raw).hexdigest() == manifest["rows_sha256"]
    waves = {r: WaveFixture(1, r) for r in (32, 128, 512)}
    result = {
        "scope": "post-run binary-order audit with rational all-phase bounds",
        "code_version": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT)
        .decode()
        .strip(),
        "files": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), DESIGN)
        },
        "primary_raw_sha256": manifest["rows_sha256"],
        "rows": [],
    }
    for n in (128, 256, 512):
        for seed in (118, 819):
            points = background_points(n, seed)
            points = points[np.argsort(points[:, 1])]
            i, j = np.triu_indices(n, 1)
            dm, dv = points[j, 1] - points[i, 1], points[j, 2] - points[i, 2]
            dy = points[j, 3:] - points[i, 3:]
            baseline = dv - np.sum(dy * dy, axis=1) / (2 * dm)
            chart_delta_t = (dm + dv) / math.sqrt(2)
            denominator = chart_delta_t + np.sqrt(
                (dm - dv) ** 2 / 2 + np.sum(dy * dy, axis=1)
            )
            for r, wave in waves.items():
                _, exact = forcing(points, i, j, denominator, wave)
                shifted = baseline + exact * denominator / (2 * dm)
                row = {
                    "n": n,
                    "seed": seed,
                    "r": r,
                    "numerical_pair_sign_changes": int(
                        np.count_nonzero((baseline > 0) != (shifted > 0))
                    ),
                }
                if r == 512:
                    values = [tuple(F(float(v)) for v in point[1:]) for point in points]
                    assert all(0 < item[0] < F(131, 120) for item in values)
                    unresolved, minimum = 0, None
                    lam = F(1, r)
                    for left in range(n):
                        for right in range(left + 1, n):
                            delta = [
                                values[right][k] - values[left][k] for k in range(4)
                            ]
                            action = (delta[2] ** 2 + delta[3] ** 2) / (2 * delta[0])
                            bound = min(2 * lam * delta[0], F(316, 49) * lam**2)
                            clearance = abs(delta[1] - action) - action * bound / (
                                delta[0] - bound
                            )
                            unresolved += clearance <= 0
                            minimum = (
                                clearance
                                if minimum is None
                                else min(minimum, clearance)
                            )
                    row["rational_bound_unresolved_pairs"] = unresolved
                    row["certified_identical_order_all_phases"] = unresolved == 0
                    row["minimum_certificate_clearance"] = float(minimum)
                result["rows"].append(row)
                print(json.dumps(row), flush=True)
    target.write_text(json.dumps(result, indent=2))


def domain_audit(folder, *, tight=False):
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit domain check before running")
    target = folder / ("tight_domain_audit.json" if tight else "domain_audit.json")
    if target.exists():
        raise FileExistsError(target)
    previous = (folder / "signed_information_audit.json").read_bytes()
    result = {
        "code_version": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT)
        .decode()
        .strip(),
        "scope": "exact-rational all-phase domain certificate at r=512",
        "signed_audit_sha256": hashlib.sha256(previous).hexdigest(),
        "files": {
            p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), DESIGN)
        },
        "rows": [],
    }
    lam = F(1, 512)
    endpoint = (
        (series(F(1, 4), 1, 128) / series(F(1, 4), 0, 128)).lo if tight else F(131, 120)
    )
    result["endpoint_lower_bound"] = str(endpoint)
    result["tight_endpoint_enclosure"] = tight
    for n in (128, 256, 512):
        for seed in (118, 819):
            slacks = []
            for point in background_points(n, seed):
                m, v, y1, y2 = map(lambda x: F(float(x)), point[1:])
                for length, available in ((m, v), (endpoint - m, 1 - v)):
                    assert length > 0
                    error = min(2 * lam * length, F(316, 49) * lam**2)
                    slacks.append(available - (y1**2 + y2**2) / (2 * (length - error)))
            row = {
                "n": n,
                "seed": seed,
                "r": 512,
                "all_points_inside_wave_diamond_all_phases": min(slacks) > 0,
                "minimum_certified_slack": float(min(slacks)),
            }
            result["rows"].append(row)
            print(json.dumps(row))
    target.write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--signed-audit", type=Path)
    parser.add_argument("--domain-audit", type=Path)
    parser.add_argument("--tight-domain-audit", type=Path)
    args = parser.parse_args()
    if args.self_check:
        self_check()
    if args.output:
        run(args.output)
    if args.signed_audit:
        signed_audit(args.signed_audit)
    if args.domain_audit:
        domain_audit(args.domain_audit)
    if args.tight_domain_audit:
        domain_audit(args.tight_domain_audit, tight=True)
