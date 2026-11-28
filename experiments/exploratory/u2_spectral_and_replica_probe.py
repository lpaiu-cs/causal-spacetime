"""Finite-sample exploratory U2 diagnostics; no asymptotic or detection verdict.

A: a chosen smoothed constraint Gram matrix, NOT the order likelihood Fisher
information. Gaussian smoothing lifts even conformal directions away from null.
B: blind energy and truth-aided, common-frame replica overlaps. The latter are
not order-only statistics and are not likelihood-ratio second moments (Route B1).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

# Avoid CPU oversubscription in independent replica workers.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from scipy import linalg, sparse
from subspacing_observer import (
    BackgroundModel,
    Calibration,
    bump_window,
    observe,
    tangent_project,
)
from subspacing_wave_sim import OrderOracle, WaveGeometry


def conformal_basis(points):
    """Tip-fixing infinitesimal conformal fields, in (u,v,y) coordinates.

    a=u+1/2, b=v; T=(a+b)/2, Z=(a-b)/2, Y=y/sqrt(2).
    Thus 2 da db-dy^2 = 2(dT^2-dZ^2-dY^2), with tips (T,x)=(0,0),(1,0).
    The /sqrt(2) longitudinal conversion suggested in the prompt would mix
    incompatible normalizations. In 4D: 3 rotations, D-K_time, 3 boosts-K_space.
    """
    points = np.asarray(points, dtype=float)
    n, dim = points.shape
    if dim not in (2, 4) or not np.all(np.isfinite(points)):
        raise ValueError("finite (N,2) or (N,4) points required")
    a, b = points[:, 0] + 0.5, points[:, 1]
    if dim == 2:
        fields = np.zeros((n, 2, 2))
        fields[:, 0, 0], fields[:, 1, 1] = a * (1 - a), b * (1 - b)
        return fields.reshape(2 * n, 2)
    t = (a + b) / 2
    x = np.column_stack(((a - b) / 2, points[:, 2:] / np.sqrt(2)))
    radius = np.sum(x * x, axis=1)
    fields = []

    def append(dt, dx):
        fields.append(
            np.column_stack((dt + dx[:, 0], dt - dx[:, 0], np.sqrt(2) * dx[:, 1:]))
        )

    for j, k in ((0, 1), (0, 2), (1, 2)):
        dx = np.zeros_like(x)
        dx[:, j], dx[:, k] = -x[:, k], x[:, j]
        append(np.zeros(n), dx)
    append(t * (1 - t) - radius, (1 - 2 * t[:, None]) * x)
    for k in range(3):
        dx = -2 * x[:, k, None] * x
        dx[:, k] += t * (1 - t) + radius
        append((1 - 2 * t) * x[:, k], dx)
    return np.stack(fields, axis=-1).reshape(4 * n, 7)


def diffeomorphism_basis(points, harmonics=8):
    """Endpoint-vanishing sin(k*pi*a), sin(k*pi*b) fields in 1+1D."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or harmonics < 1:
        raise ValueError("2D points and positive harmonics required")
    phases = (points + [0.5, 0])[:, :, None] * np.arange(1, harmonics + 1) * np.pi
    fields = np.zeros((len(points), 2, 2 * harmonics))
    fields[:, 0, :harmonics] = np.sin(phases[:, 0])
    fields[:, 1, harmonics:] = np.sin(phases[:, 1])
    return fields.reshape(2 * len(points), -1)


def pair_margin_gradient(left, right):
    """Rational null margin and derivative with respect to the right event."""
    delta = np.asarray(right) - np.asarray(left)
    if delta.shape[-1] != 4 or np.any(delta[..., 0] <= 0):
        raise ValueError("4D pairs require strictly positive delta u")
    du, dv, dy = delta[..., 0], delta[..., 1], delta[..., 2:]
    radius = np.sum(dy * dy, axis=-1)
    margin = dv - radius / (2 * du)
    gradient = np.concatenate(
        (
            (radius / (2 * du**2))[..., None],
            np.ones_like(du)[..., None],
            -dy / du[..., None],
        ),
        axis=-1,
    )
    return margin, gradient


def constraint_jacobian(points, eta, tips=True):
    """CSR rows sqrt(w)*dm, all unordered pairs, no Gaussian truncation.

    2D uses BOTH null faces, for either sign of the other coordinate.
    Crossing delta a=0 at negative delta b creates the REVERSE causal order;
    sorting by a must not discard that boundary of an unordered pair.
    Fixed-tip constraints remove unconstrained translations. Counts/density
    refer to pairs (union of the two faces in 2D), not the extra tip rows.
    """
    points = np.asarray(points, dtype=float)
    n, dim = points.shape
    if dim not in (2, 4) or n < 2 or not np.all(np.isfinite(points)):
        raise ValueError("finite points in dimension 2 or 4 required")
    if not np.isfinite(eta) or eta <= 0:
        raise ValueError("eta must be finite and positive")
    i, j = np.triu_indices(n, 1)
    swap = points[i, 0] > points[j, 0]
    i, j = np.where(swap, j, i), np.where(swap, i, j)
    rows, columns, values = [], [], []
    offset = 0

    def add(left, right, margin, gradient):
        nonlocal offset
        weight = np.exp(-0.25 * (margin / eta) ** 2) / np.sqrt(eta)
        gradient = gradient * weight[:, None]
        at = np.arange(offset, offset + len(margin))[:, None]
        for ids, sign in ((left, -1), (right, 1)):
            valid = ids >= 0
            rows.append(np.broadcast_to(at, gradient.shape)[valid].ravel())
            columns.append((ids[valid, None] * dim + np.arange(dim)).ravel())
            values.append((sign * gradient[valid]).ravel())
        offset += len(margin)

    def faces(left_points, right_points, left_ids, right_ids):
        delta = right_points - left_points
        if dim == 4:
            margin, gradient = pair_margin_gradient(left_points, right_points)
            add(left_ids, right_ids, margin, gradient)
            return abs(margin) < 3 * eta
        near = np.zeros(len(delta), dtype=bool)
        for axis in (0, 1):
            margin = delta[:, axis]
            gradient = np.zeros((len(delta), dim))
            gradient[:, axis] = 1
            add(left_ids, right_ids, margin, gradient)
            near |= abs(margin) < 3 * eta
        return near

    near = faces(points[i], points[j], i, j)
    if tips:
        lo, hi = np.zeros(dim), np.zeros(dim)
        lo[0], hi[0], hi[1] = -0.5, 0.5, 1
        ids, fixed = np.arange(n), np.full(n, -1)
        faces(np.broadcast_to(lo, points.shape), points, fixed, ids)
        faces(points, np.broadcast_to(hi, points.shape), ids, fixed)
    jac = sparse.coo_matrix(
        (np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
        shape=(offset, dim * n),
    ).tocsr()
    jac.eliminate_zeros()
    return jac, {
        "pair_count": len(i),
        "near_pair_count": int(near.sum()),
        "near_pair_density": float(near.mean()),
        "constraint_count": offset,
    }


def orthogonal_compression(matrix, basis, xp=np):
    """Householder similarity then discard symmetry coordinates, no fake zeros.

    Unlike eigsh(P H P), the returned matrix has no added projection nullspace.
    Only k Householder updates are needed, instead of a full dense QR basis.
    """
    matrix = xp.asarray(matrix).copy()
    basis = xp.asarray(basis).copy()
    k = basis.shape[1]
    if np.linalg.matrix_rank(np.asarray(basis.get() if xp is not np else basis)) != k:
        raise ValueError("symmetry basis is rank deficient")
    for col in range(k):
        v = basis[col:, col].copy()
        norm = xp.linalg.norm(v)
        v[0] += xp.where(v[0] >= 0, norm, -norm)
        v /= xp.linalg.norm(v)
        basis[col:, col:] -= 2 * xp.outer(v, v @ basis[col:, col:])
        h = xp.zeros(matrix.shape[0], dtype=matrix.dtype)
        h[col:] = v
        z = matrix @ h
        matrix += (
            -2 * xp.outer(h, z) - 2 * xp.outer(z, h) + 4 * (h @ z) * xp.outer(h, h)
        )
    return matrix[k:, k:]


def spectral_cell(n, dim, eta, seed, backend="gpu"):
    start = time.perf_counter()
    rng = np.random.default_rng(seed)
    # Fixed N is explicitly a Poisson sprinkling CONDITIONAL on its cardinality.
    points = (
        WaveGeometry(0).sprinkle(n, rng) if dim == 4 else rng.random((n, 2)) - [0.5, 0]
    )
    jac, counts = constraint_jacobian(points, eta)
    hessian = (jac.T @ jac).toarray()
    basis = conformal_basis(points)
    q = linalg.qr(basis, mode="economic")[0]
    symmetry_energy = np.sum((jac @ q) ** 2, axis=0)
    if backend == "gpu":
        import cupy as cp

        matrix = cp.asarray(hessian)
        eigen = cp.asnumpy(cp.linalg.eigvalsh(matrix))
        compressed = orthogonal_compression(matrix, basis, cp)
        projected = cp.asnumpy(cp.linalg.eigvalsh(compressed))
    else:
        eigen = linalg.eigvalsh(hessian)
        projected = linalg.eigvalsh(orthogonal_compression(hessian, basis))
    tolerance = float(np.max(abs(eigen)) * len(eigen) * np.finfo(float).eps)
    positive = eigen[eigen > tolerance]
    median = float(np.median(positive)) if len(positive) else None
    row = dict(
        N=n,
        dim=dim,
        eta=eta,
        seed=seed,
        symmetry_rank=basis.shape[1],
        lambda_1_to_16=eigen[:16].tolist(),
        lambda_1_to_12=eigen[:12].tolist(),
        lambda_min_perp=float(projected[0]),
        projected_first_16=projected[:16].tolist(),
        median_positive_eigenvalue=median,
        gap_over_median=float(projected[0] / median) if median else None,
        numerical_zero_tolerance=tolerance,
        symmetry_rayleigh=symmetry_energy.tolist(),
        symmetry_rayleigh_over_median=(symmetry_energy / median).tolist(),
        hessian_trace=float(np.trace(hessian)),
        **counts,
    )
    if dim == 2:
        modes = diffeomorphism_basis(points)
        rayleigh = np.sum((jac @ modes) ** 2, axis=0) / np.sum(modes**2, axis=0)
        row["diffeomorphism_rayleigh"] = rayleigh.tolist()
        row["diffeomorphism_rayleigh_over_median"] = (rayleigh / median).tolist()
    row["elapsed_seconds"] = time.perf_counter() - start
    return row


def transported_scores(points, moments, calibration, frequency, phase_u):
    """Truth-frame spatial symbol, with a separately specified carrier coordinate.

    Used ONLY after blind observation. Substituting u_hat here isolates phase
    transport; it does not turn this truth-aided score into a blind statistic.
    """
    model = BackgroundModel(calibration.q)
    predicted, jac = model.moments(points, calibration.early, calibration.late, True)
    shear = tangent_project(
        model.shear_response(points, calibration.early, calibration.late), jac
    )
    residual = tangent_project(moments - predicted, jac)
    weights = bump_window(points)
    numerator = weights * np.sum(np.conj(shear) * residual, axis=1)
    denominator = float(np.sum(weights * np.sum(abs(shear) ** 2, axis=1)))
    if denominator <= 1e-30:
        raise ValueError("no truth-frame window mass")
    carrier = np.exp(-1j * frequency * np.asarray(phase_u))
    return np.array([numerator @ carrier, numerator @ carrier.conj()]) / denominator


def replica_dataset(n, r, phi, seed, mode, max_probes, backend, population_check=False):
    """Generate truth; call a blind API; only then evaluate truth diagnostics."""
    start = time.perf_counter()
    geometry = WaveGeometry(4, mode, 1 / r, phi)
    points = geometry.sprinkle(n, np.random.default_rng(seed), backend)
    oracle = OrderOracle(geometry, points, backend)
    pilot = np.arange(min(16384, n // 4))
    references = np.arange(n // 2, n)
    probes = np.arange(len(pilot), n // 2)
    result = observe(
        oracle.query,
        pilot,
        references,
        probes,
        [r],
        q=4,
        max_probes=max_probes,
        count_query=oracle.counts if backend == "gpu" else None,
        diagnostics=True,
    )
    diagnostic = result.pop("diagnostics")
    exact = points[diagnostic["probe_ids"]]
    reconstructed = diagnostic["points"]
    calibration = Calibration(
        4, points[result["early_ids"]], points[result["late_ids"]], 0, 0, 0
    )
    exact_score = transported_scores(
        exact, diagnostic["moments"], calibration, r, exact[:, 0]
    )
    transported = transported_scores(
        exact, diagnostic["moments"], calibration, r, reconstructed[:, 0]
    )
    error = reconstructed[:, 0] - exact[:, 0]
    coherence = np.mean(np.exp(1j * r * error))
    # Population integration on a documented finite subset, not an exact global
    # continuum expectation. Double quadrature order checks integration error.
    subset = exact[:256]
    continuum = None
    if population_check:
        degree = max(96, int(np.ceil(r / 2)))
        scores = []
        for order in (degree, 2 * degree):
            population = geometry.population(
                subset, calibration.early, calibration.late, degree=order
            )
            scores.append(
                transported_scores(subset, population, calibration, r, subset[:, 0])
            )
        continuum = {
            "score": complex_list(scores[1]),
            "transported_score": complex_list(
                transported_scores(
                    subset, population, calibration, r, reconstructed[: len(subset), 0]
                )
            ),
            "probe_count": len(subset),
            "quadrature_degrees": [degree, 2 * degree],
            "quadrature_max_difference": float(np.max(abs(scores[1] - scores[0]))),
            "amplitude_over_lambda_squared": float(np.linalg.norm(scores[1]) * r * r),
            "scope": (
                "conditional population residual on selected finite probes; "
                "not global expectation"
            ),
        }
    return {
        "seed": seed,
        "phi": phi,
        "mode": mode,
        "N": n,
        "r": r,
        "blind_score": complex_list([result["score"][0], result["score_neg"][0]]),
        "blind_energy": float(
            (abs(result["score"][0]) ** 2 + abs(result["score_neg"][0]) ** 2) / 2
        ),
        "exact_score": complex_list(exact_score),
        "transported_score": complex_list(transported),
        "phase_cos_mean": float(coherence.real),
        "phase_coherence_abs": float(abs(coherence)),
        "u_error_rms": float(np.sqrt(np.mean(error**2))),
        "u_error_rms_over_lambda": float(r * np.sqrt(np.mean(error**2))),
        "u_error_abs_median": float(np.median(abs(error))),
        "accepted": result["accepted"],
        "active": result["active"],
        "n_ref": result["n_ref"],
        "n_probe": result["n_probe"],
        "calibration_rms": result["calibration_rms"],
        "continuum_population": continuum,
        "elapsed_seconds": time.perf_counter() - start,
    }


def complex_list(values):
    return [[float(complex(v).real), float(complex(v).imag)] for v in values]


def overlap(first, second):
    """Re sum over +/- r / 2, in the fixed truth transverse frame."""
    a, b = np.asarray(first) @ [1, 1j], np.asarray(second) @ [1, 1j]
    return float(np.real(a @ b.conj()) / 2)


def replica_pair(n, r, index, max_probes, backend):
    start = time.perf_counter()
    # Independent geometry samples for all five arms; only phase is shared.
    seeds = np.random.SeedSequence([20261128, n, r, index]).generate_state(5)
    phi = 0.37
    arms = (
        ("wave1", "circular_wave", phi),
        ("wave2", "circular_wave", phi),
        ("orthogonal", "circular_wave", phi + np.pi / 2),
        ("background1", "background", 0),
        ("background2", "background", 0),
    )
    records = {}
    for seed, (name, mode, phase) in zip(seeds, arms, strict=True):
        try:
            records[name] = replica_dataset(
                n,
                r,
                phase,
                int(seed),
                mode,
                max_probes,
                backend,
                population_check=(index == 0 and name == "wave1"),
            )
        except (ValueError, np.linalg.LinAlgError) as error:
            return {
                "status": "failed",
                "N": n,
                "r": r,
                "pair": index,
                "failed_arm": name,
                "error": str(error),
                "arms": records,
            }
    metrics = {
        "blind_energy_excess": (
            records["wave1"]["blind_energy"]
            + records["wave2"]["blind_energy"]
            - records["background1"]["blind_energy"]
            - records["background2"]["blind_energy"]
        )
        / 2
    }
    for field in ("exact_score", "transported_score"):
        baseline = overlap(records["background1"][field], records["background2"][field])
        for name, second in (("same", "wave2"), ("orthogonal", "orthogonal")):
            value = overlap(records["wave1"][field], records[second][field])
            metrics[f"{field}_Z_{name}"] = value
            metrics[f"{field}_Z_{name}_minus_background"] = value - baseline
        metrics[f"{field}_Z_background"] = baseline
        metrics[f"{field}_same_minus_orthogonal"] = (
            metrics[f"{field}_Z_same"] - metrics[f"{field}_Z_orthogonal"]
        )
    return {
        "status": "completed",
        "N": n,
        "r": r,
        "pair": index,
        "alpha": float(np.log(n) / np.log(r)),
        "arms": records,
        "metrics": metrics,
        "elapsed_seconds": time.perf_counter() - start,
    }


def aggregate_replicas(records):
    groups = {}
    for row in records:
        if row["status"] == "completed":
            groups.setdefault((row["N"], row["r"]), []).append(row)
    result = []
    for (n, r), rows in sorted(groups.items()):
        summary = {
            "N": n,
            "r": r,
            "alpha": float(np.log(n) / np.log(r)),
            "completed_pairs": len(rows),
            "metrics": {},
        }
        for key in rows[0]["metrics"]:
            values = np.array([row["metrics"][key] for row in rows])
            summary["metrics"][key] = {
                "mean": float(values.mean()),
                "sd_across_pairs": float(values.std(ddof=1))
                if len(values) > 1
                else None,
                "se_of_mean": float(values.std(ddof=1) / np.sqrt(len(values)))
                if len(values) > 1
                else None,
            }
        for key in ("phase_cos_mean", "phase_coherence_abs", "u_error_rms_over_lambda"):
            summary[key] = float(
                np.mean(
                    [
                        row["arms"][arm][key]
                        for row in rows
                        for arm in ("wave1", "wave2")
                    ]
                )
            )
        result.append(summary)
    return result


def spectral_slopes(rows):
    results = []
    for dim in (2, 4):
        for schedule in ("fixed", "quarter", "half"):
            group = [
                row for row in rows if row["dim"] == dim and row["schedule"] == schedule
            ]
            valid = [
                row
                for row in group
                if row["lambda_min_perp"] > row["numerical_zero_tolerance"]
            ]
            results.append(
                {
                    "dim": dim,
                    "schedule": schedule,
                    "points": len(valid),
                    "N": [row["N"] for row in valid],
                    "log_gap_vs_log_N_slope": float(
                        np.polyfit(
                            np.log([row["N"] for row in valid]),
                            np.log([row["lambda_min_perp"] for row in valid]),
                            1,
                        )[0]
                    )
                    if len(valid) > 1
                    else None,
                    "scope": (
                        "single seed per N, unnormalized constraint sum; "
                        "descriptive only"
                    ),
                }
            )
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("cpu", "gpu"), default="gpu")
    parser.add_argument("--budget-seconds", type=float, default=120)
    parser.add_argument("--pairs", type=int, default=4)
    parser.add_argument("--max-probes", type=int, default=16384)
    parser.add_argument(
        "--output", type=Path, default=Path("outputs/u2_and_replica_summary.json")
    )
    args = parser.parse_args()
    if args.budget_seconds <= 10 or args.pairs < 1 or args.max_probes < 8:
        parser.error("require budget > 10, pairs >= 1, max-probes >= 8")
    started = time.perf_counter()
    deadline = started + args.budget_seconds - 5
    summary = {
        "kind": "exploratory; no significance, accuracy gate or asymptotic verdict",
        "configuration": vars(args) | {"output": str(args.output)},
        "limitations": [
            "No inference about the (4/3,2) asymptotic bracket.",
            "H is a smoothed constraint Gram, not likelihood Fisher information.",
            "Fixed N means Poisson conditional on N; q=0 spectra fix the tips "
            "and use the coordinate Euclidean norm.",
            "Finite eta lifts symmetry modes; a positive gap is not rigidity proof.",
            "2D projects 2 global conformal modes; 4D projects 7. "
            "Sixteen harmonic fields are separate diagnostics.",
            "Overlap uses a truth frame; it is not a likelihood second moment.",
            "Orthogonal means phi+pi/2, not a uniform phase mixture; "
            "overlap need not vanish exactly.",
            "Transported scores replace only carrier u by blind u_hat; "
            "other factors remain truth-aided.",
            "Blind energy averages squared +/- scores; "
            "independent replica complex blind scores are not gauge-aligned.",
            "Continuum scores use population moments on 256 selected probes "
            "of wave1, pair 0 per cell; not an exact global expectation.",
            "Finite probes and noisy anchor calibration can dominate; "
            "four pairs do not certify detection or a threshold.",
            "Budget may censor jobs; omitted/terminated jobs are listed, not zeroed.",
        ],
        "spectral": [],
        "replicas": [],
    }
    source = Path(__file__).resolve().parent
    summary["source_sha256"] = {
        name: hashlib.sha256((source / name).read_bytes()).hexdigest()
        for name in (
            Path(__file__).name,
            "subspacing_observer.py",
            "subspacing_wave_sim.py",
        )
    }
    if args.backend == "gpu":
        import cupy as cp

        summary["device"] = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    # CPU futures overlap GPU spectra with independent finite-order observations.
    cells = [(2**20, r) for r in (64, 128, 256, 512, 1024)] + [
        (2**18, r) for r in (512, 1024)
    ]
    jobs = [(n, r, pair) for pair in range(args.pairs) for n, r in cells]
    executor = ProcessPoolExecutor(max_workers=4)
    pending = {
        executor.submit(replica_pair, n, r, pair, args.max_probes, args.backend): (
            n,
            r,
            pair,
        )
        for n, r, pair in jobs
    }
    try:
        for n in (128, 256, 512, 1024):
            for dim in (2, 4):
                for schedule, eta in (
                    ("fixed", 0.05),
                    ("quarter", 0.5 * n ** (-0.25)),
                    ("half", n ** (-0.5)),
                ):
                    if time.perf_counter() >= deadline - 10:
                        summary.setdefault("spectral_omitted", []).append(
                            {
                                "N": n,
                                "dim": dim,
                                "schedule": schedule,
                                "status": "not_run_budget",
                            }
                        )
                        continue
                    row = spectral_cell(n, dim, eta, 20261128 + n + dim, args.backend)
                    row["schedule"] = schedule
                    summary["spectral"].append(row)
                    print(
                        f"spectrum d={dim} N={n} {schedule}: "
                        f"gap={row['lambda_min_perp']:.5g}",
                        flush=True,
                    )
        while pending and time.perf_counter() < deadline:
            done, _ = wait(
                pending,
                timeout=min(1, max(0, deadline - time.perf_counter())),
                return_when=FIRST_COMPLETED,
            )
            for future in done:
                n, r, pair = pending.pop(future)
                try:
                    row = future.result()
                except Exception as error:
                    row = {
                        "status": "worker_failed",
                        "N": n,
                        "r": r,
                        "pair": pair,
                        "error": f"{type(error).__name__}: {error}",
                    }
                summary["replicas"].append(row)
                print(f"replica N={n} r={r} pair={pair}: {row['status']}", flush=True)
    finally:
        for future, (n, r, pair) in pending.items():
            summary["replicas"].append(
                {
                    "status": "budget_terminated"
                    if future.running()
                    else "not_run_budget",
                    "N": n,
                    "r": r,
                    "pair": pair,
                }
            )
        if pending:
            executor.terminate_workers()
        else:
            executor.shutdown()
    summary["spectral_2048"] = (
        "not_run: optional sparse <15s path not benchmarked within shared GPU budget"
    )
    summary["spectral_slopes"] = spectral_slopes(summary["spectral"])
    summary["replica_aggregates"] = aggregate_replicas(summary["replicas"])
    summary["replicas"].sort(key=lambda row: (row["N"], row["r"], row["pair"]))
    summary["elapsed_seconds"] = time.perf_counter() - started
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(f"Saved {args.output} in {summary['elapsed_seconds']:.2f}s", flush=True)


if __name__ == "__main__":
    main()
