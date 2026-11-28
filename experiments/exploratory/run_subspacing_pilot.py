"""Budgeted exploratory campaign. No frozen gate or significance claim.

Default budget is 240 seconds. Every requested configuration is recorded as
completed, failed, or not run for budget, never silently represented by a smaller
sample. N counts all independent pilot/reference/probe-pool events; reference
events are all used, while expensive probe queries are explicitly capped.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

# Tiny calibration SVDs are slower when BLAS launches a full CPU thread team.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402
from subspacing_observer import (  # noqa: E402
    FUTURE,
    PAST,
    BackgroundModel,
    calibrate,
    neyman_orthogonalize,
    observe,
    order_counts,
    spectrum,
    triangulate,
)
from subspacing_wave_sim import OrderOracle, WaveGeometry  # noqa: E402


def continuum_anchors():
    """Truth-side deterministic fixture, converted to probabilities before fitting."""
    return (
        np.array(
            [
                [-0.25, 0.25, 0.12, 0],
                [-0.18, 0.22, -0.06, 0.12],
                [-0.29, 0.34, -0.10, -0.07],
                [-0.12, 0.31, 0.08, -0.09],
                [-0.32, 0.20, -0.02, 0.06],
            ]
        ),
        np.array(
            [
                [0.24, 0.76, 0.10, 0.03],
                [0.17, 0.79, -0.12, 0.06],
                [0.31, 0.68, 0.04, -0.13],
                [0.12, 0.72, -0.07, -0.05],
                [0.28, 0.84, 0.03, 0.11],
            ]
        ),
    )


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, complex):
        return [value.real, value.imag]
    if isinstance(value, (np.integer, np.floating)):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def target_metrics(result, frequencies, target):
    """Truth-side evaluation only: the blind spectrum never receives the target."""
    frequencies = np.asarray(frequencies)
    indices = np.flatnonzero(np.isclose(frequencies, target, atol=1e-12, rtol=0))
    if len(indices) != 1:
        raise ValueError("target must occur exactly once in the frequency grid")
    amplitude = result["whitened_amplitude"]
    off = abs(frequencies - target) > 20
    floor = float(np.sqrt(np.mean(amplitude[off] ** 2))) if off.any() else float("nan")
    return dict(
        whitened_target_z=float(amplitude[indices[0]]),
        whitened_target_snr=float(amplitude[indices[0]] / max(floor, 1e-30)),
    )


def stage_one():
    early, late = continuum_anchors()
    anchors = np.vstack((early, late))
    records = []
    for q in (0.05, 0.25, 1.0, 4.0):
        started = time.perf_counter()
        model = BackgroundModel(q)
        pf = np.column_stack(
            (model.interval(PAST, anchors), model.interval(anchors, FUTURE))
        )
        fit = calibrate(
            pf,
            model.interval(early[:, None], late[None]),
            q=q,
            starts=64,
            max_nfev=2000,
        )
        background = WaveGeometry(q)
        points = background.sprinkle(16000, np.random.default_rng(2718))
        active = np.all(model.moments(points, early, late) > 1e-5, axis=1)
        points = points[active]
        exact = background.population(points, early, late)
        recovered, residual, _ = triangulate(exact, fit, initial=None)
        left, _, right = np.linalg.svd(
            anchors[:, 2:].T @ np.vstack((fit.early, fit.late))[:, 2:]
        )
        aligned = points.copy()
        aligned[:, 2:] = points[:, 2:] @ (left @ right)
        coord_error = float(np.max(np.abs(recovered - aligned)))
        if coord_error >= 1e-11:
            raise AssertionError(
                f"initial=None recovery failed at q={q}: {coord_error}"
            )
        # Sampling error diagnostic with fixed population-calibrated anchors;
        # excludes anchor-estimation uncertainty and is labelled accordingly.
        coordinate_errors = []
        n_ref = 4096
        for repetition in range(4):
            refs = background.sprinkle(n_ref, np.random.default_rng(100 + repetition))
            all_points = np.vstack((anchors, points[:32], refs))
            oracle = OrderOracle(background, all_points)
            mu = order_counts(
                oracle.query,
                np.arange(42, len(all_points)),
                np.arange(10, 42),
                np.arange(5),
                np.arange(5, 10),
            )
            noisy, _, _ = triangulate(mu, fit)
            coordinate_errors.append(noisy[:, 0] - recovered[:32, 0])
        sd_constant = float(
            np.sqrt(np.mean(np.var(coordinate_errors, axis=0, ddof=1))) * np.sqrt(n_ref)
        )
        for r in (64, 128, 256, 512):
            geometry = WaveGeometry(q, "circular_wave", 1 / r)
            mu = geometry.population(points, early, late, degree=max(96, r // 2 + 64))
            estimate, error, jac = triangulate(mu, fit)
            error, dtheta_norm = neyman_orthogonalize(estimate, error, jac, fit)
            frequencies = np.array([r - 32, r - 16, r, r + 16, r + 32], dtype=float)
            score = spectrum(estimate, error, jac, fit, frequencies)
            records.append(
                {
                    **{
                        key: score[key]
                        for key in (
                            "whitened_peak_signed_frequency",
                            "whitened_peak_z",
                            "whitened_snr",
                            "shot_sd",
                            "c_ref_estimate",
                        )
                    },
                    **target_metrics(score, frequencies, r),
                    "dtheta_norm": dtheta_norm,
                    "orthogonalize": True,
                    "q": q,
                    "r": r,
                    "init_none_coord_error_max": coord_error,
                    "init_none_points": len(points),
                    "background_residual_max": float(np.max(abs(residual))),
                    "forward_map_error_max": float(
                        np.max(abs(exact - model.moments(points, early, late)))
                    ),
                    "population_calibration_rms": fit.rms,
                    "wave_score_over_lambda_squared": float(
                        score["amplitude"][2] * r * r
                    ),
                    "sd_u_times_sqrt_n_ref": sd_constant,
                    "sd_repetitions": 4,
                    "sd_scope": "reference noise only; population-calibrated anchors",
                    "n_ref": n_ref,
                    "elapsed_q_seconds": time.perf_counter() - started,
                }
            )
        print(
            f"stage1 q={q:g}: residual={np.max(abs(residual)):.3g}, "
            f"initial=None error={coord_error:.3g} ({len(points)} points), "
            f"SD*sqrt(nR)={sd_constant:.3g}",
            flush=True,
        )
    return records


def primitive_benchmark(backend):
    geometry = WaveGeometry(4, "circular_wave", 1 / 128)
    u = np.linspace(-0.5, 0.5, 10_000_000)
    start = time.perf_counter()
    expected = geometry.primitive(u)
    result = {"points": len(u), "cpu_seconds": time.perf_counter() - start}
    if backend == "gpu":
        import cupy as cp

        device_u = cp.asarray(u)
        geometry.primitive.gpu(device_u[:16])
        cp.cuda.Stream.null.synchronize()
        start = time.perf_counter()
        actual = geometry.primitive.gpu(device_u)
        cp.cuda.Stream.null.synchronize()
        result["gpu_resident_seconds"] = time.perf_counter() - start
        start = time.perf_counter()
        host = cp.asnumpy(actual)
        result["gpu_download_seconds"] = time.perf_counter() - start
        result["gpu_cpu_max_error"] = float(np.max(abs(host - expected)))
        result["gpu_under_half_second"] = result["gpu_resident_seconds"] < 0.5
        result["device"] = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
        result["gpu_memory_pool_bytes"] = int(
            cp.get_default_memory_pool().total_bytes()
        )
        del actual, device_u
        cp.get_default_memory_pool().free_all_blocks()
    return result


def finite_dataset(
    n,
    q,
    mode,
    r,
    seed,
    backend,
    max_probes,
    known_modes=(True, False),
    frequencies=None,
    diagnostics=False,
):
    start = time.perf_counter()
    geometry = WaveGeometry(q, mode, 1 / r)
    points = geometry.sprinkle(n, np.random.default_rng(seed), backend)
    oracle = OrderOracle(geometry, points, backend)
    generation = time.perf_counter() - start
    pilot_n = min(16384, n // 4)
    probe_n = n // 2 - pilot_n
    pilot = np.arange(pilot_n)
    probes = np.arange(pilot_n, pilot_n + probe_n)
    references = np.arange(n // 2, n)
    cache = {}

    def counts(refs, at, early=(), late=()):
        key = tuple(
            hashlib.blake2b(
                np.asarray(g, dtype=np.int64).tobytes(), digest_size=16
            ).digest()
            for g in (refs, at, early, late)
        )
        if key not in cache:
            cache[key] = (
                oracle.counts(refs, at, early, late)
                if backend == "gpu"
                else order_counts(oracle.query, refs, at, early, late)
            )
        return cache[key]

    frequencies = (
        np.arange(24, 225, 4, dtype=float)
        if frequencies is None
        else np.asarray(frequencies, float)
    )
    records = []
    for known in known_modes:
        fit_started = time.perf_counter()
        row = {
            "N": n,
            "q_true": q,
            "mode": mode,
            "r": r,
            "lambda": 1 / r,
            "alpha_log_N_over_log_r": float(np.log(n) / np.log(r)),
            "sub_discreteness": bool(1 / r < n ** (-0.25)),
            "known_q": known,
            "seed": seed,
            "backend": backend,
            "generation_seconds": generation,
            "pilot_size": len(pilot),
            "reference_size": len(references),
            "probe_pool_size": len(probes),
            "max_counted_probes": max_probes,
        }
        try:
            result = observe(
                oracle.query,
                pilot,
                references,
                probes,
                frequencies,
                q=q if known else None,
                count_query=counts,
                max_probes=max_probes,
                diagnostics=diagnostics,
            )
            target = int(np.argmin(abs(frequencies - r)))
            off_target = abs(frequencies - r) > 20
            high = frequencies >= 48
            off_high = off_target & (frequencies >= 44)
            floor = float(np.sqrt(np.mean(result["amplitude"][off_target] ** 2)))
            floor_high = float(np.sqrt(np.mean(result["amplitude"][off_high] ** 2)))
            result.update(
                **target_metrics(result, frequencies, r),
                whitened_spectrum_l2=float(
                    np.sqrt(np.mean(result["whitened_amplitude"] ** 2))
                ),
                log_whitened_amplitude_vs_log_omega_slope=float(
                    np.polyfit(
                        np.log(frequencies),
                        np.log(np.maximum(result["whitened_amplitude"], 1e-30)),
                        1,
                    )[0]
                ),
                target_off_peak_rms_high=floor_high,
                target_snr_high=float(
                    result["amplitude"][target] / max(floor_high, 1e-30)
                ),
                spectrum_l2_high=float(
                    np.sqrt(np.mean(result["amplitude"][high] ** 2))
                ),
                log_amplitude_vs_log_omega_slope_high=float(
                    np.polyfit(
                        np.log(frequencies[high]),
                        np.log(np.maximum(result["amplitude"][high], 1e-30)),
                        1,
                    )[0]
                ),
                target_amplitude=float(result["amplitude"][target]),
                target_off_peak_rms=floor,
                target_snr=float(result["amplitude"][target] / max(floor, 1e-30)),
                q_absolute_error=abs(result["q_hat"] - q),
                frequencies=frequencies,
                noise_scaled_by_omega_sqrt_n_ref=(
                    result["amplitude"] * frequencies * np.sqrt(len(references))
                ),
                spectrum_l2=float(np.sqrt(np.mean(result["amplitude"] ** 2))),
                log_amplitude_vs_log_omega_slope=float(
                    np.polyfit(
                        np.log(frequencies),
                        np.log(np.maximum(result["amplitude"], 1e-30)),
                        1,
                    )[0]
                ),
            )
            row.update(status="completed", **result)
        except (ValueError, np.linalg.LinAlgError) as error:
            row.update(status="failed", error=str(error))
        row["observer_seconds"] = time.perf_counter() - fit_started
        records.append(row)
        print(
            f"N={n} q={q:g} {mode} r={r} {'known' if known else 'unknown'}: "
            f"{row['status']}, q_hat={row.get('q_hat', float('nan')):.4g}, "
            f"target SNR={row.get('target_snr', float('nan')):.3g}, "
            f"high={row.get('target_snr_high', float('nan')):.3g}, "
            f"{row['observer_seconds']:.2f}s",
            flush=True,
        )
    if backend == "gpu":
        import cupy as cp

        memory = int(cp.get_default_memory_pool().total_bytes())
        for row in records:
            row["gpu_memory_pool_bytes"] = memory
            row["gpu_memory_pool_scope"] = (
                "worker-local pool; includes reusable allocations"
            )
        del oracle
        cp.get_default_memory_pool().free_all_blocks()
    return records


def campaign_analysis(summary):
    controls = summary["stage2"]
    slopes = []
    for q in (0, 1):
        for known in (True, False):
            selected = [
                r
                for r in controls
                if r.get("status") == "completed"
                and r["q_true"] == q
                and r["known_q"] == known
            ]
            if len(selected) >= 2:
                slope = np.polyfit(
                    np.log([r["n_ref"] for r in selected]),
                    np.log([r["spectrum_l2"] for r in selected]),
                    1,
                )[0]
                slopes.append(
                    {
                        "q": q,
                        "known_q": known,
                        "log_whitened_noise_vs_log_n_ref_slope": float(
                            np.polyfit(
                                np.log([r["n_ref"] for r in selected]),
                                np.log([r["whitened_spectrum_l2"] for r in selected]),
                                1,
                            )[0]
                        ),
                        "mean_log_amplitude_vs_log_omega_slope": float(
                            np.mean(
                                [
                                    r["log_amplitude_vs_log_omega_slope"]
                                    for r in selected
                                ]
                            )
                        ),
                        "mean_log_whitened_amplitude_vs_log_omega_slope": float(
                            np.mean(
                                [
                                    r["log_whitened_amplitude_vs_log_omega_slope"]
                                    for r in selected
                                ]
                            )
                        ),
                        "log_noise_vs_log_n_ref_slope": float(slope),
                        "log_noise_high_vs_log_n_ref_slope": float(
                            np.polyfit(
                                np.log([r["n_ref"] for r in selected]),
                                np.log([r["spectrum_l2_high"] for r in selected]),
                                1,
                            )[0]
                        ),
                        "mean_log_amplitude_vs_log_omega_slope_high": float(
                            np.mean(
                                [
                                    r["log_amplitude_vs_log_omega_slope_high"]
                                    for r in selected
                                ]
                            )
                        ),
                        "omega_slope_references": {
                            "inverse_omega": -1.0,
                            "flat_shot_noise": 0.0,
                        },
                        "points": len(selected),
                        "nominal_reference": -0.5,
                    }
                )
    rows = summary["stage3"]
    for row in rows:
        if row.get("status") != "completed" or row["mode"] != "circular_wave":
            continue
        baseline = next(
            (
                b
                for b in rows
                if b.get("status") == "completed"
                and b["mode"] == "background"
                and b["N"] == row["N"]
                and b["q_true"] == row["q_true"]
                and b["known_q"] == row["known_q"]
            ),
            None,
        )
        if baseline:
            idx = int(np.argmin(abs(np.asarray(baseline["frequencies"]) - row["r"])))
            denominator = baseline["amplitude"][idx]
            row["target_amplitude_over_background"] = float(
                row["target_amplitude"] / max(denominator, 1e-30)
            )
    return {
        "noise_slopes": slopes,
        "stage1_max_init_none_coord_error": max(
            r["init_none_coord_error_max"] for r in summary["stage1"]
        ),
        "stage1_max_background_residual": max(
            r["background_residual_max"] for r in summary["stage1"]
        ),
        "stage1_max_forward_map_error": max(
            r["forward_map_error_max"] for r in summary["stage1"]
        ),
        "stage1_min_wave_response_over_lambda_squared": min(
            r["wave_score_over_lambda_squared"] for r in summary["stage1"]
        ),
        "verdict": (
            "exploratory only: three Stage 2 repetitions; "
            "single Stage 3 dataset per cell, "
            "no detection or null-error certification"
        ),
        "noise_law_limit": (
            "finite-probe shot noise need not decay as 1/omega; "
            "record scaling instead of asserting it"
        ),
        "alpha_limit": (
            "finite points with alpha>2 do not establish asymptotic detectability"
        ),
        "q_zero_limit": "flat q=0 has no identifiable preferred null axis",
        "matcher_limit": (
            "exact shear differential, high-frequency endpoint carrier; "
            "not a finite-lambda optimal likelihood"
        ),
        "initialization_limit": (
            "64 orientations and multistart polishing do not guarantee a global optimum"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget-seconds", type=float, default=240)
    parser.add_argument("--backend", choices=("cpu", "gpu"), default="gpu")
    parser.add_argument("--max-probes", type=int, default=16384)
    parser.add_argument(
        "--output", type=Path, default=Path("outputs/subspacing_pilot_summary.json")
    )
    args = parser.parse_args()
    start = time.perf_counter()
    summary = {
        "kind": "exploratory, not preregistered or confirmatory",
        "budget_seconds": args.budget_seconds,
        "concurrent_cells": 4,
        "parallel_backend": "process",
        "configuration": vars(args) | {"output": str(args.output)},
        "stage1": [],
        "stage2": [],
        "stage3": [],
    }
    if args.output.exists():
        previous = json.loads(args.output.read_text())
        history = previous.pop("previous_campaigns", [])
        summary["previous_campaigns"] = history + [previous]
    source = Path(__file__).resolve().parent
    summary["source_sha256"] = {
        name: hashlib.sha256((source / name).read_bytes()).hexdigest()
        for name in (
            "subspacing_wave_sim.py",
            "subspacing_observer.py",
            "run_subspacing_pilot.py",
        )
    }
    summary["python_version"] = platform.python_version()
    summary["numpy_version"] = np.__version__
    summary["anchor_policy"] = (
        "Four central anchors per side, p in [.0005,.04], "
        "f in [.25,.80] (dual for late); "
        "support-cubed core selection; up to ten ordered calibration anchors per side."
    )

    def save():
        summary["elapsed_seconds"] = time.perf_counter() - start
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(json_safe(summary), indent=2, allow_nan=False) + "\n"
        )

    summary["primitive_benchmark"] = primitive_benchmark(args.backend)
    summary["stage1"] = stage_one()
    save()
    jobs = [
        ("stage2", n, q, "background", 128)
        for _ in range(3)
        for n in (2**17, 2**19, 2**21)
        for q in (0.0, 1.0)
    ]
    jobs += [
        ("stage3", n, q, mode, r)
        for n in (2**20, 2**22)
        for q in (1.0, 4.0, 0.25)
        for mode, r in [("background", 128)]
        + [("circular_wave", r) for r in (64, 96, 128, 192)]
    ]
    jobs += [
        ("stage3", 2**23, q, mode, r)
        for q, mode, r in [
            (1.0, "background", 128),
            (1.0, "circular_wave", 64),
            (4.0, "background", 128),
            (4.0, "circular_wave", 64),
            (4.0, "circular_wave", 96),
        ]
    ]
    previous = 8.0
    durations = {}
    index = 0
    pending = {}
    # Independent cells overlap CPU fitting with GPU counts; seeds stay fixed.
    with ProcessPoolExecutor(max_workers=4) as executor:
        while index < len(jobs) or pending:
            while index < len(jobs) and len(pending) < 4:
                stage, n, q, mode, r = jobs[index]
                remaining = args.budget_seconds - (time.perf_counter() - start)
                estimate = durations.get(n, previous * (2 if n == 2**23 else 1))
                reserve = max(12, 1.4 * estimate)
                if remaining < reserve:
                    summary[stage].append(
                        {
                            "status": "not_run_budget",
                            "N": n,
                            "q_true": q,
                            "mode": mode,
                            "r": r,
                            "both_q_modes": True,
                            "remaining_seconds": remaining,
                            "job_index": index,
                        }
                    )
                else:
                    future = executor.submit(
                        finite_dataset,
                        n,
                        q,
                        mode,
                        r,
                        20261127 + index,
                        args.backend,
                        args.max_probes,
                    )
                    pending[future] = (stage, index, n, time.perf_counter())
                index += 1
            if not pending:
                break
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                stage, job_index, n, cell_start = pending.pop(future)
                rows = future.result()
                for row in rows:
                    row["job_index"] = job_index
                summary[stage].extend(rows)
                previous = time.perf_counter() - cell_start
                durations[n] = previous
            save()
    for stage in ("stage2", "stage3"):
        summary[stage].sort(key=lambda row: row["job_index"])
    summary["analysis"] = campaign_analysis(summary)
    save()
    print(
        f"Saved {args.output.resolve()} in {summary['elapsed_seconds']:.3f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
