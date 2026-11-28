"""Exploratory blind null/wave ensemble and Route C; never confirmatory.

N counts the entire independent event set; expensive probe counts are capped.
Targets, truth q, and mode labels are used only by generation and evaluation.
No null threshold is tested on held-out data in this small exploratory run.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402
from run_subspacing_pilot import finite_dataset, json_safe  # noqa: E402
from scipy.stats import norm, t  # noqa: E402
from subspacing_observer import (  # noqa: E402
    BackgroundModel,
    bump_window,
    local_pair_correlator,
    tangent_project,
)


def wilson(hits, total):
    """Two-sided 95% Wilson score interval, including endpoint samples."""
    if total == 0:
        return [None, None]
    z = norm.ppf(0.975)
    p = hits / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * np.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return [float(max(0, center - half)), float(min(1, center + half))]


def separation(wave, background):
    """Independent-sample Welch t interval for the difference of means."""
    wave, background = np.asarray(wave, float), np.asarray(background, float)
    if min(len(wave), len(background)) < 2:
        return {"status": "insufficient_samples"}
    variances = np.array(
        [np.var(wave, ddof=1) / len(wave), np.var(background, ddof=1) / len(background)]
    )
    variance = variances.sum()
    df = (
        variance**2
        / np.sum(variances**2 / (np.array([len(wave), len(background)]) - 1))
        if variance
        else float("inf")
    )
    delta = float(wave.mean() - background.mean())
    half = float(t.ppf(0.975, df) * np.sqrt(variance))
    return dict(
        wave_mean=float(wave.mean()),
        background_mean=float(background.mean()),
        difference=delta,
        ci95=[delta - half, delta + half],
        degrees_of_freedom=float(df),
        method="two-sided 95% Welch t; small-sample normal approximation, unadjusted",
    )


def error_statistics(values):
    values = np.asarray(values, float)
    if not len(values):
        return {"status": "no_samples"}
    # Resampling seed affects interval Monte Carlo only, never the datasets.
    draws = np.random.default_rng(10841).choice(values, (10000, len(values)))
    return dict(
        median=float(np.median(values)),
        maximum=float(values.max()),
        median_ci95=np.quantile(np.median(draws, axis=1), [0.025, 0.975]).tolist(),
        sample_maximum_bootstrap_ci95=np.quantile(
            draws.max(axis=1), [0.025, 0.975]
        ).tolist(),
        method="two-sided 95% percentile bootstrap, 10000 resamples",
        limit=(
            "Maximum interval describes resampled sample maxima; "
            "it cannot bound population extremes; small-n coverage is uncertain."
        ),
    )


def cell(job):
    part, n, r, mode, seed = job
    frequencies = (
        np.arange(24, 225, 4) if part == 1 else np.arange(24, max(225, 2 * r + 1), 8)
    )
    # Part 2 uses known q to isolate the exploratory pair statistic at low N.
    rows = finite_dataset(
        n,
        4.0,
        mode,
        r,
        seed,
        "gpu",
        8192 if part == 1 else 4096,
        known_modes=(part == 2,),
        frequencies=frequencies,
        diagnostics=part == 2,
    )
    row = rows[0]
    keys = (
        "status",
        "error",
        "N",
        "r",
        "mode",
        "seed",
        "known_q",
        "backend",
        "q_hat",
        "accepted",
        "n_ref",
        "n_probe",
        "probe_pool",
        "max_counted_probes",
        "generation_seconds",
        "observer_seconds",
        "orthogonalize",
        "dtheta_norm",
        "whitened_target_z",
        "whitened_target_snr",
        "whitened_peak_signed_frequency",
        "whitened_peak_z",
        "whitened_snr",
        "shot_sd",
        "c_ref_estimate",
    )
    output = {key: row[key] for key in keys if key in row}
    output.update(part=part, alpha=float(np.log(n) / np.log(r)))
    if row["status"] != "completed":
        return output
    output.update(
        q_err=abs(row["q_hat"] - 4),
        peak_hit=abs(abs(row["whitened_peak_signed_frequency"]) - r) <= 4,
    )
    if part == 2:
        diagnostic = row["diagnostics"]
        points, jac, calibration = (
            diagnostic["points"],
            diagnostic["jac"],
            diagnostic["calibration"],
        )
        shear = tangent_project(
            BackgroundModel(calibration.q).shear_response(
                points, calibration.early, calibration.late
            ),
            jac,
        )
        residual = tangent_project(diagnostic["residual"], jac)
        weights = bump_window(points)
        denominator = max(
            float(np.sum(weights * np.sum(abs(shear) ** 2, axis=1))), 1e-30
        )
        z = weights * np.sum(shear.conj() * residual, axis=1) / denominator
        output["local_pairs"] = []
        target = int(np.flatnonzero(frequencies == r)[0])
        off = abs(frequencies - r) > 20
        for dc in (0.02, 0.05, 0.12, np.inf):
            plus = local_pair_correlator(points[:, 0], z, frequencies, dc)
            minus = local_pair_correlator(points[:, 0], z, -frequencies, dc)
            # Average chiralities without selecting a target-favored sign.
            energy = (plus + minus) / 2
            sd = float(np.std(energy[off], ddof=1))
            mean = float(np.mean(energy[off]))
            normalized = energy / max(sd, 1e-30)
            output["local_pairs"].append(
                dict(
                    dc="inf" if np.isinf(dc) else dc,
                    target_energy=float(energy[target]),
                    off_peak_sd=sd,
                    off_peak_mean=mean,
                    whitened_target_energy=float(normalized[target]),
                    centered_target_energy=float(
                        (energy[target] - mean) / max(sd, 1e-30)
                    ),
                    frequencies=frequencies,
                    whitened_energy=normalized,
                )
            )
    return output


def aggregate(rows):
    output = {"part1": [], "part2": []}
    for part, n, r in sorted({(row["part"], row["N"], row["r"]) for row in rows}):
        selected = [
            row for row in rows if (row["part"], row["N"], row["r"]) == (part, n, r)
        ]
        arms = {
            mode: [
                row
                for row in selected
                if row["mode"] == mode and row["status"] == "completed"
            ]
            for mode in ("background", "circular_wave")
        }
        bg, wave = arms["background"], arms["circular_wave"]
        item = dict(N=n, r=r, alpha=float(np.log(n) / np.log(r)), arms={})
        for mode, group in arms.items():
            hits = sum(row["peak_hit"] for row in group)
            item["arms"][mode] = dict(
                completed=len(group),
                attempted=sum(row["mode"] == mode for row in selected),
                peak_hits=hits,
                peak_hit_rate=hits / len(group) if group else None,
                peak_hit_ci95=wilson(hits, len(group)),
                q_err=error_statistics([row["q_err"] for row in group]),
            )
        item["whitened_target_snr_separation"] = separation(
            [row["whitened_target_snr"] for row in wave],
            [row["whitened_target_snr"] for row in bg],
        )
        if part == 1 and bg:
            allowed = len(bg) // 8
            threshold = sorted(row["whitened_peak_z"] for row in bg)[-allowed - 1]
            false = sum(row["whitened_peak_z"] > threshold for row in bg)
            true = sum(row["whitened_peak_z"] > threshold for row in wave)
            item["empirical_threshold"] = dict(
                whitened_peak_z=threshold,
                comparison="strictly greater than",
                false_positives=false,
                background_n=len(bg),
                false_positive_rate=false / len(bg),
                false_positive_ci95=wilson(false, len(bg)),
                wave_detections=true,
                wave_n=len(wave),
                wave_detection_rate=true / len(wave) if wave else None,
                wave_detection_ci95=wilson(true, len(wave)),
                limit=(
                    "Threshold selected on these same 8 or 6 null seeds; "
                    "Wilson intervals do not correct threshold selection. "
                    "No held-out false-positive control."
                ),
            )
        if part == 2:
            item["local_pairs"] = []
            for i, dc in enumerate((0.02, 0.05, 0.12, "inf")):
                stats = separation(
                    [row["local_pairs"][i]["whitened_target_energy"] for row in wave],
                    [row["local_pairs"][i]["whitened_target_energy"] for row in bg],
                )
                raw = separation(
                    [row["local_pairs"][i]["target_energy"] for row in wave],
                    [row["local_pairs"][i]["target_energy"] for row in bg],
                )
                centered = separation(
                    [row["local_pairs"][i]["centered_target_energy"] for row in wave],
                    [row["local_pairs"][i]["centered_target_energy"] for row in bg],
                )
                item["local_pairs"].append(
                    dict(
                        dc=dc,
                        whitened_energy_separation=stats,
                        raw_energy_separation=raw,
                        centered_energy_separation=centered,
                        nominal_positive_ci=stats.get("ci95", [-float("inf")])[0] > 0,
                    )
                )
            scale = n / r**2
            item["route_c_window"] = dict(
                N_lambda_squared=scale,
                lower=1 / scale,
                upper=scale,
                strictly_nonempty=1 / scale < scale,
                limit=(
                    "Finite numerical bounds do not establish the asymptotic "
                    "much-less-than inequalities."
                ),
            )
        output[f"part{part}"].append(item)
    output["interpretation"] = (
        "Exploratory, not preregistered or confirmatory. Four seeds per Route C arm, "
        "correlated bins, 28 dc/cell comparisons, no multiplicity correction. "
        "A nominal positive interval is a candidate only; failure to see excess "
        "does not confirm the asymptotic pinch or prove impossibility. "
        "Off-peak SD normalization is blind to wave labels but target-conditioned "
        "evaluation excludes +/-20. whitened_peak_z and target_z are amplitude/noise "
        "scales, not Gaussian significance z scores. Shared reference noise and "
        "fitted nuisance parameters can correlate probe residuals; sum|z|^2 is the "
        "independent-probe shot proxy, not a full sampling-variance estimate."
    )
    return output


def main():
    start = time.perf_counter()
    budget = 180.0
    jobs = []
    for n, r, repetitions in ((2**22, 64, 8), (2**23, 96, 6)):
        jobs.extend(
            (1, n, r, mode)
            for mode in ("background", "circular_wave")
            for _ in range(repetitions)
        )
    for n, r in [(2**20, r) for r in (64, 128, 256, 512, 1024)] + [
        (2**18, 512),
        (2**18, 1024),
    ]:
        jobs.extend(
            (2, n, r, mode)
            for mode in ("background", "circular_wave")
            for _ in range(4)
        )
    jobs = [(*job, 740000 + i) for i, job in enumerate(jobs)]
    assert len({job[-1] for job in jobs}) == 84
    source = Path(__file__).resolve().parent
    summary = dict(
        kind="exploratory, not preregistered or confirmatory",
        budget_seconds=budget,
        seed_policy=(
            "Distinct deterministic RNG seeds for all 84 independent "
            "datasets; no shared event arrays."
        ),
        part1_q="unknown",
        part2_q="known",
        q_true=4.0,
        probe_caps={"part1": 8192, "part2": 4096},
        budget_pilot=(
            "A 32768-probe run completed 21/84 cells in 177.018s; "
            "retained as whitened_ensemble_32768_budget_incomplete.json. "
            "Probe caps reduced for runtime, not selected for detection."
        ),
        backend="gpu",
        workers=4,
        source_sha256={
            name: hashlib.sha256((source / name).read_bytes()).hexdigest()
            for name in (
                "subspacing_observer.py",
                "subspacing_wave_sim.py",
                "run_subspacing_pilot.py",
                "run_whitened_ensemble.py",
            )
        },
        rows=[],
    )
    import cupy as cp

    summary["gpu"] = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    output = Path("outputs/whitened_ensemble_summary.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    executor = ProcessPoolExecutor(max_workers=4)
    pending = {executor.submit(cell, job): job for job in jobs}
    try:
        while pending:
            remaining = budget - 3 - (time.perf_counter() - start)
            if remaining <= 0:
                break
            done, _ = wait(
                pending, timeout=min(remaining, 10), return_when=FIRST_COMPLETED
            )
            for future in done:
                job = pending.pop(future)
                try:
                    row = future.result()
                except Exception as error:
                    row = dict(
                        part=job[0],
                        N=job[1],
                        r=job[2],
                        mode=job[3],
                        seed=job[4],
                        status="failed",
                        error=f"{type(error).__name__}: {error}",
                    )
                summary["rows"].append(row)
                print(
                    f"ensemble {len(summary['rows'])}/84: {job[:4]} {row['status']}",
                    flush=True,
                )
        if pending:
            for job in pending.values():
                summary["rows"].append(
                    dict(
                        part=job[0],
                        N=job[1],
                        r=job[2],
                        mode=job[3],
                        seed=job[4],
                        status="budget_incomplete",
                    )
                )
            executor.terminate_workers()
        else:
            executor.shutdown()
    finally:
        executor.shutdown(wait=False, cancel_futures=True)
    summary["rows"].sort(key=lambda row: row["seed"])
    summary["analysis"] = aggregate(summary["rows"])
    summary["completed"] = sum(row["status"] == "completed" for row in summary["rows"])
    summary["elapsed_seconds"] = time.perf_counter() - start
    output.write_text(json.dumps(json_safe(summary), indent=2, allow_nan=False) + "\n")
    print(
        f"Saved {output}: {summary['completed']}/84 completed "
        f"in {summary['elapsed_seconds']:.3f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
