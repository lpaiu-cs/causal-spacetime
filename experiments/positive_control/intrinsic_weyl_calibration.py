"""Calibration of the nominal covariance at high event counts, on flat replicates.

Frozen design: docs/paper/paper_a/intrinsic_weyl/calibration_protocol.md.
Generation reuses the benchmark generator and landmark rule unchanged. Each
worker reads one order archive, fits the fifth- and sixth-degree jets, writes
the pure fits, and only then loads the generator truth.
"""

import argparse
import json
import os
import platform
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from replay_intrinsic_weyl_threads import ROOT, THREAD_VARIABLES, sha, write_json

PROTOCOL = "docs/paper/paper_a/intrinsic_weyl/calibration_protocol.md"
LEVELS = (1048576, 16777216)
REPLICATES, SEED = 64, 2026112100
SECONDARY, SECONDARY_SEED = 16, 2026113100
ORDERS = (5, 6)
WORKERS, THREADS = 6, 2
LOW, HIGH = 0.8, 1.25
MINIMUM_COMPLETE = 60
ORACLES = (
    "profile_plane_wave",
    "intrinsic_weyl_bombelli_benchmark",
    "intrinsic_weyl_joint_audit",
    "rotating_plane_wave",
    "window_response",
    "p14_plane_wave",
)


def assignments():
    from intrinsic_weyl_bombelli_benchmark import benchmark_amplitude

    primary = [("flat", "0", SEED + r) for r in range(REPLICATES)]
    secondary = [
        ("constant", str(benchmark_amplitude()), SECONDARY_SEED + r)
        for r in range(SECONDARY)
    ]
    return primary + secondary


def calibration(statistics, replicates=REPLICATES):
    """Frozen verdict for one jet and level from the completed truth statistics."""
    import numpy as np
    from scipy.stats import chi2, kstest

    values = np.asarray(statistics, dtype=float)
    result = {"completed": len(values), "planned": replicates}
    if len(values) < MINIMUM_COMPLETE * replicates / REPLICATES:
        result["verdict"] = "INCOMPLETE"
        return result
    dof = 10 * len(values)
    kappa = float(values.sum() / dof)
    interval = [
        float(values.sum() / chi2.ppf(0.975, dof)),
        float(values.sum() / chi2.ppf(0.025, dof)),
    ]
    rng = np.random.default_rng(1)
    means = rng.choice(values, size=(10000, len(values)), replace=True).mean(axis=1)
    result.update(
        kappa=kappa,
        interval_95=interval,
        verdict="CALIBRATED"
        if LOW <= interval[0] and interval[1] <= HIGH
        else ("NOT_CALIBRATED"),
        direction="UNDERSTATED"
        if interval[0] > 1
        else ("OVERSTATED" if interval[1] < 1 else "CONSISTENT_WITH_ONE"),
        bootstrap_interval_95=[
            float(np.percentile(means, 2.5) / 10),
            float(np.percentile(means, 97.5) / 10),
        ],
        exceedances_of_q_star=int((values > chi2.ppf(0.999, 10)).sum()),
        ks_distance_scaled=float(kstest(values / kappa, chi2(10).cdf).statistic),
        mean=float(values.mean()),
        median=float(np.median(values)),
    )
    return result


def selfcheck():
    """Fixed inputs with known verdicts."""
    checks = {
        "nominal": calibration([10.0] * 64),
        "hinted": calibration([13.4] * 64),
        "edge": calibration([11.5] * 64),
        "short": calibration([10.0] * 59),
    }
    assert checks["nominal"]["verdict"] == "CALIBRATED"
    assert checks["nominal"]["direction"] == "CONSISTENT_WITH_ONE"
    assert checks["hinted"]["verdict"] == "NOT_CALIBRATED"
    assert checks["hinted"]["direction"] == "UNDERSTATED"
    assert checks["edge"]["verdict"] == "NOT_CALIBRATED"
    assert checks["edge"]["direction"] == "UNDERSTATED"
    assert checks["short"]["verdict"] == "INCOMPLETE"
    return {name: check.get("verdict") for name, check in checks.items()}


def evaluate(fits, truth_path, n):
    """Truth statistics after the pure fits are saved."""
    import numpy as np
    from intrinsic_weyl_bombelli_benchmark import build_wave
    from intrinsic_weyl_joint import axis_signs, confidence_ball, gauge
    from intrinsic_weyl_joint_audit import aligned_tensor_error

    with np.load(truth_path, allow_pickle=False) as saved:
        points, kind = saved["points"], str(saved["kind"].item())
        amplitude = str(saved["amplitude"].item())
    value = float(amplitude)
    if value:
        wave = build_wave(kind, amplitude)
        normal = wave.normal_coordinates(points)
        expected = np.asarray(
            aligned_tensor_error(
                {
                    "coefficients": np.zeros(10),
                    "coordinates_inferred_from_counts": gauge(normal),
                },
                normal,
                value,
            )["aligned_truth_coefficients"]
        )
        candidates = (expected, axis_signs([1, 1, -1]) * expected)
    else:
        candidates = (np.zeros(10),)
    out = {"kind": kind, "amplitude": amplitude, "measurement_n": n}
    for name, fit in fits.items():
        ball = confidence_ball(fit)
        c = np.asarray(ball["coefficients_xy_positive"])
        s = np.asarray(ball["nominal_covariance_xy_positive"])
        values = [float((c - t) @ np.linalg.solve(s, c - t)) for t in candidates]
        out[name] = {
            "q": float(c @ np.linalg.solve(s, c)),
            "truth_statistic": min(values),
            "truth_statistic_candidates": values,
            "tidal_norm_estimate": ball["tidal_norm"],
            "statistical_radius_95": ball["statistical_radius"],
            "optimizer_success": fit["optimizer_success"],
            "projected_curvature_rank": fit["projected_curvature_rank"],
            "nuisance_rank": fit["nuisance_rank"],
            "nuisance_parameters": fit["nuisance_parameters"],
            "fit_seconds": fit["fit_seconds"],
        }
    return out


def worker(archive, truth, output):
    import intrinsic_weyl_joint as joint
    import numpy as np
    from intrinsic_weyl_gpu_counts import clock_second_moment
    from intrinsic_weyl_high_count import peak_working_set

    joint._clock_second_moment = clock_second_moment
    wall, cpu = time.perf_counter(), time.process_time()
    record = {
        "archive_sha256": sha(archive),
        "thread_environment": {key: os.environ.get(key) for key in THREAD_VARIABLES},
    }
    try:
        assert not any(name in sys.modules for name in ORACLES)
        counts, mask, covariance, ids, n = joint.read_order_clocks(archive)
        fits = {}
        for order in ORDERS:
            begin = time.perf_counter()
            fit = joint.fit_joint(
                counts,
                clock_mask=mask,
                clock_covariance=covariance,
                model_order=order,
                square_correction=True,
            )
            fit["fit_seconds"] = time.perf_counter() - begin
            fits[f"order_{order}"] = fit
        imported = [name for name in ORACLES if name in sys.modules]
        if imported:
            raise RuntimeError(f"Oracle imported before pure-fit boundary: {imported}")
        write_json(
            output.with_name(output.stem + "_pure.json"),
            {
                "fits": fits,
                "measurement_n": n,
                "selected_ids": ids.tolist(),
                "smallest_measurement_clock": int(counts[mask].min()),
                "archive_sha256": record["archive_sha256"],
                "oracle_modules_before_return": imported,
            },
        )
        record["evaluation"] = evaluate(fits, truth, n)
    except (ValueError, ArithmeticError, MemoryError, np.linalg.LinAlgError) as error:
        record["failure"] = f"{type(error).__name__}: {error}"
    record.update(
        wall_seconds=time.perf_counter() - wall,
        process_cpu_seconds=time.process_time() - cpu,
        peak_working_set_bytes=peak_working_set(),
    )
    write_json(output, record)


def run(destination):
    import cupy as cp
    import numpy as np
    import scipy
    from intrinsic_weyl_bombelli_benchmark import generate, lf_sha, shutil_free
    from intrinsic_weyl_gpu_counts import selfcheck as gpu_selfcheck
    from intrinsic_weyl_joint import selfcheck_jet, selfcheck_second
    from profile_plane_wave import selfcheck as wave_selfcheck

    if os.name != "nt":
        raise RuntimeError("The resource accounting of this runner uses Windows")
    if destination.exists():
        raise FileExistsError("Preserve every existing calibration run")
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip()
    if dirty:
        raise RuntimeError("Commit the protocol and runner before calibration data")
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    if shutil_free(destination.parent) < 70 * 2**30:
        raise RuntimeError("At least 70 GiB free storage is required")
    here = Path(__file__).resolve().parent
    sources = [
        here / name
        for name in (
            Path(__file__).name,
            "intrinsic_weyl_bombelli_benchmark.py",
            "profile_plane_wave.py",
            "replay_intrinsic_weyl_threads.py",
            "intrinsic_weyl_joint.py",
            "intrinsic_weyl.py",
            "intrinsic_weyl_square.py",
            "intrinsic_weyl_gpu_counts.py",
            "intrinsic_weyl_high_count.py",
            "intrinsic_weyl_joint_audit.py",
        )
    ] + [ROOT / PROTOCOL]
    hashes = {str(p.relative_to(ROOT)).replace("\\", "/"): lf_sha(p) for p in sources}
    checked = {
        "generator": wave_selfcheck(),
        "fifth_degree_jet": selfcheck_jet(),
        "sixth_degree_jet": selfcheck_second(),
        "gpu_second_moment": gpu_selfcheck(),
        "verdict": selfcheck(),
    }
    destination.mkdir(parents=True)
    started, lock = time.time(), threading.Lock()

    def status(**fields):
        with lock:
            temporary = destination / ".status.tmp"
            temporary.write_text(
                json.dumps(
                    {"pid": os.getpid(), "started_unix": started, **fields}, indent=2
                )
                + "\n"
            )
            os.replace(temporary, destination / "status.json")

    environment = os.environ.copy()
    for key in THREAD_VARIABLES:
        environment.pop(key, None)
    environment.update(OPENBLAS_NUM_THREADS=str(THREADS), OMP_NUM_THREADS=str(THREADS))
    plan = assignments()
    manifest = {
        "code_version": stamp,
        "worktree_clean_at_start": True,
        "source_sha256_lf": hashes,
        "selfcheck": checked,
        "assignments": [
            {"replicate": i, "kind": k, "amplitude": a, "seed": s}
            for i, (k, a, s) in enumerate(plan)
        ],
        "measurement_levels": LEVELS,
        "orders": ORDERS,
        "workers": WORKERS,
        "worker_thread_environment": {
            key: environment.get(key) for key in THREAD_VARIABLES
        },
        "tolerance_on_kappa": [LOW, HIGH],
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cupy": cp.__version__,
            "device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        },
        "scope": "Flat-replicate calibration of the nominal covariance; two levels",
    }
    write_json(destination / "manifest.json", manifest)

    def invoke(archive, truth, output):
        with output.with_suffix(".log").open("xb") as log:
            code = subprocess.call(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--worker",
                    str(archive),
                    "--truth",
                    str(truth),
                    "--output",
                    str(output),
                ],
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW,
            )
        if code:
            return {"failure": f"Worker exit {code}"}
        return json.loads(output.read_text())

    futures = []
    try:
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            for index, (kind, amplitude, seed) in enumerate(plan):
                generation, truth = generate(
                    destination, index, kind, amplitude, seed, LEVELS, status
                )
                for n in LEVELS:
                    archive = destination / f"case_{index}_n{n}.npz"
                    key = {"replicate": index, "kind": kind, "seed": seed, "n": n}
                    if archive.exists():
                        output = archive.with_suffix(".json")
                        job = pool.submit(invoke, archive, truth, output)
                        futures.append((key, job))
                    else:
                        futures.append(
                            (key, {"failure": generation.get("failure", "No archive")})
                        )
                status(state="RUNNING", generated=index + 1, planned=len(plan))
            rows = []
            for key, future in futures:
                result = future if isinstance(future, dict) else future.result()
                rows.append({**key, **result})
        verdicts = {}
        for kind, planned in (("flat", REPLICATES), ("constant", SECONDARY)):
            for order in ORDERS:
                for n in LEVELS:
                    values = [
                        row["evaluation"][f"order_{order}"]["truth_statistic"]
                        for row in rows
                        if row["kind"] == kind
                        and row["n"] == n
                        and "evaluation" in row
                        and "failure" not in row
                    ]
                    verdicts[f"{kind}/order_{order}/n{n}"] = calibration(
                        values, planned
                    )
        for name, digest in hashes.items():
            if lf_sha(ROOT / name) != digest:
                raise RuntimeError(f"A bound source changed during execution: {name}")
        write_json(
            destination / "summary.json",
            {
                **manifest,
                "rows": rows,
                "calibration": verdicts,
                "primary_population": "flat",
                "elapsed_seconds": time.time() - started,
            },
        )
        status(
            state="COMPLETED",
            verdicts={
                k: (v["verdict"], v.get("kappa"), v.get("interval_95"))
                for k, v in verdicts.items()
            },
        )
    except BaseException as error:
        status(state="STOPPED", error=f"{type(error).__name__}: {error}")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        print(json.dumps(selfcheck(), indent=2))
    elif args.worker:
        worker(args.worker, args.truth, args.output)
    else:
        run(args.destination)
