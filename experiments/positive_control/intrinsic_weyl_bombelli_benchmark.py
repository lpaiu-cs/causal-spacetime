"""Ibrahim-Bombelli plane-wave benchmark for the unchanged order-only Weyl inverse.

The geometry is the exact vacuum plane wave whose tidal profile along the
central observer is the linearized plus-polarized wave of arXiv:2611.09056 at
h0 = 1e-2, on a causal diamond of height lambda/2 centered on a crest. Only
the generator is new; landmark selection, clocks, covariance, the fifth-degree
inverse with curvature-square terms and the confidence construction are the
frozen ones. See docs/paper/paper_a/intrinsic_weyl/bombelli_benchmark_protocol.md.
"""

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

from replay_intrinsic_weyl_threads import ROOT, THREAD_VARIABLES, sha, write_json

PROTOCOL = "docs/paper/paper_a/intrinsic_weyl/bombelli_benchmark_protocol.md"
H0 = "0.01"
SELECTION_N, CHUNK, BATCH = 4096, 4096, 65536
LEVELS = (1048576, 4194304, 16777216, 67108864)
SEED = 2026111100
# kind, seed offset. All four are co-required at every level.
CASES = (("wave", 0), ("wave", 1), ("flat", 2), ("flat", 3))
FALSE_POSITIVE = 1e-3
# Engineering population: an audited constant-profile case, disjoint seeds.
ENGINEERING_LEVELS = (1048576, 4194304, 16777216)
ENGINEERING_CASES = (("constant", 90),)
ENGINEERING_AMPLITUDE = "0.25"
ORACLES = (
    "profile_plane_wave",
    "rotating_plane_wave",
    "rotating_plane_wave_gpu",
    "window_response",
    "p14_plane_wave",
    "intrinsic_weyl_rotating_pilot",
    "intrinsic_weyl_joint_audit",
    "intrinsic_weyl_joint_pilot",
)


def benchmark_amplitude():
    """A0 = pi^2 h0 / 2 in outer-tip units; see the protocol for the map."""
    from decimal import Decimal, localcontext

    from profile_plane_wave import PRECISION, decimal_pi

    with localcontext() as context:
        context.prec = PRECISION
        return decimal_pi() ** 2 * Decimal(H0) / 2


def build_wave(kind, amplitude):
    from profile_plane_wave import constant, half_cycle

    if kind == "wave":
        return half_cycle(amplitude)
    if kind in ("flat", "constant"):
        return constant(amplitude)
    raise ValueError("Unknown benchmark geometry")


def threshold():
    from scipy.stats import chi2

    return float(chi2.ppf(1 - FALSE_POSITIVE, 10))


def verdict(rows):
    """Frozen existence gate for one measurement level; four rows required."""
    if len(rows) != len(CASES) or [r["kind"] for r in rows] != [k for k, _ in CASES]:
        raise ValueError("All four assigned cases are required, in order")
    if any("failure" in r or "detection" not in r for r in rows):
        return "INCOMPLETE"
    limit = threshold()
    wave = [r["detection"]["q"] > limit for r in rows if r["kind"] == "wave"]
    flat = [r["detection"]["q"] > limit for r in rows if r["kind"] == "flat"]
    if any(flat):
        return "CONTROL_FAILED"
    if all(wave):
        return "DETECTED"
    return "MIXED" if any(wave) else "NOT_DETECTED"


def selfcheck():
    def row(kind, q=None, failed=False):
        value = {"kind": kind}
        if failed:
            value["failure"] = "kept"
        elif q is not None:
            value["detection"] = {"q": q}
        return value

    limit = threshold()
    assert abs(limit - 29.588298445074) < 1e-9
    high, low = limit + 1, limit - 1
    assert (
        verdict(
            [row("wave", high), row("wave", high), row("flat", low), row("flat", low)]
        )
        == "DETECTED"
    )
    assert (
        verdict(
            [row("wave", high), row("wave", low), row("flat", low), row("flat", low)]
        )
        == "MIXED"
    )
    assert (
        verdict(
            [row("wave", low), row("wave", low), row("flat", low), row("flat", low)]
        )
        == "NOT_DETECTED"
    )
    assert (
        verdict(
            [row("wave", high), row("wave", high), row("flat", high), row("flat", low)]
        )
        == "CONTROL_FAILED"
    )
    assert (
        verdict(
            [
                row("wave", high),
                row("wave", failed=True),
                row("flat", low),
                row("flat", low),
            ]
        )
        == "INCOMPLETE"
    )
    assert (
        verdict(
            [row("wave", limit), row("wave", high), row("flat", low), row("flat", low)]
        )
        == "MIXED"
    )
    try:
        verdict([row("wave", high)] * 3)
    except ValueError:
        pass
    else:
        raise AssertionError("A missing assignment was accepted")
    return {
        "threshold_chi2_10": limit,
        "strict_exceedance_required": True,
        "control_exceedance_invalidates": True,
        "failure_cannot_disappear": True,
        "all_four_required": True,
    }


def pattern(inference, direction):
    """Expected and realized signal-to-noise along one fixed tidal pattern."""
    import numpy as np
    from intrinsic_weyl_joint import TIDAL_ROOT

    target = TIDAL_ROOT @ direction
    size = float(np.linalg.norm(target))
    unit = target / size
    estimate = TIDAL_ROOT @ np.asarray(inference["coefficients_xy_positive"])
    covariance = (
        TIDAL_ROOT
        @ np.asarray(inference["nominal_covariance_xy_positive"])
        @ TIDAL_ROOT.T
    )
    sigma = float(np.sqrt(unit @ covariance @ unit))
    return {
        "truth_amplitude": size,
        "estimate": float(unit @ estimate),
        "nominal_sigma": sigma,
        "expected_snr": size / sigma,
        "realized_z": float(unit @ estimate) / sigma,
    }


def evaluate_fit(fit, clocks, truth_path, n):
    """Oracle evaluation after the pure fit has been saved."""
    import numpy as np
    from intrinsic_weyl_joint import (
        axis_signs,
        confidence_ball,
        fit_joint,
        gauge,
        tidal_distance,
    )
    from intrinsic_weyl_joint_audit import aligned_tensor_error
    from scipy.stats import chi2

    coefficients = np.asarray(fit["coefficients"])
    covariance = np.asarray(fit["nominal_component_covariance"])
    q = float(coefficients @ np.linalg.solve(covariance, coefficients))
    inference = confidence_ball(fit)
    detection = {
        "q": q,
        "threshold": threshold(),
        "nominal_p_value": float(chi2.sf(q, 10)),
        "tidal_norm_estimate": inference["tidal_norm"],
        "largest_tidal_variance": inference["largest_tidal_variance"],
        "variance_times_N": n * inference["largest_tidal_variance"],
    }
    with np.load(truth_path, allow_pickle=False) as saved:
        points, kind = saved["points"], str(saved["kind"].item())
        amplitude = str(saved["amplitude"].item())
    wave, value = build_wave(kind, amplitude), float(amplitude)
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
    counts, mask, clock_covariance, _, _ = clocks
    volume, coarse = wave.volumes(points), wave.volumes(points, 32)
    refinement = float(np.max(abs(coarse[mask] / volume[mask] - 1)))
    if refinement >= 1e-10:
        raise ArithmeticError("Interval quadrature did not stabilize")
    continuum = fit_joint(
        counts[0, -1] * volume / volume[0, -1],
        clock_mask=mask,
        clock_covariance=clock_covariance,
        model_order=5,
        square_correction=True,
    )
    mean = confidence_ball(continuum)
    estimate = inference["coefficients_xy_positive"]
    error = tidal_distance(estimate, expected)
    bias = tidal_distance(mean["coefficients_xy_positive"], expected)
    evaluation = {
        "truth_coefficients_xy_positive": expected.tolist(),
        "truth_norm": 2 * value,
        "full_tidal_error_mod_reflection": error,
        "relative_full_tidal_error": error / (2 * value) if value else None,
        "continuum_bias_mod_reflection": bias,
        "relative_continuum_bias": bias / (2 * value) if value else None,
        "continuum_tidal_norm": mean["tidal_norm"],
        "confidence_radius_95": inference["radius"],
        "nominal_ball_contains_truth": bool(error <= inference["radius"]),
        "quadrature_relative_change": refinement,
    }
    if value:
        # The residual z reflection is unobservable; keep the less favorable one.
        candidates = (expected, axis_signs([1, 1, -1]) * expected)
        electric, magnetic = (
            np.r_[np.ones(5), np.zeros(5)],
            np.r_[np.zeros(5), np.ones(5)],
        )
        rows = {}
        for name, part in (
            ("electric", electric),
            ("magnetic", magnetic),
            ("full", electric + magnetic),
        ):
            options = [pattern(inference, part * candidate) for candidate in candidates]
            rows[name] = min(options, key=lambda item: item["expected_snr"])
            rows[name]["events_for_expected_snr_one"] = (
                n / rows[name]["expected_snr"] ** 2
            )
        evaluation["patterns"] = rows
    return {
        "fit": fit,
        "inference": inference,
        "continuum_fit": continuum,
        "detection": detection,
        "evaluation": evaluation,
    }


def worker(archive, truth, output):
    """Pure order-only fit first; oracle modules load only after it is saved."""
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
    original_load, opened = np.load, []

    def guard(selected, *args, **options):
        if Path(selected).resolve() != archive.resolve():
            raise ValueError("Only the assigned order archive may be opened")
        opened.append(str(Path(selected).resolve()))
        return original_load(selected, *args, **options)

    try:
        assert not any(name in sys.modules for name in ORACLES)
        np.load = guard
        try:
            fit, clocks = joint.fit_order_archive(
                archive, with_clocks=True, model_order=5, square_correction=True
            )
        finally:
            np.load = original_load
        imported = [name for name in ORACLES if name in sys.modules]
        if imported:
            raise RuntimeError(f"Oracle imported before pure-fit boundary: {imported}")
        record["finite_fit_wall_seconds"] = time.perf_counter() - wall
        write_json(
            output.with_name(output.stem + "_pure.json"),
            {
                "fit": fit,
                "opened_archives": opened,
                "oracle_modules_before_return": imported,
                "inverse_module": joint.__file__,
                "archive_sha256": sha(archive),
            },
        )
        record["full"] = evaluate_fit(fit, clocks, truth, fit["measurement_n"])
    except (ValueError, ArithmeticError, MemoryError, np.linalg.LinAlgError) as error:
        record["failure"] = f"{type(error).__name__}: {error}"
    record.update(
        wall_seconds=time.perf_counter() - wall,
        process_cpu_seconds=time.process_time() - cpu,
        peak_working_set_bytes=peak_working_set(),
    )
    write_json(output, record)


def generate(folder, index, kind, amplitude, seed, levels, status):
    """Selection order, landmarks, and nested measurement prefixes for one case."""
    import cupy as cp
    import numpy as np
    from intrinsic_weyl_high_count import peak_working_set
    from intrinsic_weyl_joint import selection

    wave = build_wave(kind, amplitude)
    wall, cpu = time.perf_counter(), time.process_time()
    generation = {
        "selection_events_drawn": 0,
        "sampled_events": 0,
        "classified_events": 0,
    }
    truth_path = folder / f"truth_{index}.npz"
    try:
        status(state="RUNNING", phase="selection", case=index)
        rng = np.random.default_rng([seed, 0])
        points = np.vstack(
            ([0, 0.5, 0, 0], wave.sample(SELECTION_N, rng), [1, -0.5, 0, 0])
        )
        generation["selection_events_drawn"] = SELECTION_N
        # Labels are randomized independently of coordinates; only the two
        # supplied tips have distinguished IDs.
        permutation = np.r_[
            0, rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
        ]
        points = points[permutation]
        order = wave.order(points)
        try:
            ids, _ = selection(order, design="nested")
        except ValueError:
            np.savez(
                folder / f"case_{index}_selection_failed.npz",
                selection_order=order,
                selection_design="nested",
            )
            raise
        selected = {
            "selection_order": order,
            "selection_design": np.array("nested"),
            "selected_ids": ids,
        }
        points = points[ids]
        np.savez(
            truth_path,
            points=points,
            kind=np.array(kind),
            amplitude=np.array(str(amplitude)),
        )
        generation["landmarks_including_tips"] = len(points)
        top = levels[-1]
        after = np.empty((len(points), top // 8), dtype=np.uint8)
        before = np.empty_like(after)
        rng = np.random.default_rng([seed, 1])
        for start in range(0, top, BATCH):
            end = min(start + BATCH, top)
            events = np.vstack(
                [wave.sample(CHUNK, rng) for _ in range((end - start) // CHUNK)]
            )
            generation["sampled_events"] += len(events)
            future, past = wave.packed_relations(points, events, cp)
            after[:, start // 8 : end // 8], before[:, start // 8 : end // 8] = (
                future,
                past,
            )
            generation["classified_events"] = end
            if end in levels:
                with (folder / f"case_{index}_n{end}.npz").open("xb") as stream:
                    np.savez(
                        stream,
                        **selected,
                        after=after[:, : end // 8],
                        before=before[:, : end // 8],
                        measurement_n=end,
                    )
            if end % (64 * BATCH) == 0 or end == top:
                status(
                    state="RUNNING",
                    phase="generate",
                    case=index,
                    events=end,
                    target=top,
                )
    except (ValueError, ArithmeticError, MemoryError) as error:
        generation["failure"] = f"{type(error).__name__}: {error}"
    # Return pooled GPU memory before the worker process needs the device.
    cp.get_default_memory_pool().free_all_blocks()
    generation.update(
        margins_decided_in_decimal_arithmetic=wave.resolved_margins,
        wall_seconds=time.perf_counter() - wall,
        process_cpu_seconds=time.process_time() - cpu,
        peak_working_set_bytes=peak_working_set(),
    )
    if "landmarks_including_tips" in generation:
        generation["directional_predicates_classified"] = (
            2 * generation["landmarks_including_tips"] * generation["classified_events"]
        )
    write_json(folder / f"case_{index}_generation.json", generation)
    return generation, truth_path


def run(destination, engineering=False):
    import cupy as cp
    import numpy as np
    import scipy
    from intrinsic_weyl_gpu_counts import selfcheck as gpu_selfcheck
    from intrinsic_weyl_joint import selfcheck_jet
    from profile_plane_wave import selfcheck as wave_selfcheck

    if os.name != "nt":
        raise RuntimeError("The resource accounting of this runner uses Windows")
    if destination.exists():
        raise FileExistsError("Preserve every existing benchmark run")
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip()
    if dirty and not engineering:
        raise RuntimeError("Commit the protocol and runner before benchmark data")
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    here = Path(__file__).resolve().parent
    sources = [
        here / name
        for name in (
            Path(__file__).name,
            "profile_plane_wave.py",
            "replay_intrinsic_weyl_threads.py",
            "intrinsic_weyl_joint.py",
            "intrinsic_weyl.py",
            "intrinsic_weyl_square.py",
            "intrinsic_weyl_gpu_counts.py",
            "intrinsic_weyl_high_count.py",
            "intrinsic_weyl_joint_audit.py",
            "intrinsic_weyl_joint_pilot.py",
            "rotating_plane_wave.py",
            "window_response.py",
            "p14_plane_wave.py",
        )
    ]
    if not engineering:
        sources.append(ROOT / PROTOCOL)
    source_hashes = {
        str(p.relative_to(ROOT)).replace("\\", "/"): lf_sha(p) for p in sources
    }
    checked = {
        "generator": wave_selfcheck(),
        "jet_inverse": selfcheck_jet(),
        "gpu_second_moment": gpu_selfcheck(),
        "verdict": selfcheck(),
    }
    if engineering:
        cases, levels, amplitudes = (
            ENGINEERING_CASES,
            ENGINEERING_LEVELS,
            {"constant": ENGINEERING_AMPLITUDE},
        )
    else:
        cases, levels = CASES, LEVELS
        amplitudes = {"wave": str(benchmark_amplitude()), "flat": "0"}
    if shutil_free(destination.parent) < 30 * 2**30:
        raise RuntimeError("At least 30 GiB free storage is required")
    destination.mkdir(parents=True)
    started = time.time()

    def status(**fields):
        temporary = destination / ".status.tmp"
        temporary.write_text(
            json.dumps(
                {"pid": os.getpid(), "started_unix": started, **fields}, indent=2
            )
            + "\n"
        )
        os.replace(temporary, destination / "status.json")

    worker_env = os.environ.copy()
    for key in THREAD_VARIABLES:
        worker_env.pop(key, None)
    worker_env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    manifest = {
        "engineering_only": bool(engineering),
        "code_version": stamp,
        "worktree_clean_at_start": not dirty,
        "source_sha256_lf": source_hashes,
        "selfcheck": checked,
        "h0": None if engineering else H0,
        "amplitudes_outer_tip_units": amplitudes,
        "cases": [
            {"case": i, "kind": k, "seed": SEED + offset}
            for i, (k, offset) in enumerate(cases)
        ],
        "measurement_levels": levels,
        "selection_n": SELECTION_N,
        "sampling_chunk": CHUNK,
        "classification_batch": BATCH,
        "false_positive_level_per_test": FALSE_POSITIVE,
        "threshold_chi2_10": threshold(),
        "worker_thread_environment": {
            key: worker_env.get(key) for key in THREAD_VARIABLES
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cupy": cp.__version__,
            "device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        },
        "scope": (
            "Engineering replay of an audited constant-profile case; "
            "no benchmark verdict"
            if engineering
            else "Existence gate on one benchmark wave and window; no accuracy, "
            "generality or superiority claim"
        ),
    }
    write_json(destination / "manifest.json", manifest)

    def invoke(archive, truth, output):
        with output.with_suffix(".log").open("xb") as log:
            child = subprocess.Popen(
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
                env=worker_env,
                stdout=log,
                stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW,
            )
            status(
                state="RUNNING", phase="fit", output=output.name, child_pid=child.pid
            )
            code = child.wait()
        if code:
            return {
                "failure": f"Worker exit {code}; see {output.with_suffix('.log').name}"
            }
        return json.loads(output.read_text())

    rows = []
    try:
        for index, (kind, offset) in enumerate(cases):
            seed = SEED + offset
            generation, truth = generate(
                destination, index, kind, amplitudes[kind], seed, levels, status
            )
            for n in levels:
                row = {"case": index, "kind": kind, "seed": seed, "measurement_n": n}
                archive = destination / f"case_{index}_n{n}.npz"
                if not archive.exists():
                    row["failure"] = generation.get(
                        "failure", "Measurement prefix not completed"
                    )
                else:
                    result = invoke(archive, truth, archive.with_suffix(".json"))
                    row.update({k: v for k, v in result.items() if k != "full"})
                    if "full" in result:
                        row["detection"] = result["full"]["detection"]
                        row["evaluation"] = result["full"]["evaluation"]
                    row.update(
                        archive_sha256=sha(archive),
                        archive_bytes=archive.stat().st_size,
                    )
                write_json(destination / f"case_{index}_n{n}_result.json", row)
                rows.append(row)
                print(
                    json.dumps(
                        {
                            "case": index,
                            "kind": kind,
                            "n": n,
                            "q": row.get("detection", {}).get("q"),
                            "failure": row.get("failure"),
                            "seconds": row.get("wall_seconds"),
                            "peak_gib": row.get("peak_working_set_bytes", 0) / 2**30,
                        }
                    ),
                    flush=True,
                )
        summary = {
            **manifest,
            "results": rows,
            "elapsed_seconds": time.time() - started,
        }
        if not engineering:
            summary["level_verdicts"] = {
                str(n): verdict([r for r in rows if r["measurement_n"] == n])
                for n in levels
            }
        for name, digest in source_hashes.items():
            if lf_sha(ROOT / name) != digest:
                raise RuntimeError(f"A bound source changed during execution: {name}")
        write_json(destination / "summary.json", summary)
        status(
            state="COMPLETED",
            assignments=len(rows),
            level_verdicts=summary.get("level_verdicts"),
        )
    except BaseException as error:
        status(state="STOPPED", error=f"{type(error).__name__}: {error}")
        raise


def shutil_free(path):
    import shutil

    return shutil.disk_usage(path).free


def lf_sha(path):
    """Source digest under the repository's LF storage-boundary rule."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def compare(first, second):
    """Engineering replays must agree: exact order bits, 1e-8 on every fit field."""
    import numpy as np
    from replay_intrinsic_weyl_threads import differences

    report = {"first": str(first), "second": str(second), "levels": {}}
    for n in ENGINEERING_LEVELS:
        name = f"case_0_n{n}"
        with (
            np.load(first / f"{name}.npz", allow_pickle=False) as a,
            np.load(second / f"{name}.npz", allow_pickle=False) as b,
        ):
            bits = set(a.files) == set(b.files) and all(
                np.array_equal(a[key], b[key]) for key in a.files
            )
        old = json.loads((first / f"{name}.json").read_text())["full"]
        new = json.loads((second / f"{name}.json").read_text())["full"]
        # The zip container stores timestamps; the order bits are compared above.
        for record in (old, new):
            record["fit"].pop("order_archive_sha256")
        mismatch = differences(new, old)
        report["levels"][str(n)] = {
            "order_archives_bit_exact": bool(bits),
            "fit_and_evaluation_differences": mismatch,
        }
    report["reproduced"] = all(
        level["order_archives_bit_exact"]
        and not level["fit_and_evaluation_differences"]
        for level in report["levels"].values()
    )
    write_json(second / "replay_comparison.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--engineering", action="store_true")
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--compare", type=Path, nargs=2)
    args = parser.parse_args()
    if args.compare:
        outcome = compare(*args.compare)
        print(json.dumps({"reproduced": outcome["reproduced"]}))
        sys.exit(0 if outcome["reproduced"] else 1)
    elif args.check:
        print(
            json.dumps(
                {
                    "verdict": selfcheck(),
                    "benchmark_amplitude": float(benchmark_amplitude()),
                }
            )
        )
    elif args.worker:
        if args.truth is None or args.output is None:
            parser.error("truth and output are required for a worker")
        worker(args.worker, args.truth, args.output)
    else:
        if args.destination is None:
            parser.error("destination is required")
        run(args.destination, engineering=args.engineering)
