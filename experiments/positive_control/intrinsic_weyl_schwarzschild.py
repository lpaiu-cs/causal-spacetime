"""Schwarzschild detection test for the unchanged order-only Weyl inverse.

The geometry is a Schwarzschild exterior; the interval is the Alexandrov interval
of two events on the radial geodesic momentarily at rest at areal radius 6M, a
proper time 3M apart. Its tidal field is static and of Petrov type D, with the
same dimensionless strength as the plane-wave benchmark. Only the generator is
new; landmark selection, clocks, covariance, the fifth-degree inverse with
curvature-square terms and the confidence construction are the frozen ones. See
docs/paper/paper_a/intrinsic_weyl/schwarzschild_protocol.md.
"""

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

from intrinsic_weyl_bombelli_benchmark import lf_sha, pattern, shutil_free, threshold
from replay_intrinsic_weyl_threads import (
    ROOT,
    THREAD_VARIABLES,
    differences,
    sha,
    write_json,
)

PROTOCOL = "docs/paper/paper_a/intrinsic_weyl/schwarzschild_protocol.md"
MASS, APEX, PROPER_TIME = 1.0, 6.0, 3.0
SELECTION_N, CHUNK, BATCH = 4096, 4096, 65536
LEVELS = (1048576, 4194304, 16777216, 67108864)
SEED = 2026111200
# kind, seed offset. All four are co-required at every level.
CASES = (("curved", 0), ("curved", 1), ("flat", 2), ("flat", 3))
# Engineering population: a stronger field on a larger interval, disjoint seeds.
ENGINEERING_LEVELS = (1048576, 4194304, 16777216)
ENGINEERING_CASES = (("curved", 90),)
ENGINEERING_PROPER_TIME = 4.0
ORACLES = (
    "schwarzschild_diamond",
    "profile_plane_wave",
    "rotating_plane_wave",
    "rotating_plane_wave_gpu",
    "window_response",
    "p14_plane_wave",
    "intrinsic_weyl_rotating_pilot",
    "intrinsic_weyl_joint_audit",
    "intrinsic_weyl_joint_pilot",
)


def build(kind, proper_time):
    from schwarzschild_diamond import SchwarzschildDiamond

    if kind not in ("curved", "flat"):
        raise ValueError("Unknown test geometry")
    return SchwarzschildDiamond(MASS if kind == "curved" else 0.0, APEX, proper_time)


def verdict(rows):
    """Frozen existence gate for one measurement level; four rows required."""
    if len(rows) != len(CASES) or [r["kind"] for r in rows] != [k for k, _ in CASES]:
        raise ValueError("All four assigned cases are required, in order")
    if any("failure" in r or "detection" not in r for r in rows):
        return "INCOMPLETE"
    limit = threshold()
    curved = [r["detection"]["q"] > limit for r in rows if r["kind"] == "curved"]
    flat = [r["detection"]["q"] > limit for r in rows if r["kind"] == "flat"]
    if any(flat):
        return "CONTROL_FAILED"
    if all(curved):
        return "DETECTED"
    return "MIXED" if any(curved) else "NOT_DETECTED"


def selfcheck():
    def row(kind, q=None, failed=False):
        value = {"kind": kind}
        if failed:
            value["failure"] = "kept"
        elif q is not None:
            value["detection"] = {"q": q}
        return value

    limit = threshold()
    high, low = limit + 1, limit - 1

    def outcome(*values):
        kinds = ("curved", "curved", "flat", "flat")
        return verdict(
            [
                row(kind, failed=True) if value is None else row(kind, value)
                for kind, value in zip(kinds, values, strict=True)
            ]
        )

    assert outcome(high, high, low, low) == "DETECTED"
    assert outcome(high, low, low, low) == "MIXED"
    assert outcome(limit, high, low, low) == "MIXED"
    assert outcome(low, low, low, low) == "NOT_DETECTED"
    assert outcome(high, high, high, low) == "CONTROL_FAILED"
    assert outcome(high, None, low, low) == "INCOMPLETE"
    try:
        verdict([row("curved", high)] * 3)
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


def relation_selfcheck():
    """GPU relation bits against the direct solver, for both geometry classes."""
    import cupy as cp
    import numpy as np

    report = {}
    for kind in ("curved", "flat"):
        geometry = build(kind, PROPER_TIME)
        rng = np.random.default_rng(20261116)
        landmarks = np.vstack(
            (geometry.tips[0], geometry.sample(62, rng), geometry.tips[1])
        )
        events = geometry.sample(4096, rng)
        packed = geometry.packed_relations(landmarks, events, cp)
        direct = geometry.relations(landmarks[:, None], events[None])
        assert all(
            np.array_equal(np.unpackbits(p, axis=1).astype(bool), d)
            for p, d in zip(packed, direct, strict=True)
        )
        report[kind] = {
            "pairs": int(direct[0].size),
            "related": int(direct[0].sum() + direct[1].sum()),
            "resolved_by_the_direct_solver": geometry.resolved_margins,
        }
    cp.get_default_memory_pool().free_all_blocks()
    return report


def truth_coefficients(geometry, points):
    """True tidal tensor in the internal gauge of the landmarks' normal coordinates."""
    import numpy as np
    from intrinsic_weyl import WEYL

    # The inverse fixes its frame by the first two landmarks after the past tip.
    normal = geometry.normal_coordinates(points[1:3])
    first = normal[0, 1:] / np.linalg.norm(normal[0, 1:])
    second = normal[1, 1:] - (normal[1, 1:] @ first) * first
    second /= np.linalg.norm(second)
    frame = np.column_stack((first, second, np.cross(first, second)))
    electric = frame.T @ geometry.tidal() @ frame
    basis = WEYL[:5, 0, 1:, 0, 1:]
    gram = np.einsum("kij,lij->kl", basis, basis)
    return np.r_[
        np.linalg.solve(gram, np.einsum("kij,ij->k", basis, electric)), np.zeros(5)
    ]


def evaluate_fit(fit, truth_path, n):
    """Oracle evaluation after the pure fit has been saved."""
    import numpy as np
    from intrinsic_weyl import WEYL
    from intrinsic_weyl_joint import (
        TIDAL_ROOT,
        axis_signs,
        confidence_ball,
        tidal_distance,
    )
    from scipy.stats import chi2

    coefficients = np.asarray(fit["coefficients"])
    covariance = np.asarray(fit["nominal_component_covariance"])
    q = float(coefficients @ np.linalg.solve(covariance, coefficients))
    inference = confidence_ball(fit)
    estimate = np.asarray(inference["coefficients_xy_positive"])
    basis = WEYL[:5, 0, 1:, 0, 1:]
    detection = {
        "q": q,
        "threshold": threshold(),
        "nominal_p_value": float(chi2.sf(q, 10)),
        "tidal_norm_estimate": inference["tidal_norm"],
        "largest_tidal_variance": inference["largest_tidal_variance"],
        "variance_times_N": n * inference["largest_tidal_variance"],
        "electric_eigenvalues": np.linalg.eigvalsh(
            np.einsum("k,kij->ij", estimate[:5], basis)
        ).tolist(),
        "magnetic_norm": float(
            np.linalg.norm(np.einsum("k,kij->ij", estimate[5:], basis))
        ),
    }
    with np.load(truth_path, allow_pickle=False) as saved:
        points, kind = saved["points"], str(saved["kind"].item())
        proper_time = float(saved["proper_time"])
    geometry = build(kind, proper_time)
    expected = truth_coefficients(geometry, points)
    norm = float(np.linalg.norm(geometry.tidal()))
    error = tidal_distance(estimate, expected)
    evaluation = {
        "truth_coefficients_xy_positive": expected.tolist(),
        "truth_norm": norm,
        "truth_electric_eigenvalues": np.linalg.eigvalsh(geometry.tidal()).tolist(),
        "full_tidal_error_mod_reflection": error,
        "relative_full_tidal_error": error / norm if norm else None,
        "confidence_radius_95": inference["radius"],
        "nominal_ball_contains_truth": bool(error <= inference["radius"]),
        "continuum_bias": None,
        "continuum_bias_note": (
            "Exact interval volumes are not available for this geometry"
        ),
    }
    if norm:
        # The residual z reflection is unobservable. The expected figure keeps the
        # less favorable candidate; the realized one, the candidate nearer to the
        # estimate.
        candidates = (expected, axis_signs([1, 1, -1]) * expected)
        options = [pattern(inference, candidate) for candidate in candidates]
        snr = min(v["expected_snr"] for v in options)
        nearer = min(
            zip(candidates, options, strict=True),
            key=lambda item: np.linalg.norm(TIDAL_ROOT @ (estimate - item[0])),
        )[1]
        evaluation["truth_pattern"] = {
            "expected_snr": snr,
            "events_for_expected_snr_one": n / snr**2,
            "expected_snr_matched_reflection": nearer["expected_snr"],
            "realized_z_matched_reflection": nearer["realized_z"],
        }
    return {
        "fit": fit,
        "inference": inference,
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
            fit = joint.fit_order_archive(
                archive, model_order=5, square_correction=True
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
        record["full"] = evaluate_fit(fit, truth, fit["measurement_n"])
    except (ValueError, ArithmeticError, MemoryError, np.linalg.LinAlgError) as error:
        record["failure"] = f"{type(error).__name__}: {error}"
    record.update(
        wall_seconds=time.perf_counter() - wall,
        process_cpu_seconds=time.process_time() - cpu,
        peak_working_set_bytes=peak_working_set(),
    )
    write_json(output, record)


def generate(folder, index, kind, proper_time, seed, levels, status):
    """Selection order, landmarks, and nested measurement prefixes for one case."""
    import cupy as cp
    import numpy as np
    from intrinsic_weyl_high_count import peak_working_set
    from intrinsic_weyl_joint import selection

    geometry = build(kind, proper_time)
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
            (geometry.tips[0], geometry.sample(SELECTION_N, rng), geometry.tips[1])
        )
        generation["selection_events_drawn"] = SELECTION_N
        # Labels are randomized independently of coordinates; only the two
        # supplied tips have distinguished IDs.
        permutation = np.r_[
            0, rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
        ]
        points = points[permutation]
        order = geometry.order(points)
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
            proper_time=np.array(float(proper_time)),
        )
        generation["landmarks_including_tips"] = len(points)
        top = levels[-1]
        after = np.empty((len(points), top // 8), dtype=np.uint8)
        before = np.empty_like(after)
        rng = np.random.default_rng([seed, 1])
        for start in range(0, top, BATCH):
            end = min(start + BATCH, top)
            events = np.vstack(
                [geometry.sample(CHUNK, rng) for _ in range((end - start) // CHUNK)]
            )
            generation["sampled_events"] += len(events)
            future, past = geometry.packed_relations(points, events, cp)
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
        margins_resolved_by_the_direct_solver=geometry.resolved_margins,
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
    from schwarzschild_diamond import selfcheck as geometry_selfcheck

    if os.name != "nt":
        raise RuntimeError("The resource accounting of this runner uses Windows")
    if destination.exists():
        raise FileExistsError("Preserve every existing run")
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip()
    if dirty and not engineering:
        raise RuntimeError("Commit the protocol and runner before test data")
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    here = Path(__file__).resolve().parent
    sources = [
        here / name
        for name in (
            Path(__file__).name,
            "schwarzschild_diamond.py",
            "intrinsic_weyl_bombelli_benchmark.py",
            "replay_intrinsic_weyl_threads.py",
            "intrinsic_weyl_joint.py",
            "intrinsic_weyl.py",
            "intrinsic_weyl_square.py",
            "intrinsic_weyl_gpu_counts.py",
            "intrinsic_weyl_high_count.py",
        )
    ]
    if not engineering:
        sources.append(ROOT / PROTOCOL)
    source_hashes = {
        str(p.relative_to(ROOT)).replace("\\", "/"): lf_sha(p) for p in sources
    }
    checked = {
        "generator": geometry_selfcheck(),
        "gpu_relation_bits": relation_selfcheck(),
        "jet_inverse": selfcheck_jet(),
        "gpu_second_moment": gpu_selfcheck(),
        "verdict": selfcheck(),
    }
    if engineering:
        cases, levels, proper_time = (
            ENGINEERING_CASES,
            ENGINEERING_LEVELS,
            ENGINEERING_PROPER_TIME,
        )
    else:
        cases, levels, proper_time = CASES, LEVELS, PROPER_TIME
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
    reference = build("curved", proper_time)
    manifest = {
        "engineering_only": bool(engineering),
        "code_version": stamp,
        "worktree_clean_at_start": not dirty,
        "source_sha256_lf": source_hashes,
        "selfcheck": checked,
        "geometry": {
            "mass": MASS,
            "apex_areal_radius": APEX,
            "tip_proper_time": proper_time,
            "apex_isotropic_radius": reference.apex,
            "tip_coordinate_half_time": reference.half_time,
            "tip_isotropic_radius": float(reference.base[0]),
            "true_tidal_eigenvalues": np.linalg.eigvalsh(reference.tidal()).tolist(),
            "true_full_tidal_norm": float(np.linalg.norm(reference.tidal())),
        },
        "cases": [
            {"case": i, "kind": k, "seed": SEED + offset}
            for i, (k, offset) in enumerate(cases)
        ],
        "measurement_levels": levels,
        "selection_n": SELECTION_N,
        "sampling_chunk": CHUNK,
        "classification_batch": BATCH,
        "false_positive_level_per_test": 1e-3,
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
            "Engineering replay on a stronger field; no test verdict"
            if engineering
            else "Existence gate on one Schwarzschild interval; no accuracy, "
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
                destination, index, kind, proper_time, seed, levels, status
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
                            "relative_error": row.get("evaluation", {}).get(
                                "relative_full_tidal_error"
                            ),
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


def compare(first, second):
    """Engineering replays must agree: exact order bits, 1e-8 on every fit field."""
    import numpy as np

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
        report["levels"][str(n)] = {
            "order_archives_bit_exact": bool(bits),
            "fit_and_evaluation_differences": differences(new, old),
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
        print(json.dumps({"verdict": selfcheck()}))
    elif args.worker:
        if args.truth is None or args.output is None:
            parser.error("truth and output are required for a worker")
        worker(args.worker, args.truth, args.output)
    else:
        if args.destination is None:
            parser.error("destination is required")
        run(args.destination, engineering=args.engineering)
