"""Frozen independent precision test; unchanged inverse, powered event counts."""

import argparse
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

from replay_intrinsic_weyl_threads import (
    ROOT,
    THREAD_VARIABLES,
    differences,
    sha,
    write_json,
)

PROTOCOL = "docs/paper/paper_a/intrinsic_weyl/weak_curvature_accuracy_protocol.md"
ENGINEERING_SHA = "baa7a69e68f2becb98eb41e776618a491618375fdde64483f5386408112a97c9"
ALPHA, BETA, DELTA, ETA = 0.05, 0.1, 0.02, 1e-10
CHUNK, BATCH, SEED = 4096, 32768, 2026112000
# A, w, retained pilot case (None means fresh selection), frozen V cap.
CASES = (
    (0.0, 0.0, 0, 4000),
    (0.0, 0.5, 1, 4000),
    (0.25, 0.0, 2, 5500),
    (0.25, 0.5, 4, 5700),
    (0.25, -0.5, 6, 2900),
    (0.25, -0.5, None, 8000),
    (0.25, 0.0, None, 8000),
    (0.25, 0.5, None, 8000),
)
EXPECTED_N = (
    25108480,
    25108480,
    49311744,
    51101696,
    26001408,
    71720960,
    71720960,
    71720960,
)


def read(path):
    return json.loads(path.read_text())


def plan():
    from scipy.stats import chi2

    factor = (
        math.sqrt(chi2.ppf(1 - ALPHA / 8, 10)) + math.sqrt(chi2.ppf(1 - BETA / 8, 10))
    ) ** 2
    rows = []
    for i, (amplitude, rotation, retained, cap) in enumerate(CASES):
        margin = (
            ((1 - DELTA) * 0.25 - 2 * DELTA) * (2 * amplitude)
            if amplitude
            else (1 - DELTA) * 0.125
        )
        n = CHUNK * math.ceil(cap * factor / (margin - 2 * ETA) ** 2 / CHUNK)
        assert n == EXPECTED_N[i]
        rows.append(
            {
                "case": i,
                "amplitude": amplitude,
                "rotation_rate": rotation,
                "retained_pilot_case": retained,
                "variance_cap": cap,
                "measurement_n": n,
                "seed": SEED + i,
            }
        )
    return rows


def verdict(rows):
    if len(rows) != 8 or {r["case"] for r in rows} != set(range(8)):
        raise ValueError("All eight distinct assignments are required")
    return (
        "SUPPORTED"
        if all(
            "failure" not in r and r.get("precision", {}).get("passed", False)
            for r in rows
        )
        else "NOT_SUPPORTED"
    )


def precision(full, n, cap):
    from intrinsic_weyl_joint import confidence_ball

    ball = confidence_ball(full["fit"], coverage=1 - ALPHA / 8, relative_bias=DELTA)
    confidence_ball(full["continuum_fit"], coverage=1 - ALPHA / 8, relative_bias=DELTA)
    evaluation = full["evaluation"]
    norm = evaluation["truth_norm"]
    error = evaluation["full_tidal_error_mod_reflection"]
    radius = ball["radius"] + ETA / (1 - DELTA)
    threshold = 0.25 * norm if norm else 0.125
    price = n * ball["largest_tidal_variance"]
    conditions = {
        "bias": evaluation["continuum_bias_mod_reflection"] <= DELTA * norm + ETA,
        "quadrature": 0 <= evaluation["quadrature_relative_change"] < 1e-10,
        "contains_truth": error <= radius,
        "whole_ball_within_accuracy_target": error + radius <= threshold,
        "variance_cap": price <= cap,
    }
    return {
        "nominal_case_coverage": 1 - ALPHA / 8,
        "family_coverage": 1 - ALPHA,
        "target_joint_power": 1 - BETA,
        "radius": radius,
        "statistical_radius": ball["statistical_radius"],
        "error": error,
        "truth_norm": norm,
        "accuracy_threshold": threshold,
        "error_plus_radius": error + radius,
        "relative_error_plus_radius": (error + radius) / norm if norm else None,
        "variance_times_N": price,
        "variance_cap": cap,
        "conditions": conditions,
        "passed": all(conditions.values()),
    }


def selfcheck(reference=None):
    assert sum(r["measurement_n"] for r in plan()) == 391794688
    rows = [{"case": i, "precision": {"passed": True}} for i in range(8)]
    assert verdict(rows) == "SUPPORTED"
    rows[7]["failure"] = "kept even beside a favorable precision flag"
    assert verdict(rows) == "NOT_SUPPORTED"
    try:
        verdict(rows[:-1])
    except ValueError:
        pass
    else:
        raise AssertionError("A missing assignment was accepted")
    checks = {
        "all_power_counts_match": True,
        "all_eight_required": True,
        "failure_cannot_disappear": True,
    }
    if reference is not None:
        old = read(reference / "case_2_n4194304.json")["full"]
        assert old["evaluation"]["relative_full_tidal_error"] < 0.25
        tested = precision(old, 4194304, 5500)
        assert not tested["conditions"]["whole_ball_within_accuracy_target"]
        assert not tested["passed"]
        checks["old_good_point_with_wide_interval_rejected"] = True
    return checks


def worker(archive, truth, output):
    import intrinsic_weyl_joint as joint
    import numpy as np
    from intrinsic_weyl_gpu_counts import clock_second_moment
    from intrinsic_weyl_high_count import worker as unchanged_worker

    joint._clock_second_moment = clock_second_moment
    original_fit = joint.fit_order_archive
    forbidden = (
        "rotating_plane_wave",
        "rotating_plane_wave_gpu",
        "window_response",
        "p14_plane_wave",
        "intrinsic_weyl_rotating_pilot",
        "intrinsic_weyl_joint_audit",
        "intrinsic_weyl_joint_pilot",
    )

    def blind_fit(path, **kwargs):
        original_load, opened = np.load, []

        def guard(selected, *args, **options):
            if Path(selected).resolve() != archive.resolve():
                raise ValueError("Only the assigned order archive may be opened")
            opened.append(str(Path(selected).resolve()))
            return original_load(selected, *args, **options)

        assert not any(name in sys.modules for name in forbidden)
        np.load = guard
        try:
            result = original_fit(path, **kwargs)
            assert not any(name in sys.modules for name in forbidden)
            write_json(
                output.with_name(output.stem + "_blind_boundary.json"),
                {
                    "opened_archives": opened,
                    "oracle_modules_loaded": [],
                    "archive_sha256": sha(archive),
                    "gpu_counts_module": str(
                        Path(__file__).with_name("intrinsic_weyl_gpu_counts.py")
                    ),
                },
            )
            return result
        finally:
            np.load = original_load

    joint.fit_order_archive = blind_fit
    unchanged_worker(archive, truth, output)


def generate(folder, pilot, case, status):
    import numpy as np
    from intrinsic_weyl_high_count import peak_working_set
    from intrinsic_weyl_joint import selection
    from rotating_plane_wave import RotatingWave
    from rotating_plane_wave_gpu import RotatingWaveGPU

    i, n = case["case"], case["measurement_n"]
    amplitude, rotation = case["amplitude"], case["rotation_rate"]
    wave = RotatingWave(amplitude, rotation)
    wall, cpu = time.perf_counter(), time.process_time()
    generation = {
        "selection_events_drawn": 0,
        "sampled_events": 0,
        "classified_events": 0,
        "measurement_target": n,
    }
    archive = folder / f"case_{i}.npz"
    truth_path = folder / f"truth_{i}.npz"
    try:
        retained = case["retained_pilot_case"]
        if retained is not None:
            with np.load(
                pilot / f"case_{retained}_nested.npz", allow_pickle=False
            ) as old:
                selected = {
                    k: old[k]
                    for k in ("selection_order", "selection_design", "selected_ids")
                }
            with np.load(
                pilot / f"truth_{retained}_nested.npz", allow_pickle=False
            ) as truth:
                points = truth["points"]
                assert (float(truth["amplitude"]), float(truth["rotation_rate"])) == (
                    amplitude,
                    rotation,
                )
            shutil.copyfile(pilot / f"truth_{retained}_nested.npz", truth_path)
        else:
            status(state="RUNNING", phase="selection", case=i)
            rng = np.random.default_rng([case["seed"], 0])
            points = np.vstack(
                ([0, 0.5, 0, 0], wave.sample(4096, rng), [1, -0.5, 0, 0])
            )
            generation["selection_events_drawn"] = 4096
            permutation = np.r_[
                0, rng.permutation(np.arange(1, len(points) - 1)), len(points) - 1
            ]
            points = points[permutation]
            order = wave.order(points)
            try:
                ids, _ = selection(order, design="nested")
            except ValueError:
                np.savez(
                    folder / f"case_{i}_selection_failed.npz",
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
                truth_path, points=points, amplitude=amplitude, rotation_rate=rotation
            )
        expected, _ = selection(selected["selection_order"], design="nested")
        assert np.array_equal(selected["selected_ids"], expected)
        generation["landmarks_including_tips"] = len(points)
        gpu = RotatingWaveGPU(amplitude, rotation, wave.coefficients)
        after = np.empty((len(points), n // 8), dtype=np.uint8)
        before = np.empty_like(after)
        rng = np.random.default_rng([case["seed"], 1])
        try:
            for start in range(0, n, BATCH):
                end = min(start + BATCH, n)
                events = []
                for _ in range((end - start) // CHUNK):
                    events.append(wave.sample(CHUNK, rng))
                    generation["sampled_events"] += CHUNK
                future, past = gpu.packed_relations(points, np.vstack(events))
                after[:, start // 8 : end // 8] = future
                before[:, start // 8 : end // 8] = past
                generation["classified_events"] = end
                if end % (32 * BATCH) == 0 or end == n:
                    status(
                        state="RUNNING", phase="generate", case=i, events=end, target=n
                    )
        except (ValueError, ArithmeticError, MemoryError):
            end = generation["classified_events"]
            if end:
                with (folder / f"case_{i}_partial.npz").open("xb") as partial:
                    np.savez(
                        partial,
                        **selected,
                        after=after[:, : end // 8],
                        before=before[:, : end // 8],
                        measurement_n=end,
                    )
            raise
        status(state="RUNNING", phase="archive", case=i, events=n)
        with archive.open("xb") as stream:
            np.savez(stream, **selected, after=after, before=before, measurement_n=n)
    except (ValueError, ArithmeticError, MemoryError) as error:
        generation["failure"] = f"{type(error).__name__}: {error}"
    generation.update(
        wall_seconds=time.perf_counter() - wall,
        process_cpu_seconds=time.process_time() - cpu,
        parent_lifetime_peak_working_set_bytes=peak_working_set(),
    )
    if "landmarks_including_tips" in generation:
        generation["directional_predicates_classified"] = (
            2 * generation["landmarks_including_tips"] * generation["classified_events"]
        )
    write_json(folder / f"case_{i}_generation.json", generation)
    return generation, archive, truth_path


def run(destination, pilot, reference, engineering):
    import cupy as cp
    import numpy as np
    import scipy

    if destination.exists():
        raise FileExistsError("Preserve every existing accuracy run")
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT):
        raise RuntimeError("Commit the protocol and runner before data")
    if shutil.disk_usage(destination.parent).free < 30 * 2**30:
        raise RuntimeError("At least 30 GiB free storage is required")
    assert sha(engineering) == ENGINEERING_SHA
    engineering_record = read(engineering)
    assert engineering_record["verdict"] == "EQUIVALENT_ON_FROZEN_ENGINEERING_CHECKS"
    assert len(engineering_record["cases"]) == 8
    assert all(c["packed_relations_bit_exact"] for c in engineering_record["cases"])
    assert engineering_record["pure_estimator_replay"]["all_clock_arrays_bit_exact"]
    assert engineering_record["pure_estimator_replay"]["all_fit_fields_exact"]
    old_manifest = read(reference / "manifest.json")
    assert {key: os.environ.get(key) for key in THREAD_VARIABLES} == old_manifest[
        "thread_environment"
    ]
    assert os.environ.get("CUDA_VISIBLE_DEVICES") == "0"
    assert old_manifest["runtime"] == {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    assert cp.__version__ == engineering_record["runtime"]["cupy"] == "14.2.0"
    assert (
        cp.cuda.runtime.runtimeGetVersion()
        == engineering_record["runtime"]["cuda_runtime"]
    )
    assert cp.cuda.runtime.driverGetVersion() == engineering_record["runtime"]["driver"]
    assert (
        cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
        == engineering_record["runtime"]["device"]
    )
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    for ancestor in (old_manifest["code_version"], engineering_record["code_version"]):
        subprocess.run(
            ["git", "merge-base", "--is-ancestor", ancestor, stamp],
            cwd=ROOT,
            check=True,
        )
    sources = dict(old_manifest["source_sha256"])
    sources.update(engineering_record["source_sha256"])
    sources.update(
        {
            str(Path(__file__).relative_to(ROOT)).replace("\\", "/"): sha(
                Path(__file__)
            ),
            PROTOCOL: sha(ROOT / PROTOCOL),
            "pyproject.toml": sha(ROOT / "pyproject.toml"),
        }
    )
    planning = {
        str(
            reference / "summary.json"
        ): "e85c39648777b23ab638dbd32b435aca9089392b97d6b75e75d837d96c48486c",
        str(
            reference.parent
            / "intrinsic_weyl_nested_coverage_v1_continued/summary.json"
        ): "238ac2400deda3d6553e7da628bc0c27fdb277aa6bb221d924a2d21d6ed86d27",
        str(engineering): ENGINEERING_SHA,
    }

    def verify_inputs():
        for name, digest in sources.items():
            assert sha(ROOT / name) == digest, name
        for name, digest in old_manifest["pilot_input_sha256"].items():
            assert sha(pilot / name) == digest, name
        for name, digest in planning.items():
            assert sha(Path(name)) == digest, name

    verify_inputs()
    assignments = plan()
    checked = selfcheck(reference)
    destination.mkdir()
    start = time.time()

    def status(**fields):
        temporary = destination / ".status.tmp"
        temporary.write_text(
            json.dumps({"pid": os.getpid(), "started_unix": start, **fields}, indent=2)
            + "\n"
        )
        os.replace(temporary, destination / "status.json")

    manifest = {
        "code_version": stamp,
        "source_sha256": sources,
        "pilot_input_sha256": old_manifest["pilot_input_sha256"],
        "planning_evidence_sha256": planning,
        "engineering_receipt_sha256": ENGINEERING_SHA,
        "cases": assignments,
        "maximum_measurement_events": sum(c["measurement_n"] for c in assignments),
        "maximum_new_selection_events": 3 * 4096,
        "supplied_new_tips": 6,
        "sampling_chunk": CHUNK,
        "classification_batch": BATCH,
        "alpha_family": ALPHA,
        "beta_joint_power": BETA,
        "relative_bias_allowance": DELTA,
        "absolute_bias_slack": ETA,
        "weak_relative_accuracy_target": 0.25,
        "flat_absolute_accuracy_target": 0.125,
        "thread_environment": old_manifest["thread_environment"],
        "runtime": {**old_manifest["runtime"], **engineering_record["runtime"]},
        "gpu_used": True,
        "selfcheck": checked,
        "scope": (
            "Prospective simultaneous interval accuracy on six weak cases, "
            "including three independent selections, with two flat controls. "
            "Gaussian/variance/bias planning assumptions remain explicit; "
            "no universal weak-field claim."
        ),
    }
    write_json(destination / "manifest.json", manifest)
    status(state="RUNNING", phase="starting", total_cases=8)

    def invoke(archive, truth, output, case):
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
                stdout=log,
                stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW,
            )
            status(state="RUNNING", phase="fit", case=case, child_pid=child.pid)
            exit_code = child.wait()
        if exit_code:
            return {
                "failure": (
                    f"Worker exit {exit_code}; "
                    f"see {output.with_suffix('.log').name}"
                )
            }
        return read(output)

    rows = []
    try:
        preflight = invoke(
            reference / "case_4_n262144.npz",
            pilot / "truth_4_nested.npz",
            destination / "preflight_worker.json",
            "preflight",
        )
        original = read(reference / "case_4_n262144.json")["full"]
        mismatch = differences(preflight.get("full"), original)
        write_json(
            destination / "preflight_comparison.json",
            {
                "case": 4,
                "measurement_n": 262144,
                "differences": mismatch,
                "within_frozen_tolerance": not mismatch,
                "all_full_fields_exact": preflight.get("full") == original,
                "scientific_fresh_events": 0,
            },
        )
        if mismatch:
            raise RuntimeError("Integrated worker replay failed before new data")
        for case in assignments:
            i = case["case"]
            status(state="RUNNING", phase="prepare", case=i)
            generation, archive, truth = generate(destination, pilot, case, status)
            row = {**case, "generation": generation}
            if "failure" in generation:
                row["failure"] = generation["failure"]
            else:
                output = destination / f"case_{i}_full.json"
                result = invoke(archive, truth, output, i)
                row.update({k: v for k, v in result.items() if k != "full"})
                if "full" in result:
                    full = result["full"]
                    assert full["evaluation"]["truth_norm"] == 2 * case["amplitude"]
                    row.update(full["evaluation"])
                    try:
                        row["precision"] = precision(
                            full, case["measurement_n"], case["variance_cap"]
                        )
                    except (ValueError, ArithmeticError) as error:
                        row["failure"] = (
                            f"Precision assessment: {type(error).__name__}: {error}"
                        )
                row["archive_bytes"] = archive.stat().st_size
            paths = sorted(destination.glob(f"case_{i}_*")) + [archive, truth]
            row["artifact_sha256"] = {
                path.name: sha(path) for path in paths if path.is_file()
            }
            write_json(destination / f"case_{i}_result.json", row)
            rows.append(row)
            print(
                json.dumps(
                    {
                        "case": i,
                        "precision": row.get("precision"),
                        "failure": row.get("failure"),
                        "completed_assignments": len(rows),
                    }
                ),
                flush=True,
            )
        verify_inputs()
        for row in rows:
            for name, digest in row["artifact_sha256"].items():
                assert sha(destination / name) == digest, name
        summary = {
            **manifest,
            "results": rows,
            "completed_assignments": len(rows),
            "precision_gate": verdict(rows),
            "sampled_measurement_events": sum(
                r["generation"]["sampled_events"] for r in rows
            ),
            "classified_measurement_events": sum(
                r["generation"]["classified_events"] for r in rows
            ),
            "new_selection_events": sum(
                r["generation"]["selection_events_drawn"] for r in rows
            ),
            "elapsed_seconds": time.time() - start,
            "audit_status": "PENDING_INDEPENDENT_POST_COMPLETION_AUDIT",
        }
        write_json(destination / "summary.json", summary)
        status(
            state="COMPLETED",
            assignments=len(rows),
            precision_gate=summary["precision_gate"],
        )
    except BaseException as error:
        status(state="STOPPED", error=f"{type(error).__name__}: {error}")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--pilot", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--engineering", type=Path)
    args = parser.parse_args()
    if args.check:
        print(
            json.dumps({"selfcheck": selfcheck(args.reference), "assignments": plan()})
        )
    elif args.worker:
        if args.truth is None or args.output is None:
            parser.error("truth and output are required for a worker")
        worker(args.worker, args.truth, args.output)
    else:
        if any(
            v is None
            for v in (args.destination, args.pilot, args.reference, args.engineering)
        ):
            parser.error("destination, pilot, reference and engineering are required")
        run(args.destination, args.pilot, args.reference, args.engineering)
