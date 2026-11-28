"""Fixed independent measurement ladder; no accuracy or coverage verdict."""

import argparse
import ctypes
import json
import os
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from replay_intrinsic_weyl_threads import (
    ROOT,
    THREAD_VARIABLES,
    differences,
    sha,
    write_json,
)

LEVELS = (262144, 1048576, 4194304)
CASES = (
    (0, 0),
    (0, 0.5),
    (0.25, 0),
    (0.75, 0),
    (0.25, 0.5),
    (0.75, 0.5),
    (0.25, -0.5),
    (0.75, -0.5),
)
SEED, CHUNK = 2026103000, 4096


def peak_working_set():
    """Windows process peak, including native NumPy allocations."""
    from ctypes import wintypes

    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("faults", wintypes.DWORD)] + [
            (name, ctypes.c_size_t)
            for name in (
                "peak",
                "working",
                "peak_paged",
                "paged",
                "peak_nonpaged",
                "nonpaged",
                "pagefile",
                "peak_pagefile",
            )
        ]

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    api = ctypes.WinDLL("psapi", use_last_error=True).GetProcessMemoryInfo
    api.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
    api.restype = wintypes.BOOL
    counter = Counters()
    counter.cb = ctypes.sizeof(counter)
    if not api(kernel.GetCurrentProcess(), ctypes.byref(counter), counter.cb):
        raise ctypes.WinError(ctypes.get_last_error())
    return counter.peak


def worker(archive, truth, output):
    # No generator, truth evaluator or profile is imported before this returns.
    import intrinsic_weyl_joint as joint

    wall, cpu = time.perf_counter(), time.process_time()
    record = {
        "archive_sha256": sha(archive),
        "thread_environment": {key: os.environ.get(key) for key in THREAD_VARIABLES},
    }
    try:
        fit, clocks = joint.fit_order_archive(
            archive, with_clocks=True, model_order=5, square_correction=True
        )
        record["finite_fit_wall_seconds"] = time.perf_counter() - wall
        record["finite_fit_cpu_seconds"] = time.process_time() - cpu
        forbidden = (
            "rotating_plane_wave",
            "window_response",
            "p14_plane_wave",
            "intrinsic_weyl_joint_audit",
            "intrinsic_weyl_joint_pilot",
            "intrinsic_weyl_rotating_pilot",
        )
        imported = [name for name in forbidden if name in sys.modules]
        if imported:
            raise RuntimeError(f"Oracle imported before pure-fit boundary: {imported}")
        write_json(
            output.with_name(output.stem + "_pure.json"),
            {
                "fit": fit,
                "oracle_modules_before_return": imported,
                "inverse_module": joint.__file__,
                "archive_sha256": sha(archive),
            },
        )
        from intrinsic_weyl_rotating_pilot import evaluate_fit

        record["full"] = evaluate_fit(
            fit, clocks, truth, model_order=5, square_correction=True
        )
    except (ValueError, ArithmeticError, MemoryError) as error:
        record["failure"] = f"{type(error).__name__}: {error}"
    record.update(
        wall_seconds=time.perf_counter() - wall,
        process_cpu_seconds=time.process_time() - cpu,
        peak_working_set_bytes=peak_working_set(),
    )
    write_json(output, record)


def write_prefix(path, selection_data, after, before, n):
    import numpy as np

    if n <= 0 or n % 8 or after.shape != before.shape or n > 8 * after.shape[1]:
        raise ValueError("Invalid packed measurement prefix")
    with path.open("xb") as output:
        np.savez_compressed(
            output,
            **selection_data,
            after=after[:, : n // 8],
            before=before[:, : n // 8],
            measurement_n=n,
        )


def selfcheck():
    import numpy as np

    selection_data = {
        "selection_order": np.triu(np.ones((3, 3), bool), 1),
        "selection_design": np.array("nested"),
        "selected_ids": np.arange(3),
    }
    after = np.arange(12, dtype=np.uint8).reshape(3, 4)
    before = 255 - after
    with tempfile.TemporaryDirectory() as folder:
        for n in (8, 24, 32):
            path = Path(folder) / f"{n}.npz"
            write_prefix(path, selection_data, after, before, n)
            with np.load(path, allow_pickle=False) as saved:
                assert set(saved.files) == set(selection_data) | {
                    "after",
                    "before",
                    "measurement_n",
                }
                for key, value in selection_data.items():
                    assert np.array_equal(saved[key], value)
                assert int(saved["measurement_n"]) == n
                assert np.array_equal(saved["after"], after[:, : n // 8])
                assert np.array_equal(saved["before"], before[:, : n // 8])
        try:
            write_prefix(Path(folder) / "bad.npz", selection_data, after, before, 40)
        except ValueError:
            pass
        else:
            raise AssertionError("Oversized prefix was accepted")
    assert peak_working_set() > 0
    return {"prefix_fields_and_bits_exact": True, "oversized_prefix_rejected": True}


def run(pilot, destination):
    if os.name != "nt":
        raise RuntimeError("This fixed resource measurement uses Windows")
    expected_env = {
        key: "1" if key in THREAD_VARIABLES[:2] else None for key in THREAD_VARIABLES
    }
    if {key: os.environ.get(key) for key in THREAD_VARIABLES} != expected_env:
        raise RuntimeError("Set the two fixed thread limits before starting Python")
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit the protocol and implementation before new data")
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", stamp, "HEAD"], cwd=ROOT, check=True
    )
    sources = [
        Path(__file__).with_name(name)
        for name in (
            Path(__file__).name,
            "replay_intrinsic_weyl_threads.py",
            "intrinsic_weyl_rotating_pilot.py",
            "intrinsic_weyl_joint.py",
            "intrinsic_weyl.py",
            "intrinsic_weyl_square.py",
            "intrinsic_weyl_joint_audit.py",
            "intrinsic_weyl_joint_pilot.py",
            "rotating_plane_wave.py",
            "window_response.py",
            "p14_plane_wave.py",
        )
    ]
    sources.append(ROOT / "docs/paper/paper_a/intrinsic_weyl/high_count_diagnostic.md")
    inputs = (
        [pilot / "manifest.json"]
        + [
            pilot / f"{kind}_{i}_nested.npz"
            for i in range(8)
            for kind in ("case", "truth")
        ]
        + [pilot / f"case_{i}_nested_full.json" for i in range(8)]
    )
    source_hashes = {
        str(p.relative_to(ROOT)).replace("\\", "/"): sha(p) for p in sources
    }
    input_hashes = {p.name: sha(p) for p in inputs}
    old_manifest = json.loads((pilot / "manifest.json").read_text())
    if old_manifest["cases_fixed_before_data"] != [list(case) for case in CASES]:
        raise ValueError("Unexpected old pilot population")

    def verify_bound_files():
        if any(sha(ROOT / name) != value for name, value in source_hashes.items()):
            raise RuntimeError("A frozen source changed during execution")
        if any(sha(pilot / name) != value for name, value in input_hashes.items()):
            raise RuntimeError("An old pilot input changed during execution")

    checked = selfcheck()
    destination.mkdir(parents=True, exist_ok=False)
    started = time.time()

    def status(**fields):
        (destination / "status.json").write_text(
            json.dumps(
                {"pid": os.getpid(), "started_unix": started, **fields}, indent=2
            )
            + "\n"
        )

    def invoke(archive, truth, output, env):
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
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
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

    import numpy as np
    import scipy

    manifest = {
        "code_version": stamp,
        "source_sha256": source_hashes,
        "pilot_input_sha256": input_hashes,
        "selfcheck": checked,
        "cases": CASES,
        "measurement_levels": LEVELS,
        "measurement_seed_base": SEED,
        "chunk": CHUNK,
        "maximum_fresh_events": 8 * LEVELS[-1],
        "fresh_selection_events": 0,
        "thread_environment": expected_env,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "scope": "Fixed-shape descriptive diagnostic; no accuracy or coverage verdict",
    }
    write_json(destination / "manifest.json", manifest)
    default_env = os.environ.copy()
    for key in THREAD_VARIABLES:
        default_env.pop(key, None)
    replay = invoke(
        pilot / "case_5_nested.npz",
        pilot / "truth_5_nested.npz",
        destination / "preflight.json",
        default_env,
    )
    original = json.loads((pilot / "case_5_nested_full.json").read_text())
    mismatch = differences(replay.get("full"), original)
    write_json(
        destination / "preflight_comparison.json",
        {
            "differences": mismatch,
            "fresh_events": 0,
            "within_frozen_tolerance": not mismatch,
            "original_sha256": sha(pilot / "case_5_nested_full.json"),
        },
    )
    if mismatch:
        status(state="STOPPED_BEFORE_DATA", phase="preflight_failed")
        raise RuntimeError("Old-data evaluation replay failed; preserve the comparison")

    from rotating_plane_wave import RotatingWave

    rows, generated = [], []
    for case, (amplitude, rotation) in enumerate(CASES):
        verify_bound_files()
        with np.load(pilot / f"case_{case}_nested.npz", allow_pickle=False) as saved:
            selected = {
                key: saved[key]
                for key in ("selection_order", "selection_design", "selected_ids")
            }
        with np.load(pilot / f"truth_{case}_nested.npz", allow_pickle=False) as truth:
            if (float(truth["amplitude"]), float(truth["rotation_rate"])) != (
                amplitude,
                rotation,
            ):
                raise ValueError("Truth metadata differs from the fixed population")
            points = truth["points"]
        after = np.empty((len(points), LEVELS[-1] // 8), dtype=np.uint8)
        before = np.empty_like(after)
        wave, rng = (
            RotatingWave(amplitude, rotation),
            np.random.default_rng([SEED + case, 1]),
        )
        generation = {
            "case": case,
            "sampled_events": 0,
            "classified_events": 0,
            "prefix_costs": [],
        }
        wall, cpu = time.perf_counter(), time.process_time()
        status(state="RUNNING", phase="generate", case=case, events=0)
        try:
            for start in range(0, LEVELS[-1], CHUNK):
                events = wave.sample(CHUNK, rng)
                generation["sampled_events"] += CHUNK
                future, past = wave.relations(points[:, None], events[None])
                end = start + CHUNK
                after[:, start // 8 : end // 8] = np.packbits(future, axis=1)
                before[:, start // 8 : end // 8] = np.packbits(past, axis=1)
                generation["classified_events"] = end
                if end in LEVELS:
                    write_prefix(
                        destination / f"case_{case}_n{end}.npz",
                        selected,
                        after,
                        before,
                        end,
                    )
                    generation["prefix_costs"].append(
                        {
                            "measurement_n": end,
                            "wall_seconds": time.perf_counter() - wall,
                            "process_cpu_seconds": time.process_time() - cpu,
                        }
                    )
                if end % (64 * CHUNK) == 0:
                    status(state="RUNNING", phase="generate", case=case, events=end)
        except (ValueError, ArithmeticError, MemoryError) as error:
            generation["failure"] = f"{type(error).__name__}: {error}"
        generation.update(
            wall_seconds=time.perf_counter() - wall,
            process_cpu_seconds=time.process_time() - cpu,
            parent_lifetime_peak_working_set_bytes=peak_working_set(),
            directional_predicates_classified=2
            * len(points)
            * generation["classified_events"],
            directional_predicates_attempted=2
            * len(points)
            * generation["sampled_events"],
        )
        write_json(destination / f"case_{case}_generation.json", generation)
        generated.append(generation)
        del after, before
        for n in LEVELS:
            archive = destination / f"case_{case}_n{n}.npz"
            row = {
                "case": case,
                "amplitude": amplitude,
                "rotation_rate": rotation,
                "seed": SEED + case,
                "measurement_n": n,
            }
            if not archive.exists():
                row["failure"] = generation.get(
                    "failure", "Measurement prefix not completed"
                )
            else:
                output = archive.with_suffix(".json")
                result = invoke(
                    archive,
                    pilot / f"truth_{case}_nested.npz",
                    output,
                    os.environ.copy(),
                )
                row.update(
                    {key: value for key, value in result.items() if key != "full"}
                )
                if "full" in result:
                    row.update(result["full"]["evaluation"])
                    row["variance_price_V"] = n * row["largest_tidal_variance"]
                row.update(
                    archive_sha256=sha(archive), archive_bytes=archive.stat().st_size
                )
                if output.exists():
                    row["full_record_sha256"] = sha(output)
            write_json(destination / f"case_{case}_n{n}_result.json", row)
            rows.append(row)
            print(json.dumps(row), flush=True)
    verify_bound_files()
    write_json(
        destination / "summary.json",
        {
            **manifest,
            "results": rows,
            "generation": generated,
            "fresh_events_sampled": sum(g["sampled_events"] for g in generated),
            "completed_assignments": len(rows),
        },
    )
    status(state="COMPLETED", assignments=len(rows))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--pilot", type=Path, default=ROOT / "outputs/intrinsic_weyl_nested_pilot_v1"
    )
    parser.add_argument(
        "--destination",
        type=Path,
        default=ROOT / "outputs/intrinsic_weyl_high_count_v1",
    )
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.check:
        print(json.dumps(selfcheck()))
    elif args.worker:
        worker(args.worker, args.truth, args.output)
    else:
        run(args.pilot, args.destination)
