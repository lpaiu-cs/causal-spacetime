"""Fixed order-only replay of thread settings; no new scientific samples."""

import argparse
import cProfile
import hashlib
import json
import math
import os
import pstats
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_DEFAULT_NUM_THREADS",
    "GOTO_NUM_THREADS",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as output:
        output.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def differences(actual, expected, path=""):
    """Return every discrepancy under the frozen numeric/discrete rule."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or actual.keys() != expected.keys():
            return {path: "dictionary keys differ"}
        result = {}
        for key in expected:
            result.update(differences(actual[key], expected[key], f"{path}/{key}"))
        return result
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            return {path: "array lengths differ"}
        result = {}
        for i, (left, right) in enumerate(zip(actual, expected, strict=True)):
            result.update(differences(left, right, f"{path}/{i}"))
        return result
    if isinstance(expected, float):
        same = isinstance(actual, (float, int)) and math.isclose(
            actual, expected, rel_tol=1e-8, abs_tol=1e-8
        )
    else:
        same = type(actual) is type(expected) and actual == expected
    return {} if same else {path: {"actual": actual, "expected": expected}}


def worker(archive, output):
    import intrinsic_weyl_joint as joint

    profiler = cProfile.Profile()
    wall, cpu = time.perf_counter(), time.process_time()
    fit = profiler.runcall(
        joint.fit_order_archive, archive, model_order=5, square_correction=True
    )
    inference = joint.confidence_ball(fit)
    wall, cpu = time.perf_counter() - wall, time.process_time() - cpu
    forbidden = ("rotating_plane_wave", "window_response", "p14_plane_wave")
    if any(name in sys.modules for name in forbidden) or any(
        name.startswith(("intrinsic_weyl_joint_audit", "intrinsic_weyl_joint_pilot"))
        for name in sys.modules
    ):
        raise RuntimeError("Oracle or generator imported in a pure replay")
    stats = pstats.Stats(profiler).stats
    profile = {
        name: sum(value[3] for key, value in stats.items() if key[2] == name)
        for name in ("read_order_clocks", "_clock_second_moment", "fit_joint")
    }
    write_json(
        output,
        {
            "fit": fit,
            "inference": inference,
            "wall_seconds": wall,
            "process_cpu_seconds": cpu,
            "profile_cumulative_seconds": profile,
            "thread_environment": {
                key: os.environ.get(key) for key in THREAD_VARIABLES
            },
            "oracle_or_generator_imported": False,
        },
    )


def run(folder, destination):
    assert not differences([1.0], [1.0 + 1e-9])
    assert differences([1.0], [1.0 + 1e-6])
    assert differences({"rank": 10}, {"rank": 9})
    assert differences(True, 1)
    assert differences(float("nan"), 0.0)
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip():
        raise RuntimeError("Commit the fixed replay plan and code before execution")
    base_environment = os.environ.copy()
    if any(base_environment.get(key) is not None for key in THREAD_VARIABLES):
        raise RuntimeError(
            "The fixed default comparison expects unset thread variables"
        )
    destination.mkdir(parents=True, exist_ok=False)
    sources = [Path(__file__)] + [
        Path(__file__).with_name(name)
        for name in (
            "intrinsic_weyl_joint.py",
            "intrinsic_weyl.py",
            "intrinsic_weyl_square.py",
        )
    ]
    sources.append(ROOT / "docs/paper/paper_a/intrinsic_weyl/nested_pilot_results.md")
    source_hashes = {p.relative_to(ROOT).as_posix(): sha(p) for p in sources}
    manifest = {
        "code_version": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256": source_hashes,
        "pilot_manifest_sha256": sha(folder / "manifest.json"),
        "scope": "Numerical replay and descriptive profiled cost at N=262144",
        "fresh_events": 0,
        "continuous_field_tolerance": {"absolute": 1e-8, "relative": 1e-8},
        "runtime_counter_reported_separately": "function_evaluations",
        "concurrent_campaign_settings_changed": False,
    }
    write_json(destination / "manifest.json", manifest)
    schedule = [("default_before", 5)] + [("single", case) for case in range(8)]
    schedule.append(("default_after", 5))
    rows = []
    for mode, case in schedule:
        archive = folder / f"case_{case}_nested.npz"
        saved = folder / f"case_{case}_nested_full.json"
        output = destination / f"{mode}_{case}.json"
        row = {
            "mode": mode,
            "case": case,
            "archive_sha256": sha(archive),
            "original_wrapper_sha256": sha(saved),
        }
        environment = base_environment.copy()
        if mode == "single":
            environment.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
        child = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--worker",
                str(archive),
                "--output",
                str(output),
            ],
            cwd=ROOT,
            env=environment,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
        if child.returncode:
            row.update(returncode=child.returncode, failure=child.stderr)
        else:
            returned = json.loads(output.read_text())
            # Saved evaluation/truth fields never enter the pure subprocess.
            original = json.loads(saved.read_text())
            fit, expected = dict(returned["fit"]), dict(original["fit"])
            row["function_evaluations"] = {
                "replay": fit.pop("function_evaluations"),
                "original": expected.pop("function_evaluations"),
            }
            row["fit_differences"] = differences(fit, expected)
            row["confidence_differences"] = differences(
                returned["inference"], original["inference"]
            )
            row["all_original_fields_exact"] = (
                returned["fit"] == original["fit"]
                and returned["inference"] == original["inference"]
            )
            row["within_frozen_tolerance"] = (
                not row["fit_differences"] and not row["confidence_differences"]
            )
            for key in (
                "wall_seconds",
                "process_cpu_seconds",
                "profile_cumulative_seconds",
                "thread_environment",
                "oracle_or_generator_imported",
            ):
                row[key] = returned[key]
            row["replay_file_sha256"] = sha(output)
        if (
            sha(archive) != row["archive_sha256"]
            or sha(saved) != row["original_wrapper_sha256"]
        ):
            raise RuntimeError("An original pilot input changed during replay")
        rows.append(row)
        write_json(destination / f"{mode}_{case}_comparison.json", row)
        print(
            json.dumps(
                {
                    key: row.get(key)
                    for key in (
                        "mode",
                        "case",
                        "wall_seconds",
                        "process_cpu_seconds",
                        "within_frozen_tolerance",
                        "failure",
                    )
                }
            ),
            flush=True,
        )
    # All pure subprocesses have returned before the parent opens truth values.
    from intrinsic_weyl_joint import tidal_distance

    for row in rows:
        if "failure" in row:
            continue
        returned = json.loads(
            (destination / f"{row['mode']}_{row['case']}.json").read_text()
        )
        original = json.loads(
            (folder / f"case_{row['case']}_nested_full.json").read_text()
        )
        error = tidal_distance(
            returned["inference"]["coefficients_xy_positive"],
            original["evaluation"]["truth_coefficients_xy_positive"],
        )
        row["ball_contains_truth"] = error <= returned["inference"]["radius"]
        row["old_ball_flag_preserved"] = (
            row["ball_contains_truth"]
            == original["evaluation"]["nominal_ball_contains_truth"]
        )
        row["full_tidal_error"] = error
    if any(sha(ROOT / name) != value for name, value in source_hashes.items()):
        raise RuntimeError("Replay sources changed during execution")
    manifest.update(cases=rows, completed_replays=len(rows))
    manifest["all_single_replays_match"] = all(
        row.get("within_frozen_tolerance", False)
        and row.get("old_ball_flag_preserved", False)
        for row in rows
        if row["mode"] == "single"
    )
    write_json(destination / "summary.json", manifest)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--folder", type=Path, default=ROOT / "outputs/intrinsic_weyl_nested_pilot_v1"
    )
    parser.add_argument(
        "--destination",
        type=Path,
        default=ROOT / "outputs/intrinsic_weyl_thread_replay_v1",
    )
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.output)
    else:
        run(args.folder, args.destination)
