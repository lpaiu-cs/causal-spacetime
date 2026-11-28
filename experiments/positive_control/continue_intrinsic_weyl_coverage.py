"""Continue the interrupted, fixed 512-case dataset without replacing any data."""

import argparse
import hashlib
import importlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import scipy

MANIFEST_SHA = "df14768daefa8f92a9ef8b9fae070cccb996d565d6228679cf30c5e54fc28480"
COMPLETED = range(226)
SAVED = range(226, 230)


def sha(path, *, lf=False):
    content = path.read_bytes()
    return hashlib.sha256(
        content.replace(b"\r\n", b"\n") if lf else content
    ).hexdigest()


def write_json(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def git(root, *arguments):
    return subprocess.check_output(["git", *arguments], cwd=root, text=True).strip()


def load_frozen(root, name):
    source = root / "experiments/positive_control"
    sys.path.insert(0, str(source))
    module = importlib.import_module(name)
    if Path(module.__file__).resolve().parent != source.resolve():
        raise RuntimeError(f"Wrong source checkout for {name}")
    return module


def preflight(root, source):
    manifest_path = source / "manifest.json"
    if sha(manifest_path) != MANIFEST_SHA:
        raise ValueError("This continuation is only for the recorded frozen manifest")
    manifest = json.loads(manifest_path.read_text())
    expected = {"manifest.json"}
    for case in range(230):
        expected.update((f"case_{case}_nested.npz", f"truth_{case}_nested.npz"))
        if case in COMPLETED:
            expected.update((f"case_{case}_full.json", f"case_{case}_result.json"))
    if {p.name for p in source.iterdir()} != expected:
        raise ValueError("Interrupted checkpoint differs from the recorded 226+4 state")
    if git(root, "status", "--porcelain"):
        raise RuntimeError("Frozen source checkout must remain clean")
    git(root, "merge-base", "--is-ancestor", manifest["code_version"], "HEAD")
    for name, digest in manifest["source_sha256_lf"].items():
        assert sha(root / name, lf=True) == digest, name
    runtime = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    assert runtime == manifest["runtime"]
    assert {
        k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")
    } == manifest["worker_thread_environment"]
    rng = np.random.default_rng(2026100950)
    population = np.column_stack(
        (rng.choice([0.0, 0.25, 0.75], 512), rng.choice([-0.5, 0.0, 0.5], 512))
    ).tolist()
    assert manifest["population"] == population
    for case in COMPLETED:
        row = json.loads((source / f"case_{case}_result.json").read_text())
        assert (row["case"], row["seed"]) == (case, 2026101000 + case)
        assert [row["amplitude"], row["rotation_rate"]] == population[case]
        for name, digest in row["artifact_sha256"].items():
            assert sha(source / name) == digest, name
    for case in SAVED:
        with np.load(source / f"case_{case}_nested.npz", allow_pickle=False) as saved:
            assert set(saved.files) == {
                "selection_order",
                "selection_design",
                "selected_ids",
                "after",
                "before",
                "measurement_n",
            }
            # Read every array now, so a truncated ZIP member fails before work starts.
            arrays = {name: saved[name] for name in saved.files}
            assert int(arrays["measurement_n"]) == 262144
            assert str(arrays["selection_design"]) == "nested"
            assert arrays["selection_order"].shape == (4098, 4098)
            assert (
                arrays["after"].shape
                == arrays["before"].shape
                == (len(arrays["selected_ids"]), 32768)
            )
        with np.load(source / f"truth_{case}_nested.npz", allow_pickle=False) as saved:
            assert saved["points"].shape == (len(arrays["selected_ids"]), 4)
            assert [float(saved["amplitude"]), float(saved["rotation_rate"])] == (
                population[case]
            )
    return manifest, {name: sha(source / name) for name in sorted(expected)}


def pure_replay(arguments):
    root, source = arguments
    joint = load_frozen(root, "intrinsic_weyl_joint")
    fit = joint.fit_order_archive(
        source / "case_225_nested.npz", model_order=5, square_correction=True
    )
    assert "rotating_plane_wave" not in sys.modules
    assert "intrinsic_weyl_rotating_pilot" not in sys.modules
    return {"fit": fit, "inference": joint.confidence_ball(fit)}


def continue_case(arguments):
    root, destination, case, amplitude, rotation, seed = arguments
    frozen = load_frozen(root, "intrinsic_weyl_jet_coverage")
    path = destination / f"case_{case}_nested.npz"
    truth = destination / f"truth_{case}_nested.npz"
    full = destination / f"case_{case}_full.json"
    result_path = destination / f"case_{case}_result.json"
    if full.exists() or result_path.exists():
        raise FileExistsError("Never re-evaluate a returned fit or replace a result")
    if case not in SAVED:
        if path.exists() or truth.exists():
            raise FileExistsError("Never replace previously generated inputs")
        return frozen.run_case((destination, case, amplitude, rotation, seed), True)
    # Only these four saved input pairs lack fits. Reuse them byte-for-byte.
    started = time.perf_counter()
    with np.load(path, allow_pickle=False) as saved:
        landmarks = len(saved["selected_ids"])
    row = {
        "case": case,
        "amplitude": amplitude,
        "rotation_rate": rotation,
        "seed": seed,
        "generation": {
            "selection_n": 4096,
            "measurement_n_generated": 262144,
            "union_landmarks_including_tips": landmarks,
            "landmarks": landmarks - 2,
        },
        "reused_saved_input": True,
        "seconds_scope": (
            "Continuation evaluation only; original generation time unavailable"
        ),
    }
    try:
        result = frozen.evaluate(path, truth, model_order=5, square_correction=True)
        write_json(full, result)
        row.update(result["evaluation"])
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        row["failure"] = f"{type(error).__name__}: {error}"
    row["seconds"] = time.perf_counter() - started
    row["artifact_sha256"] = {p.name: sha(p) for p in (path, truth, full) if p.exists()}
    write_json(result_path, row)
    return row


def run(root, source, destination):
    coordinator_root = Path(__file__).resolve().parents[2]
    if git(coordinator_root, "status", "--porcelain"):
        raise RuntimeError("Commit the continuation and its addendum before execution")
    manifest, preserved = preflight(root, source)
    if destination.exists():
        raise FileExistsError("Keep every existing continuation intact")
    with ProcessPoolExecutor(max_workers=1) as pool:
        replay = pool.submit(pure_replay, (root, source)).result()
    original = json.loads((source / "case_225_full.json").read_text())
    assert replay["fit"] == original["fit"]
    assert replay["inference"] == original["inference"]
    print(
        json.dumps({"preflight": "passed", "pure_case_225_replay": "all fields exact"}),
        flush=True,
    )
    destination.mkdir(parents=True)
    for name, digest in preserved.items():
        shutil.copy2(source / name, destination / name)
        assert sha(destination / name) == digest, name
    receipt = {
        "original_directory": str(source),
        "original_status": "INCOMPLETE_AFTER_INTERRUPTION",
        "source_code_version": manifest["code_version"],
        "coordinator_code_version": git(coordinator_root, "rev-parse", "HEAD"),
        "coordinator_source": Path(__file__).relative_to(coordinator_root).as_posix(),
        "coordinator_sha256_lf": sha(Path(__file__), lf=True),
        "started_utc": datetime.now(UTC).isoformat(),
        "preserved_sha256": preserved,
        "completed_cases_reused": list(COMPLETED),
        "saved_inputs_reused": list(SAVED),
        "new_cases": list(range(230, 512)),
        "pure_replay": {"case": 225, "all_fit_and_inference_fields_exact": True},
        "workers": 4,
    }
    write_json(destination / "continuation.json", receipt)
    arguments = [
        (root, destination, i, *manifest["population"][i], 2026101000 + i)
        for i in range(226, 512)
    ]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for row in pool.map(continue_case, arguments):
            print(
                json.dumps(
                    {
                        k: row[k]
                        for k in (
                            "case",
                            "failure",
                            "bias_condition_met",
                            "nominal_ball_contains_truth",
                            "seconds",
                        )
                        if k in row
                    }
                ),
                flush=True,
            )
    # Both original and copied artifacts must survive; all failures stay included.
    _, final_preserved = preflight(root, source)
    assert preserved == final_preserved
    assert sha(Path(__file__), lf=True) == receipt["coordinator_sha256_lf"]
    for name, digest in preserved.items():
        assert sha(destination / name) == digest, name
    rows = [
        json.loads((destination / f"case_{i}_result.json").read_text())
        for i in range(512)
    ]
    frozen = load_frozen(root, "intrinsic_weyl_jet_coverage")
    verdict = frozen.verdict(rows, 512)
    write_json(
        destination / "summary.json",
        manifest
        | verdict
        | {
            "cases": rows,
            "continuation": receipt,
            "scope": (
                "Fixed 512-case population completed across an interruption; "
                "not accuracy"
            ),
        },
    )
    print(json.dumps(verdict), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--frozen-root", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        _, files = preflight(args.frozen_root.resolve(), args.source.resolve())
        print(
            json.dumps(
                {"preserved_files": len(files), "saved_inputs": 4, "new_cases": 282}
            )
        )
    else:
        status_path = args.output.with_suffix(".status.json")
        status = {
            "pid": os.getpid(),
            "started_utc": datetime.now(UTC).isoformat(),
            "state": "RUNNING",
        }
        write_json(status_path, status)
        try:
            run(
                args.frozen_root.resolve(), args.source.resolve(), args.output.resolve()
            )
        except BaseException as error:
            status.update(state="FAILED", error=f"{type(error).__name__}: {error}")
            raise
        else:
            status["state"] = "COMPLETE"
        finally:
            status["ended_utc"] = datetime.now(UTC).isoformat()
            status_path.write_text(
                json.dumps(status, indent=2) + "\n", encoding="utf-8", newline="\n"
            )
