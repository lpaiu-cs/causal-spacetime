"""Old-data equivalence gate before a new powered precision experiment."""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = "docs/paper/paper_a/intrinsic_weyl/gpu_equivalence_protocol.md"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def replay(data):
    import intrinsic_weyl_joint as joint
    import numpy as np
    from intrinsic_weyl_gpu_counts import clock_second_moment

    archive = data / "case_4_n262144.npz"
    load = np.load
    opened = []

    def guard(path, *args, **kwargs):
        assert Path(path).resolve() == archive.resolve()
        opened.append(str(Path(path).resolve()))
        return load(path, *args, **kwargs)

    np.load = guard
    forbidden = (
        "rotating_plane_wave",
        "rotating_plane_wave_gpu",
        "window_response",
        "p14_plane_wave",
        "intrinsic_weyl_rotating_pilot",
        "intrinsic_weyl_joint_audit",
        "intrinsic_weyl_joint_pilot",
    )
    assert not any(name in sys.modules for name in forbidden)
    start = time.perf_counter()
    cpu = joint.read_order_clocks(archive)
    cpu_seconds = time.perf_counter() - start
    joint._clock_second_moment = clock_second_moment
    start = time.perf_counter()
    gpu = joint.read_order_clocks(archive)
    gpu_seconds = time.perf_counter() - start
    assert all(np.array_equal(a, b) for a, b in zip(cpu, gpu, strict=True))
    start = time.perf_counter()
    fit = joint.fit_order_archive(archive, model_order=5, square_correction=True)
    fit_seconds = time.perf_counter() - start
    assert not any(name in sys.modules for name in forbidden)
    saved = read(data / "case_4_n262144_pure.json")
    assert fit == saved["fit"], "GPU count path changed the full saved fit"
    return {
        "case": 4,
        "measurement_n": 262144,
        "archive_sha256": sha(archive),
        "all_clock_arrays_bit_exact": True,
        "all_fit_fields_exact": True,
        "cpu_clock_seconds": cpu_seconds,
        "gpu_clock_seconds": gpu_seconds,
        "gpu_path_pure_fit_seconds": fit_seconds,
        "opened_archives": opened,
        "oracle_modules_loaded": [],
    }


def check(data, pilot, output):
    import ast
    import statistics

    import cupy as cp
    import numpy as np
    import scipy
    from intrinsic_weyl_gpu_counts import selfcheck
    from rotating_plane_wave import RotatingWave
    from rotating_plane_wave_gpu import RotatingWaveGPU

    assert not output.exists()
    assert not subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT)
    stamp = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    manifest = read(data / "manifest.json")
    for name, digest in manifest["source_sha256"].items():
        assert sha(ROOT / name) == digest, name
    for name, digest in manifest["pilot_input_sha256"].items():
        assert sha(pilot / name) == digest, name
    assert {k: os.environ.get(k) for k in manifest["thread_environment"]} == manifest[
        "thread_environment"
    ]
    assert np.__version__ == manifest["runtime"]["numpy"]
    assert scipy.__version__ == manifest["runtime"]["scipy"]
    names = (
        "rotating_plane_wave_gpu.py",
        "intrinsic_weyl_gpu_counts.py",
        Path(__file__).name,
    )
    sources = {
        "experiments/positive_control/" + name: sha(
            ROOT / "experiments/positive_control" / name
        )
        for name in names
    }
    sources[PROTOCOL] = sha(ROOT / PROTOCOL)
    methods = {"blocks", "inverse", "action", "relations"}

    def method_trees(path):
        tree = ast.parse(path.read_text())
        cls = next(n for n in tree.body if isinstance(n, ast.ClassDef))
        return {
            n.name: ast.dump(n)
            for n in cls.body
            if isinstance(n, ast.FunctionDef) and n.name in methods
        }

    assert method_trees(
        ROOT / "experiments/positive_control/rotating_plane_wave.py"
    ) == method_trees(ROOT / "experiments/positive_control/rotating_plane_wave_gpu.py")
    start = time.perf_counter()
    result = {
        "code_version": stamp,
        "source_sha256": sources,
        "moment_selfcheck": selfcheck(),
        "cases": [],
        "benchmark": {},
        "scientific_fresh_events": 0,
    }
    for case, (amplitude, rotation) in enumerate(manifest["cases"]):
        with np.load(pilot / f"truth_{case}_nested.npz", allow_pickle=False) as truth:
            points = truth["points"]
            assert (float(truth["amplitude"]), float(truth["rotation_rate"])) == (
                amplitude,
                rotation,
            )
        with np.load(data / f"case_{case}_n262144.npz", allow_pickle=False) as old:
            expected_after, expected_before = old["after"], old["before"]
        wave = RotatingWave(amplitude, rotation)
        accelerated = RotatingWaveGPU(amplitude, rotation, wave.coefficients)
        lengths = np.array([1e-9, 1e-5, 0.1, 0.75, 1.0])
        left = np.tile([0.05, -0.08], (len(lengths), 1))
        right = np.tile([-0.03, 0.02], (len(lengths), 1))
        cost, bound = wave.action(lengths, left, right)
        g_cost, g_bound = accelerated.action(
            cp.asarray(lengths), cp.asarray(left), cp.asarray(right)
        )
        g_cost, g_bound = cp.asnumpy(g_cost), cp.asnumpy(g_bound)
        assert np.all(abs(cost - g_cost) < np.minimum(bound, g_bound))
        p = np.column_stack((np.zeros_like(lengths), np.zeros_like(lengths), left))
        q = np.column_stack((lengths, -cost, right))
        for backend, lp, rq in (
            (wave, p, q),
            (accelerated, cp.asarray(p), cp.asarray(q)),
        ):
            try:
                backend.relations(lp, rq)
            except ArithmeticError:
                pass
            else:
                raise AssertionError("A constructed null-margin query was accepted")
        for sign in (-1, 1):
            shifted = q.copy()
            shifted[:, 1] += sign * 100 * np.maximum(bound, g_bound)
            expected = wave.relations(p, shifted)
            observed = accelerated.relations(cp.asarray(p), cp.asarray(shifted))
            assert all(
                np.array_equal(a, cp.asnumpy(b))
                for a, b in zip(expected, observed, strict=True)
            )
        rng = np.random.default_rng([2026103000 + case, 1])
        for start_event in range(0, 262144, 32768):
            events = np.vstack([wave.sample(4096, rng) for _ in range(8)])
            after, before = accelerated.packed_relations(points, events)
            interval = slice(start_event // 8, (start_event + len(events)) // 8)
            assert np.array_equal(after, expected_after[:, interval])
            assert np.array_equal(before, expected_before[:, interval])
            if case == 2 and start_event == 0:
                cpu_times, gpu_times = [], []
                for _ in range(3):
                    clock = time.perf_counter()
                    f, b = wave.relations(points[:, None], events[None])
                    cpu_packed = np.packbits(f, axis=1), np.packbits(b, axis=1)
                    cpu_times.append(time.perf_counter() - clock)
                    clock = time.perf_counter()
                    gpu_packed = accelerated.packed_relations(points, events)
                    gpu_times.append(time.perf_counter() - clock)
                    assert all(
                        np.array_equal(a, b)
                        for a, b in zip(cpu_packed, gpu_packed, strict=True)
                    )
                result["benchmark"] = {
                    "case": case,
                    "events": len(events),
                    "cpu_seconds": cpu_times,
                    "gpu_seconds_including_transfer": gpu_times,
                    "median_ratio_cpu_over_gpu": statistics.median(cpu_times)
                    / statistics.median(gpu_times),
                }
        result["cases"].append(
            {
                "case": case,
                "events": 262144,
                "landmarks": len(points),
                "packed_relations_bit_exact": True,
                "action_max_difference": float(np.max(abs(cost - g_cost))),
                "null_guard_rejects": True,
                "resolved_margin_signs_match": True,
            }
        )
        print(json.dumps(result["cases"][-1]), flush=True)
    child = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve()),
            "--data",
            str(data),
            "--replay",
        ],
        capture_output=True,
        text=True,
        creationflags=subprocess.CREATE_NO_WINDOW,
    )
    if child.returncode:
        raise RuntimeError(child.stdout + child.stderr)
    result["pure_estimator_replay"] = json.loads(child.stdout)
    for name, digest in sources.items():
        assert sha(ROOT / name) == digest, name
    result.update(
        {
            "wall_seconds": time.perf_counter() - start,
            "verdict": "EQUIVALENT_ON_FROZEN_ENGINEERING_CHECKS",
            "runtime": {
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "cupy": cp.__version__,
                "cuda_runtime": cp.cuda.runtime.runtimeGetVersion(),
                "driver": cp.cuda.runtime.driverGetVersion(),
                "device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
            },
            "gpu_reserved_bytes_at_end": cp.get_default_memory_pool().total_bytes(),
        }
    )
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "output": str(output),
                "sha256": sha(output),
                "benchmark": result["benchmark"],
                "replay": result["pure_estimator_replay"],
                "wall_seconds": result["wall_seconds"],
            }
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--pilot", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    if args.replay:
        print(json.dumps(replay(args.data)))
    else:
        if args.pilot is None or args.output is None:
            parser.error("pilot and output are required")
        check(args.data, args.pilot, args.output)
