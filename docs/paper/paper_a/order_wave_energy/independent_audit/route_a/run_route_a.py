"""Stamped runner for the Route A flip-realisability probe. CPU only.

Usage, from a clean worktree, with an interpreter that has NumPy and SciPy:

    python run_route_a.py --stage control
    python run_route_a.py --stage pilot
    python run_route_a.py --stage main

Each stage refuses to start on a dirty tree (untracked files count), refuses
to overwrite its output folder, and records the commit it ran from. The main
stage reads its number of orders from evidence/pilot/sizing.json, which
analyse_route_a.py writes from the pilot rows and which must be committed
before the main stage runs.
"""

from __future__ import annotations

import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import argparse  # noqa: E402
import concurrent.futures  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import pathlib  # noqa: E402
import platform  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402

import flip_probe as fp  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
EVENTS = 8
STAGES = {
    # geometry, sizes, orders per size, checkpoints, sweeps per checkpoint
    "control": ("flat_1p1", (32, 64, 128), 8, 2, 5),
    "pilot": ("background_q1", (96,), 8, 8, 250),
    "main": ("background_q1", (64, 128, 256, 512), None, 8, 250),
}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=HERE).decode().strip()


def work(task):
    name, n, seed, checkpoints, sweeps = task
    started = time.time()
    meta, rows = fp.probe_order(
        fp.GEOMETRIES[name], n, seed, EVENTS, checkpoints, sweeps, True
    )
    meta["seconds"] = round(time.time() - started, 1)
    return meta, rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=sorted(STAGES), required=True)
    parser.add_argument("--workers", type=int, default=14)
    args = parser.parse_args()

    dirty = git("status", "--porcelain")
    if dirty:
        raise SystemExit("refusing to run: worktree is not clean\n" + dirty)
    output = HERE / "evidence" / args.stage
    if output.exists():
        raise SystemExit(f"refusing to overwrite {output}")
    name, sizes, orders, checkpoints, sweeps = STAGES[args.stage]
    if orders is None:
        sizing = json.loads((HERE / "evidence/pilot/sizing.json").read_text())
        orders = int(sizing["orders_per_size"])
    stamp = git("rev-parse", "HEAD")
    output.mkdir(parents=True)
    sources = ("flip_probe.py", "run_route_a.py", "analyse_route_a.py", "design.md")
    manifest = {
        "stage": args.stage,
        "status": "in_progress",
        "code_version": stamp,
        "worktree_clean_at_start": True,
        "utc_started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "device": "cpu",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "geometry": name,
        "sizes": list(sizes),
        "orders_per_size": orders,
        "events_per_order": EVENTS,
        "checkpoints": checkpoints,
        "sweeps_per_checkpoint": sweeps,
        "files": {
            s: hashlib.sha256((HERE / s).read_bytes()).hexdigest() for s in sources
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    tasks = [
        (name, n, seed, checkpoints, sweeps)
        for n in sorted(sizes, reverse=True)
        for seed in range(orders)
    ]
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for meta, rows in pool.map(work, tasks):
            results.append((meta, rows))
            print(meta["n"], meta["seed"], len(rows), meta["seconds"], flush=True)
    results.sort(key=lambda item: (item[0]["n"], item[0]["seed"]))
    with (output / "orders.jsonl").open("w", encoding="utf-8", newline="\n") as stream:
        for meta, _ in results:
            stream.write(json.dumps(meta) + "\n")
    with (output / "rows.jsonl").open("w", encoding="utf-8", newline="\n") as stream:
        for _, rows in results:
            for row in rows:
                stream.write(json.dumps(row) + "\n")
    manifest["status"] = "complete"
    manifest["utc_finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    manifest["rows_sha256"] = hashlib.sha256(
        (output / "rows.jsonl").read_bytes()
    ).hexdigest()
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
