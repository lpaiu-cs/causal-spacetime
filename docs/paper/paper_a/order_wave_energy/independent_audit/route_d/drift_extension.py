"""Route D, first measurement, extension to N = 512. See extension_design.md.

Same estimator as drift_probe.py, whose functions are called unchanged. The
only difference is the schedule: chains of N^2 sweeps, lags N 2^j up to 64 N
and N^2 / 2. No verdict is computed anywhere in this file.

    python drift_extension.py --stage control
    python drift_extension.py --stage main

Twelve workers, one compute thread each, priority class below normal.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import pathlib
import platform
import time

import drift_probe as dp
import numpy as np
import scipy

HERE = pathlib.Path(__file__).resolve().parent
N_EXTENSION, ORDERS, FIT_CONFIGS = 512, 12, 24
fp, hv = dp.fp, dp.hv


def schedule(n):
    """Lags and stored sweep counts of the extension."""
    lags = sorted({n * 2**j for j in range(7)} | {n * n // 2})
    marks = sorted({t for lag in lags for t in (lag, 2 * lag)})
    return lags, marks


def extension_task(task):
    """drift_probe.main_task with the extension schedule; nothing else differs."""
    n, seed, folder = task
    started = time.time()
    priority = hv.lower_priority()
    geom = fp.BACKGROUND
    rng = np.random.default_rng(np.random.SeedSequence([dp.SEED_TAG, 4, n, seed]))
    xa = fp.sample(geom, n, rng)
    rel = fp.order_matrix(xa)
    central = dp.central_mask(xa[:n])
    related = rel[:n, :n]
    outer = ~central & related.any(axis=0) & related.any(axis=1)
    lags, marks = schedule(n)
    states, skipped = dp.run_chain(
        xa, rel, geom, rng, marks + [t + 1 for t in marks] + [1]
    )
    one_sweep = [
        np.mean(np.sum((states[t + 1] - states[t])[central] ** 2, axis=1))
        for t in [0, *marks]
    ]
    hessian = dp.family_hessian(xa[:n])
    values, vectors = np.linalg.eigh(hessian)
    rows = []
    for lag in lags:
        for start in (0, lag):
            row = dp.pair_stats(states[start], states[start + lag], central, outer)
            projected = vectors.T @ np.array(row["params"])
            row["family_relaxed_stiff"] = float(
                np.mean(projected[-4:] ** 2 * values[-4:] / 2)
            )
            row.update(n=n, seed=seed, lag=lag, start=start)
            rows.append(row)
    meta = {
        "n": n,
        "seed": seed,
        "central_events": int(central.sum()),
        "outer_events": int(outer.sum()),
        "s1_squared": float(np.mean(one_sweep)),
        "hessian_eigenvalues_over_n": (values / n).tolist(),
        "skipped_moves": skipped,
        "sweeps": max(marks) + 1,
        "priority_class": priority,
        "seconds": round(time.time() - started, 1),
    }
    target = pathlib.Path(folder) / "tasks" / f"main_{n}_{seed}.json"
    target.write_text(json.dumps({"meta": meta, "rows": rows}), newline="\n")
    return "main", n, seed, meta["seconds"], priority


def dispatch(task):
    kind, *rest = task
    worker = extension_task if kind == "main" else dp.fit_control_task
    return worker(tuple(rest))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("control", "main"), required=True)
    args = parser.parse_args()
    dirty = dp.git("status", "--porcelain")
    if dirty:
        raise SystemExit("refusing to run: worktree is not clean\n" + dirty)
    stamp = dp.git("rev-parse", "HEAD")
    output = HERE / "evidence_extension" / args.stage
    if output.exists():
        raise SystemExit(f"refusing to overwrite {output}")
    if args.stage == "main":
        control = json.loads(
            (HERE / "evidence_extension/control/summary.json").read_text()
        )
        if not control["fit_pass"]:
            raise SystemExit("refusing to run: the control did not pass")
    (output / "tasks").mkdir(parents=True)
    folder = str(output)
    if args.stage == "control":
        tasks = [("fit", N_EXTENSION, seed, folder) for seed in range(FIT_CONFIGS)]
    else:
        tasks = [("main", N_EXTENSION, seed, folder) for seed in range(ORDERS)]
    sources = (
        "drift_extension.py",
        "analyse_extension.py",
        "drift_probe.py",
        "analyse_drift.py",
        "extension_design.md",
    )
    manifest = {
        "stage": args.stage,
        "status": "in_progress",
        "code_version": stamp,
        "worktree_clean_at_start": True,
        "utc_started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "device": "cpu",
        "workers": dp.WORKERS,
        "threads_per_worker": 1,
        "runner_priority_class": hv.lower_priority(),
        "logical_processors": os.cpu_count(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "tasks": len(tasks),
        "files": {
            s: hashlib.sha256((HERE / s).read_bytes()).hexdigest() for s in sources
        },
    }

    def save():
        (output / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", newline="\n"
        )

    save()
    priorities = set()
    with concurrent.futures.ProcessPoolExecutor(max_workers=dp.WORKERS) as pool:
        for kind, n, seed, seconds, priority in pool.map(dispatch, tasks):
            priorities.add(priority)
            print(kind, n, seed, seconds, priority, flush=True)
    records = []
    for path in sorted((output / "tasks").glob("*.json")):
        records.append({"file": path.stem, **json.loads(path.read_text())})
        path.unlink()
    (output / "tasks").rmdir()
    with (output / "records.jsonl").open("w", encoding="utf-8", newline="\n") as out:
        out.writelines(json.dumps(r) + "\n" for r in records)
    manifest.update(
        status="complete",
        utc_finished=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        worker_priority_classes=sorted(p for p in priorities if p is not None),
        records_sha256=hashlib.sha256(
            (output / "records.jsonl").read_bytes()
        ).hexdigest(),
    )
    save()


if __name__ == "__main__":
    main()
