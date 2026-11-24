"""Route A, third instrument: wall-hit participation, stratified. See design_v3.md.

Generator coordinates; a diagnostic of the background order law, not an
observer. The chain law is that of flip_probe.chain and hit_probe.run_block.
Changes against hit_probe.py, each from the referee report:

- tracked events carry a stratum (central / extremal / other);
- chains are longer than the estimator window where affordable;
- the candidate check uses the exact orientation of the crossed pair;
- a move whose segment end is not finite is skipped and counted.

    python hit_probe_v3.py --stage control
    python hit_probe_v3.py --stage main

Twelve workers, one compute thread each, priority class below normal.
"""

from __future__ import annotations

import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import argparse  # noqa: E402
import concurrent.futures  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import pathlib  # noqa: E402
import platform  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402

import flip_probe as fp  # noqa: E402
import numpy as np  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
WORKERS = 12
BELOW_NORMAL = 0x00004000
WINDOW = 8  # blocks used by the gated estimator
CENTRAL_SCALE = 0.6
MAX_CENTRAL, MAX_OTHER, CONTROL_TRACKED = 64, 32, 64
# stage -> geometry, {N: number of blocks of 2N sweeps}, orders per N
STAGES = {
    "control": ("flat_1p1", {64: 8, 128: 8, 256: 8}, 8),
    "main": (
        "background_q1",
        {64: 64, 128: 64, 256: 64, 512: 32, 1024: 16, 2048: 8},
        12,
    ),
}


def _kernel():
    import ctypes

    kernel = ctypes.windll.kernel32
    kernel.GetCurrentProcess.restype = ctypes.c_void_p
    kernel.SetPriorityClass.argtypes = [ctypes.c_void_p, ctypes.c_ulong]
    kernel.GetPriorityClass.argtypes = [ctypes.c_void_p]
    kernel.GetPriorityClass.restype = ctypes.c_ulong
    return kernel


def lower_priority():
    """Set this process to below-normal priority; return the class in force."""
    if os.name != "nt":
        return None
    kernel = _kernel()
    handle = kernel.GetCurrentProcess()
    kernel.SetPriorityClass(handle, BELOW_NORMAL)
    return int(kernel.GetPriorityClass(handle))


def run_block(x, geom, rng, sweeps, row_of, counts, rel, link, crit, tally):
    """Advance the chain in place; count segment-end partners of tracked events."""
    total = len(x)
    n = total - 2
    for _ in range(sweeps):
        for k in rng.permutation(n):
            direction = rng.normal(size=x.shape[1])
            direction /= np.linalg.norm(direction)
            delta = x - x[k]
            a = direction[0] ** 2 - np.sum(direction[1:] ** 2)
            if abs(a) < 1e-9:
                tally["skipped"] += 1
                continue
            b = delta[:, 0] * direction[0] - delta[:, 1:] @ direction[1:]
            c = delta[:, 0] ** 2 - np.sum(delta[:, 1:] ** 2, axis=1)
            disc = b * b - a * c
            disc[k] = -1.0
            root = np.sqrt(np.where(disc >= 0, disc, np.nan))
            roots = np.concatenate(((b - root) / a, (b + root) / a))
            above = np.where(roots > 0, roots, np.inf)
            below = np.where(roots < 0, roots, -np.inf)
            hi, lo = int(np.argmin(above)), int(np.argmax(below))
            upper, lower = above[hi], below[lo]
            if not (math.isfinite(upper) and math.isfinite(lower)):
                tally["skipped"] += 1
                continue
            row = row_of[k]
            if row >= 0:
                for index, end in ((hi, upper), (lo, lower)):
                    partner = index % total
                    counts[row, partner] += 1
                    if partner >= n:
                        continue
                    if rel[k, partner]:
                        fine = link[k, partner]
                    elif rel[partner, k]:
                        fine = link[partner, k]
                    elif x[k, 0] + end * direction[0] < x[partner, 0]:
                        fine = crit[k, partner]
                    else:
                        fine = crit[partner, k]
                    if not fine:
                        tally["outside"] += 1
            t = (lower + (upper - lower) * rng.random()) * (1 - 1e-9)
            proposal = x[k] + t * direction
            if math.log(rng.random()) < geom.log_weight(proposal) - geom.log_weight(
                x[k]
            ):
                x[k] = proposal


def pair_mean(cross, pairs):
    return float(np.mean([cross[b, c] for b, c in pairs]))


def probe(task):
    name, n, seed, blocks, folder = task
    started = time.time()
    priority = lower_priority()
    geom = fp.GEOMETRIES[name]
    rng = np.random.default_rng(np.random.SeedSequence([20261126, geom.dim, n, seed]))
    xa = fp.sample(geom, n, rng)
    rel = fp.order_matrix(xa)
    link, crit = fp.links_and_critical(rel)
    partners = link | link.T | crit | crit.T
    extremal = ~rel[:n, :n].any(axis=0) | ~rel[:n, :n].any(axis=1)
    if geom.dim == 2:
        central = np.zeros(n, dtype=bool)
        tracked = np.sort(rng.choice(n, size=min(CONTROL_TRACKED, n), replace=False))
        rank_u = np.argsort(np.argsort(xa[:n, 0] + xa[:n, 1]))
        rank_v = np.argsort(np.argsort(xa[:n, 0] - xa[:n, 1]))
        at_u, at_v = np.argsort(rank_u), np.argsort(rank_v)
    else:
        p, q = np.array(geom.tip_p), np.array(geom.tip_q)
        low = (p + q) / 2 - CENTRAL_SCALE * (q - p) / 2
        high = (p + q) / 2 + CENTRAL_SCALE * (q - p) / 2
        inner = np.vstack((xa[:n], low, high))
        central = (fp.margins(inner, n, np.arange(n))[0] > 0) & (
            fp.margins(inner, np.arange(n), n + 1)[0] > 0
        )
        inside, rest = np.nonzero(central)[0], np.nonzero(~central)[0]
        tracked = np.sort(
            np.concatenate(
                (
                    rng.choice(inside, min(MAX_CENTRAL, len(inside)), replace=False),
                    rng.choice(rest, min(MAX_OTHER, len(rest)), replace=False),
                )
            )
        )
    row_of = -np.ones(n, dtype=np.int64)
    row_of[tracked] = np.arange(len(tracked))
    x = xa.copy()
    counts = np.zeros((blocks, len(tracked), n + 2), dtype=np.int64)
    tally = {"skipped": 0, "outside": 0}
    drift = []
    for block in range(blocks):
        run_block(x, geom, rng, 2 * n, row_of, counts[block], rel, link, crit, tally)
        if not np.array_equal(fp.order_matrix(x), rel):
            raise RuntimeError("chain left the order cell")
        drift.append(float(np.sqrt(np.mean(np.sum((x - xa)[:n] ** 2, axis=1)))))
    window = [(b, c) for b in range(WINDOW) for c in range(b + 4, WINDOW)]
    near = [(b, c) for b in range(WINDOW) for c in range(b + 1, min(b + 3, WINDOW))]
    far = [(b, c) for b in range(blocks) for c in range(b + blocks // 2, blocks)]
    events = []
    for row, event in enumerate(tracked):
        hits = counts[:, row, :n].astype(float)
        totals = hits.sum(axis=1)
        if np.any(totals == 0):
            raise RuntimeError("a block has no hit on another event")
        freq = hits / totals[:, None]
        cross = freq @ freq.T
        record = {
            "n": n,
            "seed": seed,
            "event": int(event),
            "central": bool(central[event]),
            "extremal": bool(extremal[event]),
            "F": int(partners[event, :n].sum()),
            "tip_hit_fraction": float(counts[:, row, n:].sum() / counts[:, row].sum()),
            "partners_hit": int(np.count_nonzero(hits.sum(axis=0))),
            "q_window": pair_mean(cross, window),
            "q_near": pair_mean(cross, near),
            "q_long": pair_mean(cross, far),
            "overlap_by_lag": [
                float(np.mean(np.diagonal(cross, offset=lag)))
                for lag in range(1, blocks)
            ],
        }
        if geom.dim == 2:
            u, v = rank_u[event], rank_v[event]
            if 0 < u < n - 1 and 0 < v < n - 1:
                around = {at_u[u - 1], at_u[u + 1], at_v[v - 1], at_v[v + 1]}
                record["klass"] = "four" if len(around) == 4 else "shared"
            else:
                record["klass"] = "edge"
        events.append(record)
    meta = {
        "n": n,
        "seed": seed,
        "geometry": name,
        "blocks": blocks,
        "sweeps": 2 * n * blocks,
        "central_events": int(central.sum()),
        "extremal_events": int(extremal.sum()),
        "outside_candidates": tally["outside"],
        "skipped_moves": tally["skipped"],
        "chain_rms_displacement": drift,
        "priority_class": priority,
        "seconds": round(time.time() - started, 1),
    }
    target = pathlib.Path(folder) / "tasks" / f"{n}_{seed}.json"
    target.write_text(json.dumps({"meta": meta, "events": events}), newline="\n")
    return n, seed, meta["seconds"], priority


def git(*args):
    return subprocess.check_output(["git", *args], cwd=HERE).decode().strip()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=sorted(STAGES), required=True)
    args = parser.parse_args()
    dirty = git("status", "--porcelain")
    if dirty:
        raise SystemExit("refusing to run: worktree is not clean\n" + dirty)
    stamp = git("rev-parse", "HEAD")
    output = HERE / "evidence_v3" / args.stage
    if output.exists():
        raise SystemExit(f"refusing to overwrite {output}")
    name, sizes, orders = STAGES[args.stage]
    (output / "tasks").mkdir(parents=True)
    sources = ("flip_probe.py", "hit_probe_v3.py", "analyse_hits_v3.py", "design_v3.md")
    manifest = {
        "stage": args.stage,
        "status": "in_progress",
        "code_version": stamp,
        "worktree_clean_at_start": True,
        "utc_started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "device": "cpu",
        "workers": WORKERS,
        "threads_per_worker": 1,
        "runner_priority_class": lower_priority(),
        "below_normal_priority_class": BELOW_NORMAL,
        "logical_processors": os.cpu_count(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "geometry": name,
        "blocks_by_size": sizes,
        "orders_per_size": orders,
        "files": {
            s: hashlib.sha256((HERE / s).read_bytes()).hexdigest() for s in sources
        },
    }

    def save():
        (output / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", newline="\n"
        )

    save()
    tasks = sorted(
        (
            (name, n, seed, blocks, str(output))
            for n, blocks in sizes.items()
            for seed in range(orders)
        ),
        key=lambda t: -(t[1] ** 2) * t[3],
    )
    priorities = set()
    with concurrent.futures.ProcessPoolExecutor(max_workers=WORKERS) as pool:
        for n, seed, seconds, priority in pool.map(probe, tasks):
            priorities.add(priority)
            print(n, seed, seconds, priority, flush=True)
    metas, events = [], []
    for _, n, seed, _, _ in sorted(tasks, key=lambda t: (t[1], t[2])):
        path = output / "tasks" / f"{n}_{seed}.json"
        loaded = json.loads(path.read_text())
        metas.append(loaded["meta"])
        events.extend(loaded["events"])
        path.unlink()
    (output / "tasks").rmdir()
    with (output / "orders.jsonl").open("w", encoding="utf-8", newline="\n") as stream:
        stream.writelines(json.dumps(m) + "\n" for m in metas)
    with (output / "events.jsonl").open("w", encoding="utf-8", newline="\n") as stream:
        stream.writelines(json.dumps(e) + "\n" for e in events)
    manifest.update(
        status="complete",
        utc_finished=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        worker_priority_classes=sorted(p for p in priorities if p is not None),
        events_sha256=hashlib.sha256(
            (output / "events.jsonl").read_bytes()
        ).hexdigest(),
    )
    save()


if __name__ == "__main__":
    main()
