"""Route A, second instrument: wall participation from an order-preserving chain.

See design_v2.md. Generator coordinates; a diagnostic of the background order
law, not an observer.

The chain of flip_probe.chain is run from the sampled configuration, which is
an exact draw from the conditional law given its order C. Each move takes one
event along a uniformly random line and finds the segment on which the order
is unchanged. Each end of that segment lies on the cone of one other event:
just beyond it exactly one relation has changed. That partner is a wall of
the moved event, found with exact geometry and no linear programme.

For a tracked event the distribution of the partner over chain time is its
wall-hit distribution. Its participation is the inverse of the sum of squared
probabilities, estimated from products of frequencies in different blocks of
the chain so that counting noise does not bias it.

    python hit_probe.py --stage control
    python hit_probe.py --stage main_a
    python hit_probe.py --stage main_b
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
BLOCKS = 8
SWEEPS_PER_EVENT = 16  # total sweeps = 16 N
TRACKED = 64
STAGES = {
    "control": ("flat_1p1", (64, 128, 256), 8),
    "main_a": ("background_q1", (128, 256, 512, 1024), 12),
    "main_b": ("background_q1", (2048,), 12),
}


def run_block(x, geom, rng, sweeps, row_of, counts):
    """Advance the chain in place; count chord-end partners of tracked events."""
    total = len(x)
    n = total - 2
    for _ in range(sweeps):
        for k in rng.permutation(n):
            direction = rng.normal(size=x.shape[1])
            direction /= np.linalg.norm(direction)
            delta = x - x[k]
            a = direction[0] ** 2 - np.sum(direction[1:] ** 2)
            if abs(a) < 1e-9:
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
            row = row_of[k]
            if row >= 0:
                counts[row, hi % total] += 1
                counts[row, lo % total] += 1
            t = (lower + (upper - lower) * rng.random()) * (1 - 1e-9)
            proposal = x[k] + t * direction
            if math.log(rng.random()) < geom.log_weight(proposal) - geom.log_weight(
                x[k]
            ):
                x[k] = proposal


def probe(task):
    name, n, seed = task
    started = time.time()
    geom = fp.GEOMETRIES[name]
    rng = np.random.default_rng(np.random.SeedSequence([20261125, geom.dim, n, seed]))
    xa = fp.sample(geom, n, rng)
    rel = fp.order_matrix(xa)
    link, crit = fp.links_and_critical(rel)
    candidate = link | link.T | crit | crit.T
    tracked = np.sort(rng.choice(n, size=min(TRACKED, n), replace=False))
    row_of = -np.ones(n, dtype=np.int64)
    row_of[tracked] = np.arange(len(tracked))
    x = xa.copy()
    sweeps = SWEEPS_PER_EVENT * n // BLOCKS
    counts = np.zeros((BLOCKS, len(tracked), n + 2), dtype=np.int64)
    drift = []
    for block in range(BLOCKS):
        run_block(x, geom, rng, sweeps, row_of, counts[block])
        if not np.array_equal(fp.order_matrix(x), rel):
            raise RuntimeError("chain left the order cell")
        drift.append(float(np.sqrt(np.mean(np.sum((x - xa)[:n] ** 2, axis=1)))))
    events = []
    if geom.dim == 2:
        rank_u = np.argsort(np.argsort(xa[:n, 0] + xa[:n, 1]))
        rank_v = np.argsort(np.argsort(xa[:n, 0] - xa[:n, 1]))
    for row, event in enumerate(tracked):
        inner = counts[:, row, :n].astype(float)
        hits = inner.sum(axis=1)
        freq = inner / hits[:, None]
        cross = freq @ freq.T
        by_lag = [
            float(np.mean(np.diagonal(cross, offset=lag))) for lag in range(1, BLOCKS)
        ]
        pooled = inner.sum(axis=0)
        outside = float(pooled[~candidate[event, :n]].sum())
        record = {
            "n": n,
            "seed": seed,
            "event": int(event),
            "F": int(candidate[event, :n].sum()),
            "hits_per_block": float(hits.mean()),
            "tip_hit_fraction": float(
                counts[:, row, n:].sum() / counts[:, row, :].sum()
            ),
            "hits_outside_candidates": outside,
            "partners_hit": int(np.count_nonzero(pooled)),
            "overlap_by_lag": by_lag,
            "same_block": float(np.mean(np.diagonal(cross))),
        }
        if geom.dim == 2:
            record["interior"] = bool(
                0 < rank_u[event] < n - 1 and 0 < rank_v[event] < n - 1
            )
        events.append(record)
    meta = {
        "n": n,
        "seed": seed,
        "geometry": name,
        "sweeps": sweeps * BLOCKS,
        "relations": int(rel[:n, :n].sum()),
        "chain_rms_displacement": drift,
        "seconds": round(time.time() - started, 1),
    }
    return meta, events


def git(*args):
    return subprocess.check_output(["git", *args], cwd=HERE).decode().strip()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=sorted(STAGES), required=True)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    dirty = git("status", "--porcelain")
    if dirty:
        raise SystemExit("refusing to run: worktree is not clean\n" + dirty)
    stamp = git("rev-parse", "HEAD")
    output = HERE / "evidence_v2" / args.stage
    if output.exists():
        raise SystemExit(f"refusing to overwrite {output}")
    name, sizes, orders = STAGES[args.stage]
    output.mkdir(parents=True)
    sources = ("flip_probe.py", "hit_probe.py", "analyse_hits.py", "design_v2.md")
    manifest = {
        "stage": args.stage,
        "status": "in_progress",
        "code_version": stamp,
        "worktree_clean_at_start": True,
        "utc_started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "device": "cpu",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "geometry": name,
        "sizes": list(sizes),
        "orders_per_size": orders,
        "tracked_events_per_order": TRACKED,
        "blocks": BLOCKS,
        "sweeps": "16 N",
        "files": {
            s: hashlib.sha256((HERE / s).read_bytes()).hexdigest() for s in sources
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", newline="\n"
    )
    tasks = [
        (name, n, seed) for n in sorted(sizes, reverse=True) for seed in range(orders)
    ]
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for meta, events in pool.map(probe, tasks):
            results.append((meta, events))
            print(meta["n"], meta["seed"], meta["seconds"], flush=True)
    results.sort(key=lambda item: (item[0]["n"], item[0]["seed"]))
    with (output / "orders.jsonl").open("w", encoding="utf-8", newline="\n") as stream:
        for meta, _ in results:
            stream.write(json.dumps(meta) + "\n")
    with (output / "events.jsonl").open("w", encoding="utf-8", newline="\n") as stream:
        for _, events in results:
            for event in events:
                stream.write(json.dumps(event) + "\n")
    manifest["status"] = "complete"
    manifest["utc_finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    manifest["events_sha256"] = hashlib.sha256(
        (output / "events.jsonl").read_bytes()
    ).hexdigest()
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", newline="\n"
    )


if __name__ == "__main__":
    main()
