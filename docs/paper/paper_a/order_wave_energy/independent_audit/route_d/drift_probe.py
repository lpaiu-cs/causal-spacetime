"""Route D, first measurement: how far central events move at fixed order.

See design.md. Generator coordinates; a diagnostic of the background order
law, not an observer. No verdict is computed anywhere in this file.

The chain is the Route A chain (hit_probe_v3.run_block with no tracked
event). States are stored on a fixed schedule; for pairs of stored states
the displacement of the central events is split into the best element of
the seven-parameter order-preserving conformal group of the diamond and a
residual.

    python drift_probe.py --stage control
    python drift_probe.py --stage main

Twelve workers, one compute thread each, priority class below normal.
"""

from __future__ import annotations

import os
import pathlib
import sys

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "route_a"))

import argparse  # noqa: E402
import concurrent.futures  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import platform  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402

import flip_probe as fp  # noqa: E402
import hit_probe_v3 as hv  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402
from scipy.linalg import expm  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

WORKERS = 12
SEED_TAG = 20261127
CENTRAL_SCALE = 0.6
BINS = (0.25, 0.5)
MAIN_SIZES, MAIN_ORDERS = (128, 192, 256), 12
FIT_SIZES, FIT_CONFIGS = (128, 192, 256), 24
CHAIN_CONTROL_N, CHAIN_CONTROL_ORDERS = 64, 96
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
GENERATORS = []
for _j in (1, 2, 3):  # boosts
    _m = np.zeros((4, 4))
    _m[0, _j] = _m[_j, 0] = 1.0
    GENERATORS.append(_m)
for _j, _k in ((1, 2), (1, 3), (2, 3)):  # rotations
    _m = np.zeros((4, 4))
    _m[_j, _k], _m[_k, _j] = -1.0, 1.0
    GENERATORS.append(_m)
TIP_P = np.array(fp.BACKGROUND.tip_p)
TIP_Q = np.array(fp.BACKGROUND.tip_q)


def msq(y):
    """Minkowski square, signature (-, +, +, +)."""
    return -(y[..., 0] ** 2) + np.sum(y[..., 1:] ** 2, axis=-1)


W0 = (TIP_P - TIP_Q) / msq(TIP_P - TIP_Q)


def group_map(params, x, with_log_jacobian=False):
    """g = phi^-1 o (lambda Lambda) o phi, phi(x) = I(x - q) - I(p - q).

    params = (log lambda, three boosts, three rotations). Fixes both tips and
    preserves the causal order of events inside the diamond.
    """
    lam = math.exp(params[0])
    lorentz = expm(sum(c * g for c, g in zip(params[1:], GENERATORS, strict=True)))
    y = x - TIP_Q
    u = y / msq(y)[:, None] - W0
    w = lam * (u @ lorentz.T) + W0
    square = msq(w)
    out = TIP_Q + w / square[:, None]
    if with_log_jacobian:
        # conformal factors: inversion |y.y|^-4, dilation lambda^4, inversion
        return out, -4 * np.log(np.abs(msq(y))) + 4 * params[0] - 4 * np.log(
            np.abs(square)
        )
    return out


def fit_group(xa, xb):
    """Least-squares group element taking xa to xb. Returns params, residuals."""
    result = least_squares(
        lambda p: (group_map(p, xa) - xb).ravel(),
        np.zeros(7),
        method="lm",
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=4000,
    )
    return result.x, group_map(result.x, xa) - xb, bool(result.success)


def central_mask(x):
    low = (TIP_P + TIP_Q) / 2 - CENTRAL_SCALE * (TIP_Q - TIP_P) / 2
    high = (TIP_P + TIP_Q) / 2 + CENTRAL_SCALE * (TIP_Q - TIP_P) / 2

    def inside(d):
        return d[:, 0] > np.sqrt(np.sum(d[:, 1:] ** 2, axis=1))

    return inside(x - low) & inside(high - x)


def family_hessian(x):
    """Minus the second derivative of the log conditional density along the
    seven group parameters, at the identity, by central differences."""

    def log_density(p):
        moved, log_jac = group_map(p, x, with_log_jacobian=True)
        m = (moved[:, 0] + moved[:, 1]) / math.sqrt(2.0)
        return float(np.sum(-2.0 * np.log1p((m - fp.TAN_HALF) ** 2) + log_jac))

    h = 1e-3
    out = np.zeros((7, 7))
    for i in range(7):
        for j in range(i, 7):
            total = 0.0
            for si in (1, -1):
                for sj in (1, -1):
                    p = np.zeros(7)
                    p[i] += si * h
                    p[j] += sj * h
                    total += si * sj * log_density(p)
            out[i, j] = out[j, i] = -total / (4 * h * h)
    return out


def run_chain(xa, rel, geom, rng, times):
    """Advance the chain; return {sweeps: positions} at the requested times."""
    n = len(xa) - 2
    x = xa.copy()
    idle = -np.ones(n, dtype=np.int64)
    tally = {"skipped": 0, "outside": 0}
    states, now = {0: xa[:n].copy()}, 0
    for target in sorted(set(times)):
        hv.run_block(x, geom, rng, target - now, idle, None, None, None, None, tally)
        now = target
        if not np.array_equal(fp.order_matrix(x), rel):
            raise RuntimeError("chain left the order cell")
        states[target] = x[:n].copy()
    return states, tally["skipped"]


def pair_stats(xa, xb, central, outer):
    a, b = xa[central], xb[central]
    count = len(a)
    params, residual, success = fit_group(a, b)
    factor = 4 * count / (4 * count - 7)
    wave = (residual[:, 0] + residual[:, 1]) / math.sqrt(2.0)
    distance = np.sqrt(np.sum((a[:, None, :] - a[None, :, :]) ** 2, axis=2))
    apart = np.sum((residual[:, None, :] - residual[None, :, :]) ** 2, axis=2)
    upper = np.triu(np.ones((count, count), dtype=bool), 1)
    edges = (0.0, *BINS, np.inf)
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        pick = upper & (distance >= lo) & (distance < hi)
        bins.append(float(apart[pick].mean()) if pick.any() else None)
    out = group_map(params, xa[outer]) - xb[outer]
    return {
        "fit_success": success,
        "params": params.tolist(),
        "central_total": float(np.mean(np.sum((b - a) ** 2, axis=1))),
        "central_residual": float(np.mean(np.sum(residual**2, axis=1)) * factor),
        "central_residual_wave": float(np.mean(wave**2) * factor),
        "outer_total": float(np.mean(np.sum((xb[outer] - xa[outer]) ** 2, axis=1))),
        "outer_residual": float(np.mean(np.sum(out**2, axis=1))),
        "pair_difference_bins": bins,
    }


def main_task(task):
    n, seed, folder = task
    started = time.time()
    priority = hv.lower_priority()
    geom = fp.BACKGROUND
    rng = np.random.default_rng(np.random.SeedSequence([SEED_TAG, 4, n, seed]))
    xa = fp.sample(geom, n, rng)
    rel = fp.order_matrix(xa)
    central = central_mask(xa[:n])
    related = rel[:n, :n]
    outer = ~central & related.any(axis=0) & related.any(axis=1)
    short = [n * 2**j for j in range(7)]
    long = [n * n // 2, n * n, 2 * n * n]
    lags = sorted(set(short + long))
    marks = sorted({t for lag in lags for t in (lag, 2 * lag)})
    states, skipped = run_chain(
        xa, rel, geom, rng, marks + [t + 1 for t in marks] + [1]
    )
    one_sweep = [
        np.mean(np.sum((states[t + 1] - states[t])[central] ** 2, axis=1))
        for t in [0, *marks]
    ]
    hessian = family_hessian(xa[:n])
    values, vectors = np.linalg.eigh(hessian)
    rows = []
    for lag in lags:
        for start in (0, lag):
            row = pair_stats(states[start], states[start + lag], central, outer)
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


def fit_control_task(task):
    n, seed, folder = task
    priority = hv.lower_priority()
    rng = np.random.default_rng(np.random.SeedSequence([SEED_TAG, 9, n, seed]))
    xa = fp.sample(fp.BACKGROUND, n, rng)
    x = xa[:n]
    central = central_mask(x)
    direction = rng.normal(size=7)
    direction /= np.linalg.norm(direction)
    scale = 0.1
    for _ in range(6):  # scale the element to the target central displacement
        moved = group_map(scale * direction, x)
        rms = math.sqrt(np.mean(np.sum((moved - x)[central] ** 2, axis=1)))
        scale *= (3.2 / math.sqrt(n)) / rms
    truth = scale * direction
    image = group_map(truth, x)
    same_order = bool(
        np.array_equal(fp.order_matrix(np.vstack((image, xa[n:]))), fp.order_matrix(xa))
    )
    sigma = 10.0 / n / 2.0  # per coordinate; rms 10/N per event
    start = x + sigma * rng.normal(size=x.shape)
    end = image + sigma * rng.normal(size=x.shape)
    a, b = start[central], end[central]
    injected = np.mean(np.sum((group_map(truth, a) - b) ** 2, axis=1))
    _, residual, success = fit_group(a, b)
    count = len(a)
    recovered = np.mean(np.sum(residual**2, axis=1)) * 4 * count / (4 * count - 7)
    record = {
        "n": n,
        "seed": seed,
        "central_events": count,
        "same_order": same_order,
        "fit_success": success,
        "ratio": float(math.sqrt(recovered / injected)),
        "priority_class": priority,
    }
    target = pathlib.Path(folder) / "tasks" / f"fit_{n}_{seed}.json"
    target.write_text(json.dumps(record), newline="\n")
    return "fit", n, seed, 0.0, priority


def chain_control_task(task):
    n, seed, folder = task
    started = time.time()
    priority = hv.lower_priority()
    geom = fp.FLAT2
    rng = np.random.default_rng(np.random.SeedSequence([SEED_TAG, 2, n, seed]))
    xa = fp.sample(geom, n, rng)
    rel = fp.order_matrix(xa)
    lags = [n * n // 4, n * n // 2, n * n, 2 * n * n, 4 * n * n]
    states, skipped = run_chain(xa, rel, geom, rng, [*lags, 8 * n * n])
    exact = 1.0 / (3.0 * (n + 1))
    record = {
        "n": n,
        "seed": seed,
        "skipped_moves": skipped,
        "ratio_by_lag": {
            str(lag): [
                float(np.mean(np.sum((states[s + lag] - states[s]) ** 2, axis=1)))
                / exact
                for s in (0, lag)
            ]
            for lag in lags
        },
        "priority_class": priority,
        "seconds": round(time.time() - started, 1),
    }
    target = pathlib.Path(folder) / "tasks" / f"chain_{n}_{seed}.json"
    target.write_text(json.dumps(record), newline="\n")
    return "chain", n, seed, record["seconds"], priority


def dispatch(task):
    kind, *rest = task
    return {"main": main_task, "fit": fit_control_task, "chain": chain_control_task}[
        kind
    ](tuple(rest))


def git(*args):
    return subprocess.check_output(["git", *args], cwd=HERE).decode().strip()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("control", "main"), required=True)
    args = parser.parse_args()
    dirty = git("status", "--porcelain")
    if dirty:
        raise SystemExit("refusing to run: worktree is not clean\n" + dirty)
    stamp = git("rev-parse", "HEAD")
    output = HERE / "evidence" / args.stage
    if output.exists():
        raise SystemExit(f"refusing to overwrite {output}")
    if args.stage == "main":
        control = json.loads((HERE / "evidence/control/summary.json").read_text())
        if not control["controls_pass"]:
            raise SystemExit("refusing to run: a control did not pass")
    (output / "tasks").mkdir(parents=True)
    folder = str(output)
    if args.stage == "control":
        tasks = [
            ("chain", CHAIN_CONTROL_N, seed, folder)
            for seed in range(CHAIN_CONTROL_ORDERS)
        ] + [("fit", n, seed, folder) for n in FIT_SIZES for seed in range(FIT_CONFIGS)]
    else:
        tasks = [
            ("main", n, seed, folder)
            for n in sorted(MAIN_SIZES, reverse=True)
            for seed in range(MAIN_ORDERS)
        ]
    sources = ("drift_probe.py", "analyse_drift.py", "design.md")
    manifest = {
        "stage": args.stage,
        "status": "in_progress",
        "code_version": stamp,
        "worktree_clean_at_start": True,
        "utc_started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "device": "cpu",
        "workers": WORKERS,
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
    with concurrent.futures.ProcessPoolExecutor(max_workers=WORKERS) as pool:
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
