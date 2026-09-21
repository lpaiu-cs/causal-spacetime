"""Prospective development stages; confirmation is deliberately a later freeze."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
import time
from pathlib import Path

import numpy as np
import window_response as wr
from scipy.stats import norm, t
from window_response_g1 import ROOT, clean_stamp

AMPLITUDES = (0.0, 1 / 64, 1 / 16, 1 / 4, 1.0, 2.0)
PILOT_N = (200, 800)
DOCUMENTS = ROOT / "docs/paper/paper_a"


def summarize(raw, region, amplitude):
    raw = np.asarray(raw)
    reps = len(raw)
    p, both = raw.mean(axis=0)
    p_estimate = float(p)
    p_se = raw[:, 0].std(ddof=1) / math.sqrt(reps)
    zeta = both - p * p + raw[:, 0].var(ddof=1) / reps
    zeta_se = (raw[:, 1] - 2 * p * raw[:, 0]).std(ddof=1) / math.sqrt(reps)
    if amplitude == 0:
        p = 0.1 if region == "D" else wr.box_coefficients(region)[0]
        p_se = 0.0
        zeta = both - p * p
        zeta_se = raw[:, 1].std(ddof=1) / math.sqrt(reps)
        if region == "D":
            zeta, zeta_se = 0.01, 0.0
    critical = t.ppf(1 - 0.05 / 128, reps - 1)
    return {
        "p": float(p),
        "zeta": float(zeta),
        "p_interval": [
            max(0, float(p - critical * p_se)),
            min(1, float(p + critical * p_se)),
        ],
        "zeta_interval": [
            max(0, float(zeta - critical * zeta_se)),
            float(zeta + critical * zeta_se),
        ],
        "raw_p_estimate": p_estimate,
        "reps": reps,
        "numerical_interval_method": "approximate Bonferroni Student t, alpha=.05/64",
    }


def variance(moment, n):
    p, zeta = moment["p"], moment["zeta"]
    return 2 * p * (1 - p) / (n * (n - 1)) + 4 * (n - 2) * zeta / (n * (n - 1))


def threshold(flat, n):
    half = norm.ppf(0.995) * math.sqrt(variance(flat, n))
    return [flat["p"] - half, flat["p"] + half]


def predict(moment, n, bounds):
    sd = math.sqrt(variance(moment, n))
    return float(
        norm.cdf((bounds[0] - moment["p"]) / sd)
        + norm.sf((bounds[1] - moment["p"]) / sd)
    )


def grid_worker(task):
    ri, ai = task
    region, amplitude = wr.REGIONS[ri], AMPLITUDES[ai]
    started = time.process_time()
    seed = [20260921, 6201, ri, ai]
    raw = wr.quadrature_replicates(region, amplitude, 17, 16, seed)
    row = {
        "region": region,
        "A": amplitude,
        "seed_prefix": seed,
        "power": 17,
        "replicates_p_both": raw.tolist(),
        "moments": summarize(raw, region, amplitude),
    }
    if region == "D":
        # Same input coordinates isolate inverse-CDF discretization.
        coarse = wr.quadrature_replicates(
            region, amplitude, 15, 8, [20260921, 6201, ri, ai, 999], 4097
        )
        fine = wr.quadrature_replicates(
            region, amplitude, 15, 8, [20260921, 6201, ri, ai, 999], 8193
        )
        difference = fine - coarse
        row["cdf_grid_comparison"] = {
            "grids": [4097, 8193],
            "replicate_differences": difference.tolist(),
            "max_abs_replicate_difference": float(abs(difference).max()),
        }
    row["cpu_seconds"] = time.process_time() - started
    return row


def pilot_worker(task):
    ri, ni, ai, flat = task
    region, n, amplitude = wr.REGIONS[ri], PILOT_N[ni], AMPLITUDES[ai]
    seed = [20260921, 6301, ri, ni, ai]
    started = time.process_time()
    rows = []
    for rep in range(80):
        rng = np.random.default_rng([*seed, rep])
        points = wr.sample_points(region, amplitude, n, rng)
        rows.append(wr.statistics(points, amplitude))
    return {
        "region": region,
        "N": n,
        "A": amplitude,
        "seed_prefix": seed,
        "threshold": threshold(flat, n),
        "raw_statistics": rows,
        "cpu_seconds": time.process_time() - started,
    }


def run_stage(stage):
    stamp = clean_stamp()
    output = DOCUMENTS / f"window_response_{stage}_results.json"
    if output.exists():
        raise FileExistsError(output)
    started = time.perf_counter()
    if stage == "grid":
        worker = grid_worker
        tasks = [(ri, ai) for ri in range(4) for ai in range(6)]
    else:
        grid_path = DOCUMENTS / "window_response_grid_results.json"
        grid = json.loads(grid_path.read_text(encoding="utf-8"))
        flats = {
            row["region"]: row["moments"] for row in grid["cells"] if row["A"] == 0
        }
        tasks = [
            (ri, ni, ai, flats[wr.REGIONS[ri]])
            for ri in range(4)
            for ni in range(2)
            for ai in range(6)
        ]
        worker = pilot_worker
    cells = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        pending = [pool.submit(worker, task) for task in tasks]
        for future in concurrent.futures.as_completed(pending):
            row = future.result()
            cells.append(row)
            print(
                json.dumps(
                    {
                        "stage": stage,
                        "done": len(cells),
                        "total": len(tasks),
                        "region": row["region"],
                        "A": row["A"],
                        "cpu_seconds": round(row["cpu_seconds"], 3),
                    }
                ),
                flush=True,
            )
            if sum(r["cpu_seconds"] for r in cells) > 8 * 3600:
                raise RuntimeError("development CPU cap reached")
    cells.sort(
        key=lambda row: (wr.REGIONS.index(row["region"]), row.get("N", 0), row["A"])
    )
    result = {
        "stage": stage,
        "code_version": stamp,
        "workers": 4,
        "logical_cpus": os.cpu_count(),
        "wall_seconds": time.perf_counter() - started,
        "cpu_seconds": sum(row["cpu_seconds"] for row in cells),
        "cells": cells,
    }
    if stage == "grid":
        for region in wr.REGIONS:
            rows = [r for r in cells if r["region"] == region]
            flat = next(r["moments"] for r in rows if r["A"] == 0)
            for row in rows:
                row["N300_prediction"] = predict(
                    row["moments"], 300, threshold(flat, 300)
                )
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "result": str(output),
                "wall_seconds": result["wall_seconds"],
                "cpu_seconds": result["cpu_seconds"],
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("grid", "pilot"))
    run_stage(parser.parse_args().stage)
