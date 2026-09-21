"""Freeze forecasts and power before running independent confirmation posets."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import time

import numpy as np
import window_response as wr
from scipy.stats import beta, binom, norm, t
from window_response_g1 import clean_stamp
from window_response_run import DOCUMENTS, predict, summarize, threshold, variance

ALPHA = 0.05 / 16
FREEZE = DOCUMENTS / "window_response_confirmation_freeze.json"


def read_stage(stage):
    path = DOCUMENTS / f"window_response_{stage}_results.json"
    return json.loads(path.read_text(encoding="utf-8")), hashlib.sha256(
        path.read_bytes()
    ).hexdigest()


def cp(hits, count, alpha=ALPHA):
    return [
        0.0 if hits == 0 else float(beta.ppf(alpha / 2, hits, count - hits + 1)),
        1.0
        if hits == count
        else float(beta.ppf(1 - alpha / 2, hits + 1, count - hits)),
    ]


def prediction_interval(moment, n, bounds):
    pl, ph = moment["p_interval"]
    zl, zh = moment["zeta_interval"]
    ps = [pl, ph] + ([0.5] if pl < 0.5 < ph else [])
    variances = [variance({"p": p, "zeta": z}, n) for p in ps for z in (zl, zh)]
    sd_lo, sd_hi = math.sqrt(min(variances)), math.sqrt(max(variances))
    lower_tails = [
        norm.cdf((bounds[0] - p) / sd) for p in (pl, ph) for sd in (sd_lo, sd_hi)
    ]
    upper_tails = [
        norm.sf((bounds[1] - p) / sd) for p in (pl, ph) for sd in (sd_lo, sd_hi)
    ]
    return [
        float(min(lower_tails) + min(upper_tails)),
        min(1.0, float(max(lower_tails) + max(upper_tails))),
    ]


def forecast_worker(task):
    cell, region, amplitude, power = task
    start = time.process_time()
    seed = [20260921, 6202, cell, power]
    raw = wr.quadrature_replicates(region, amplitude, power, 32, seed)
    row = {
        "region": region,
        "A": amplitude,
        "power": power,
        "seed_prefix": seed,
        "replicates_p_both": raw.tolist(),
        "moments": summarize(raw, region, amplitude),
        "cpu_seconds": time.process_time() - start,
    }
    if region == "D":
        row["independent_volume_integrals"] = wr.diamond_moments(amplitude)
    return row


def pilot_review(pilot, grid):
    lookup = {(r["region"], r["A"]): r["moments"] for r in grid["cells"]}
    incidents, diagnostics, controls = [], [], []
    for cell in pilot["cells"]:
        f = np.array([row["f"] for row in cell["raw_statistics"]])
        lo, hi = cell["threshold"]
        hits = int(((f < lo) | (f > hi)).sum())
        interval = cp(hits, len(f), 0.05 / 48)
        mo = lookup[cell["region"], cell["A"]]
        forecast = predict(mo, cell["N"], [lo, hi])
        scenario = (
            [0.01, 0.03]
            if cell["A"] == 0
            else [max(0, forecast - 0.04), min(1, forecast + 0.04)]
        )
        if interval[1] < scenario[0] or interval[0] > scenario[1]:
            incidents.append(
                [cell["region"], cell["N"], cell["A"], "tail scenario excluded"]
            )
        diagnostics.append(
            {
                "region": cell["region"],
                "N": cell["N"],
                "A": cell["A"],
                "hits": hits,
                "predicted": forecast,
                "pilot_CI": interval,
                "mean_z_diagnostic": float(
                    (f.mean() - mo["p"]) / (f.std(ddof=1) / math.sqrt(len(f)))
                ),
                "SD_ratio": float(f.std(ddof=1) / math.sqrt(variance(mo, cell["N"]))),
            }
        )
    for n in (200, 800):
        arrays = [
            np.array([v["f"] for v in c["raw_statistics"]])
            for a in (0, 1)
            for c in pilot["cells"]
            if c["region"] == "R+" and c["N"] == n and c["A"] == a
        ]
        flat, curved = arrays
        v0, v1 = flat.var(ddof=1) / 80, curved.var(ddof=1) / 80
        df = (v0 + v1) ** 2 / (v0 * v0 / 79 + v1 * v1 / 79)
        difference = float(curved.mean() - flat.mean())
        half = float(t.ppf(0.975, df) * math.sqrt(v0 + v1))
        controls.append(
            {
                "N": n,
                "difference": difference,
                "CI95": [difference - half, difference + half],
            }
        )
        if difference - half <= 0:
            incidents.append(["R+", n, 1, "positive control failed"])
    return {
        "incidents": incidents,
        "positive_controls": controls,
        "diagnostics": diagnostics,
        "scope": "Pilot compatibility does not establish the power scenarios.",
    }


def acceptance(n, low, high):
    k = np.arange(n + 1)
    lower_ok = (
        np.ones(n + 1, dtype=bool) if low <= 0 else binom.sf(k - 1, n, low) <= ALPHA / 2
    )
    upper_ok = (
        np.ones(n + 1, dtype=bool) if high >= 1 else binom.cdf(k, n, high) <= ALPHA / 2
    )
    accepted = k[lower_ok & upper_ok]
    return (int(accepted[0]), int(accepted[-1])) if len(accepted) else None


def acceptance_power(n, probability, interval):
    if interval is None:
        return 0.0
    lo, hi = interval
    return float(
        max(0, binom.cdf(hi, n, probability) - binom.cdf(lo - 1, n, probability))
    )


def power_table(cells, pilot):
    prices = {}
    for region in wr.REGIONS:
        prices[region] = max(
            c["cpu_seconds"] / 80 / math.comb(c["N"], 2)
            for c in pilot["cells"]
            if c["region"] == region
        )
    results = []
    for m in (1200, 1600, 2000, 2400, 3200, 4800, 6400, 9600):
        flat_accept = acceptance(m, 0, 0.05)
        curve_accept = [
            acceptance(
                m,
                c["prediction"] - 0.1 + c["numeric_error"],
                c["prediction"] + 0.1 - c["numeric_error"],
            )
            for c in cells
        ]
        scenarios = []
        for fpr in (0.01, 0.02, 0.03):
            for bias in (0, 0.02, 0.04):
                powers = [acceptance_power(m, fpr, flat_accept)] * 8
                for c, interval in zip(cells, curve_accept, strict=True):
                    powers.append(
                        min(
                            acceptance_power(m, p, interval)
                            for p in (
                                max(0, c["prediction"] - bias),
                                min(1, c["prediction"] + bias),
                            )
                        )
                    )
                scenarios.append(
                    {
                        "FPR": fpr,
                        "absolute_bias": bias,
                        "joint_power_lower_bound": max(
                            0.0, 1 - sum(1 - p for p in powers)
                        ),
                    }
                )
        price = 2 * sum(
            2 * m * math.comb(c["N"], 2) * prices[c["region"]] for c in cells
        )
        results.append({"m": m, "scenarios": scenarios, "cpu_seconds_2x": price})
    return results


def freeze_forecasts():
    stamp = clean_stamp()
    if FREEZE.exists():
        raise FileExistsError(FREEZE)
    grid, grid_hash = read_stage("grid")
    pilot, pilot_hash = read_stage("pilot")
    review = pilot_review(pilot, grid)
    if review["incidents"]:
        raise RuntimeError(f"G2 unresolved: {review['incidents']}")
    selected = {}
    for region in wr.REGIONS:
        rows = [c for c in grid["cells"] if c["region"] == region and c["A"] > 0]
        crossings = [
            (a, b)
            for a, b in zip(rows[:-1], rows[1:], strict=True)
            if (a["N300_prediction"] - 0.5) * (b["N300_prediction"] - 0.5) <= 0
        ]
        if crossings:
            a, b = crossings[0]
            selected[region] = {
                "A": math.sqrt(a["A"] * b["A"]),
                "transition_located": True,
                "bracket": [a["A"], b["A"]],
            }
        else:
            row = min(
                (rows[0], rows[-1]), key=lambda c: abs(c["N300_prediction"] - 0.5)
            )
            selected[region] = {
                "A": row["A"],
                "transition_located": False,
                "warning": "Boundary selected; half-power transition not found.",
            }
    all_levels = []
    for power in (19, 20):
        tasks = [
            (2 * ri + arm, region, 0 if arm == 0 else selected[region]["A"], power)
            for ri, region in enumerate(wr.REGIONS)
            for arm in (0, 1)
        ]
        with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
            rows = list(pool.map(forecast_worker, tasks))
        cells = []
        for ri, region in enumerate(wr.REGIONS):
            flat, curved = rows[2 * ri]["moments"], rows[2 * ri + 1]["moments"]
            if region == "D":
                # Keep the predeclared QMC mean/error; the independent integrals
                # test the instrument and do not manufacture zero uncertainty.
                exact = rows[2 * ri + 1]["independent_volume_integrals"]
                for key in ("p", "zeta"):
                    interval = curved[key + "_interval"]
                    if not interval[0] <= exact[key] <= interval[1]:
                        raise ArithmeticError(
                            "Diamond QMC excludes its independent integral"
                        )
                allowance = 10 * max(
                    c.get("cdf_grid_comparison", {}).get(
                        "max_abs_replicate_difference", 0
                    )
                    for c in grid["cells"]
                )
                curved = dict(curved)
                for key in ("p", "zeta"):
                    lo, hi = curved[key + "_interval"]
                    curved[key + "_interval"] = [max(0, lo - allowance), hi + allowance]
            for n in (300, 1200):
                bounds = threshold(flat, n)
                prediction = predict(curved, n, bounds)
                ci = prediction_interval(curved, n, bounds)
                cells.append(
                    {
                        "region": region,
                        "N": n,
                        "A": selected[region]["A"],
                        "threshold": bounds,
                        "prediction": prediction,
                        "prediction_interval": ci,
                        "numeric_error": max(prediction - ci[0], ci[1] - prediction),
                        "flat_moments": flat,
                        "curved_moments": curved,
                    }
                )
        maximum = max(c["numeric_error"] for c in cells)
        all_levels.append(
            {"power": power, "cells": rows, "maximum_beta_error": maximum}
        )
        print(
            json.dumps({"precision_level": power, "maximum_beta_error": maximum}),
            flush=True,
        )
        if maximum <= 0.01:
            break
    table = power_table(cells, pilot)
    feasible = [
        r
        for r in table
        if min(s["joint_power_lower_bound"] for s in r["scenarios"]) >= 0.9
        and r["cpu_seconds_2x"] <= 48 * 3600
    ]
    development_cpu = (
        grid["cpu_seconds"]
        + pilot["cpu_seconds"]
        + sum(sum(c["cpu_seconds"] for c in level["cells"]) for level in all_levels)
    )
    passed = bool(feasible) and maximum <= 0.01 and development_cpu <= 8 * 3600
    result = {
        "code_version": stamp,
        "G2": "FEASIBLE" if passed else "INFEASIBLE",
        "grid_sha256": grid_hash,
        "pilot_sha256": pilot_hash,
        "pilot_review": review,
        "selected": selected,
        "cells": cells,
        "numerical_levels": all_levels,
        "power_sensitivity": table,
        "m_per_cell": feasible[0]["m"] if passed else None,
        "development_cpu_seconds": development_cpu,
        "confirmation_seed_prefix": [20260921, 6401],
        "family_alpha": 0.05,
        "per_interval_alpha": ALPHA,
        "accuracy_tolerance": 0.1,
        "FPR_upper": 0.05,
        "forecast_method": "Normal tails; U-statistic variance; numerical p,zeta",
        "power_scope": (
            "Conditional on declared FPR<=.03 and forecast bias<=.04 scenarios; "
            "not a guarantee of actual power."
        ),
    }
    FREEZE.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "G2": result["G2"],
                "m": result["m_per_cell"],
                "selected": selected,
                "development_cpu_seconds": development_cpu,
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("freeze",))
    parser.parse_args()
    freeze_forecasts()
