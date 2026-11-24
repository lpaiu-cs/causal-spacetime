"""Summaries and frozen gates for the Route A probe. Reads rows, runs nothing.

    python analyse_route_a.py --stage control
    python analyse_route_a.py --stage pilot      # also writes sizing.json
    python analyse_route_a.py --stage main

All thresholds below are the ones frozen in design.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
from collections import defaultdict

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
Z95, Z90 = 1.959964, 1.281552
BIAS_ALLOWANCE = 0.05
DESIGN_SLOPE_F = 0.40  # combinatorial slope seen while timing, used for sizing
GRID = (64, 128, 256, 512)
COMPLETENESS_FLOOR = 0.90
CONTROL_FLOOR = 0.99
FAR_LAG = 4  # checkpoints; 1000 sweeps
BOOT = 10000


def wilson(successes, trials):
    if trials == 0:
        return None
    p = successes / trials
    denom = 1 + Z95**2 / trials
    centre = (p + Z95**2 / (2 * trials)) / denom
    half = Z95 * math.sqrt(p * (1 - p) / trials + Z95**2 / (4 * trials**2)) / denom
    return [centre - half, centre + half]


def certified(row):
    """Certificates by source: start, chain checkpoints, second start."""
    start = row["single0"] is not None or bool(row.get("joint0"))
    chain = any(hit is not None for hit in row["chain_single"])
    second = bool(row.get("joint1"))
    return start, chain or second


def per_event(rows):
    events = defaultdict(list)
    for row in rows:
        events[(row["n"], row["seed"], row["event"])].append(row)
    out = []
    for (n, seed, event), group in sorted(events.items()):
        states = [[row["single0"] is not None for row in group]] + [
            [row["chain_single"][t] is not None for row in group]
            for t in range(len(group[0]["chain_single"]))
        ]
        walls = np.array(states, dtype=float)
        lag_overlap = {
            lag: float(np.mean(np.sum(walls[:-lag] * walls[lag:], axis=1)))
            for lag in range(1, len(walls))
        }
        start = [certified(row)[0] for row in group]
        other = [certified(row)[1] for row in group]
        out.append(
            {
                "n": n,
                "seed": seed,
                "event": event,
                "F": len(group),
                "S": int(sum(row["single0"] is not None for row in group)),
                "R": int(sum(start)),
                "U": int(sum(a or b for a, b in zip(start, other, strict=True))),
                "independent": int(sum(other)),
                "independent_found_at_start": int(
                    sum(a and b for a, b in zip(start, other, strict=True))
                ),
                "walls_mean": float(walls.sum(axis=1).mean()),
                "lag_overlap": lag_overlap,
            }
        )
    return out


def slope(sizes, values):
    x = np.log(np.asarray(sizes, dtype=float))
    y = np.log(np.asarray(values, dtype=float))
    x = x - x.mean()
    return float(np.sum(x * (y - y.mean())) / np.sum(x * x))


def order_means(events, sizes, key):
    """For each size, the per-order mean of an event-level quantity."""
    table = {n: defaultdict(list) for n in sizes}
    for e in events:
        table[e["n"]][e["seed"]].append(e[key])
    return {
        n: np.array([np.mean(v) for _, v in sorted(table[n].items())]) for n in sizes
    }


def participation(events, n):
    """Mean wall count over the fraction of walls shared at far lags."""
    group = [e for e in events if e["n"] == n]
    walls = np.mean([e["walls_mean"] for e in group])
    lags = sorted(group[0]["lag_overlap"])
    far = [lag for lag in lags if lag >= FAR_LAG]
    by_lag = {
        lag: float(np.mean([e["lag_overlap"][lag] for e in group]) / walls)
        for lag in lags
    }
    if not far:  # the control stage has two checkpoints only
        return None, by_lag
    shared = np.mean([np.mean([e["lag_overlap"][lag] for lag in far]) for e in group])
    return float(walls**2 / shared) if shared > 0 else None, by_lag


def summarise(stage):
    folder = HERE / "evidence" / stage
    rows = [json.loads(line) for line in (folder / "rows.jsonl").open(encoding="utf-8")]
    orders = [
        json.loads(line) for line in (folder / "orders.jsonl").open(encoding="utf-8")
    ]
    manifest = json.loads((folder / "manifest.json").read_text())
    events = per_event(rows)
    sizes = sorted({e["n"] for e in events})
    summary = {
        "stage": stage,
        "code_version": manifest["code_version"],
        "rows_sha256": hashlib.sha256((folder / "rows.jsonl").read_bytes()).hexdigest(),
        "candidates": len(rows),
        "events": len(events),
        "by_size": {},
    }
    for n in sizes:
        group = [e for e in events if e["n"] == n]
        mine = [r for r in rows if r["n"] == n]
        cert = [r for r in mine if any(certified(r))]
        scaled = np.array([n * abs(r["margin"]) for r in mine])
        scaled_cert = np.array([n * abs(r["margin"]) for r in cert])
        joint_only = sum(
            1 for r in mine if r["single0"] is None and bool(r.get("joint0"))
        )
        found = sum(e["independent_found_at_start"] for e in group)
        indep = sum(e["independent"] for e in group)
        k_hat, by_lag = participation(events, n)
        drift = np.array(
            [o["chain_rms_displacement"] for o in orders if o["n"] == n], dtype=float
        )
        summary["by_size"][n] = {
            "orders": len({e["seed"] for e in group}),
            "events": len(group),
            "mean_F": float(np.mean([e["F"] for e in group])),
            "mean_S": float(np.mean([e["S"] for e in group])),
            "mean_R": float(np.mean([e["R"] for e in group])),
            "mean_U": float(np.mean([e["U"] for e in group])),
            "certified_by_joint_only_at_start": joint_only,
            "completeness": {
                "independent_certificates": indep,
                "also_found_at_start": found,
                "wilson95": wilson(found, indep),
            },
            "participation_from_far_lags": k_hat,
            "wall_overlap_fraction_by_lag": by_lag,
            "scaled_margin_quantiles_all": np.quantile(
                scaled, [0.25, 0.5, 0.75, 0.9]
            ).tolist(),
            "scaled_margin_quantiles_certified": (
                np.quantile(scaled_cert, [0.25, 0.5, 0.75, 0.9]).tolist()
                if len(scaled_cert)
                else None
            ),
            "chain_rms_displacement_times_n": (drift.mean(axis=0) * n).tolist(),
            "mean_seconds_per_order": float(
                np.mean([o["seconds"] for o in orders if o["n"] == n])
            ),
        }
    if stage == "control":
        known = [r for r in rows if r["known_1p1"]]
        hit = sum(1 for r in known if certified(r)[0])
        outside = sum(1 for r in rows if not r["known_1p1"] and any(certified(r)))
        interval = wilson(hit, len(known))
        summary["control"] = {
            "known_realisable": len(known),
            "found_at_start": hit,
            "wilson95": interval,
            "certified_outside_the_rule": outside,
            "gate_sensitivity_pass": bool(interval[0] >= CONTROL_FLOOR),
        }
    if stage == "pilot":
        n = sizes[0]
        means = order_means(events, sizes, "U")[n]
        relative = float(means.std(ddof=1) / means.mean())
        sxx = float(np.sum((np.log(GRID) - np.mean(np.log(GRID))) ** 2))
        target_d = (DESIGN_SLOPE_F / 2 - BIAS_ALLOWANCE) / (Z95 + Z90)
        target_log = target_d * math.sqrt(sxx)
        orders_needed = int(min(48, max(12, math.ceil((relative / target_log) ** 2))))
        sizing = {
            "pilot_size": n,
            "pilot_orders": len(means),
            "order_mean_U": means.tolist(),
            "relative_sd_between_orders": relative,
            "target_se_of_D": target_d,
            "target_se_of_log_mean_per_size": target_log,
            "orders_per_size": orders_needed,
            "rule": "ceil((relative_sd / target_se_log)^2), clipped to [12, 48]",
            "pilot_code_version": manifest["code_version"],
        }
        (folder / "sizing.json").write_text(json.dumps(sizing, indent=2) + "\n")
        summary["sizing"] = sizing
    if stage == "main":
        rng = np.random.default_rng(20261124)
        tables = {key: order_means(events, sizes, key) for key in "FSRU"}
        point = {
            key: slope(sizes, [tables[key][n].mean() for n in sizes]) for key in "FSRU"
        }
        draws = {key: [] for key in ("F", "S", "R", "U", "D_U", "D_R")}
        for _ in range(BOOT):
            picks = {
                n: rng.integers(0, len(tables["F"][n]), len(tables["F"][n]))
                for n in sizes
            }
            k = {
                key: slope(sizes, [tables[key][n][picks[n]].mean() for n in sizes])
                for key in "FSRU"
            }
            for key in "FSRU":
                draws[key].append(k[key])
            draws["D_U"].append(k["U"] - k["F"] / 2)
            draws["D_R"].append(k["R"] - k["F"] / 2)
        point["D_U"] = point["U"] - point["F"] / 2
        point["D_R"] = point["R"] - point["F"] / 2
        stats = {
            key: {
                "estimate": point[key],
                "se": float(np.std(draws[key], ddof=1)),
                "ci95": np.quantile(draws[key], [0.025, 0.975]).tolist(),
            }
            for key in draws
        }
        found = sum(e["independent_found_at_start"] for e in events)
        indep = sum(e["independent"] for e in events)
        completeness = wilson(found, indep)
        required_se = (point["F"] / 2 - BIAS_ALLOWANCE) / (Z95 + Z90)
        gates = {
            "F_grows": bool(stats["F"]["ci95"][0] > 0.25),
            "completeness_wilson95": completeness,
            "completeness_pass": bool(
                completeness is not None and completeness[0] >= COMPLETENESS_FLOOR
            ),
            "required_se_of_D": required_se,
            "powered": bool(stats["D_U"]["se"] <= required_se),
        }
        low, high = stats["D_U"]["ci95"]
        if not (gates["F_grows"] and gates["completeness_pass"]):
            verdict = "unresolved: instrument gate failed"
        elif not gates["powered"]:
            verdict = "unresolved: underpowered"
        elif high < 0:
            verdict = "cell-limited realisability"
        elif low > 0:
            verdict = "count-limited realisability"
        else:
            verdict = "undecided"
        summary["slopes"] = stats
        summary["gates"] = gates
        summary["verdict"] = verdict
    (folder / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("control", "pilot", "main"), required=True)
    print(json.dumps(summarise(parser.parse_args().stage), indent=2))
