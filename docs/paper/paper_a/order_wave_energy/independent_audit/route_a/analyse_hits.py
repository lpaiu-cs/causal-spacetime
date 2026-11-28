"""Summaries and frozen gates for the wall-hit instrument (design_v2.md).

python analyse_hits.py --stage control
python analyse_hits.py --stage main        # main_a, plus main_b if present
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
Z95, Z90 = 1.959964, 1.281552
BIAS_ALLOWANCE = 0.05
GATE_SIZES = (256, 512, 1024, 2048)
NEAR, FAR = slice(0, 2), slice(3, 7)  # block lags 1-2 and 4-7
MIXING_CEILING = 1.10
CONTROL_TARGET, CONTROL_TOLERANCE = 4.0, 0.3
BOOT = 10000


def load(stage):
    folder = HERE / "evidence_v2" / stage
    if not (folder / "events.jsonl").exists():
        return [], [], None
    events = [json.loads(s) for s in (folder / "events.jsonl").open(encoding="utf-8")]
    orders = [json.loads(s) for s in (folder / "orders.jsonl").open(encoding="utf-8")]
    return events, orders, json.loads((folder / "manifest.json").read_text())


def order_table(events, key):
    """{n: array over orders of the per-order mean of key(event)}."""
    table = {}
    for e in events:
        table.setdefault(e["n"], {}).setdefault(e["seed"], []).append(key(e))
    return {
        n: np.array([np.mean(v) for _, v in sorted(seeds.items())])
        for n, seeds in sorted(table.items())
    }


def slope(sizes, values):
    x = np.log(np.asarray(sizes, dtype=float))
    y = np.log(np.asarray(values, dtype=float))
    x = x - x.mean()
    return float(np.sum(x * (y - y.mean())) / np.sum(x * x))


def interval(draws):
    return np.quantile(draws, [0.025, 0.975]).tolist()


def summarise(stage):
    rng = np.random.default_rng(20261125)
    if stage == "control":
        events, orders, manifest = load("control")
        stamps = {"control": manifest["code_version"]}
    else:
        events, orders, stamps = [], [], {}
        for part in ("main_a", "main_b"):
            e, o, manifest = load(part)
            events += e
            orders += o
            if manifest is not None:
                stamps[part] = manifest["code_version"]
    q_far = order_table(events, lambda e: float(np.mean(e["overlap_by_lag"][FAR])))
    q_near = order_table(events, lambda e: float(np.mean(e["overlap_by_lag"][NEAR])))
    f_mean = order_table(events, lambda e: e["F"])
    sizes = sorted(q_far)
    summary = {"stage": stage, "code_versions": stamps, "by_size": {}}
    summary["hits_outside_candidates"] = float(
        sum(e["hits_outside_candidates"] for e in events)
    )
    for n in sizes:
        mine = [e for e in events if e["n"] == n]
        picks = rng.integers(0, len(q_far[n]), (BOOT, len(q_far[n])))
        ratio = q_near[n][picks].mean(axis=1) / q_far[n][picks].mean(axis=1)
        drift = np.array(
            [o["chain_rms_displacement"] for o in orders if o["n"] == n]
        ).mean(axis=0)
        summary["by_size"][n] = {
            "orders": len(q_far[n]),
            "events": len(mine),
            "mean_F": float(f_mean[n].mean()),
            "participation_far": float(1 / q_far[n].mean()),
            "participation_far_ci95": interval(1 / q_far[n][picks].mean(axis=1)),
            "participation_near": float(1 / q_near[n].mean()),
            "participation_by_lag": [
                float(1 / np.mean([e["overlap_by_lag"][lag] for e in mine]))
                for lag in range(7)
            ],
            "mixing_ratio": float(q_near[n].mean() / q_far[n].mean()),
            "mixing_ratio_ci95": interval(ratio),
            "between_order_relative_sd": float(
                (1 / q_far[n]).std(ddof=1) / (1 / q_far[n]).mean()
            ),
            "mean_partners_hit": float(np.mean([e["partners_hit"] for e in mine])),
            "tip_hit_fraction": float(np.mean([e["tip_hit_fraction"] for e in mine])),
            "hits_per_block": float(np.mean([e["hits_per_block"] for e in mine])),
            "chain_rms_displacement_times_n": (drift * n).tolist(),
            "chain_rms_displacement_times_sqrt_n": (drift * np.sqrt(n)).tolist(),
            "mean_seconds_per_order": float(
                np.mean([o["seconds"] for o in orders if o["n"] == n])
            ),
        }
    if stage == "control":
        inner = [e for e in events if e["interior"]]
        table = order_table(inner, lambda e: float(np.mean(e["overlap_by_lag"][FAR])))
        pooled = np.concatenate([table[n] for n in sorted(table)])
        picks = rng.integers(0, len(pooled), (BOOT, len(pooled)))
        ci = interval(1 / pooled[picks].mean(axis=1))
        summary["control"] = {
            "interior_events": len(inner),
            "participation": float(1 / pooled.mean()),
            "ci95": ci,
            "by_size": {n: float(1 / table[n].mean()) for n in sorted(table)},
            "gate_pass": bool(
                ci[0] >= CONTROL_TARGET - CONTROL_TOLERANCE
                and ci[1] <= CONTROL_TARGET + CONTROL_TOLERANCE
                and summary["hits_outside_candidates"] == 0
            ),
        }
    else:
        gate = [n for n in GATE_SIZES if n in sizes]
        summary["gate_sizes_present"] = gate
        if len(gate) >= 3:
            draws = {"K": [], "F": [], "D": []}
            for _ in range(BOOT):
                pick = {n: rng.integers(0, len(q_far[n]), len(q_far[n])) for n in gate}
                k = slope(gate, [1 / q_far[n][pick[n]].mean() for n in gate])
                f = slope(gate, [f_mean[n][pick[n]].mean() for n in gate])
                draws["K"].append(k)
                draws["F"].append(f)
                draws["D"].append(k - f / 2)
            point = {
                "K": slope(gate, [1 / q_far[n].mean() for n in gate]),
                "F": slope(gate, [f_mean[n].mean() for n in gate]),
            }
            point["D"] = point["K"] - point["F"] / 2
            stats = {
                key: {
                    "estimate": point[key],
                    "se": float(np.std(draws[key], ddof=1)),
                    "ci95": interval(draws[key]),
                }
                for key in draws
            }
            required = (point["F"] / 2 - BIAS_ALLOWANCE) / (Z95 + Z90)
            mixed = all(
                summary["by_size"][n]["mixing_ratio_ci95"][1] <= MIXING_CEILING
                for n in gate
            )
            gates = {
                "complete_grid": gate == list(GATE_SIZES),
                "no_hits_outside_candidates": summary["hits_outside_candidates"] == 0,
                "F_grows": bool(stats["F"]["ci95"][0] > 0.25),
                "required_se_of_D": required,
                "powered": bool(stats["D"]["se"] <= required),
                "mixed": bool(mixed),
            }
            low, high = stats["D"]["ci95"]
            if not (gates["no_hits_outside_candidates"] and gates["F_grows"]):
                verdict = "unresolved: instrument gate failed"
            elif not gates["powered"]:
                verdict = "unresolved: underpowered"
            elif low > 0:
                verdict = "count-limited participation"
            elif high < 0:
                verdict = (
                    "cell-limited participation"
                    if mixed
                    else "unresolved: flat but not shown mixed"
                )
            else:
                verdict = "undecided"
            if not gates["complete_grid"]:
                verdict = "provisional (grid incomplete): " + verdict
            summary.update(slopes=stats, gates=gates, verdict=verdict)
    name = "summary_control.json" if stage == "control" else "summary_main.json"
    (HERE / "evidence_v2" / name).write_text(
        json.dumps(summary, indent=2) + "\n", newline="\n"
    )
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("control", "main"), required=True)
    print(json.dumps(summarise(parser.parse_args().stage), indent=2))
