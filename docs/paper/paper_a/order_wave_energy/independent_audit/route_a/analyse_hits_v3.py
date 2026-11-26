"""Summaries and frozen gates for the third Route A design (design_v3.md).

python analyse_hits_v3.py --stage control
python analyse_hits_v3.py --stage main
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
MIXING_SIZES = (256, 512, 1024)
MIXING_CEILING, MIXING_TREND = 1.10, 0.05
CONTROL_TARGET, CONTROL_TOLERANCE = 4.0, 0.1
BOOT = 10000


def per_order(events, key):
    """{n: array over orders of the per-order mean of key}."""
    table = {}
    for e in events:
        table.setdefault(e["n"], {}).setdefault(e["seed"], []).append(e[key])
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


def describe(events, orders, rng):
    """Per-size summary of one stratum, and its per-order tables."""
    out = {}
    tables = {k: per_order(events, k) for k in ("q_window", "q_near", "q_long", "F")}
    for n in sorted(tables["F"]):
        mine = [e for e in events if e["n"] == n]
        w, near, far = (tables[k][n] for k in ("q_window", "q_near", "q_long"))
        picks = rng.integers(0, len(w), (BOOT, len(w)))
        ratio = w[picks].mean(axis=1) / far[picks].mean(axis=1)
        drift = np.array(
            [o["chain_rms_displacement"] for o in orders if o["n"] == n]
        ).mean(axis=0)
        lags = np.mean([e["overlap_by_lag"] for e in mine], axis=0)
        marks = (0, 7, len(drift) - 1)
        out[n] = {
            "orders": len(w),
            "events": len(mine),
            "mean_F": float(tables["F"][n].mean()),
            "K_window": float(1 / w.mean()),
            "K_window_ci95": interval(1 / w[picks].mean(axis=1)),
            "K_near": float(1 / near.mean()),
            "K_long": float(1 / far.mean()),
            "K_long_over_K_window": float(w.mean() / far.mean()),
            "K_long_over_K_window_ci95": interval(ratio),
            "median_event_participation": float(
                np.median([1 / e["q_window"] for e in mine])
            ),
            "between_order_relative_sd": float((1 / w).std(ddof=1) / (1 / w).mean()),
            "K_by_block_lag": {
                str(lag): float(1 / lags[lag - 1])
                for lag in (1, 2, 4, 7, 15, 31, 63)
                if lag <= len(lags)
            },
            "mean_partners_hit": float(np.mean([e["partners_hit"] for e in mine])),
            "tip_hit_fraction": float(np.mean([e["tip_hit_fraction"] for e in mine])),
            "chain_rms_displacement_times_n": [float(drift[i] * n) for i in marks],
            "chain_rms_displacement_times_sqrt_n": [
                float(drift[i] * np.sqrt(n)) for i in marks
            ],
        }
    return out, tables


def summarise(stage):
    rng = np.random.default_rng(20261126)
    folder = HERE / "evidence_v3" / stage
    events = [json.loads(s) for s in (folder / "events.jsonl").open(encoding="utf-8")]
    orders = [json.loads(s) for s in (folder / "orders.jsonl").open(encoding="utf-8")]
    manifest = json.loads((folder / "manifest.json").read_text())
    summary = {
        "stage": stage,
        "code_version": manifest["code_version"],
        "workers": manifest["workers"],
        "runner_priority_class": manifest["runner_priority_class"],
        "worker_priority_classes": manifest.get("worker_priority_classes"),
        "hits_outside_candidates": int(sum(o["outside_candidates"] for o in orders)),
        "skipped_moves": int(sum(o["skipped_moves"] for o in orders)),
        "seconds_per_order": {
            str(n): float(np.mean([o["seconds"] for o in orders if o["n"] == n]))
            for n in sorted({o["n"] for o in orders})
        },
    }
    if stage == "control":
        classes = {}
        for klass in ("four", "shared", "edge"):
            mine = [e for e in events if e["klass"] == klass]
            table = per_order(mine, "q_window")
            pooled = np.concatenate([table[n] for n in sorted(table)])
            picks = rng.integers(0, len(pooled), (BOOT, len(pooled)))
            classes[klass] = {
                "events": len(mine),
                "participation": float(1 / np.mean([e["q_window"] for e in mine])),
                "ci95": interval(1 / pooled[picks].mean(axis=1)),
            }
        low, high = classes["four"]["ci95"]
        summary["classes"] = classes
        summary["gate_pass"] = bool(
            low >= CONTROL_TARGET - CONTROL_TOLERANCE
            and high <= CONTROL_TARGET + CONTROL_TOLERANCE
            and summary["hits_outside_candidates"] == 0
        )
    else:
        strata = {
            "central": [e for e in events if e["central"]],
            "other_extremal": [e for e in events if not e["central"] and e["extremal"]],
            "other": [e for e in events if not e["central"] and not e["extremal"]],
        }
        summary["strata"] = {}
        for label, group in strata.items():
            summary["strata"][label], tables = describe(group, orders, rng)
            if label == "central":
                central = tables
        present = [n for n in GATE_SIZES if n in central["F"]]
        w, f = central["q_window"], central["F"]
        draws = {"K": [], "F": [], "D": []}
        for _ in range(BOOT):
            pick = {n: rng.integers(0, len(w[n]), len(w[n])) for n in present}
            k = slope(present, [1 / w[n][pick[n]].mean() for n in present])
            g = slope(present, [f[n][pick[n]].mean() for n in present])
            draws["K"].append(k)
            draws["F"].append(g)
            draws["D"].append(k - g / 2)
        point = {
            "K": slope(present, [1 / w[n].mean() for n in present]),
            "F": slope(present, [f[n].mean() for n in present]),
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
        table = summary["strata"]["central"]
        ratios = [table[n]["K_long_over_K_window"] for n in MIXING_SIZES]
        uppers = [table[n]["K_long_over_K_window_ci95"][1] for n in MIXING_SIZES]
        trend = max(ratios[1] - ratios[0], ratios[2] - ratios[0], ratios[2] - ratios[1])
        required = (point["F"] / 2 - BIAS_ALLOWANCE) / (Z95 + Z90)
        gates = {
            "complete_grid": present == list(GATE_SIZES),
            "no_hits_outside_candidates": summary["hits_outside_candidates"] == 0,
            "F_grows": bool(stats["F"]["ci95"][0] > 0.25),
            "mixing_upper_ends": uppers,
            "mixing_trend": float(trend),
            "mixed": bool(max(uppers) <= MIXING_CEILING and trend <= MIXING_TREND),
            "required_se_of_D": required,
            "powered": bool(stats["D"]["se"] <= required),
        }
        low, high = stats["D"]["ci95"]
        if not (
            gates["complete_grid"]
            and gates["no_hits_outside_candidates"]
            and gates["F_grows"]
        ):
            verdict = "unresolved: instrument gate failed"
        elif not gates["mixed"]:
            verdict = "unresolved: mixing not shown"
        elif not gates["powered"]:
            verdict = "unresolved: underpowered"
        elif low > 0:
            verdict = "count-limited participation"
        elif high < 0:
            verdict = "no growth seen by this instrument"
        else:
            verdict = "undecided"
        summary.update(slopes=stats, gates=gates, verdict=verdict)
    (folder / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", newline="\n"
    )
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("control", "main"), required=True)
    print(json.dumps(summarise(parser.parse_args().stage), indent=2))
