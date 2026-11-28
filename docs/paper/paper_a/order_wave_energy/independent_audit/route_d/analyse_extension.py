"""Summaries for the N = 512 extension (extension_design.md). No verdict.

Reads the existing records of the first measurement (N = 128, 192, 256) and
the extension records, and writes evidence_extension/summary.json.

    python analyse_extension.py --stage control
    python analyse_extension.py --stage main
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np
from analyse_drift import t_interval
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent
EXTENSION = HERE / "evidence_extension"
FIT_TOLERANCE = 0.05
PLATEAU = 1.10
KEYS = {"full": "central_residual", "wave": "central_residual_wave"}


def read(path):
    return [json.loads(s) for s in path.open(encoding="utf-8")]


def control():
    folder = EXTENSION / "control"
    records = read(folder / "records.jsonl")
    manifest = json.loads((folder / "manifest.json").read_text())
    summary = t_interval([r["ratio"] for r in records])
    summary["n"] = records[0]["n"]
    summary["all_orders_identical"] = all(r["same_order"] for r in records)
    summary["all_fits_succeeded"] = all(r["fit_success"] for r in records)
    low, high = summary["ci95"]
    out = {
        "stage": "control",
        "code_version": manifest["code_version"],
        "fit": summary,
        "fit_pass": bool(
            low >= 1 - FIT_TOLERANCE
            and high <= 1 + FIT_TOLERANCE
            and summary["all_orders_identical"]
            and summary["all_fits_succeeded"]
        ),
        "workers": manifest["workers"],
        "priority_classes": [
            manifest["runner_priority_class"],
            manifest.get("worker_priority_classes"),
        ],
    }
    (folder / "summary.json").write_text(json.dumps(out, indent=2) + "\n", newline="\n")
    return out


def per_order(records, key, lag):
    """sqrt(mean over the two pairs of key / s1^2), one value per order."""
    return np.array(
        [
            math.sqrt(
                np.mean([row[key] for row in r["rows"] if row["lag"] == lag])
                / r["meta"]["s1_squared"]
            )
            for r in records
        ]
    )


def describe(ratio):
    if ratio["ci95"][1] <= PLATEAU:
        return "plateau"
    if ratio["ci95"][0] > PLATEAU:
        return "creep"
    return "neither"


def main_stage():
    old = read(HERE / "evidence" / "main" / "records.jsonl")
    new = read(EXTENSION / "main" / "records.jsonl")
    manifest = json.loads((EXTENSION / "main" / "manifest.json").read_text())
    records = old + new
    sizes = sorted({r["meta"]["n"] for r in records})
    extension_size = new[0]["meta"]["n"]
    out = {
        "stage": "main",
        "code_version": manifest["code_version"],
        "workers": manifest["workers"],
        "priority_classes": [
            manifest["runner_priority_class"],
            manifest.get("worker_priority_classes"),
        ],
        "sizes": sizes,
        "extension_size": extension_size,
        "by_size": {},
        "size_dependence": {},
    }
    points = {(k, lag): [] for k in KEYS for lag in ("16N", "half")}
    scale = []
    for n in sizes:
        mine = [r for r in records if r["meta"]["n"] == n]
        s1 = np.sqrt([r["meta"]["s1_squared"] for r in mine])
        scale += list(s1 * n)
        entry = {
            "orders": len(mine),
            "central_events": float(
                np.mean([r["meta"]["central_events"] for r in mine])
            ),
            "s1_times_n": t_interval(s1 * n),
        }
        for label, key in KEYS.items():
            short = per_order(mine, key, 16 * n)
            half = per_order(mine, key, n * n // 2)
            ratio = t_interval(half / short)
            entry[label] = {
                "R_16N": t_interval(short),
                "R_half_N2": t_interval(half),
                "ratio_half_N2_over_16N": ratio,
                "description": describe(ratio),
                "chart_residual_times_n_16N": t_interval(short * s1 * n),
                "chart_residual_times_n_half_N2": t_interval(half * s1 * n),
            }
            points[(label, "16N")] += [(n, v**2) for v in short]
            points[(label, "half")] += [(n, v**2) for v in half]
        if n == extension_size:
            lags = sorted({row["lag"] for row in mine[0]["rows"]})
            entry["curves"] = {
                str(lag): {
                    "lag_over_n": lag / n,
                    "lag_over_n_squared": lag / n**2,
                    "R": t_interval(per_order(mine, "central_residual", lag)),
                    "R_wave": t_interval(per_order(mine, "central_residual_wave", lag)),
                    "R_outer": t_interval(per_order(mine, "outer_residual", lag)),
                    "total": t_interval(per_order(mine, "central_total", lag)),
                    "outer_total": t_interval(per_order(mine, "outer_total", lag)),
                    "family_relaxed_stiff": t_interval(
                        [
                            np.mean(
                                [
                                    row["family_relaxed_stiff"]
                                    for row in r["rows"]
                                    if row["lag"] == lag
                                ]
                            )
                            for r in mine
                        ]
                    ),
                }
                for lag in lags
            }
            entry["outer_events"] = float(
                np.mean([r["meta"]["outer_events"] for r in mine])
            )
            entry["hessian_eigenvalues_over_n"] = np.mean(
                [r["meta"]["hessian_eigenvalues_over_n"] for r in mine], axis=0
            ).tolist()
            entry["fit_failures"] = int(
                sum(not row["fit_success"] for r in mine for row in r["rows"])
            )
            entry["skipped_moves"] = int(sum(r["meta"]["skipped_moves"] for r in mine))
            entry["seconds_per_order"] = float(
                np.mean([r["meta"]["seconds"] for r in mine])
            )
        out["by_size"][n] = entry
    sigma = float(np.mean(scale))
    out["size_dependence"]["s1_times_n_pooled"] = sigma
    for (label, lag), pairs in points.items():
        x = np.array([p[0] for p in pairs], dtype=float)
        y = np.array([p[1] for p in pairs], dtype=float)
        fit = stats.linregress(x, y)
        t = float(stats.t.ppf(0.975, len(x) - 2))
        low, high = fit.slope - t * fit.stderr, fit.slope + t * fit.stderr
        out["size_dependence"][f"{label}_{lag}"] = {
            "points": len(x),
            "A": float(fit.intercept),
            "B": float(fit.slope),
            "B_ci95": [float(low), float(high)],
            "description": "count-scale term found" if low > 0 else "none found",
            "c": float(sigma * math.sqrt(max(fit.slope, 0.0))),
            "c_upper": float(sigma * math.sqrt(max(high, 0.0))),
            "smallest_c_distinguishable_from_zero": float(
                sigma * math.sqrt(t * fit.stderr)
            ),
        }
    (EXTENSION / "summary.json").write_text(
        json.dumps(out, indent=2) + "\n", newline="\n"
    )
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("control", "main"), required=True)
    chosen = parser.parse_args().stage
    print(json.dumps(control() if chosen == "control" else main_stage(), indent=2))
