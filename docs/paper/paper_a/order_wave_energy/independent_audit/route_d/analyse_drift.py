"""Summaries for the first Route D measurement (design.md). No verdict.

python analyse_drift.py --stage control
python analyse_drift.py --stage main
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent
FIT_TOLERANCE, CHAIN_TOLERANCE = 0.05, 0.2
CHAIN_REQUIRED_SE = 0.2 / 3.2415
PLATEAU = 1.10


def t_interval(values):
    """Mean, standard error and two-sided 95% t interval over orders."""
    v = np.asarray(values, dtype=float)
    mean, se = float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v)))
    half = float(stats.t.ppf(0.975, len(v) - 1)) * se
    return {
        "mean": mean,
        "se": se,
        "ci95": [mean - half, mean + half],
        "orders": len(v),
    }


def load(stage):
    folder = HERE / "evidence" / stage
    records = [json.loads(s) for s in (folder / "records.jsonl").open(encoding="utf-8")]
    return folder, records, json.loads((folder / "manifest.json").read_text())


def control():
    folder, records, manifest = load("control")
    fits = [r for r in records if r["file"].startswith("fit_")]
    chains = [r for r in records if r["file"].startswith("chain_")]
    out = {"stage": "control", "code_version": manifest["code_version"], "fit": {}}
    fit_pass = True
    for n in sorted({r["n"] for r in fits}):
        mine = [r for r in fits if r["n"] == n]
        summary = t_interval([r["ratio"] for r in mine])
        summary["all_orders_identical"] = all(r["same_order"] for r in mine)
        summary["all_fits_succeeded"] = all(r["fit_success"] for r in mine)
        low, high = summary["ci95"]
        fit_pass &= (
            low >= 1 - FIT_TOLERANCE
            and high <= 1 + FIT_TOLERANCE
            and summary["all_orders_identical"]
            and summary["all_fits_succeeded"]
        )
        out["fit"][n] = summary
    lags = sorted(int(k) for k in chains[0]["ratio_by_lag"])
    n = chains[0]["n"]
    out["chain"] = {
        "n": n,
        "by_lag_over_n_squared": {
            str(lag / n**2): t_interval(
                [np.mean(r["ratio_by_lag"][str(lag)]) for r in chains]
            )
            for lag in lags
        },
        "skipped_moves": int(sum(r["skipped_moves"] for r in chains)),
    }
    gated = out["chain"]["by_lag_over_n_squared"]["4.0"]
    chain_pass = (
        gated["ci95"][0] >= 1 - CHAIN_TOLERANCE
        and gated["ci95"][1] <= 1 + CHAIN_TOLERANCE
    )
    out["chain"]["powered"] = bool(gated["se"] <= CHAIN_REQUIRED_SE)
    out["fit_pass"] = bool(fit_pass)
    out["chain_pass"] = bool(chain_pass and out["chain"]["powered"])
    out["controls_pass"] = bool(out["fit_pass"] and out["chain_pass"])
    out["workers"] = manifest["workers"]
    out["priority_classes"] = [
        manifest["runner_priority_class"],
        manifest.get("worker_priority_classes"),
    ]
    (folder / "summary.json").write_text(json.dumps(out, indent=2) + "\n", newline="\n")
    return out


def main_stage():
    folder, records, manifest = load("main")
    out = {"stage": "main", "code_version": manifest["code_version"], "by_size": {}}
    out["workers"] = manifest["workers"]
    out["priority_classes"] = [
        manifest["runner_priority_class"],
        manifest.get("worker_priority_classes"),
    ]
    points = {"16N": [], "half": []}
    scale = []
    for n in sorted({r["meta"]["n"] for r in records}):
        mine = [r for r in records if r["meta"]["n"] == n]
        s1 = np.array([r["meta"]["s1_squared"] for r in mine])
        lags = sorted({row["lag"] for row in mine[0]["rows"]})

        def per_order(key, lag, mine=mine, s1=s1):
            """sqrt(mean over the two pairs of key / s1^2), one value per order."""
            return np.array(
                [
                    math.sqrt(
                        np.mean([row[key] for row in r["rows"] if row["lag"] == lag])
                        / s1[i]
                    )
                    for i, r in enumerate(mine)
                ]
            )

        curves = {}
        for lag in lags:
            near_far = []
            for r in mine:
                pairs = [
                    row["pair_difference_bins"]
                    for row in r["rows"]
                    if row["lag"] == lag
                ]
                near = np.mean(
                    [p[0] for p in pairs if p[0] is not None and p[2] is not None]
                )
                far = np.mean(
                    [p[2] for p in pairs if p[0] is not None and p[2] is not None]
                )
                near_far.append(math.sqrt(far / near))
            curves[str(lag)] = {
                "lag_over_n": lag / n,
                "lag_over_n_squared": lag / n**2,
                "R": t_interval(per_order("central_residual", lag)),
                "R_wave": t_interval(per_order("central_residual_wave", lag)),
                "R_outer": t_interval(per_order("outer_residual", lag)),
                "total": t_interval(per_order("central_total", lag)),
                "outer_total": t_interval(per_order("outer_total", lag)),
                "far_over_near_pair_difference": t_interval(
                    [v for v in near_far if math.isfinite(v)]
                ),
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
        ratio = t_interval(
            per_order("central_residual", 2 * n * n)
            / per_order("central_residual", 16 * n)
        )
        if ratio["ci95"][1] <= PLATEAU:
            description = "plateau"
        elif ratio["ci95"][0] > PLATEAU:
            description = "creep"
        else:
            description = "neither"
        for label, lag in (("16N", 16 * n), ("half", n * n // 2)):
            points[label] += [
                (n, value**2) for value in per_order("central_residual", lag)
            ]
        scale += list(np.sqrt(s1) * n)
        out["by_size"][n] = {
            "orders": len(mine),
            "central_events": float(
                np.mean([r["meta"]["central_events"] for r in mine])
            ),
            "outer_events": float(np.mean([r["meta"]["outer_events"] for r in mine])),
            "s1_times_n": t_interval(np.sqrt(s1) * n),
            "hessian_eigenvalues_over_n": np.mean(
                [r["meta"]["hessian_eigenvalues_over_n"] for r in mine], axis=0
            ).tolist(),
            "fit_failures": int(
                sum(not row["fit_success"] for r in mine for row in r["rows"])
            ),
            "skipped_moves": int(sum(r["meta"]["skipped_moves"] for r in mine)),
            "seconds_per_order": float(np.mean([r["meta"]["seconds"] for r in mine])),
            "R_2N2_over_R_16N": ratio,
            "description": description,
            "curves": curves,
        }
    sigma = float(np.mean(scale))
    out["size_dependence"] = {"s1_times_n_pooled": sigma}
    for label, pairs in points.items():
        x = np.array([p[0] for p in pairs], dtype=float)
        y = np.array([p[1] for p in pairs], dtype=float)
        fit = stats.linregress(x, y)
        t = float(stats.t.ppf(0.975, len(x) - 2))
        b_low, b_high = fit.slope - t * fit.stderr, fit.slope + t * fit.stderr
        out["size_dependence"][label] = {
            "A": float(fit.intercept),
            "A_ci95": [
                float(fit.intercept - t * fit.intercept_stderr),
                float(fit.intercept + t * fit.intercept_stderr),
            ],
            "B": float(fit.slope),
            "B_ci95": [float(b_low), float(b_high)],
            "c": float(sigma * math.sqrt(max(fit.slope, 0.0))),
            "c_ci95": [
                float(sigma * math.sqrt(max(b_low, 0.0))),
                float(sigma * math.sqrt(max(b_high, 0.0))),
            ],
            "smallest_c_distinguishable_from_zero": float(
                sigma * math.sqrt(t * fit.stderr)
            ),
        }
    (folder / "summary.json").write_text(json.dumps(out, indent=2) + "\n", newline="\n")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("control", "main"), required=True)
    chosen = parser.parse_args().stage
    print(json.dumps(control() if chosen == "control" else main_stage(), indent=2))
