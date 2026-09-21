"""Audit saved confirmation, then produce descriptive comparisons and figures.

No posets are added to confirmation and no primary decision is fitted here.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import window_response as wr
from scipy.stats import binomtest, mannwhitneyu, t
from window_response_g1 import ROOT, clean_stamp
from window_response_run import DOCUMENTS

DEST = DOCUMENTS / "window_response"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((DOCUMENTS / f"window_response_{name}.json").read_text())


def audit(frozen, result, raw):
    """Independent public binomial-test API; do not import the producer scorer."""
    assert (
        digest(DOCUMENTS / "window_response_confirmation_raw.npz")
        == result["raw_sha256"]
    )
    assert (
        digest(DOCUMENTS / "window_response_confirmation_freeze.json")
        == result["freeze_sha256"]
    )
    assert raw.shape == (8, 2, 3200, 3) and np.isfinite(raw).all()
    producer = result["code_version"]
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", producer, "HEAD"], cwd=ROOT, check=True
    )
    committed = subprocess.check_output(
        [
            "git",
            "show",
            f"{producer}:docs/paper/paper_a/window_response_confirmation_freeze.json",
        ],
        cwd=ROOT,
    )
    # Git stores LF; Python's Windows text writer stored CRLF. The recorded raw
    # SHA above checks the actual input bytes; this separately binds its content
    # to the pre-data commit without treating newline conversion as new content.
    working = (DOCUMENTS / "window_response_confirmation_freeze.json").read_bytes()
    assert committed.replace(b"\r\n", b"\n") == working.replace(b"\r\n", b"\n")
    verified = []
    for i, (design, saved) in enumerate(
        zip(frozen["cells"], result["cells"], strict=True)
    ):
        n = design["N"]
        for column, denominator in ((0, math.comb(n, 2)), (1, math.comb(n, 3))):
            counts = raw[i, :, :, column] * denominator
            assert np.all((counts >= 0) & (counts <= denominator))
            np.testing.assert_allclose(counts, counts.round(), atol=1e-5, rtol=0)
        intervals = []
        for arm, label in enumerate(("flat", "curved")):
            values = raw[i, arm, :, 0]
            low, high = design["threshold"]
            k = int(np.count_nonzero((values < low) | (values > high)))
            ci = binomtest(k, len(values)).proportion_ci(
                confidence_level=1 - frozen["per_interval_alpha"], method="exact"
            )
            assert k == saved[label]["hits"]
            np.testing.assert_allclose(
                ci, saved[label]["CP_interval"], atol=2e-12, rtol=0
            )
            assert abs(values.mean() - saved[label]["f_mean"]) < 1e-15
            assert abs(values.std(ddof=1) - saved[label]["f_SD"]) < 1e-15
            intervals.append(ci)
        difference = [
            intervals[1].low - design["prediction"] - design["numeric_error"],
            intervals[1].high - design["prediction"] + design["numeric_error"],
        ]
        g3 = intervals[0].high <= frozen["FPR_upper"]
        g4 = min(difference) >= -0.1 and max(difference) <= 0.1
        assert g3 == saved["G3"] and g4 == saved["G4"]
        assert (g3 and intervals[1].low >= 0.8) == saved["detectable"]
        np.testing.assert_allclose(
            difference, saved["difference_interval"], atol=2e-12, rtol=0
        )
        verified.append({"region": design["region"], "N": n, "G3": g3, "G4": g4})
    assert result["G3_all"] == all(c["G3"] for c in verified)
    assert result["G4_all"] == all(c["G4"] for c in verified)
    assert not np.any(raw[..., 2])
    return {
        "hashes_and_predata_forecast_verified": True,
        "independent_recomputed_cells": verified,
    }


def auc(flat, curved):
    return float(
        mannwhitneyu(curved, flat, alternative="greater").statistic
        / (len(flat) * len(curved))
    )


def secondary(frozen, pilot, raw):
    rows = []
    for i, cell in enumerate(frozen["cells"]):
        n, region = cell["N"], cell["region"]
        nearest = min(
            (c for c in pilot["cells"] if c["region"] == region and c["A"] > 0),
            key=lambda c: abs(math.log(c["A"] / cell["A"])) + abs(math.log(c["N"] / n)),
        )
        baseline = next(
            c
            for c in pilot["cells"]
            if c["region"] == region and c["N"] == nearest["N"] and c["A"] == 0
        )
        directions, scores = [], []
        for j, key in enumerate(("f", "chain3")):
            delta = np.mean([v[key] for v in nearest["raw_statistics"]]) - np.mean(
                [v[key] for v in baseline["raw_statistics"]]
            )
            direction = 1 if delta >= 0 else -1
            directions.append(direction)
            scores.append(auc(direction * raw[i, 0, :, j], direction * raw[i, 1, :, j]))
        transfer = next(
            c for c in frozen["cells"] if c["region"] == "R+" and c["N"] == n
        )
        lo, hi = transfer["threshold"]
        rates = [
            float(np.mean((raw[i, arm, :, 0] < lo) | (raw[i, arm, :, 0] > hi)))
            for arm in (0, 1)
        ]
        # Diagnostic only: one balanced label permutation, never a new gate.
        values = np.abs(raw[i, :, :, 0].reshape(-1) - cell["flat_moments"]["p"])
        permuted = np.random.default_rng([20260921, 6501, i]).permutation(values)
        m = len(values) // 2
        rows.append(
            {
                "region": region,
                "N": n,
                "A": cell["A"],
                "pilot_orientation_source": {"A": nearest["A"], "N": nearest["N"]},
                "directions_f_chain3": directions,
                "AUC_f": scores[0],
                "AUC_chain3": scores[1],
                "delta_AUC_chain3_minus_f": scores[1] - scores[0],
                "Rplus_transfer_FPR_TPR": rates,
                "random_label_AUC": auc(permuted[:m], permuted[m:]),
            }
        )
    return rows


def timing():
    rows = []
    for n in (300, 1200):
        points = wr.sample_points(
            "R0", 1, n, np.random.default_rng([20260921, 6502, n])
        )
        samples = []
        for _ in range(9):
            start = time.perf_counter_ns()
            matrix, _ = wr.order_matrix(points, 1)
            middle = time.perf_counter_ns()
            f = matrix.sum() / math.comb(n, 2)
            second = time.perf_counter_ns()
            chain = np.dot(matrix.sum(axis=0), matrix.sum(axis=1)) / math.comb(n, 3)
            end = time.perf_counter_ns()
            assert 0 <= chain <= f <= 1
            samples.append(
                [(middle - start) / 1e6, (second - middle) / 1e6, (end - second) / 1e6]
            )
        rows.append(
            {
                "N": n,
                "region": "R0",
                "A": 1,
                "median_ms_order_f_chain3": np.median(samples, axis=0).tolist(),
                "repeats": 9,
            }
        )
    return rows


def figures(frozen, results, pilot, grid):
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 180,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 5.8), constrained_layout=True)
    for ax, region in zip(axes.flat, wr.REGIONS, strict=True):
        cells = [c for c in grid["cells"] if c["region"] == region]
        p0 = 0.1 if region == "D" else wr.box_coefficients(region)[0]
        x = [c["A"] for c in cells]
        y = [c["moments"]["p"] - p0 for c in cells]
        ax.plot(
            x,
            y,
            "o-",
            color="#164e78",
            label="Full integral (development)",
            lw=1.4,
            ms=4,
        )
        observations = [
            c for c in pilot["cells"] if c["region"] == region and c["N"] == 800
        ]
        means, errors = [], []
        for c in observations:
            values = np.array([r["f"] for r in c["raw_statistics"]])
            means.append(values.mean() - p0)
            errors.append(t.ppf(0.975, 79) * values.std(ddof=1) / math.sqrt(80))
        ax.errorbar(
            [c["A"] for c in observations],
            means,
            yerr=errors,
            fmt="x",
            color="#b55b13",
            label="Pilot N=800 (pointwise 95% CI)",
            ms=5,
        )
        ax.axhline(0, color="#777777", lw=0.5)
        ax.set(title=region, xlabel="A", ylabel=r"$p_R(A)-p_R(0)$")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, fontsize=8)
    fig.savefig(DEST / "response.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.4), constrained_layout=True)
    positions = np.arange(8)
    obs = np.array([r["curved"]["rate"] for r in results["cells"]]) * 100
    limits = np.array([r["curved"]["CP_interval"] for r in results["cells"]]) * 100
    predicted = np.array([c["prediction"] for c in frozen["cells"]]) * 100
    numerr = np.array([c["numeric_error"] for c in frozen["cells"]]) * 100
    ax.errorbar(
        positions - 0.08,
        predicted,
        yerr=numerr,
        fmt="D",
        color="#164e78",
        label="Frozen prediction + numerical guard",
        ms=5,
    )
    ax.errorbar(
        positions + 0.08,
        obs,
        yerr=np.stack((obs - limits[:, 0], limits[:, 1] - obs)),
        fmt="o",
        color="#b55b13",
        label="Observed + simultaneous sampling intervals",
        ms=5,
        capsize=3,
    )
    ax.set_xticks(positions, [f"{c['region']}\nN={c['N']}" for c in frozen["cells"]])
    ax.set(ylabel="Detection probability (%)", ylim=(-4, 105))
    ax.axhline(80, color="#999999", lw=0.7, ls="--")
    ax.legend(loc="center right", fontsize=8)
    fig.savefig(DEST / "prediction.png")
    plt.close(fig)


def main():
    stamp = clean_stamp()
    frozen, result, pilot, grid = [
        read(n)
        for n in (
            "confirmation_freeze",
            "confirmation_results",
            "pilot_results",
            "grid_results",
        )
    ]
    raw = np.load(
        DOCUMENTS / "window_response_confirmation_raw.npz", allow_pickle=False
    )["values"]
    checks = audit(frozen, result, raw)
    DEST.mkdir(exist_ok=True)
    output = {
        "code_version": stamp,
        "primary_audit": checks,
        "secondary_status": (
            "Descriptive post-confirmation analysis; no primary gate is changed."
        ),
        "secondary": secondary(frozen, pilot, raw),
        "incremental_timing": timing(),
    }
    figures(frozen, result, pilot, grid)
    (DEST / "analysis.json").write_text(
        json.dumps(output, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "primary_audit": "PASS",
                "secondary": output["secondary"],
                "timing": output["incremental_timing"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
