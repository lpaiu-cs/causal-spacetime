"""Regenerate Paper A figures from committed experiment summary CSVs.

Reads only the version-controlled CSVs under figures/data/ (copied from the
cited experiment scripts' outputs) so figures are reproducible. Palette:
Okabe-Ito (colorblind-safe; validated). One y-axis per axes; multi-measure
figures use small multiples.

Usage: python docs/paper/paper_a/figures/make_figures.py
"""

from __future__ import annotations

import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = Path("docs/paper/paper_a/figures/data")
OUT = Path("docs/paper/paper_a/figures")
BLUE = "#0072B2"
VERM = "#D55E00"
GREEN = "#009E73"
INK = "#222222"
MUTED = "#888888"
GRID = "#DDDDDD"


def _rows(name: str) -> list[dict]:
    return list(csv.DictReader(open(DATA / name)))


def _f(r: dict, k: str) -> float:
    return float(r[k])


def _style(ax) -> None:
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=INK, labelsize=8)
    ax.grid(True, color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)


# --------------------------------------------------------------------------
# Figure 1: the reconstruction ladder (schematic)
# --------------------------------------------------------------------------
LADDER = [
    ("R0", "order only", "dimension; timelike shape", "no absolute scale"),
    ("R1", "+ measure", "timelike proper time", "needs a density"),
    ("R2", "+ observer", "radar time; unsigned distance", "sign undetermined"),
    ("R3", "+ orientation", "signed coords; Lorentz map", "supplied orientation"),
    ("R4", "+ atlas", "Poincare transition maps", "supplied charts"),
    ("R5", "+ conformal weights", "volume; coarse-grain stability", "factor supplied"),
]


def figure_ladder() -> None:
    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    ax.axis("off")
    n = len(LADDER)
    for i, (rung, ingredient, recovered, bound) in enumerate(LADDER):
        y = n - 1 - i
        ax.add_patch(
            plt.Rectangle((0.0, y + 0.08), 9.7, 0.84, facecolor="#F4F4F2",
                          edgecolor=MUTED, linewidth=0.8, zorder=1)
        )
        ax.text(0.28, y + 0.5, rung, fontsize=12, fontweight="bold", color=BLUE,
                va="center", zorder=2)
        # top line: ingredient -> reconstructed
        ax.text(1.05, y + 0.62, ingredient, fontsize=9.5, color=INK, va="center",
                zorder=2)
        ax.annotate("", xy=(4.5, y + 0.62), xytext=(3.95, y + 0.62),
                    arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.2))
        ax.text(4.7, y + 0.62, recovered, fontsize=9.5, color=INK, va="center",
                fontweight="medium", zorder=2)
        # bottom line: bounded-by note (own line, no collision)
        ax.text(1.05, y + 0.28, f"bounded by: {bound}", fontsize=8, color=MUTED,
                va="center", zorder=2, style="italic")
    ax.text(1.05, n + 0.02, "supplied ingredient  ->  reconstructed quantity",
            fontsize=8.5, color=MUTED)
    ax.set_xlim(-0.1, 9.8)
    ax.set_ylim(0, n + 0.4)
    ax.set_title(
        "The reconstruction ladder: causal order + one supplied ingredient at a "
        "time",
        fontsize=11, color=INK, loc="left",
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig1_ladder.png", dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 2: grounded convergence panel (2 x 2)
# --------------------------------------------------------------------------
def _mean_by(rows, key, value, where=None):
    acc = defaultdict(list)
    for r in rows:
        if where and not where(r):
            continue
        acc[_f(r, key)].append(_f(r, value))
    xs = sorted(acc)
    return xs, [st.mean(acc[x]) for x in xs]


def figure_convergence() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.2))

    n_ticks = [300, 600, 1200, 2400]

    def _n_axis(ax):
        ax.set_xscale("log")
        ax.set_xticks(n_ticks)
        ax.set_xticklabels([str(v) for v in n_ticks], fontsize=8)
        ax.minorticks_off()

    # (a) dimension vs N
    ax = axes[0, 0]
    _style(ax)
    dim_rows = _rows("dimension_reconstruction_summary.csv")
    for d, color in zip((2, 3, 4), (BLUE, VERM, GREEN), strict=True):
        sub = [r for r in dim_rows if int(_f(r, "spacetime_dim")) == d]
        sub.sort(key=lambda r: _f(r, "N"))
        xs = [_f(r, "N") for r in sub]
        ys = [_f(r, "mean_estimated_dim") for r in sub]
        ax.plot(xs, ys, marker="o", ms=4, color=color, lw=1.8)
        ax.axhline(d, color=color, ls=":", lw=0.8, alpha=0.6)
        ax.annotate(f"D = {d}", (xs[-1], ys[-1] + 0.06), fontsize=8, color=color,
                    ha="right")
    _n_axis(ax)
    ax.set_xlabel("N (events)", fontsize=9, color=INK)
    ax.set_ylabel("estimated dimension", fontsize=9, color=INK)
    ax.set_ylim(1.8, 4.25)
    ax.set_title("(a) Myrheim-Meyer dimension (exp10)", fontsize=9.5, color=INK,
                 loc="left")

    # (b) proper-time relative RMSE vs N
    ax = axes[0, 1]
    _style(ax)
    tp = _rows("timelike_pair_reconstruction_summary.csv")
    tp.sort(key=lambda r: _f(r, "N"))
    ns = [_f(r, "N") for r in tp]
    ax.plot(ns, [_f(r, "tau_volume_relative_rmse") for r in tp], marker="o", ms=4,
            color=BLUE, lw=1.8, label="volume estimator")
    ax.plot(ns, [_f(r, "tau_chain_relative_rmse") for r in tp], marker="s", ms=4,
            color=VERM, lw=1.8, label="chain estimator")
    _n_axis(ax)
    ax.set_xlabel("N (events)", fontsize=9, color=INK)
    ax.set_ylabel("proper-time relative RMSE", fontsize=9, color=INK)
    ax.set_title("(b) timelike proper time (exp07)", fontsize=9.5, color=INK,
                 loc="left")
    ax.legend(fontsize=8, frameon=False)

    # (c) radar RMSE vs ticks (mean over N)
    ax = axes[1, 0]
    _style(ax)
    rr = _rows("discrete_radar_reconstruction_summary.csv")
    xt, dist = _mean_by(rr, "tick_count", "radar_distance_rmse")
    _, tim = _mean_by(rr, "tick_count", "radar_time_rmse")
    ax.plot(xt, dist, marker="o", ms=4, color=BLUE, lw=1.8, label="radar distance")
    ax.plot(xt, tim, marker="s", ms=4, color=VERM, lw=1.8, label="radar time")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("observer ticks", fontsize=9, color=INK)
    ax.set_ylabel("RMSE (mean over N)", fontsize=9, color=INK)
    ax.set_title("(c) observer radar (exp11)", fontsize=9.5, color=INK, loc="left")
    ax.legend(fontsize=8, frameon=False)

    # (d) Lorentz beta RMSE vs ticks (mean over N), per beta
    ax = axes[1, 1]
    _style(ax)
    lz = _rows("oriented_radar_lorentz_summary.csv")
    for beta, color in zip((0.3, 0.6), (BLUE, VERM), strict=True):
        xt, br = _mean_by(lz, "tick_count", "fitted_beta_rmse",
                          where=lambda r, b=beta: abs(_f(r, "beta") - b) < 1e-9)
        ax.plot(xt, br, marker="o", ms=4, color=color, lw=1.8, label=f"beta = {beta}")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("observer ticks", fontsize=9, color=INK)
    ax.set_ylabel("fitted-beta RMSE (mean over N)", fontsize=9, color=INK)
    ax.set_title("(d) Lorentz-map recovery (exp13)", fontsize=9.5, color=INK,
                 loc="left")
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle(
        "Reconstruction accuracy improves with sampling / clock resolution "
        "(controlled 1+1D, higher-D for dimension)",
        fontsize=10.5, color=INK, x=0.01, ha="left",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / "fig2_convergence.png", dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 3: measure dependence (R5) -- weighted vs unweighted volume
# --------------------------------------------------------------------------
def figure_measure() -> None:
    rows = [r for r in _rows("weighted_conformal_volume_summary.csv")
            if r["profile"] == "constant_1.5"]
    rows.sort(key=lambda r: _f(r, "N"))
    ns = [_f(r, "N") for r in rows]
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    _style(ax)
    ax.plot(ns, [_f(r, "unweighted_relative_rmse") for r in rows], marker="s", ms=5,
            color=VERM, lw=2, label="unweighted (no measure supplied)")
    ax.plot(ns, [_f(r, "weighted_relative_rmse") for r in rows], marker="o", ms=5,
            color=BLUE, lw=2, label="weighted (measure supplied)")
    ax.set_xscale("log")
    ax.set_xticks([600, 1200, 2400])
    ax.set_xticklabels(["600", "1200", "2400"], fontsize=9)
    ax.minorticks_off()
    ax.set_xlabel("N (events)", fontsize=10, color=INK)
    ax.set_ylabel("volume relative RMSE", fontsize=10, color=INK)
    ax.set_ylim(0, None)
    ax.legend(fontsize=9, frameon=False, loc="center right")
    ax.set_title(
        "R5: volume reconstruction needs a supplied measure\n"
        "(unweighted is biased and does not converge; weighted does)",
        fontsize=10.5, color=INK, loc="left",
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig3_measure.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    figure_ladder()
    figure_convergence()
    figure_measure()
    for name in ("fig1_ladder", "fig2_convergence", "fig3_measure"):
        print("wrote", OUT / f"{name}.png")
