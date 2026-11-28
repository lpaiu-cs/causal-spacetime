"""Fig. 1 of the Weyl note: benchmark statistic and tidal error against event count.

Reads only the committed run summary. Sized for a 119 mm text width with 8 pt
lettering and no title inside the artwork. Wave cases use filled symbols and
flat controls open ones; circles and squares mark the independent cases, so
nothing depends on color alone.
Usage: make_figure.py [preview.png]
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
SUMMARY = HERE.parent / "paper_a" / "intrinsic_weyl" / "bombelli_benchmark_summary.json"
BLUE, LIGHT, GRAY = "#2a78d6", "#86b6ef", "#898781"
INK, SECONDARY, AXIS, GRID = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9"
MARKERS = {0: "o", 1: "s", 2: "o", 3: "s"}
MARK = {"markersize": 5.5, "markeredgewidth": 1.1, "linewidth": 1.5}
DASH = (0, (4, 3))


def style(ax, letter, ylabel, levels):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.text(-0.13, 1.0, letter, transform=ax.transAxes, fontsize=10, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=8, color=INK)
    ax.set_xticks(levels)
    ax.set_xticklabels([f"{n / 1e6:.3g}" for n in levels])
    ax.minorticks_off()
    ax.set_xlim(levels[0] / 1.5, levels[-1] * 1.5)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=SECONDARY, labelsize=8, length=3, width=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)


def legend(ax, location):
    ax.legend(
        fontsize=8,
        loc=location,
        ncol=2,
        handlelength=1.8,
        columnspacing=1.2,
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=1,
        labelcolor=INK,
    )


def main(preview=None):
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    levels = summary["measurement_levels"]
    rows = {(r["case"], r["measurement_n"]): r for r in summary["results"]}
    plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42})
    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(4.685, 4.9), sharex=True, constrained_layout=True
    )

    style(top, "a", "detection statistic $Q$", levels)
    for case, kind in enumerate(c["kind"] for c in summary["cases"]):
        wave = kind == "wave"
        top.plot(
            levels,
            [rows[case, n]["detection"]["q"] for n in levels],
            marker=MARKERS[case],
            color=BLUE if wave else GRAY,
            markerfacecolor=BLUE if wave else "white",
            markeredgecolor="white" if wave else GRAY,
            label=f"{'wave' if wave else 'flat control'}, case {case}",
            **MARK,
        )
    limit = summary["threshold_chi2_10"]
    top.axhline(limit, color=INK, linewidth=0.8, linestyle=DASH)
    top.text(
        levels[-1] * 1.4,
        limit * 1.12,
        "threshold $q^*$",
        fontsize=8,
        color=INK,
        ha="right",
    )
    top.set_ylim(3, 900)
    top.set_yticks([3, 10, 30, 100, 300])
    top.set_yticklabels(["3", "10", "30", "100", "300"])
    legend(top, "upper left")

    style(bottom, "b", "relative to crest tidal norm", levels)
    for case in (0, 1):
        scale = rows[case, levels[0]]["evaluation"]["truth_norm"]
        values = [rows[case, n]["evaluation"] for n in levels]
        bottom.plot(
            levels,
            [v["confidence_radius_95"] / scale for v in values],
            marker=MARKERS[case],
            color=LIGHT,
            markerfacecolor="white",
            markeredgecolor=LIGHT,
            label=f"95% radius, case {case}",
            **MARK,
        )
        bottom.plot(
            levels,
            [v["full_tidal_error_mod_reflection"] / scale for v in values],
            marker=MARKERS[case],
            color=BLUE,
            markeredgecolor="white",
            label=f"error, case {case}",
            **MARK,
        )
    bias = rows[0, levels[-1]]["evaluation"]["relative_continuum_bias"]
    bottom.axhline(bias, color=INK, linewidth=0.8, linestyle=DASH)
    bottom.text(
        levels[0] / 1.4,
        bias * 0.78,
        "continuum bias of the model",
        fontsize=8,
        color=INK,
    )
    bottom.set_ylim(0.07, 9)
    bottom.set_yticks([0.1, 0.3, 1, 3])
    bottom.set_yticklabels(["10%", "30%", "100%", "300%"])
    bottom.set_xlabel("measurement events $N$ (millions)", fontsize=8, color=INK)
    legend(bottom, "upper right")

    fig.savefig(HERE / "Fig1.pdf", metadata={"CreationDate": None, "ModDate": None})
    if preview:
        fig.savefig(preview, dpi=220)


if __name__ == "__main__":
    main(*sys.argv[1:2])
