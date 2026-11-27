"""Generate publication-quality vector figures for docs/paper/subspacing_wave.

Reads only committed/saved JSON summaries. Sized to 119 mm text width with
8 pt sans-serif lettering, panel labels a/b, and colorblind-safe shape+stroke encodings.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / "outputs"

BLUE, TEAL, AMBER, GRAY = "#1f6feb", "#0f766e", "#d97706", "#78716c"
INK, AXIS, GRID = "#111827", "#cbd5e1", "#e2e8f0"


def _style(ax, letter, xlabel, ylabel, xlog=True, ylog=True):
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.text(
        -0.12,
        1.04,
        letter,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="bottom",
    )
    ax.set_xlabel(xlabel, fontsize=8, color=INK)
    ax.set_ylabel(ylabel, fontsize=8, color=INK)
    ax.grid(True, which="major", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK, labelsize=8, length=3, width=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)


def make_fig1():
    """Fig. 1: U2 conformal spectral scaling and two-replica phase decoherence."""
    p = OUT / "u2_and_replica_summary.json"
    if not p.exists():
        raise FileNotFoundError(p)
    d = json.loads(p.read_text(encoding="utf-8"))
    plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42})
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(4.685, 2.35), constrained_layout=True
    )

    _style(
        ax_a,
        "a",
        "sprinkled events $N$",
        "gap scaling exponent $s_{\\lambda_1}$",
        xlog=False,
        ylog=False,
    )
    # Plot spectral slopes across schedules for 1+1D vs 3+1D
    scheds = ["fixed", "quarter", "half"]
    labels = ["$\\eta=0.05$", "$\\eta=N^{-1/4}$", "$\\eta=N^{-1/2}$"]
    slope_map = {
        (r["dim"], r["schedule"]): r["log_gap_vs_log_N_slope"]
        for r in d["spectral_slopes"]
    }
    x = np.arange(len(scheds))
    ax_a.plot(
        x,
        [slope_map[(2, s)] for s in scheds],
        marker="o",
        color=BLUE,
        label="1+1D (2 conf.)",
        linewidth=1.5,
    )
    ax_a.plot(
        x,
        [slope_map[(4, s)] for s in scheds],
        marker="s",
        color=AMBER,
        markerfacecolor="white",
        label="3+1D (7 conf.)",
        linewidth=1.5,
    )
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(labels, fontsize=7.5)
    ax_a.set_xlabel("cone smoothing schedule", fontsize=8)
    ax_a.axhline(0.0, color=GRAY, linestyle="--", linewidth=0.7)
    ax_a.legend(fontsize=7.2, frameon=False, loc="upper right")

    _style(
        ax_b,
        "b",
        "scaling exponent $\\alpha = \\log N / \\log r$",
        "phase coherence $|\\langle e^{i\\Delta u/\\lambda}\\rangle|$",
        xlog=False,
        ylog=False,
    )
    rep20 = sorted(
        [r for r in d["replica_aggregates"] if r["N"] == 2**20],
        key=lambda x: x["alpha"],
    )
    alphas = [r["alpha"] for r in rep20]
    cohs = [r["phase_coherence_abs"] for r in rep20]
    jitters = [r["u_error_rms_over_lambda"] for r in rep20]
    ax_b.plot(
        alphas, cohs, marker="o", color=TEAL, linewidth=1.6, label="phase coherence"
    )
    for al, ch, jit in zip(alphas, cohs, jitters, strict=True):
        if abs(al - 2.0) < 0.05:
            xytext, va = (8, 5), "bottom"
        elif abs(al - 2.22) < 0.05:
            xytext, va = (10, -9), "top"
        else:
            xytext, va = (0, 5), "bottom"
        ax_b.annotate(
            f"{jit:.1f}",
            (al, ch),
            textcoords="offset points",
            xytext=xytext,
            ha="center",
            va=va,
            fontsize=6.5,
            color=GRAY,
        )
    ax_b.axvline(2.0, color=GRAY, linestyle=":", linewidth=0.9)
    ax_b.set_ylim(-0.05, 0.98)
    ax_b.set_xlim(1.9, 3.5)
    ax_b.text(2.04, 0.85, "$\\alpha=2$", fontsize=7.2, color=GRAY)

    fig.savefig(HERE / "Fig1.pdf")
    plt.close(fig)


def make_fig2_continuum():
    """Fig. 2 preview (continuum Weyl rolloff prior to finite crossover integration)."""
    p = OUT / "exploratory_step3_crossover.json"
    if not p.exists():
        return
    d = json.loads(p.read_text(encoding="utf-8"))
    cont = d.get("continuum", {})
    if not cont:
        return
    plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42})
    fig, ax = plt.subplots(1, 1, figsize=(4.685, 2.5), constrained_layout=True)
    _style(
        ax,
        "a",
        "dimensionless wavenumber $r = T/\\lambda$",
        "Weyl ratio $\\sqrt{\\widehat E^2}/(\\sqrt{2}C_{\\mathrm{peak}})$",
        xlog=True,
        ylog=False,
    )

    for label, cond, marker, color in (
        (
            "$q=0.01$ fixed energy",
            [
                x
                for x in cont["fixed_q"]
                if x["q"] == 0.01 and x["mode"] == "circular_wave"
            ],
            "o",
            BLUE,
        ),
        (
            "$q=0.25$ fixed energy",
            [
                x
                for x in cont["fixed_q"]
                if x["q"] == 0.25 and x["mode"] == "circular_wave"
            ],
            "s",
            AMBER,
        ),
        ("$C_0=0.5$ fixed crest Weyl", cont["fixed_crest_c0"], "^", TEAL),
    ):
        rs = [x["r"] for x in cond]
        ratios = [x["electric_ratio"] for x in cond]
        ax.plot(
            rs,
            ratios,
            marker=marker,
            color=color,
            linewidth=1.4,
            markersize=4.5,
            label=label,
        )
    ax.axhline(1.0, color=GRAY, linestyle="--", linewidth=0.8)
    ax.set_ylim(-0.05, 1.12)
    ax.legend(fontsize=7.2, frameon=False, loc="upper right")
    fig.savefig(HERE / "Fig2_continuum_preview.pdf")
    plt.close(fig)


if __name__ == "__main__":
    make_fig1()
    make_fig2_continuum()
    print("Generated preliminary figures in", HERE)
