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

S1_RES = ROOT / "docs/prereg/frozen/subspacing_confirmatory/results.json"
S2_RES = ROOT / "docs/prereg/frozen/subspacing_s2/results.json"

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
        "gap scaling exponent $s_{\\lambda_1^\\perp}$",
        xlog=False,
        ylog=False,
    )
    # Plot spectral slopes across schedules for 1+1D vs 3+1D
    scheds = ["fixed", "quarter", "half"]
    labels = ["$\\eta=0.05$", "$\\eta=0.5N^{-1/4}$", "$\\eta=N^{-1/2}$"]
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
        "phase coherence $|\\langle e^{i(\\hat{u}-u)/\\lambda}\\rangle|$",
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


def make_fig2():
    """Fig. 2: Continuum Weyl rolloff and finite-N Weyl-to-spectral crossover."""
    p = OUT / "exploratory_step3_crossover.json"
    if not p.exists():
        raise FileNotFoundError(p)
    if not S2_RES.exists():
        raise FileNotFoundError(S2_RES)
    d = json.loads(p.read_text(encoding="utf-8"))
    s2 = json.loads(S2_RES.read_text(encoding="utf-8"))
    cont = d["continuum"]
    fin = d["finite"]
    arms = s2["summary"]["arm_summaries"]

    plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42})
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(4.685, 2.55), constrained_layout=True
    )

    # Panel a: continuum fifth-degree Weyl electric ratio
    _style(
        ax_a,
        "a",
        "dimensionless wavenumber $r = T/\\lambda$",
        "Weyl ratio $\\sqrt{\\widehat E^2}/(\\sqrt{2}C_{\\mathrm{peak}})$",
        xlog=True,
        ylog=False,
    )
    for label, cond, marker, mfc, color in (
        (
            "$q=0.01$ fixed energy",
            [
                x
                for x in cont["fixed_q"]
                if x["q"] == 0.01 and x["mode"] == "circular_wave"
            ],
            "o",
            BLUE,
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
            "white",
            AMBER,
        ),
        ("$C_0=0.5$ fixed crest", cont["fixed_crest_c0"], "^", TEAL, TEAL),
    ):
        rs = [x["r"] for x in cond]
        ratios = [x["electric_ratio"] for x in cond]
        ax_a.plot(
            rs,
            ratios,
            marker=marker,
            color=color,
            markerfacecolor=mfc,
            linewidth=1.35,
            markersize=4.0,
            label=label,
        )
    ax_a.axvspan(16.0, 48.0, color="#f1f5f9", zorder=0)
    ax_a.axhline(1.0, color=GRAY, linestyle="--", linewidth=0.75)
    ax_a.set_ylim(-0.05, 1.12)
    ax_a.legend(fontsize=6.5, frameon=False, loc="lower left")

    # Panel b: finite-N Weyl jet collapse vs. blind spectral turn-on
    _style(
        ax_b,
        "b",
        "dimensionless wavenumber $r = T/\\lambda$",
        "Weyl ratio (5th-deg. jet)",
        xlog=True,
        ylog=False,
    )
    ax_b.axvspan(16.0, 48.0, color="#f1f5f9", zorder=0)

    w_fin = sorted(
        [x for x in fin["weyl_finite"] if x["mode"] == "circular_wave"],
        key=lambda x: x["r"],
    )
    ax_b.plot(
        [x["r"] for x in w_fin],
        [x["electric_ratio"] for x in w_fin],
        marker="s",
        markerfacecolor="white",
        color=AMBER,
        linestyle="--",
        linewidth=1.15,
        markersize=3.8,
        label="Weyl ($2^{21}$)",
    )
    for r_w, arm_w in ((4.0, "X_weyl_r4"), (32.0, "X_weyl_r32")):
        er = arms[arm_w]["electric_ratio"]
        yerr = [[er["mean"] - er["ci95_lo"]], [er["ci95_hi"] - er["mean"]]]
        ax_b.errorbar(
            [r_w],
            [er["mean"]],
            yerr=yerr,
            fmt="D",
            color=AMBER,
            markersize=4.3,
            capsize=2.2,
            linewidth=1.2,
            label="Weyl S2 ($2^{21}$)" if r_w == 4.0 else None,
        )
    ax_b.set_ylim(-0.06, 1.42)
    ax_b.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_b.yaxis.label.set_color(AMBER)
    ax_b.tick_params(axis="y", colors=AMBER, labelsize=7.5)

    ax_r = ax_b.twinx()
    ax_r.set_ylabel("spectral peak $Z_{\\max}$", fontsize=8, color=TEAL)
    ax_r.spines["top"].set_visible(False)
    ax_r.spines["left"].set_visible(False)
    ax_r.spines["bottom"].set_visible(False)
    ax_r.spines["right"].set_color(TEAL)
    ax_r.spines["right"].set_linewidth(0.8)
    ax_r.tick_params(axis="y", colors=TEAL, labelsize=7.5, length=3, width=0.6)

    sub_q1 = sorted(
        [
            x
            for x in fin["subspacing_finite"]
            if x["q"] == 1.0 and x["mode"] == "circular_wave"
        ],
        key=lambda x: x["r"],
    )
    ax_r.plot(
        [x["r"] for x in sub_q1],
        [x["wide_peak_z"] for x in sub_q1],
        marker="o",
        markerfacecolor="white",
        color="#64748b",
        linestyle=":",
        linewidth=1.0,
        markersize=3.2,
        label="scout ($q=1$)",
    )
    for r_s, arm_s, marker, col, lbl in (
        (32.0, "X_sub_q4_r32", "^", TEAL, "S2 ($q=4$)"),
        (64.0, "W1", "^", TEAL, None),
        (96.0, "W2", "^", TEAL, None),
        (128.0, "W3", "^", TEAL, None),
        (48.0, "X_sub_q1_r48", "o", BLUE, "S2 ($q=1$)"),
    ):
        zm = arms[arm_s]["z_max"]
        yerr = [[zm["mean"] - zm["ci95_lo"]], [zm["ci95_hi"] - zm["mean"]]]
        ax_r.errorbar(
            [r_s],
            [zm["mean"]],
            yerr=yerr,
            fmt=marker,
            color=col,
            markersize=4.5,
            capsize=2.2,
            linewidth=1.3,
            label=lbl,
        )
    ax_r.axhline(2.50, color=BLUE, linestyle="--", linewidth=0.8)
    ax_r.text(1.12, 2.58, "$Z_{\\mathrm{th}}=2.50$", fontsize=6.4, color=BLUE)
    ax_r.set_ylim(1.2, 6.7)
    ax_r.set_yticks([1.5, 2.5, 3.5, 4.5, 5.5])

    h1, l1 = ax_b.get_legend_handles_labels()
    h2, l2 = ax_r.get_legend_handles_labels()
    ax_r.legend(
        h1 + h2,
        l1 + l2,
        fontsize=5.6,
        frameon=False,
        loc="upper right",
        ncol=2,
        columnspacing=0.6,
        handletextpad=0.35,
    )

    fig.savefig(HERE / "Fig2.pdf")
    plt.close(fig)


def make_fig3():
    """Fig. 3: S1-to-S2 null whitening correction and S2 confirmatory accuracy."""
    if not S1_RES.exists() or not S2_RES.exists():
        raise FileNotFoundError("Missing frozen confirmatory results.json")
    s1 = json.loads(S1_RES.read_text(encoding="utf-8"))
    s2 = json.loads(S2_RES.read_text(encoding="utf-8"))
    a1 = s1["summary"]["arm_summaries"]
    a2 = s2["summary"]["arm_summaries"]

    plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42})
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(4.685, 2.75), constrained_layout=True
    )

    # Panel a: S1 vs S2 whitened peak Z_max means & 95% CIs across arms
    _style(
        ax_a,
        "a",
        "",
        "whitened peak $Z_{\\max}$ (mean $\\pm$ 95% CI)",
        xlog=False,
        ylog=False,
    )
    arm_order = [
        ("N22", "N22"),
        ("N23", "N23"),
        ("X_null_q1", "Nq1"),
        ("X_sub_q4_r32", "X4,32"),
        ("X_sub_q1_r48", "X1,48"),
        ("W1", "W1"),
        ("W1p", "W1p"),
        ("W2", "W2"),
        ("W3", "W3"),
    ]
    xs = np.arange(len(arm_order))
    for i, (arm, _) in enumerate(arm_order):
        if arm in a1 and "z_max" in a1[arm]:
            z1 = a1[arm]["z_max"]
            ax_a.errorbar(
                [i - 0.14],
                [z1["mean"]],
                yerr=[[z1["mean"] - z1["ci95_lo"]], [z1["ci95_hi"] - z1["mean"]]],
                fmt="s",
                color=GRAY,
                markerfacecolor="white",
                markersize=3.8,
                capsize=1.8,
                linewidth=1.0,
                label="S1 ($1+c_1/\\omega^2$)" if i == 0 else None,
            )
        z2 = a2[arm]["z_max"]
        is_null = arm in ("N22", "N23", "X_null_q1")
        col = TEAL if is_null else BLUE
        mk = "D" if is_null else "o"
        lbl = (
            "S2 null ($\\omega^{-2}+\\omega^{-4}$)"
            if arm == "N22"
            else ("S2 wave ($\\omega^{-2}+\\omega^{-4}$)" if arm == "W1" else None)
        )
        ax_a.errorbar(
            [i + 0.14],
            [z2["mean"]],
            yerr=[[z2["mean"] - z2["ci95_lo"]], [z2["ci95_hi"] - z2["mean"]]],
            fmt=mk,
            color=col,
            markersize=4.0,
            capsize=1.8,
            linewidth=1.2,
            label=lbl,
        )
    ax_a.axhline(2.50, color=BLUE, linestyle="--", linewidth=0.8)
    ax_a.axhline(3.30, color=GRAY, linestyle=":", linewidth=0.8)
    ax_a.text(
        -0.25,
        3.52,
        "N22+N23 peak $\\leq 32$:\n28/48 $\\rightarrow$ 3/56",
        fontsize=6.0,
        color=INK,
        va="bottom",
    )
    ax_a.set_xticks(xs)
    ax_a.set_xticklabels(
        [lbl for _, lbl in arm_order], fontsize=6.6, rotation=32, ha="right"
    )
    ax_a.set_ylim(1.4, 5.85)
    ax_a.legend(fontsize=5.8, frameon=False, loc="upper left")

    # Panel b: S2 blind focusing error |q_hat - q| and detection/lock counts
    _style(
        ax_b,
        "b",
        "blind focusing error $|\\hat{q} - q|$",
        "",
        xlog=False,
        ylog=False,
    )
    ys = np.arange(len(arm_order))[::-1]
    for y_pos, (arm, _short_lbl) in zip(ys, arm_order, strict=True):
        info = a2[arm]
        qe = info["q_abs_err"]
        is_null = arm in ("N22", "N23", "X_null_q1")
        col = TEAL if is_null else BLUE
        mk = "D" if is_null else "o"
        ax_b.errorbar(
            [qe["mean"]],
            [y_pos],
            xerr=[[qe["mean"] - qe["ci95_lo"]], [qe["ci95_hi"] - qe["mean"]]],
            fmt=mk,
            color=col,
            markersize=3.8,
            capsize=1.8,
            linewidth=1.15,
        )
        n_arm = info["n"]
        det_k = info["detected_k"]
        lock_k = info["freq_locked_k"]
        tag = (
            f"FP {det_k}/{n_arm}" if is_null else f"{det_k}/{n_arm}, L {lock_k}/{n_arm}"
        )
        ax_b.text(0.0365, y_pos, tag, va="center", ha="left", fontsize=6.2, color=INK)
    ax_b.axvline(0.025, color=GRAY, linestyle=":", linewidth=0.75)
    ax_b.axvline(0.035, color=GRAY, linestyle="--", linewidth=0.75)
    ax_b.set_yticks(ys)
    ax_b.set_yticklabels([lbl for _, lbl in arm_order], fontsize=6.8)
    ax_b.set_xlim(0.0, 0.056)
    ax_b.set_ylim(-0.6, len(arm_order) - 0.3)

    fig.savefig(HERE / "Fig3.pdf")
    plt.close(fig)


if __name__ == "__main__":
    make_fig1()
    make_fig2()
    make_fig3()
    print("Generated Fig1.pdf, Fig2.pdf, Fig3.pdf in", HERE)
