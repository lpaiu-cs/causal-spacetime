"""Regenerate the journal-format (vector PDF) figures for the Paper A draft.

Companion to `docs/paper/paper_a/figures/make_figures.py` /
`make_ladder_figure.py`, which own the manuscript.md PNG figures and stay
untouched. This script reads the SAME committed inputs -- the 19 summary
CSVs under `docs/paper/paper_a/figures/data/` and the preregistered
artifacts under `docs/prereg/` -- and emits the six PDFs the LaTeX draft
includes. Nothing is typed in: every plotted number is parsed from a
committed file, and the two new result figures (fig4, fig5) draw the raw
per-reading arrays stored inside the frozen artifacts.

Journal styling: text set in Latin Modern (the iopart preprint form is
Computer Modern; the document mimic uses its fontspec successor, so the
figures match the body), Computer Modern mathtext, Okabe-Ito palette
(colorblind-safe; per IOP figure guidance, color is never the sole
channel -- series also differ by marker or hatch), no in-figure titles --
captions carry the framing.

Usage: python docs/paper/paper_a/latex/figures/make_journal_figures.py
"""

from __future__ import annotations

import csv
import json
import statistics as st
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

DATA = Path("docs/paper/paper_a/figures/data")
PREREG = Path("docs/prereg")
OUT = Path("docs/paper/paper_a/latex/figures")

BLUE = "#0072B2"
VERM = "#D55E00"
GREEN = "#009E73"
ORANGE = "#E69F00"
SKY = "#56B4E9"
INK = "#222222"
GREY = "#555555"
MUTED = "#888888"
LIGHT = "#CCCCCC"
GRID = "#DDDDDD"

LM = Path("/usr/local/texlive/2025/texmf-dist/fonts/opentype/public/lm")


def _setup_fonts() -> None:
    for stem in ("lmroman10-regular", "lmroman10-bold",
                 "lmroman10-italic", "lmroman10-bolditalic"):
        f = LM / f"{stem}.otf"
        if f.exists():
            font_manager.fontManager.addfont(str(f))
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Latin Modern Roman", "CMU Serif", "STIXGeneral"],
        "mathtext.fontset": "cm",
        "font.size": 8.5,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "pdf.fonttype": 42,
    })


def _rows(name: str) -> list[dict]:
    with (DATA / name).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _f(r: dict, k: str) -> float:
    return float(r[k])


def _style(ax) -> None:
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=INK)
    ax.grid(True, color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def _panel(ax, tag: str) -> None:
    ax.set_title(tag, loc="left", fontsize=9, color=INK, pad=6)


# --------------------------------------------------------------------------
# Figure 0 (new): setup illustration -- sprinkled diamond, radar protocol,
# Rindler wedge. An ILLUSTRATION, not a result: fixed seed, drawn with the
# repository's own foundation modules, no measured quantity displayed.
# --------------------------------------------------------------------------
def figure_setup() -> None:
    import sys

    sys.path.insert(0, "src")
    import numpy as np

    from causal_spacetime_lab.causal import causal_matrix_1p1
    from causal_spacetime_lab.chains import longest_chain_indices
    from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(6.1, 2.3))

    # -- (a) Poisson sprinkling into a causal diamond (interval endpoints
    #    (-T/2, 0) and (T/2, 0)), Hasse links, and one longest chain
    #    between the diamond tips.
    T = 2.0
    half = T / 2
    ev = sprinkle_1p1_causal_diamond(110, T=T, seed=7)
    ev = np.vstack([[-half, 0.0], ev, [half, 0.0]])
    C = causal_matrix_1p1(ev)
    hasse = C & ~(C @ C)
    chain = longest_chain_indices(C, 0, len(ev) - 1, event_times=ev[:, 0])
    for i, j in zip(*np.nonzero(hasse), strict=True):
        ax1.plot([ev[i, 1], ev[j, 1]], [ev[i, 0], ev[j, 0]], color=LIGHT,
                 lw=0.4, zorder=1)
    ax1.plot(ev[:, 1], ev[:, 0], ".", ms=2.6, color=GREY, zorder=2)
    ax1.plot(ev[chain, 1], ev[chain, 0], "-o", ms=3.2, lw=1.3, color=BLUE,
             zorder=3)
    ax1.plot([0, half, 0, -half, 0], [-half, 0, half, 0, -half], ls="--",
             color=INK, lw=0.7, zorder=2)
    ax1.set_xlim(-1.18, 1.18)
    ax1.set_ylim(-half - 0.14, half + 0.14)
    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$t$")
    _panel(ax1, "(a) sprinkled diamond")

    # -- (b) radar protocol: observer chain with tick labels, one target
    #    event, the bounding null signals.
    t_e, x_e = 1.05, 0.62
    ax2.plot([0, 0], [-0.1, 2.1], color=INK, lw=1.2)
    for tk in np.arange(0.0, 2.01, 0.25):
        ax2.plot([-0.03, 0.03], [tk, tk], color=INK, lw=0.8)
    tau_m, tau_p = t_e - x_e, t_e + x_e
    ax2.plot([0, x_e], [tau_m, t_e], color=GREY, lw=0.9, ls=":")
    ax2.plot([x_e, 0], [t_e, tau_p], color=GREY, lw=0.9, ls=":")
    ax2.plot([x_e], [t_e], "o", ms=4, color=BLUE)
    ax2.annotate("$e$", (x_e, t_e), xytext=(5, 0), textcoords="offset points",
                 fontsize=8, color=BLUE)
    ax2.plot([0, 0], [tau_m, tau_p], ls="none")
    ax2.plot(0, tau_m, "o", ms=3, color=INK)
    ax2.plot(0, tau_p, "o", ms=3, color=INK)
    ax2.annotate(r"$\tau_-$", (0, tau_m), xytext=(-14, -3),
                 textcoords="offset points", fontsize=8, color=INK)
    ax2.annotate(r"$\tau_+$", (0, tau_p), xytext=(-14, -2),
                 textcoords="offset points", fontsize=8, color=INK)
    ax2.set_xlim(-0.35, 1.05)
    ax2.set_ylim(-0.15, 2.25)
    ax2.set_xlabel("$x$")
    _panel(ax2, "(b) radar protocol")
    ax2.set_yticks([])

    # -- (c) Rindler observer: hyperbolic worldline, wedge, horizon.
    a = 1.0
    tt = np.linspace(-1.6, 1.6, 200)
    xx = np.sqrt(1.0 / a**2 + tt**2)
    lim = 1.9
    wedge_t = np.linspace(-lim, lim, 100)
    ax3.fill_betweenx(wedge_t, np.abs(wedge_t), lim, color="#EFEFEC",
                      zorder=0)
    ax3.plot([0, lim], [0, lim], color=INK, ls="--", lw=0.8)
    ax3.plot([0, lim], [0, -lim], color=INK, ls="--", lw=0.8)
    ax3.plot(xx, tt, color=BLUE, lw=1.4)
    ax3.annotate("horizon", (0.9 * lim, 0.9 * lim), xytext=(-4, -12),
                 textcoords="offset points", fontsize=7.5, color=INK,
                 ha="right", rotation=45)
    ax3.annotate("wedge", (1.45, -1.05), fontsize=8, color=GREY)
    ax3.set_xlim(-0.4, lim)
    ax3.set_ylim(-lim, lim)
    ax3.set_xlabel("$x$")
    _panel(ax3, "(c) Rindler observer")
    ax3.set_yticks([])

    for ax in (ax1, ax2, ax3):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.grid(False)

    fig.tight_layout(pad=0.5)
    fig.savefig(OUT / "fig0_setup.pdf")
    plt.close(fig)
    print(f"  fig0: sprinkle N={len(ev)} (seed 7), longest chain "
          f"length {len(chain)} (illustration only)")


# --------------------------------------------------------------------------
# Figure 1: the reconstruction ladder (schematic; mirrors Table 1 / Sec. 2)
# --------------------------------------------------------------------------
LADDER = [
    ("R0", "order only",
     "flat-sprinkling dimension; raw chain statistics",
     "model assumptions; no metric scale"),
    ("R1", "order + global density", "timelike proper time",
     "needs a density"),
    ("R2", "order + observer clock", "radar time; unsigned distance",
     "sign undetermined"),
    ("R3", "order + observer + orientation", "signed coords; Lorentz map",
     "calibrated separation"),
    ("R4", "order + observer + oriented atlas", "Poincaré transition maps",
     "calibrated charts"),
    ("R5", "order + global/local measure", "volume; coarse-grain stability",
     "global/local measure supplied"),
]


def figure_ladder() -> None:
    fig, ax = plt.subplots(figsize=(6.1, 3.7))
    ax.axis("off")
    n = len(LADDER)
    for i, (rung, ingredient, recovered, bound) in enumerate(LADDER):
        y = n - 1 - i
        ax.add_patch(plt.Rectangle((0.0, y + 0.07), 9.7, 0.86,
                                   facecolor="#F5F5F3", edgecolor=MUTED,
                                   linewidth=0.6, zorder=1))
        ax.text(0.25, y + 0.5, rung, fontsize=10.5, fontweight="bold",
                color=BLUE, va="center", zorder=2)
        ax.text(1.0, y + 0.63, ingredient, fontsize=8.5, color=INK,
                va="center", zorder=2)
        ax.annotate("", xy=(4.72, y + 0.63), xytext=(4.42, y + 0.63),
                    arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.0))
        ax.text(4.9, y + 0.63, recovered, fontsize=8.5, color=INK,
                va="center", zorder=2)
        ax.text(1.0, y + 0.28, f"bounded by: {bound}", fontsize=7.5,
                color=MUTED, va="center", zorder=2, style="italic")
    ax.text(1.0, n + 0.06,
            "required supplied structure $\\rightarrow$ reconstructed quantity",
            fontsize=8, color=MUTED)
    ax.set_xlim(-0.1, 9.8)
    ax.set_ylim(0, n + 0.42)
    fig.tight_layout(pad=0.4)
    fig.savefig(OUT / "fig1_ladder.pdf")
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
    fig, axes = plt.subplots(2, 2, figsize=(6.1, 4.6))
    n_ticks = [300, 600, 1200, 2400]

    def _n_axis(ax):
        ax.set_xscale("log")
        ax.set_xticks(n_ticks)
        ax.set_xticklabels([str(v) for v in n_ticks])
        ax.minorticks_off()

    ax = axes[0, 0]
    _style(ax)
    dim_rows = _rows("dimension_reconstruction_summary.csv")
    for d, color in zip((2, 3, 4), (BLUE, VERM, GREEN), strict=True):
        sub = [r for r in dim_rows if int(_f(r, "spacetime_dim")) == d]
        sub.sort(key=lambda r: _f(r, "N"))
        xs = [_f(r, "N") for r in sub]
        ys = [_f(r, "mean_estimated_dim") for r in sub]
        ax.plot(xs, ys, marker="o", ms=3.5, color=color, lw=1.4)
        ax.axhline(d, color=color, ls=":", lw=0.7, alpha=0.6)
        ax.annotate(f"$D = {d}$", (xs[-1], ys[-1] + 0.07), fontsize=8,
                    color=color, ha="right")
    _n_axis(ax)
    ax.set_xlabel("$N$ (events)")
    ax.set_ylabel("estimated dimension")
    ax.set_ylim(1.8, 4.3)
    _panel(ax, "(a) dimension")

    ax = axes[0, 1]
    _style(ax)
    tp = _rows("timelike_pair_reconstruction_summary.csv")
    tp.sort(key=lambda r: _f(r, "N"))
    ns = [_f(r, "N") for r in tp]
    ax.plot(ns, [_f(r, "tau_volume_relative_rmse") for r in tp], marker="o",
            ms=3.5, color=BLUE, lw=1.4, label="volume estimator")
    ax.plot(ns, [_f(r, "tau_chain_relative_rmse") for r in tp], marker="s",
            ms=3.5, color=VERM, lw=1.4, label="chain estimator")
    _n_axis(ax)
    ax.set_xlabel("$N$ (events)")
    ax.set_ylabel("proper-time relative RMSE")
    _panel(ax, "(b) timelike proper time")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    _style(ax)
    rr = _rows("discrete_radar_reconstruction_summary.csv")
    xt, dist = _mean_by(rr, "tick_count", "radar_distance_rmse")
    _, tim = _mean_by(rr, "tick_count", "radar_time_rmse")
    ax.plot(xt, dist, marker="o", ms=3.5, color=BLUE, lw=1.4,
            label="radar distance")
    ax.plot(xt, tim, marker="s", ms=3.5, color=VERM, lw=1.4,
            label="radar time")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("observer ticks")
    ax.set_ylabel("RMSE (mean over $N$)")
    _panel(ax, "(c) observer radar")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    _style(ax)
    lz = _rows("oriented_radar_lorentz_summary.csv")
    for beta, color in zip((0.3, 0.6), (BLUE, VERM), strict=True):
        xt, br = _mean_by(lz, "tick_count", "fitted_beta_rmse",
                          where=lambda r, b=beta: abs(_f(r, "beta") - b) < 1e-9)
        ax.plot(xt, br, marker="o", ms=3.5, color=color, lw=1.4,
                label=f"$\\beta = {beta}$")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("observer ticks")
    ax.set_ylabel("fitted-$\\beta$ RMSE (mean over $N$)")
    _panel(ax, "(d) Lorentz-map recovery")
    ax.legend(frameon=False)

    fig.tight_layout(pad=0.5)
    fig.savefig(OUT / "fig2_convergence.pdf")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 3: measure dependence (R5)
# --------------------------------------------------------------------------
def figure_measure() -> None:
    """R5 re-pointed to the position-dependent profile; the constant
    arm is not plotted (its weighted relative RMSE is the flat
    profile's identically, so it demonstrates M, not W)."""

    def profile(name):
        rows = [r for r in _rows("weighted_conformal_volume_summary.csv")
                if r["profile"] == name]
        rows.sort(key=lambda r: _f(r, "N"))
        return rows

    sin, flat = profile("sinusoidal_0.3"), profile("flat")
    ns = [_f(r, "N") for r in sin]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.1, 2.6))

    _style(ax1)
    ax1.plot(ns, [_f(r, "unweighted_relative_rmse") for r in sin], marker="s",
             ms=3.5, color=VERM, lw=1.4, label="unweighted (density only)")
    ax1.plot(ns, [_f(r, "weighted_relative_rmse") for r in sin], marker="o",
             ms=3.5, color=BLUE, lw=1.4, label="local-measure weighted")
    ax1.plot(ns, [_f(r, "weighted_relative_rmse") for r in flat], ls=":",
             color=MUTED, lw=1.3, label="flat-profile baseline")
    ax1.set_xscale("log")
    ax1.set_xticks([600, 1200, 2400])
    ax1.set_xticklabels(["600", "1200", "2400"])
    ax1.minorticks_off()
    ax1.set_xlabel("$N$ (events)")
    ax1.set_ylabel("volume relative RMSE")
    ax1.set_ylim(0, None)
    ax1.legend(frameon=False, loc="center right", fontsize=7)
    _panel(ax1, "(a) floor vs recovered scaling")

    _style(ax2)
    ax2.plot(ns, [_f(r, "unweighted_volume_rmse") for r in sin], marker="s",
             ms=3.5, color=VERM, lw=1.4, label="unweighted (density only)")
    ax2.plot(ns, [_f(r, "weighted_volume_rmse") for r in sin], marker="o",
             ms=3.5, color=BLUE, lw=1.4, label="local-measure weighted")
    ax2.plot(ns, [_f(r, "weighted_volume_rmse") for r in flat], ls=":",
             color=MUTED, lw=1.3, label="flat-profile baseline")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_yticks([0.015, 0.02, 0.03, 0.04, 0.06])
    ax2.set_yticklabels(["0.015", "0.02", "0.03", "0.04", "0.06"])
    ax2.set_xticks([600, 1200, 2400])
    ax2.set_xticklabels(["600", "1200", "2400"])
    ax2.minorticks_off()
    ax2.set_xlabel("$N$ (events)")
    ax2.set_ylabel("volume RMSE")
    ax2.legend(frameon=False, loc="lower left", fontsize=7)
    _panel(ax2, "(b) logarithmic axes")

    fig.tight_layout(pad=0.5)
    fig.savefig(OUT / "fig3_measure.pdf")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 4 (new): plane-wave capstone -- C1 paired shift and C2 separation
# --------------------------------------------------------------------------
def figure_capstone() -> None:
    art = json.loads((PREREG / "p14_prereg_results.json").read_text(
        encoding="utf-8"))
    assert art["stage_positive"] is True and art["c1"]["verdict"] == "confirmed" \
        and art["c2"]["verdict"] == "confirmed", "unexpected verdicts"
    delta = art["c1"]["raw"]["delta"]
    mean = art["c1"]["metrics"]["mean"]
    eps = art["eps_delta"]
    f_c = art["c2"]["raw"]["f_curved"]
    f_f = art["c2"]["raw"]["f_flat"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.1, 2.6))

    _style(ax1)
    ax1.hist(delta, bins=48, color=BLUE, alpha=0.85, edgecolor="white",
             linewidth=0.2)
    ax1.axvline(0.0, color=GREY, lw=0.8)
    ax1.axvline(eps, color=VERM, lw=1.2, ls="--")
    ax1.axvline(mean, color=INK, lw=1.0)
    ax1.annotate("frozen margin $\\varepsilon_\\Delta$", (eps, 1.0),
                 xycoords=("data", "axes fraction"), xytext=(6, -12),
                 textcoords="offset points", fontsize=7.5, color=VERM,
                 va="top")
    ax1.annotate("mean", (mean, 0.86), xycoords=("data", "axes fraction"),
                 xytext=(4, 0), textcoords="offset points", fontsize=7.5,
                 color=INK)
    ax1.set_xlabel("paired difference $f_A - f_0$")
    ax1.set_ylabel("sprinklings")
    _panel(ax1, f"(a) C1 paired shift ($n = {len(delta)}$)")

    _style(ax2)
    bins = 48
    ax2.hist(f_f, bins=bins, color=SKY, alpha=0.85, edgecolor="white",
             linewidth=0.2, label="flat arm")
    ax2.hist(f_c, bins=bins, color=ORANGE, alpha=0.85, edgecolor=GREY,
             linewidth=0.4, hatch="///", label="curved arm")
    ax2.set_xlabel("relation fraction $f$ (single poset)")
    ax2.set_ylabel("sprinklings")
    _panel(ax2, f"(b) C2 arms ($n = {len(f_c)}$ per arm)")
    ax2.legend(frameon=False)

    fig.tight_layout(pad=0.5)
    fig.savefig(OUT / "fig4_capstone.pdf")
    plt.close(fig)
    print(f"  fig4: C1 mean {mean:.7f} vs eps_delta {eps:.7g}; "
          f"C2 separation min(curved)={min(f_c):.6f} > max(flat)={max(f_f):.6f}: "
          f"{min(f_c) > max(f_f)}")


# --------------------------------------------------------------------------
# Figure 5 (new): Schwarzschild extension -- S4 paired shift, S5 arms
# --------------------------------------------------------------------------
def figure_schwarzschild() -> None:
    s4 = json.loads((PREREG / "p14_s4_results.json").read_text(encoding="utf-8"))
    s5 = json.loads((PREREG / "p14_s5_results.json").read_text(encoding="utf-8"))
    assert s4["verdict"] == "CONFIRMED" and s5["outcome"] == "DETECTED", \
        "unexpected verdicts"
    delta = s4["delta_lower"]["per_reading"]
    ci = s4["identified_ci95"]
    eps_det = -s4["margins"]["eps_det"]
    f_c = s5["curved_f_lower"]["per_reading"]
    f_f = s5["flat_f"]["per_reading"]
    auc = s5["auc"]["auc_lower_series"][0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.1, 2.6))

    _style(ax1)
    ax1.hist(delta, bins=32, color=BLUE, alpha=0.85, edgecolor="white",
             linewidth=0.2)
    ax1.axvline(0.0, color=GREY, lw=0.8)
    ax1.axvline(eps_det, color=VERM, lw=1.2, ls="--")
    ax1.annotate("frozen gate $-\\varepsilon_{\\mathrm{det}}$",
                 (eps_det, 1.0), xycoords=("data", "axes fraction"),
                 xytext=(-6, -12), textcoords="offset points", fontsize=7.5,
                 color=VERM, ha="right", va="top")
    ax1.axvspan(ci[0], ci[1], color=INK, alpha=0.25, lw=0)
    ax1.annotate("identified CI95",
                 xy=((ci[0] + ci[1]) / 2, 24), xytext=(-0.024, 27),
                 fontsize=7.5, color=INK, ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", color=INK, lw=0.7,
                                 shrinkA=2, shrinkB=2))
    ax1.set_xlabel("paired shift $\\Delta f$ per reading")
    ax1.set_ylabel("readings")
    _panel(ax1, f"(a) S4 paired shift ($n = {len(delta)}$)")

    _style(ax2)
    bins = 26
    ax2.hist(f_f, bins=bins, color=SKY, alpha=0.85, edgecolor="white",
             linewidth=0.2, label="flat arm")
    ax2.hist(f_c, bins=bins, color=ORANGE, alpha=0.75, edgecolor=GREY,
             linewidth=0.4, hatch="///", label="Schwarzschild arm")
    ax2.annotate(f"AUC $= {auc:.4f}$", (0.03, 0.92), xycoords="axes fraction",
                 fontsize=8, color=INK)
    ax2.set_xlabel("relation fraction $f$ (single poset)")
    ax2.set_ylabel("readings")
    _panel(ax2, f"(b) S5 arms ($n = {len(f_c)}$ per arm)")
    ax2.legend(frameon=False, loc="upper right")

    fig.tight_layout(pad=0.5)
    fig.savefig(OUT / "fig5_schwarzschild.pdf")
    plt.close(fig)
    print(f"  fig5: S4 CI95 [{ci[0]:.6f}, {ci[1]:.6f}] vs gate {eps_det}; "
          f"S5 AUC {auc:.4f}, overlap exists: {min(f_c) <= max(f_f)}")


# --------------------------------------------------------------------------
# Figure 6: the executed mass ladder (journal restyle of fig4_ladder_count)
# --------------------------------------------------------------------------
RUNGS = (
    ("0.1333", "p14_o5_count.json"),
    ("0.1867", "p14_s6_m14_count.json"),
    ("0.2400", "p14_s6_m18_count.json"),
    ("0.4000", "p14_s6_m30_count.json"),
)


def _load_rungs() -> list[dict]:
    rows = []
    for mu, name in RUNGS:
        a = json.loads((PREREG / name).read_text(encoding="utf-8"))
        d, fz, scan = a["decision"], a["frozen_config"], a["scan"]
        rows.append({
            "mu": float(mu),
            "v_lo": fz["v_lo"], "v_hi": fz["v_hi"],
            "c_lo": d["c_lo_volume"], "c_hi": d["c_hi_volume"],
            "d_lo": d["d_lo"], "d_hi": d["d_hi"], "band": d["band"],
            "k": scan["k_certain"], "u": scan["u_ambiguous"],
            "verdict": d["verdict"],
        })
    return rows


def figure_count() -> None:
    rows = _load_rungs()
    assert all(r["verdict"] == "CONCORDANT" for r in rows), "unexpected verdict"

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(6.1, 2.7), gridspec_kw={"width_ratios": [1.05, 1]})

    ax1.plot([r["mu"] for r in rows],
             [(r["v_lo"] + r["v_hi"]) / 2 for r in rows],
             color=LIGHT, lw=1.0, zorder=0)
    for i, r in enumerate(rows):
        x = r["mu"]
        ax1.plot([x - 0.005] * 2, [r["v_lo"], r["v_hi"]], color=BLUE, lw=5,
                 solid_capstyle="butt",
                 label="certified $V$ (oracle)" if i == 0 else None)
        ax1.plot([x + 0.005] * 2, [r["c_lo"], r["c_hi"]], color=ORANGE, lw=5,
                 solid_capstyle="butt",
                 label="$C$ = count / intensity" if i == 0 else None)
        ax1.annotate(f'{(r["v_lo"] + r["v_hi"]) / 2:.1f}', (x, r["v_hi"]),
                     textcoords="offset points", xytext=(0, 5),
                     ha="center", fontsize=7, color=GREY)
    ax1.text(0.435, 62, "intervals overlap at this scale;\nresidual in (b)",
             fontsize=7.5, color=GREY, style="italic", va="bottom",
             ha="right")
    ax1.set_xlabel("compactness $\\mu = 2M/r_c$")
    ax1.set_ylabel("certified 4-volume of the rung's diamond")
    _panel(ax1, "(a) certified volume vs count")
    ax1.set_xticks([r["mu"] for r in rows])
    ax1.set_xticklabels([f'{r["mu"]:.4f}' for r in rows], fontsize=7)
    ax1.set_xlim(0.10, 0.44)
    ax1.grid(axis="y", color=LIGHT, lw=0.5)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=7.5, loc="upper left",
               bbox_to_anchor=(0.0, 0.97))

    ys = list(range(len(rows)))[::-1]
    ax2.axvspan(-1, 1, color=LIGHT, alpha=0.45, lw=0)
    ax2.axvline(0, color=GREY, lw=0.7)
    for y, r in zip(ys, rows, strict=True):
        lo, hi = r["d_lo"] / r["band"], r["d_hi"] / r["band"]
        ax2.plot([lo, hi], [y, y], color=BLUE, lw=5, solid_capstyle="butt")
        ax2.plot([lo, hi], [y, y], "|", color=BLUE, ms=7, mew=1.3)
    ax2.axvline(-1, color=GREY, ls="--", lw=0.8)
    ax2.axvline(1, color=GREY, ls="--", lw=0.8)
    ax2.set_yticks(ys)
    ax2.set_yticklabels(
        [f'$\\mu$ = {r["mu"]:.4f}\nK = {r["k"]:,}, U = {r["u"]}'
         for r in rows], fontsize=7)
    # The band's definition lives in the caption; spelling it out here
    # ran the label past the figure edge and clipped the subscript.
    ax2.set_xlabel("identified discrepancy $D$, in units of the band $B$")
    _panel(ax2, "(b) equivalence gate ($\\tau$ = 2.5%)")
    ax2.set_xlim(-1.35, 1.35)
    ax2.set_ylim(-0.6, len(rows) - 0.4)
    ax2.text(1.02, len(rows) - 0.55, " gate", fontsize=7, color=GREY,
             va="center")
    ax2.grid(axis="x", color=LIGHT, lw=0.5)
    ax2.set_axisbelow(True)

    for ax in (ax1, ax2):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    fig.tight_layout(pad=0.5)
    fig.savefig(OUT / "fig6_ladder_count.pdf")
    plt.close(fig)
    for r in rows:
        print(f'  fig6: mu={r["mu"]:.4f} V=[{r["v_lo"]:.4f}, {r["v_hi"]:.4f}] '
              f'D/B=[{r["d_lo"] / r["band"]:+.3f}, {r["d_hi"] / r["band"]:+.3f}] '
              f'{r["verdict"]}')


if __name__ == "__main__":
    _setup_fonts()
    OUT.mkdir(parents=True, exist_ok=True)
    figure_setup()
    figure_ladder()
    figure_convergence()
    figure_measure()
    figure_capstone()
    figure_schwarzschild()
    figure_count()
    for name in ("fig0_setup", "fig1_ladder", "fig2_convergence",
                 "fig3_measure", "fig4_capstone", "fig5_schwarzschild",
                 "fig6_ladder_count"):
        print("wrote", OUT / f"{name}.pdf")
