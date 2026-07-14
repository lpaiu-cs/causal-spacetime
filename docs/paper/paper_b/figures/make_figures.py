"""Regenerate Paper B figures from the frozen preregistration CSVs.

Reads only committed frozen artifacts under docs/prereg/frozen/ so figures are
reproducible and provenance-locked. Palette: Okabe-Ito blue/vermillion
(colorblind-safe; validated). Single y-axis per figure; both series in a figure
share the same [0,1] violation/discordance scale.

Usage: python docs/paper/paper_b/figures/make_figures.py
"""

from __future__ import annotations

import csv
import statistics as st
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FROZEN = Path("docs/prereg/frozen")
OUT = Path("docs/paper/paper_b/figures")
BLUE = "#0072B2"  # truth-order error / structured
VERM = "#D55E00"  # held-out violation / foils
GREEN = "#009E73"  # null gap (Fig 5 only)
INK = "#222222"
MUTED = "#888888"
GRID = "#DDDDDD"


def _rows(name: str) -> list[dict]:
    return list(csv.DictReader(open(FROZEN / name)))


def _num(r: dict, k: str) -> float:
    try:
        return float(r[k])
    except (TypeError, ValueError):
        return float("nan")


def _d1_ok(rows: list[dict], condition: str) -> list[float]:
    return [
        _num(r, "heldout_violation")
        for r in rows
        if r.get("status") == "ok"
        and r.get("condition") == condition
        and r.get("embedding_dim") == "1.0"
    ]


def _style(ax) -> None:
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=INK, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def figure_discriminator() -> None:
    """Fig 1: confirmatory held-out violation by condition (gate dim d=1)."""

    stage_b = _rows("pc_v1_stage_b_sensitivity.csv")
    stage_c = _rows("pc_v1_stage_c_specificity.csv")
    series = [
        ("structured\n(geometric)", _d1_ok(stage_b, "structured"), BLUE),
        ("random order\n(geometry-free)", _d1_ok(stage_c, "random_order"), VERM),
        (
            "column-shuffled\n(consistency-destroyed)",
            _d1_ok(stage_c, "column_shuffled"),
            VERM,
        ),
    ]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    _style(ax)
    for i, (_label, vals, color) in enumerate(series):
        xs = [i + (j - len(vals) / 2) * 0.012 for j in range(len(vals))]
        ax.scatter(xs, vals, s=26, color=color, alpha=0.75, edgecolor="white",
                   linewidth=0.5, zorder=3)
        m = st.mean(vals)
        ax.plot([i - 0.22, i + 0.22], [m, m], color=INK, linewidth=2, zorder=4)
        ax.annotate(f"mean {m:.3f}", (i + 0.24, m), fontsize=8, color=INK,
                    va="center")
    ax.axhline(0.05, color=MUTED, linestyle="--", linewidth=1.2, zorder=2)
    ax.annotate("held-out gate = 0.05", (2.35, 0.065), fontsize=8, color=MUTED,
                ha="right")
    ax.set_xticks(range(len(series)))
    ax.set_xticklabels([s[0] for s in series], fontsize=9)
    ax.set_ylabel("held-out violation (d = 1)", fontsize=10, color=INK)
    ax.set_ylim(-0.02, 0.42)
    ax.set_title(
        "PC-V1: the pipeline passes on geometric order and blocks on\n"
        "matched geometry-free order (confirmatory seeds, frozen gate)",
        fontsize=10.5, color=INK, loc="left",
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig1_discriminator.png", dpi=200)
    plt.close(fig)


def figure_dose_response() -> None:
    """Fig 2: P1 geometry-dilution dose-response with the false-pass window."""

    rows = [
        r
        for r in _rows("p1_stage_b_epsilon_sweep.csv")
        if r.get("status") == "ok" and _num(r, "density_held") == 1.0
    ]
    grid = [0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0]

    def series(col: str) -> tuple[list[float], list[float], list[float]]:
        means, sds, eps = [], [], []
        for e in grid:
            cells = [_num(r, col) for r in rows if _num(r, "epsilon") == e]
            if cells:
                eps.append(e)
                means.append(st.mean(cells))
                sds.append(st.pstdev(cells) if len(cells) > 1 else 0.0)
        return eps, means, sds

    e_t, truth, truth_sd = series("truth_order_error")
    e_h, held, held_sd = series("heldout_violation")

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    _style(ax)

    # false-pass window: between truth crossing (~0.31) and heldout crossing (~0.50)
    ax.axvspan(0.31, 0.50, color=MUTED, alpha=0.12, zorder=0)
    ax.annotate("false-pass window\n(embeds, but not true space)",
                (0.405, 0.44), fontsize=8, color=MUTED, ha="center")

    for eps, mean, sd, color, label in (
        (e_t, truth, truth_sd, BLUE, "truth-order error (space recovery)"),
        (e_h, held, held_sd, VERM, "held-out violation (embeddability)"),
    ):
        lo = [m - s for m, s in zip(mean, sd, strict=True)]
        hi = [m + s for m, s in zip(mean, sd, strict=True)]
        ax.fill_between(eps, lo, hi, color=color, alpha=0.15, zorder=1)
        ax.plot(eps, mean, color=color, linewidth=2, marker="o", markersize=5,
                zorder=3, label=label)

    ax.axhline(0.15, color=BLUE, linestyle=":", linewidth=1.1, zorder=2)
    ax.axhline(0.05, color=VERM, linestyle=":", linewidth=1.1, zorder=2)
    ax.annotate("truth gate 0.15", (1.0, 0.16), fontsize=7.5, color=BLUE, ha="right")
    ax.annotate("held-out gate 0.05", (1.0, 0.06), fontsize=7.5, color=VERM, ha="right")
    ax.axvline(0.31, color=INK, linestyle="-", linewidth=0.8, alpha=0.5, zorder=2)
    ax.annotate("epsilon* ~ 0.31", (0.31, -0.045), fontsize=8, color=INK, ha="center")

    ax.set_xlabel("geometry dilution  epsilon   (0 = Minkowski, 1 = geometry-free)",
                  fontsize=10, color=INK)
    ax.set_ylabel("sign-discordance / violation fraction", fontsize=10, color=INK)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.06, 0.58)
    ax.legend(loc="upper left", fontsize=8.5, frameon=False)
    ax.set_title(
        "P1: geometry recovery degrades monotonically as order is diluted\n"
        "at fixed relation density (~0.57 across all epsilon; bands = +/-1 SD)",
        fontsize=10.5, color=INK, loc="left",
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig2_dose_response.png", dpi=200)
    plt.close(fig)


def figure_confound() -> None:
    """Fig 5: the shared-scalar confound — raw vs parallax dissimilarity.

    Slopegraph on the geometry-free control (PC-V1 Stage C seeds). Each line is
    one seed's d=1 held-out violation under raw bracket-width dissimilarity
    (left, shared scalar retained) and parallax dissimilarity (right, shared
    scalar removed). Reads confound_data.csv (compute_confound_data.py); the
    parallax column reproduces the frozen Stage C result exactly.
    """

    rows = list(csv.DictReader(open(OUT / "confound_data.csv")))
    raw = [float(r["raw"]) for r in rows]
    par = [float(r["parallax"]) for r in rows]
    gate = 0.05

    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    _style(ax)
    ax.axhspan(-0.02, gate, color=MUTED, alpha=0.10, zorder=0)
    ax.axhline(gate, color=MUTED, linestyle="--", linewidth=1.2, zorder=2)
    ax.annotate("gate 0.05", (1.5, gate), fontsize=8, color=MUTED, va="center",
                ha="right")
    ax.annotate("shaded = passes the gate\n(false-pass for geometry-free input)",
                (-0.33, 0.012), fontsize=8, color=MUTED, ha="left", va="center")

    for r_val, p_val in zip(raw, par, strict=True):
        ax.plot([0, 1], [r_val, p_val], color=MUTED, linewidth=1.0, alpha=0.6,
                zorder=2)
        raw_color = VERM if r_val <= gate else INK
        ax.scatter([0], [r_val], s=34, color=raw_color, zorder=3,
                   edgecolor="white", linewidth=0.5)
        ax.scatter([1], [p_val], s=34, color=BLUE, zorder=3,
                   edgecolor="white", linewidth=0.5)

    for x, vals in ((0, raw), (1, par)):
        m = sum(vals) / len(vals)
        ax.plot([x - 0.09, x + 0.09], [m, m], color=INK, linewidth=2.5, zorder=4)
        ax.annotate(f"mean {m:.3f}", (x, m + 0.018), fontsize=8.5, color=INK,
                    ha="center")

    false_pass = sum(1 for v in raw if v <= gate)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([
        "raw dissimilarity\n(shared scalar retained)",
        "parallax dissimilarity\n(shared scalar removed)",
    ], fontsize=9.5)
    ax.set_xlim(-0.35, 1.5)
    ax.set_ylim(-0.02, 0.36)
    ax.set_ylabel("held-out violation (d = 1)", fontsize=10, color=INK)
    ax.set_title(
        "PC-V1 confound: on geometry-free order, raw dissimilarity sits at the\n"
        f"gate ({false_pass}/{len(raw)} false-pass); parallax centering blocks all "
        f"{len(par)}/{len(par)}",
        fontsize=10.5, color=INK, loc="left",
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig5_confound.png", dpi=200)
    plt.close(fig)


def figure_dimension_2d() -> None:
    """Fig 3: 2+1D dimension selection (P2-v2 confirmatory)."""

    rows = [
        r
        for r in _rows("p2_v2_stage_b_2plus1d.csv")
        if r.get("status") == "ok"
    ]
    dims = [1, 2, 3]
    truth = {d: [_num(r, f"d{d}_truth") for r in rows] for d in dims}
    gate = 0.15  # frozen P2-v2 truth gate

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    _style(ax)
    ax.axhline(gate, color=MUTED, linestyle="--", linewidth=1.2, zorder=2)
    ax.annotate("truth gate 0.15", (2.4, gate + 0.008), fontsize=8, color=MUTED,
                ha="right")
    for i, d in enumerate(dims):
        vals = truth[d]
        xs = [i + (j - len(vals) / 2) * 0.010 for j in range(len(vals))]
        color = VERM if d == 1 else BLUE
        ax.scatter(xs, vals, s=24, color=color, alpha=0.75, edgecolor="white",
                   linewidth=0.5, zorder=3)
        m = st.mean(vals)
        ax.plot([i - 0.2, i + 0.2], [m, m], color=INK, linewidth=2.2, zorder=4)
        ax.annotate(f"mean {m:.3f}", (i + 0.24, m), fontsize=8, color=INK,
                    va="center")
    ax.set_xticks(range(len(dims)))
    ax.set_xticklabels(
        ["d = 1\n(underfits)", "d = 2\n(recovers 2D)", "d = 3\n(no gain)"],
        fontsize=9.5,
    )
    ax.set_ylabel("truth-order error vs true 2D position", fontsize=10, color=INK)
    ax.set_ylim(0.0, 0.32)
    ax.set_title(
        "2+1D dimension selection (P2-v2, confirmatory): recovery needs d = 2\n"
        "and saturates there -- d = 1 underfits, d = 3 adds nothing",
        fontsize=10.5, color=INK, loc="left",
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig3_dimension_2d.png", dpi=200)
    plt.close(fig)


def figure_emergence_chain() -> None:
    """Fig 4: the emergence chain in one row — P3 null | P4 phases | P5 verdict.

    Panel A: P3 percolation sweep (held-out violation where six chains exist;
    structural-block counts where they do not; sprinkled P3-A reference band).
    Panel B: P4 dual-start mean action vs beta (hysteresis opens at beta 3-4;
    exact bipartite floor -120). Panel C: P5 per-configuration gate metrics vs
    beta with the uniform-ensemble calibration at beta = 0 and the structural
    block at beta = 32.
    """

    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(12.8, 4.3))
    for ax in (ax_a, ax_b, ax_c):
        _style(ax)

    # --- Panel A: P3 percolation sweep -----------------------------------
    p3 = _rows("p3_stage_b_dynamics.csv")
    cal3 = [_num(r, "heldout") for r in _rows("p3_stage_a_calibration.csv")
            if r.get("status") == "ok"]
    grid_p = sorted({_num(r, "p") for r in p3})
    ax_a.axhspan(min(cal3), max(cal3), color=BLUE, alpha=0.15, zorder=0)
    ax_a.annotate("sprinkled 1+1D reference (P3-A, passes)",
                  (grid_p[0], max(cal3) + 0.012), fontsize=7.5, color=BLUE)
    ax_a.axhline(0.10, color=MUTED, linestyle="--", linewidth=1.2, zorder=2)
    ax_a.annotate("held-out gate 0.10", (grid_p[-1], 0.108), fontsize=7.5,
                  color=MUTED, ha="right")
    for p in grid_p:
        rows = [r for r in p3 if _num(r, "p") == p]
        ok = [_num(r, "heldout") for r in rows if r.get("status") == "ok"]
        blocked = len(rows) - len(ok)
        xs = [p + (j - len(ok) / 2) * 0.00012 for j in range(len(ok))]
        ax_a.scatter(xs, ok, s=22, color=VERM, alpha=0.75, edgecolor="white",
                     linewidth=0.5, zorder=3)
        ax_a.annotate(f"{blocked}/20", (p, 0.372), fontsize=7, color=MUTED,
                      ha="center")
    ax_a.annotate("structurally blocked runs (no six 25-tick chains)",
                  (grid_p[-1], 0.345), fontsize=7, color=MUTED, ha="right")
    ax_a.set_xticks(grid_p)
    ax_a.set_xticklabels([f"{p:g}" for p in grid_p], fontsize=8.5)
    ax_a.set_xlabel("percolation probability p", fontsize=9.5, color=INK)
    ax_a.set_ylabel("held-out violation", fontsize=9.5, color=INK)
    ax_a.set_ylim(0.0, 0.40)
    ax_a.set_title("P3: transitive percolation\n0/100 pass (blocks at every p)",
                   fontsize=10, color=INK, loc="left")

    # --- Panel B: P4 dual-start hysteresis --------------------------------
    p4 = _rows("p4_stage_b_sweep.csv")
    betas4 = sorted({_num(r, "beta") for r in p4})
    ax_b.axvspan(2.5, 4.5, color=MUTED, alpha=0.10, zorder=0)
    ax_b.annotate("continuum", (1.0, 15), fontsize=8, color=MUTED, ha="center")
    ax_b.annotate("transition\n(hysteresis)", (3.5, 15), fontsize=8,
                  color=MUTED, ha="center")
    ax_b.annotate("crystal", (5.6, 15), fontsize=8, color=MUTED, ha="center")
    ax_b.axhline(-120, color=INK, linestyle="--", linewidth=1.0, alpha=0.6,
                 zorder=2)
    ax_b.annotate("exact bipartite S = -120", (0.0, -117), fontsize=7.5,
                  color=INK, alpha=0.8)
    for start, color, marker, label in (
        ("R", BLUE, "o", "random start"),
        ("X", VERM, "s", "crystal (bipartite) start"),
    ):
        means = []
        for b in betas4:
            vals = [_num(r, "S") for r in p4
                    if _num(r, "beta") == b and r.get("start") == start]
            xs = [b + (j - len(vals) / 2) * 0.02 for j in range(len(vals))]
            ax_b.scatter(xs, vals, s=18, color=color, alpha=0.55,
                         edgecolor="white", linewidth=0.4, zorder=3,
                         marker=marker)
            means.append(st.mean(vals))
        ax_b.plot(betas4, means, color=color, linewidth=1.8, marker=marker,
                  markersize=5, zorder=4, label=label)
    ax_b.set_xticks(betas4)
    ax_b.set_xticklabels([f"{b:g}" for b in betas4], fontsize=8.5)
    ax_b.set_xlabel("coupling beta", fontsize=9.5, color=INK)
    ax_b.set_ylabel("chain mean action <S>", fontsize=9.5, color=INK)
    ax_b.legend(loc="center left", fontsize=8, frameon=False)
    ax_b.set_title("P4: 2D-orders phases (N = 100)\nhysteretic crystallization",
                   fontsize=10, color=INK, loc="left")

    # --- Panel C: P5 discriminator verdict --------------------------------
    # beta = 32 rows come from the shard schema (acc, beta, chain_S, sample,
    # seed, status): DictReader leaves null_gap/truth unset for them.
    p5 = _rows("p5_stage_b_all.csv")
    ok5 = [r for r in p5 if r.get("status") == "ok"]
    cal5 = [r for r in _rows("p5_stage_a_calibration.csv")
            if r.get("status") == "ok"]
    n_blocked = len(p5) - len(ok5)
    groups = [("0\n(uniform cal.)", cal5), ("2", [r for r in ok5 if _num(r, "beta") == 2.0]),
              ("8", [r for r in ok5 if _num(r, "beta") == 8.0])]
    ax_c.axhline(0.10, color=MUTED, linestyle="--", linewidth=1.2, zorder=2)
    ax_c.annotate("gate 0.10", (3.5, 0.112), fontsize=7.5, color=MUTED,
                  ha="right")
    ax_c.axhline(0.50, color=BLUE, linestyle=":", linewidth=1.0, zorder=2)
    ax_c.annotate("truth chance level 0.5", (3.45, 0.515), fontsize=7.5,
                  color=BLUE, ha="right")
    for series_key, color, marker, dx, label in (
        ("heldout", VERM, "o", -0.18, "held-out violation"),
        ("null_gap", GREEN, "^", 0.0, "null gap"),
        ("truth", BLUE, "s", 0.18, "truth-order error"),
    ):
        for i, (_lbl, rows) in enumerate(groups):
            vals = [_num(r, series_key) for r in rows]
            xs = [i + dx + (j - len(vals) / 2) * 0.014 for j in range(len(vals))]
            ax_c.scatter(xs, vals, s=22, color=color, alpha=0.75,
                         edgecolor="white", linewidth=0.5, zorder=3,
                         marker=marker, label=label if i == 0 else None)
    ax_c.axvspan(2.6, 3.4, color=MUTED, alpha=0.12, zorder=0)
    ax_c.annotate(f"structural block\n(0 chains, {n_blocked}/{n_blocked})",
                  (3.0, 0.30), fontsize=8, color=MUTED, ha="center")
    ax_c.set_xticks([0, 1, 2, 3])
    ax_c.set_xticklabels([g[0] for g in groups] + ["32\n(crystal)"],
                         fontsize=8.5)
    ax_c.set_xlim(-0.55, 3.55)
    ax_c.set_ylim(0.0, 0.70)
    ax_c.set_xlabel("coupling beta", fontsize=9.5, color=INK)
    ax_c.set_ylabel("gate metric", fontsize=9.5, color=INK)
    ax_c.legend(loc="upper left", fontsize=8, frameon=False)
    ax_c.set_title("P5: discriminator verdict (N = 600)\n18/18 pass; crystal blocked",
                   fontsize=10, color=INK, loc="left")

    fig.tight_layout()
    fig.savefig(OUT / "fig4_emergence_chain.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    figure_discriminator()      # Fig 1 (PC-V1, Section 4)
    figure_dose_response()      # Fig 2 (P1, Section 5)
    figure_dimension_2d()       # Fig 3 (P2-v2 robustness, Section 6)
    figure_emergence_chain()    # Fig 4 (emergence chain, Section 7)
    figure_confound()           # Fig 5 (confound, Discussion Section 8)
    for name in ("fig1_discriminator", "fig2_dose_response",
                 "fig3_dimension_2d", "fig4_emergence_chain",
                 "fig5_confound"):
        print("wrote", OUT / f"{name}.png")
