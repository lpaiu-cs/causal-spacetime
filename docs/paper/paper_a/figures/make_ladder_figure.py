"""Regenerate the Section 6.8 mass-ladder figure from the committed artifacts.

Reads the four per-rung count artifacts directly -- the same
`docs/prereg/p14_*_count.json` files whose figures
`tests/test_paper_a_count_integration.py` re-derives against the manuscript
table -- so the figure cannot disagree with the text. Nothing is typed in
here: every plotted number is parsed from an artifact.

Palette: Okabe-Ito (colorblind-safe), matching make_figures.py.

Usage: python docs/paper/paper_a/figures/make_ladder_figure.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PREREG = Path("docs/prereg")
OUT = Path("docs/paper/paper_a/figures")

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREY = "#555555"
LIGHT = "#CCCCCC"

#: (mu as printed in the manuscript, artifact) -- the executed ladder
RUNGS = (
    ("0.1333", "p14_o5_count.json"),
    ("0.1867", "p14_s6_m14_count.json"),
    ("0.2400", "p14_s6_m18_count.json"),
    ("0.4000", "p14_s6_m30_count.json"),
)


def load() -> list[dict]:
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


def plotted_digest(rows: list[dict]) -> str:
    """A canonical digest of exactly the values the figure draws.

    Stamped into the PNG and re-derived by the contract test, so a stale
    PNG -- one generated before a rung was added or an artifact changed --
    is detectable from the displayed asset itself, not merely from the
    generator source.
    """

    payload = json.dumps(
        [[r["mu"], r["v_lo"], r["v_hi"], r["c_lo"], r["c_hi"],
          r["d_lo"], r["d_hi"], r["band"], r["k"], r["u"], r["verdict"]]
         for r in rows],
        separators=(",", ":"), sort_keys=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def main() -> None:
    rows = load()
    assert all(r["verdict"] == "CONCORDANT" for r in rows), "unexpected verdict"

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(9.4, 3.9), gridspec_kw={"width_ratios": [1.05, 1]})

    # ---- panel A: certified volume vs the count-derived volume interval.
    # At this scale the two are indistinguishable -- which IS the result;
    # panel B resolves the residual. The trend line carries the geometry.
    ax1.plot([r["mu"] for r in rows],
             [(r["v_lo"] + r["v_hi"]) / 2 for r in rows],
             color=LIGHT, lw=1.2, zorder=0)
    for i, r in enumerate(rows):
        x = r["mu"]
        ax1.plot([x - 0.005] * 2, [r["v_lo"], r["v_hi"]], color=BLUE, lw=7,
                 solid_capstyle="butt",
                 label="certified $V$ (oracle)" if i == 0 else None)
        ax1.plot([x + 0.005] * 2, [r["c_lo"], r["c_hi"]], color=ORANGE, lw=7,
                 solid_capstyle="butt",
                 label="$C$ = count / intensity" if i == 0 else None)
        ax1.annotate(f'{(r["v_lo"]+r["v_hi"])/2:.1f}', (x, r["v_hi"]),
                     textcoords="offset points", xytext=(0, 7),
                     ha="center", fontsize=8, color=GREY)
    ax1.text(0.135, 112,
             "the two overlap at this scale:\nthe residual is panel B",
             fontsize=8, color=GREY, style="italic", va="top")
    ax1.set_xlabel(r"compactness  $\mu = 2M/r_c$")
    ax1.set_ylabel(r"4-volume of the fixed diamond")
    ax1.set_title("A. the count tracks the certified volume over a 3x span in $\mu$",
                  fontsize=10, loc="left")
    ax1.set_xticks([r["mu"] for r in rows])
    ax1.set_xticklabels([f'{r["mu"]:.4f}' for r in rows])
    ax1.set_xlim(0.10, 0.44)
    ax1.grid(axis="y", color=LIGHT, lw=.6)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=8, loc="upper left")

    # ---- panel B: identified discrepancy, normalized by each rung's band
    ys = list(range(len(rows)))[::-1]
    ax2.axvspan(-1, 1, color=LIGHT, alpha=.45, lw=0)
    ax2.axvline(0, color=GREY, lw=.8)
    for y, r in zip(ys, rows, strict=True):
        lo, hi = r["d_lo"] / r["band"], r["d_hi"] / r["band"]
        ax2.plot([lo, hi], [y, y], color=BLUE, lw=7, solid_capstyle="butt")
        ax2.plot([lo, hi], [y, y], "|", color=BLUE, ms=9, mew=1.6)
    ax2.axvline(-1, color=GREY, ls="--", lw=1)
    ax2.axvline(1, color=GREY, ls="--", lw=1)
    ax2.set_yticks(ys)
    ax2.set_yticklabels(
        [rf'$\mu$ = {r["mu"]:.4f}' + "\n" + f'K = {r["k"]:,}, U = {r["u"]}'
         for r in rows], fontsize=8)
    ax2.set_xlabel(r"identified discrepancy $D$, in units of the "
                   r"band $B=\tau V_{\rm ref}$")
    ax2.set_title(r"B. equivalence gate ($\tau$ = 2.5%): all four contained",
                  fontsize=10, loc="left")
    ax2.set_xlim(-1.35, 1.35)
    ax2.set_ylim(-0.6, len(rows) - 0.4)
    ax2.text(1.0, len(rows) - 0.52, "  gate", fontsize=8, color=GREY,
             va="center")
    ax2.grid(axis="x", color=LIGHT, lw=.6)
    ax2.set_axisbelow(True)

    for ax in (ax1, ax2):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    fig.tight_layout()
    dest = OUT / "fig4_ladder_count.png"
    stamp = plotted_digest(rows)
    fig.savefig(dest, dpi=220, metadata={"LadderDigest": stamp})
    print(f"wrote {dest}")
    print(f"  LadderDigest {stamp}")
    for r in rows:
        print(f'  mu={r["mu"]:.4f}  V=[{r["v_lo"]:.4f}, {r["v_hi"]:.4f}]  '
              f'D/B=[{r["d_lo"]/r["band"]:+.3f}, {r["d_hi"]/r["band"]:+.3f}]  '
              f'{r["verdict"]}')


if __name__ == "__main__":
    main()
