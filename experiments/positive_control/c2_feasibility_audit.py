"""C2-unpaired Schwarzschild feasibility audit from STORED raw arrays
(S3 official + S4 campaign) -- zero solver calls, exploratory grade.

Marginal f_M / f_0 samples per reading are iid across readings (each
reading is an independent sprinkle); pairing only correlates the two
WITHIN a reading, so the marginal means/SDs are valid unpaired-arm
estimates and the within-reading correlation itself quantifies the
paired design's advantage. Empirical AUC uses off-diagonal cross
pairs (the diagonal is the paired probability, excluded); bootstrap
resamples READINGS jointly to respect the pairing present in this
estimation data. BA here is the descriptive midpoint-threshold value;
a preregistered C2 design would freeze a train/test split instead.

Bootstrap stream [781, 4841]: documented design lineage, never a
campaign stream (the in-session audit of 2026-08-10 used [781, 4850];
this committed rerun draws its own stream so neither replays the
other's draws -- point estimates are data-determined and identical,
bootstrap CIs agree to resampling noise).

Output: docs/prereg/p14_c2_feasibility_audit.json -- the sizing
anchor a C2 preregistration would cite. Per the review ruling, the
prereg must fix its own gates, CP-certified branch powers, the AUC CI
method, a BA train/test protocol, and anchor on the CONSERVATIVE
block (S3, d = 2.598), not the pooled estimate.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[2]
_OUT = _REPO / "docs" / "prereg" / "p14_c2_feasibility_audit.json"

BOOT_SEED = (781, 4841)
BOOT_B = 4000
COST_US_PER_PAIR = 768.0     # S1 priced operating point
PAIRS_PER_READING = 300 * 299 / 2


def _auc_offdiag(fm: np.ndarray, f0: np.ndarray) -> float:
    u = ((fm[:, None] < f0[None, :]).astype(float)
         + (fm[:, None] == f0[None, :]) * 0.5)
    np.fill_diagonal(u, np.nan)
    return float(np.nanmean(u))


def _ba_midpoint(fm: np.ndarray, f0: np.ndarray) -> float:
    thr = 0.5 * (fm.mean() + f0.mean())
    return 0.5 * (float((fm < thr).mean()) + float((f0 >= thr).mean()))


def _block_stats(fm: np.ndarray, f0: np.ndarray,
                 rng: np.random.Generator) -> dict:
    n = len(fm)
    pooled = math.sqrt(0.5 * (fm.std(ddof=1) ** 2 + f0.std(ddof=1) ** 2))
    boots = np.empty((BOOT_B, 3))
    for b in range(BOOT_B):
        idx = rng.integers(0, n, n)
        bm, b0 = fm[idx], f0[idx]
        bp = math.sqrt(0.5 * (bm.std(ddof=1) ** 2 + b0.std(ddof=1) ** 2))
        boots[b] = ((b0.mean() - bm.mean()) / bp, _auc_offdiag(bm, b0),
                    _ba_midpoint(bm, b0))
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
    return {
        "n_readings": n,
        "f_M": {"mean": float(fm.mean()), "sd": float(fm.std(ddof=1))},
        "f_0": {"mean": float(f0.mean()), "sd": float(f0.std(ddof=1))},
        "within_reading_corr": float(np.corrcoef(fm, f0)[0, 1]),
        "unpaired_d": {"point": float((f0.mean() - fm.mean()) / pooled),
                       "ci95_bootstrap": [float(lo[0]), float(hi[0])]},
        "auc_offdiag": {"point": _auc_offdiag(fm, f0),
                        "ci95_bootstrap": [float(lo[1]), float(hi[1])]},
        "ba_midpoint_descriptive": {"point": _ba_midpoint(fm, f0),
                                    "ci95_bootstrap": [float(lo[2]),
                                                       float(hi[2])]},
    }


def main() -> None:
    blocks = {}
    for name, path in (("S3", "p14_s3_probe_results.json"),
                       ("S4", "p14_s4_results.json")):
        r = json.loads((_REPO / "docs" / "prereg" / path)
                       .read_text(encoding="utf-8"))
        blocks[name] = (
            np.array(r["f_schwarzschild_lower"]["per_reading"]),
            np.array(r["f_flat"]["per_reading"]))

    rng = np.random.default_rng([*BOOT_SEED])
    stats = {name: _block_stats(fm, f0, rng)
             for name, (fm, f0) in blocks.items()}

    fm = np.concatenate([blocks["S3"][0], blocks["S4"][0]])
    f0 = np.concatenate([blocks["S3"][1], blocks["S4"][1]])
    pooled_sd = math.sqrt(0.5 * (fm.std(ddof=1) ** 2
                                 + f0.std(ddof=1) ** 2))
    a = _auc_offdiag(fm, f0)
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)

    def auc_se(n: int) -> float:
        return math.sqrt((a * (1 - a) + (n - 1) * (q1 - a * a)
                          + (n - 1) * (q2 - a * a)) / (n * n))

    d = float((f0.mean() - fm.mean()) / pooled_sd)
    sizing = []
    for n in (300, 600, 1200, 2400, 4800):
        sizing.append({
            "n_per_arm": n,
            "auc_ci95_lower_hanley_mcneil": a - 1.96 * auc_se(n),
            "d_ci95_half_width": 1.96 * math.sqrt(2.0 / n
                                                  + d * d / (4 * n)),
            "curved_arm_hours": n * PAIRS_PER_READING
            * COST_US_PER_PAIR / 1e6 / 3600,
        })

    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=_REPO,
                         capture_output=True, text=True,
                         check=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=_REPO,
                           capture_output=True, text=True,
                           check=True).stdout.strip()

    result = {
        "audit": "C2-unpaired Schwarzschild feasibility (stored-data, "
                 "exploratory; zero solver calls)",
        "inputs": ["docs/prereg/p14_s3_probe_results.json",
                   "docs/prereg/p14_s4_results.json"],
        "blocks": stats,
        "pooled": {"d": d, "auc_offdiag": a,
                   "note": "pooled shown for context; a prereg anchors "
                           "on the conservative block (S3)"},
        "conservative_anchor_block": "S3",
        "sizing_table": sizing,
        "verdict": ("computationally feasible; separation is strong "
                    "but imperfect (unlike the plane-wave C2's "
                    "complete separation) -- a limitation statement, "
                    "not a claim, absent a dedicated AUC-upper gate"),
        "params": {"bootstrap_seed": list(BOOT_SEED), "bootstrap_b": BOOT_B,
                   "cost_us_per_pair": COST_US_PER_PAIR},
        "code": {"rev": rev, "dirty": bool(dirty)},
    }
    _OUT.write_text(json.dumps(result, indent=2) + "\n",
                    encoding="utf-8", newline="\n")
    print(f"wrote {_OUT}")
    for name in ("S3", "S4"):
        s = stats[name]
        print(f"{name}: d {s['unpaired_d']['point']:.3f}  "
              f"AUC {s['auc_offdiag']['point']:.4f}  "
              f"corr {s['within_reading_corr']:.4f}")
    print(f"pooled d {d:.3f}  AUC {a:.4f}")


if __name__ == "__main__":
    main()
