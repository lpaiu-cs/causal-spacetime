# ruff: noqa
# Referee's scratch script for review.md, archived verbatim apart from these three lines.
# Written in this folder by the fresh-context AI referee of local_null_challenge_r2.md.
"""Referee check r5: exact counts behind N9's printed fractions.

Section 6 of the revision says: "on W_lam, all sampled limit boxes are inside
the window only at lambda = 2^-17 (q = 1) and 2^-20 (q = 0.25)". N9 prints
the fractions with three decimals and gates only lambda = 2^-20, so a
printed 1.000 at 2^-17 could be 1999/2000. This script re-runs N9's own
seeded sampling (imports checks.py; same seeds, same order of draws) and
prints the exact counts. It also prints, per pair, the distance of the limit
c-range from the window edges, to show where the misses are.

  E1  information: exact inside-counts at each lambda, both q.
Run: python r5_n9_exact_counts.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import checks as C  # noqa: E402


def main():
    for q in (0.25, 1.0):
        geo = C.Geo(q)
        p0q = float(C.p0_fast(q))
        rng = np.random.default_rng(90 + int(4 * q))
        for lam in (2.0**-11, 2.0**-14, 2.0**-17, 2.0**-20):
            pairs = C.sample_pairs(geo, lam, 0.37, 2000, rng, False)
            X = np.array([C.theta0(geo, u1, vt1, z1) for u1, vt1, z1, *_ in pairs])
            Y = np.array([C.theta0(geo, u2, vt2, z2) for *_, u2, vt2, z2 in pairs])
            amin, amax = np.minimum(X[:, 2], Y[:, 2]), np.maximum(X[:, 2], Y[:, 2])
            pmin, pmax = np.minimum(X[:, 0], Y[:, 0]), np.maximum(X[:, 0], Y[:, 0])
            lo, hi = amin / pmax**2, amax / pmin**2
            inside = (lo > p0q) & (hi < 0.1)
            worst = np.min(np.minimum(lo - p0q, 0.1 - hi))
            print(
                f"[INFO] E1 q={q} lam=2^{int(round(np.log2(lam)))}: inside {int(inside.sum())}/2000; "
                f"least signed margin {worst:+.2e}",
                flush=True,
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
