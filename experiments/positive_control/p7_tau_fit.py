"""Reproduce the P7 tunneling-exponent fit from the checked-in dataset.

Reads `docs/p7_tunneling_data/runs_combined.csv` -- the complete record of
every tau measurement across all four run batches -- and prints the exponent,
its bootstrap interval, and the N=600 cost projections quoted in
`docs/p7_enhanced_sampling.md`. Deterministic; no simulation.

Deduplication is the load-bearing step. All batches derived their RNG seed
from the same formula (base + 1000*rep + n), so a (n, rep) pair names ONE
random-walk trajectory: batches that re-ran it merely observed the same walk
to a different stopping point (batch B's N=30 rep=0 at 25 round trips and
batch C's at 15 are prefixes of one trajectory). An earlier headline fit
treated all 25 clean rows as independent -- 9 of them were re-measurements,
which biased the point estimate and made the bootstrap interval artificially
narrow. Here each walk contributes exactly one clean measurement: the one
with the most completed round trips.

Censored rows (budget exhausted before the target) never enter the fit; the
run-level bound reported for them is moves/target, since moves/completed
overstates the bound by including the unfinished traversal.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "docs" / "p7_tunneling_data" / "runs_combined.csv"

HALVINGS_TO_CONVERGE = 17
NUMBA_US_PER_MOVE = 135.0
BOOTSTRAP_DRAWS = 10_000
BOOTSTRAP_SEED = 0


def load_walks() -> tuple[list[tuple[int, float]], list[dict]]:
    """One clean tau per independent walk, plus the censored rows.

    Returns (clean, censored) where clean is a list of (n, tau) with tau from
    the most-round-trips clean measurement of each (n, rep) walk, and censored
    holds walks with no clean measurement at all.
    """

    with open(DATA, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    walks: dict[tuple[int, int], list[dict]] = {}
    for row in rows:
        key = (int(row["n_elements"]), int(row["rep"]))
        walks.setdefault(key, []).append(row)

    clean: list[tuple[int, float]] = []
    censored: list[dict] = []
    for (n, _rep), group in sorted(walks.items()):
        candidates = [g for g in group if g["censored"] == "False"]
        if candidates:
            best = max(candidates, key=lambda g: int(g["round_trips"]))
            clean.append((n, int(best["moves"]) / int(best["round_trips"])))
        else:
            best = max(group, key=lambda g: int(g["moves"]))
            censored.append(best)
    return clean, censored


def main() -> None:
    clean, censored = load_walks()
    print(f"independent walks: {len(clean)} clean, {len(censored)} censored")

    by_size: dict[int, list[float]] = {}
    for n, tau in clean:
        by_size.setdefault(n, []).append(tau)
    print("\nper-N tau (moves per round trip), one measurement per walk:")
    for n in sorted(by_size):
        taus = np.array(by_size[n])
        print(
            f"  N={n:4d}  walks={taus.size}  mean={taus.mean():,.0f}  "
            f"min={taus.min():,.0f}  max={taus.max():,.0f}"
        )
    for row in censored:
        bound = int(row["moves"]) / int(row["target_round_trips"])
        print(
            f"  N={row['n_elements']:>4}  rep={row['rep']} censored: "
            f"tau > {bound:,.0f} "
            f"({row['round_trips']}/{row['target_round_trips']} trips in "
            f"{int(row['moves']):,} moves)"
        )

    log_n = np.log([n for n, _ in clean])
    log_tau = np.log([tau for _, tau in clean])
    exponent, _ = (float(v) for v in np.polyfit(log_n, log_tau, 1))

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    boots = []
    for _ in range(BOOTSTRAP_DRAWS):
        pick = rng.integers(0, log_n.size, log_n.size)
        if np.unique(log_n[pick]).size < 2:
            continue
        boots.append(np.polyfit(log_n[pick], log_tau[pick], 1)[0])
    low, high = (float(v) for v in np.percentile(boots, [5, 95]))

    print(f"\ntau ~ N^{exponent:.2f}")
    print(f"90% bootstrap interval on the exponent: [{low:.2f}, {high:.2f}]")

    print(f"\nN=600 converged ln_g ({HALVINGS_TO_CONVERGE} round trips), "
          f"intercept refit per alpha:")
    for label, alpha in (("low", low), ("point", exponent), ("high", high)):
        # Least-squares intercept given the fixed slope; carrying the
        # point-fit intercept to the interval endpoints misstates costs.
        refit_intercept = float(np.mean(log_tau - alpha * log_n))
        moves = np.exp(refit_intercept) * 600.0**alpha * HALVINGS_TO_CONVERGE
        hours = moves * NUMBA_US_PER_MOVE * 1e-6 / 3600
        print(
            f"  alpha={alpha:.2f} ({label:>5}): {moves:.2e} moves = "
            f"{hours:,.0f} core-hours (Numba kernel, "
            f"{NUMBA_US_PER_MOVE:.0f} us/move)"
        )

    print(
        "\nPlanning numbers, not results: the lever arm is N=30..80 and the "
        "extrapolation multiplies any exponent error by log(600/80) in the "
        "log-cost. Censored rows (all high values at large N) are excluded, "
        "so the fit errs in the method's favour; window edge-touch flags were "
        "not recorded for these batches (limitation noted in the verdict doc)."
    )


if __name__ == "__main__":
    main()
