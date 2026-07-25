"""P10 Stage A: the calibrated continuum ladder (characterization only).

See docs/prereg/p10_continuum_limit.md Sections 6b-6c. Two arms through
the identical frozen P3 discriminator at every size:

    arm E   smeared-action ensemble at beta = 2, eps N = 12, deep in the
            continuum phase; two chains per N (random + bipartite start),
            burn 100k steps, 10 retained samples per chain at 50k-step
            spacing (both margins set at twice what Stage 0 measured)
    arm S   uniform 2D orders == Poisson sprinklings under the
            null-coordinate dictionary; 20 direct samples per N

No gate and no hypothesis verdict anywhere in this file: Stage A is
characterization, and the mixing screen is applied by the aggregator as
a reporting rule, never as a silent filter.

Usage (each chain is one process, so the six chains parallelize):
    python p10_continuum_ladder.py --arm e --n 900 --start bipartite
    python p10_continuum_ladder.py --arm s
    python p10_continuum_ladder.py --aggregate
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from p3_dynamics import analyze_order
from p5_two_orders_emergence import order_inputs
from pc_common import DEFAULT_OUTPUT_DIR, git_describe, write_rows_csv

from causal_spacetime_lab.positive_control.two_orders import (
    bipartite_perm,
    mcmc_2d_order_fast,
)

BETA = 2.0
EPS_TIMES_N = 12.0
BURN_STEPS = 100_000
SAMPLE_SPACING = 50_000
SAMPLES_PER_CHAIN = 10
TOTAL_STEPS = BURN_STEPS + SAMPLES_PER_CHAIN * SAMPLE_SPACING  # 600k
LADDER = (600, 900, 1200)
STARTS = ("random", "bipartite")

#: Frozen seed table (prereg 6c). Chain seeds also feed the discriminator
#: (seed + sample index) so every fit is reproducible row by row.
CHAIN_SEEDS = {
    (600, "random"): 820, (600, "bipartite"): 830,
    (900, "random"): 840, (900, "bipartite"): 850,
    (1200, "random"): 860, (1200, "bipartite"): 870,
}
ARM_S_SEED_BASE = {600: 900, 900: 920, 1200: 940}
ARM_S_SAMPLES = 20

#: Post-hoc mixing screen (prereg 6c): reported by the aggregator,
#: never silently applied.
MIN_ESS = 10.0


def evaluate_perm(pi: np.ndarray, seed: int) -> dict:
    """One permutation through the frozen discriminator, with P7's G."""

    causal, times, coords = order_inputs(pi)
    row = analyze_order(causal, times, coords, seed=seed, want_truth=True)
    if row.get("status") == "ok":
        m_h = (0.10 - row["heldout"]) / 0.10
        m_n = (row["null_gap"] - 0.10) / 0.10
        m_t = (0.40 - row["truth"]) / 0.40
        row["G"] = float(np.clip(0.5 + min(m_h, m_n, m_t), 0.0, 1.0))
    else:
        row["G"] = 0.0
    return row


def run_chain(n: int, start: str, output_dir: Path) -> None:
    seed = CHAIN_SEEDS[(n, start)]
    eps = EPS_TIMES_N / n
    if start == "bipartite":
        pi0 = bipartite_perm(n)
    else:
        pi0 = np.random.default_rng(seed + 1).permutation(n)

    samples, acceptance, perms = mcmc_2d_order_fast(
        pi0, beta=BETA, eps=eps, steps=TOTAL_STEPS, seed=seed,
        sample_every=SAMPLE_SPACING, burn_frac=BURN_STEPS / TOTAL_STEPS,
        collect_perms=True,
    )
    rows = []
    for k, (obs, pi) in enumerate(zip(samples, perms, strict=True)):
        row = {
            "arm": "E", "n": float(n), "start": start,
            "chain_seed": float(seed), "sample_index": float(k),
            "step": float(BURN_STEPS + k * SAMPLE_SPACING),
            "n0": float(obs["n0"]), "action": float(obs["S"]),
            "acceptance": float(acceptance),
            "code_version": git_describe(),
            **evaluate_perm(pi, seed=seed + 7 * k),
        }
        rows.append(row)
        flag = row.get("status", "?")
        print(f"E n={n} {start} sample {k}: {flag} "
              f"G={row.get('G'):.3f} n0={row['n0']:.0f}", flush=True)
    write_rows_csv(output_dir / f"p10_e_n{n}_{start}.csv", rows)


def run_arm_s(output_dir: Path) -> None:
    rows = []
    for n in LADDER:
        base = ARM_S_SEED_BASE[n]
        for k in range(ARM_S_SAMPLES):
            rng = np.random.default_rng(base + k)
            pi = rng.permutation(n)
            row = {
                "arm": "S", "n": float(n), "start": "uniform",
                "chain_seed": float(base + k), "sample_index": float(k),
                "code_version": git_describe(),
                **evaluate_perm(pi, seed=base + k),
            }
            rows.append(row)
        done = sum(1 for r in rows if r["n"] == n and r.get("status") == "ok")
        print(f"S n={n}: {done}/{ARM_S_SAMPLES} ok", flush=True)
    write_rows_csv(output_dir / "p10_s.csv", rows)


def _ess(values: list[float]) -> float:
    """Geyer initial-positive-pair ESS, matching P7's convention."""

    x = np.asarray(values, dtype=float)
    x = x - x.mean()
    m = x.size
    if m < 4 or float(np.dot(x, x)) == 0.0:
        return float(m)
    denominator = float(np.dot(x, x))
    autocorr = [
        float(np.dot(x[: m - lag], x[lag:]) / denominator)
        for lag in range(m)
    ]
    total = 0.0
    for k in range(1, m // 2):
        pair = autocorr[2 * k - 1] + autocorr[2 * k]
        if pair <= 0.0:
            break
        total += pair
    tau = 1.0 + 2.0 * total
    return float(m / tau)


def _median_and_bootstrap(values, seed, draws=4000):
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        return float("nan"), (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    medians = np.median(
        x[rng.integers(0, x.size, size=(draws, x.size))], axis=1
    )
    return float(np.median(x)), (
        float(np.percentile(medians, 2.5)),
        float(np.percentile(medians, 97.5)),
    )


def aggregate(output_dir: Path) -> None:
    import csv

    rows = []
    for path in sorted(output_dir.glob("p10_e_n*.csv")) + [output_dir / "p10_s.csv"]:
        if path.exists():
            rows.extend(csv.DictReader(path.open(encoding="utf-8")))
    summary: dict = {"code_version": git_describe(), "per_n": {}, "chains": []}

    # mixing screen, reported per chain (prereg 6c)
    for n in LADDER:
        chains = {}
        for r in rows:
            if r["arm"] == "E" and int(float(r["n"])) == n:
                chains.setdefault(r["start"], []).append(r)
        if not chains:
            continue
        random_tail = [float(r["n0"]) for r in chains.get("random", [])[5:]]
        band = (min(random_tail), max(random_tail)) if random_tail else None
        for start, chain_rows in chains.items():
            chain_rows.sort(key=lambda r: float(r["sample_index"]))
            n0s = [float(r["n0"]) for r in chain_rows]
            ess = _ess(n0s)
            melted = (
                band is not None
                and band[0] * 0.9 <= n0s[0] <= band[1] * 1.1
            ) if start == "bipartite" else True
            summary["chains"].append({
                "n": n, "start": start, "ess_n0": round(ess, 1),
                "ess_pass": bool(ess >= MIN_ESS),
                "first_sample_in_random_band": bool(melted),
                "acceptance": float(chain_rows[0]["acceptance"]),
            })

    for n in LADDER:
        per_arm = {}
        for arm in ("E", "S"):
            ok = [
                r for r in rows
                if r["arm"] == arm and int(float(r["n"])) == n
                and r.get("status") == "ok"
            ]
            entry = {"n_ok": len(ok), "n_total": sum(
                1 for r in rows if r["arm"] == arm and int(float(r["n"])) == n
            )}
            for key in ("heldout", "null_gap", "truth", "G"):
                med, ci = _median_and_bootstrap(
                    [float(r[key]) for r in ok], seed=hash((n, arm, key)) % 2**31
                )
                entry[f"median_{key}"] = med
                entry[f"ci_{key}"] = ci
            per_arm[arm] = entry
        summary["per_n"][str(n)] = per_arm

    (output_dir / "p10_stage_a_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", choices=["e", "s"], default=None)
    parser.add_argument("--n", type=int, choices=LADDER, default=None)
    parser.add_argument("--start", choices=STARTS, default=None)
    parser.add_argument("--aggregate", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.aggregate:
        aggregate(args.output_dir)
    elif args.arm == "s":
        run_arm_s(args.output_dir)
    elif args.arm == "e":
        if args.n is None or args.start is None:
            raise SystemExit("--arm e requires --n and --start")
        run_chain(args.n, args.start, args.output_dir)
    else:
        raise SystemExit("choose --arm e/--arm s/--aggregate")


if __name__ == "__main__":
    main()
