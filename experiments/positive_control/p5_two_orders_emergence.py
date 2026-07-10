"""P5: does the ACTION-WEIGHTED continuum phase pass the frozen geometry
discriminator? The first instrument-based judgment of an equilibrium quantum-
gravity-style ensemble in this program.

Design: 2D orders at N=600 (the scale where the frozen P3 protocol -- 6
disjoint chains >= 25 ticks, >= 20 bracketed targets -- is structurally
feasible: 6/6 uniform seeds), smeared action at eps = 0.02 (eps*N = 12, the
same dimensionless combination as P4 at N=100, so the crystal's action
advantage is held fixed while entropy grows).

Stages:
  a      calibration on uniform 2D orders (== sprinkled diamonds): the
         order-intrinsic discriminator with TRUE lightcone coordinates
         (t = i + pi_i, x = i - pi_i) must pass and recover x; column-shuffle
         must block. Also the N=600 profile reference for recon phase gates.
  recon  fast-sampler beta sweep to locate the crystallization scale beta_c
         at N=600 (profile-based, pre-freeze).
  b      confirmatory: equilibrium samples at frozen continuum betas and a
         crystal control beta; each sampled configuration is judged by the
         frozen discriminator gates.
  verdict aggregate stage-b rows into the decision registry.

Reuses the frozen P3 analysis pipeline (analyze_order) verbatim.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from p3_dynamics import analyze_order, cohens_d, _median  # frozen P3 pipeline
from causal_spacetime_lab.positive_control.two_orders import (
    bipartite_perm,
    chain_observables,
    mcmc_2d_order_fast,
    perm_to_causal_matrix,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "positive_control"
FROZEN = ROOT / "docs" / "prereg" / "frozen" / "p5_test_constants.json"

N_ELEMENTS = 600
EPS = 0.02


def order_inputs(pi: np.ndarray):
    """Causal matrix plus TRUE lightcone time/space coordinates of a 2D order."""
    idx = np.arange(pi.size, dtype=float)
    return (
        perm_to_causal_matrix(pi),
        idx + pi,          # t (topological time)
        idx - pi,          # x (true spatial coordinate, order-level truth)
    )


def stage_a() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for seed in range(10):
        rng = np.random.default_rng(seed)
        pi = rng.permutation(N_ELEMENTS)
        causal, times, coords = order_inputs(pi)
        row = analyze_order(causal, times, coords, seed=seed, want_truth=True)
        row["seed"] = seed
        rows.append(row)
        print(f"seed {seed}: {row}", flush=True)
    ok = [r for r in rows if r["status"] == "ok"]
    reference = {
        key: float(np.mean([chain_observables(
            perm_to_causal_matrix(np.random.default_rng(1000 + s).permutation(N_ELEMENTS))
        )[key] for s in range(20)]))
        for key in ("n0", "n1", "n2", "height")
    }
    summary = {
        "n_elements": N_ELEMENTS,
        "eps": EPS,
        "n_valid": len(ok),
        "median_heldout": _median([r["heldout"] for r in ok]),
        "max_heldout": max((r["heldout"] for r in ok), default=None),
        "median_null_gap": _median([r["null_gap"] for r in ok]),
        "min_null_gap": min((r["null_gap"] for r in ok), default=None),
        "median_truth": _median([r["truth"] for r in ok]),
        "profile_reference": reference,
    }
    with open(OUT / "p5_stage_a_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    keys = sorted({k for r in rows for k in r})
    with open(OUT / "p5_stage_a.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary, indent=2))


def stage_recon(betas: list[float], steps: int) -> None:
    ref = json.loads((OUT / "p5_stage_a_summary.json").read_text())["profile_reference"]
    n12_ref = ref["n1"] + ref["n2"]
    print(f"reference: {ref}")
    print("beta  start  <S>        n0      n12     height  acc     phase", flush=True)
    for beta in betas:
        for tag, pi0 in (
            ("R", np.random.default_rng(50 + int(beta * 10)).permutation(N_ELEMENTS)),
            ("X", bipartite_perm(N_ELEMENTS)),
        ):
            samples, acc, _ = mcmc_2d_order_fast(
                pi0, beta=beta, eps=EPS, steps=steps,
                seed=60 + int(beta * 10) * 2 + (tag == "X"),
                sample_every=max(steps // 40, 1),
            )
            m = {k: float(np.mean([s[k] for s in samples]))
                 for k in ("S", "n0", "n1", "n2", "height")}
            n12 = m["n1"] + m["n2"]
            if abs(m["n0"] / ref["n0"] - 1.0) <= 0.5 and m["height"] >= 0.7 * ref["height"]:
                phase = "continuum"
            elif n12 <= 0.5 * n12_ref and m["height"] <= 0.25 * ref["height"]:
                phase = "crystal"
            else:
                phase = "other"
            print(
                f"{beta:5.2f}  {tag}   {m['S']:+9.1f}  {m['n0']:6.0f}  {n12:6.0f}"
                f"  {m['height']:6.1f}  {acc:.3f}  {phase}",
                flush=True,
            )


def stage_b(betas: list[float], seeds: list[int]) -> None:
    frozen = json.loads(FROZEN.read_text())
    rows = []
    tag = "-".join(f"{b:g}" for b in betas) + f"_s{seeds[0]}-{seeds[-1]}"
    for beta in betas:
        crystal = beta == frozen["crystal_beta"]
        steps = frozen["steps_crystal" if crystal else "steps_continuum"]
        sample_every = frozen[
            "sample_every_crystal" if crystal else "sample_every_continuum"
        ]
        for seed in seeds:
            if crystal:
                pi0 = bipartite_perm(N_ELEMENTS)
            else:
                rng = np.random.default_rng(seed * 31 + int(round(beta * 10)))
                pi0 = rng.permutation(N_ELEMENTS)
            samples, acc, perms = mcmc_2d_order_fast(
                pi0, beta=beta, eps=EPS, steps=steps,
                seed=seed * 1000 + int(round(beta * 10)),
                sample_every=sample_every,
                burn_frac=frozen["burn_frac"], collect_perms=True,
            )
            mean_s = float(np.mean([s["S"] for s in samples]))
            print(f"beta={beta:g} seed={seed}: chain done <S>={mean_s:+.1f} "
                  f"acc={acc:.3f}, judging {len(perms)} configs", flush=True)
            for k, pi in enumerate(perms):
                causal, times, coords = order_inputs(pi)
                result = analyze_order(
                    causal, times, coords,
                    seed=seed * 100 + k, want_truth=True,
                )
                row = {"beta": beta, "seed": seed, "sample": k,
                       "chain_S": mean_s, "acc": float(acc), **result}
                rows.append(row)
                print(f"  config {k}: {result}", flush=True)
    keys = sorted({k for r in rows for k in r})
    with open(OUT / f"p5_stage_b_{tag}.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote p5_stage_b_{tag}.csv")


def verdict() -> None:
    frozen = json.loads(FROZEN.read_text())
    rows = []
    for path in sorted(OUT.glob("p5_stage_b_*.csv")):
        with open(path) as fh:
            rows.extend(list(csv.DictReader(fh)))
    registry = {"frozen_constants": FROZEN.name, "verdict_by_beta": {}}
    for beta in sorted({float(r["beta"]) for r in rows}):
        sub = [r for r in rows if float(r["beta"]) == beta]
        ok = [r for r in sub if r["status"] == "ok"]
        passes = [
            r for r in ok
            if float(r["heldout"]) <= frozen["gate_heldout"]
            and float(r["null_gap"]) >= frozen["gate_nullgap"]
            and float(r["truth"]) <= frozen["gate_truth"]
        ]
        expected = frozen["expectations"][f"{beta:g}"]
        n_pass, n_total = len(passes), len(sub)
        if expected == "pass":
            met = n_pass >= frozen["pass_min_fraction"] * n_total
        else:
            met = n_pass == 0
        registry["verdict_by_beta"][f"{beta:g}"] = {
            "n_configs": n_total,
            "n_structural_block": n_total - len(ok),
            "n_gate_pass": n_pass,
            "median_heldout": _median([float(r["heldout"]) for r in ok]),
            "median_null_gap": _median([float(r["null_gap"]) for r in ok]),
            "median_truth": _median([float(r["truth"]) for r in ok]),
            "expected": expected,
            "expectation_met": bool(met),
        }
        print(f"beta={beta:g}: {registry['verdict_by_beta'][f'{beta:g}']}")
    registry["all_expectations_met"] = all(
        v["expectation_met"] for v in registry["verdict_by_beta"].values()
    )
    with open(OUT / "p5_stage_b_decision_registry.json", "w") as fh:
        json.dump(registry, fh, indent=2)
    print(f"all expectations met: {registry['all_expectations_met']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True,
                        choices=["a", "recon", "b", "verdict"])
    parser.add_argument("--betas", default="")
    parser.add_argument("--seeds", default="100-102")
    parser.add_argument("--steps", type=int, default=800_000)
    args = parser.parse_args()
    if args.stage == "a":
        stage_a()
    elif args.stage == "recon":
        betas = [float(b) for b in (args.betas or "1,2,4,8,16").split(",")]
        stage_recon(betas, args.steps)
    elif args.stage == "b":
        betas = [float(b) for b in args.betas.split(",")]
        lo, hi = args.seeds.split("-")
        stage_b(betas, list(range(int(lo), int(hi) + 1)))
    else:
        verdict()
