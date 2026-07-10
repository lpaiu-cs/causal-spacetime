"""P4: action-weighted emergence in the restricted 2D-orders ensemble.

Stage A (calibration, pre-freeze): sprinkled-diamond reference profile at the
experiment scale plus sampler validation gates (exact Gibbs at N=6; ensemble
identity uniform-permutations == sprinkled diamond). Constants frozen into
docs/prereg/frozen/p4_test_constants.json.

Stage B (confirmatory, post-freeze): multi-seed dual-start beta sweep with the
frozen phase-classification gates. Run shards via --betas.

Usage:
  python p4_action_emergence.py --stage a
  python p4_action_emergence.py --stage b --betas 0,1 --seeds 100-104
"""

from __future__ import annotations

import argparse
import csv
import json
from itertools import permutations
from pathlib import Path

import numpy as np

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.positive_control.action import smeared_action_2d
from causal_spacetime_lab.positive_control.two_orders import (
    bipartite_perm,
    chain_observables,
    mcmc_2d_order,
    perm_to_causal_matrix,
)
from causal_spacetime_lab.sprinkling import sprinkle_minkowski_causal_diamond

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "positive_control"
FROZEN = ROOT / "docs" / "prereg" / "frozen" / "p4_test_constants.json"

N_ELEMENTS = 100
EPS = 0.12
STEPS = 400_000
SAMPLE_EVERY = 1000


def sprinkled_diamond(n: int, seed: int) -> np.ndarray:
    events = sprinkle_minkowski_causal_diamond(n, spacetime_dim=2, T=2.0, seed=seed)
    C = np.array(causal_matrix_1p1(events), dtype=bool)
    np.fill_diagonal(C, False)
    return C


def stage_a() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    # Reference profile: sprinkled 2D diamond at the experiment scale.
    rows = []
    for seed in range(30):
        obs = chain_observables(sprinkled_diamond(N_ELEMENTS, seed))
        obs["seed"] = seed
        rows.append(obs)
    ref = {
        key: float(np.mean([r[key] for r in rows]))
        for key in ("R", "n0", "n1", "n2", "height")
    }

    # Gate V1: exact Gibbs at N=6 (720 permutations).
    n, beta, eps = 6, 0.4, 0.3
    perms = [np.array(p) for p in permutations(range(n))]
    actions = np.array([smeared_action_2d(perm_to_causal_matrix(p), eps) for p in perms])
    weights = np.exp(-beta * (actions - actions.min()))
    target = weights / weights.sum()
    index = {p.tobytes(): k for k, p in enumerate(perms)}
    rng = np.random.default_rng(1)
    pi = perms[0].copy()
    action = smeared_action_2d(perm_to_causal_matrix(pi), eps)
    steps = 600_000
    counts = np.zeros(len(perms))
    for t in range(steps):
        i, j = rng.integers(0, n, 2)
        if i == j:
            continue
        pi[i], pi[j] = pi[j], pi[i]
        proposed = smeared_action_2d(perm_to_causal_matrix(pi), eps)
        if np.log(rng.uniform()) < -beta * (proposed - action):
            action = proposed
        else:
            pi[i], pi[j] = pi[j], pi[i]
        if t >= steps // 10:
            counts[index[pi.tobytes()]] += 1
    empirical = counts / counts.sum()
    gibbs_tv = float(0.5 * np.abs(target - empirical).sum())
    gibbs_visited = int((counts > 0).sum())

    # Gate V2: ensemble identity at the experiment scale (beta = 0).
    rng = np.random.default_rng(2)
    samples, acceptance = mcmc_2d_order(
        rng.permutation(N_ELEMENTS), beta=0.0, eps=EPS, steps=120_000, seed=3,
        sample_every=500,
    )
    mcmc_n0 = float(np.mean([s["n0"] for s in samples]))
    identity_dev = abs(mcmc_n0 - ref["n0"]) / ref["n0"]

    summary = {
        "n_elements": N_ELEMENTS,
        "eps": EPS,
        "reference": ref,
        "gibbs_tv": gibbs_tv,
        "gibbs_visited": gibbs_visited,
        "gibbs_states": len(perms),
        "beta0_mcmc_n0": mcmc_n0,
        "beta0_identity_rel_dev": float(identity_dev),
        "beta0_acceptance": float(acceptance),
    }
    with open(OUT / "p4_stage_a_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    with open(OUT / "p4_stage_a.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary, indent=2))
    ok = gibbs_tv < 0.05 and gibbs_visited == len(perms) and identity_dev < 0.05
    print(f"stage A gates: {'PASS' if ok else 'FAIL'}")


def classify(mean_obs: dict[str, float], gates: dict) -> str:
    ref = gates["reference"]
    n12_ref = ref["n1"] + ref["n2"]
    continuum = (
        abs(mean_obs["n0"] / ref["n0"] - 1.0) <= gates["continuum_n0_rel_tol"]
        and mean_obs["height"] >= gates["continuum_height_min"]
    )
    crystal = (
        (mean_obs["n1"] + mean_obs["n2"]) <= gates["crystal_n12_frac"] * n12_ref
        and mean_obs["height"] <= gates["crystal_height_max"]
    )
    if continuum and not crystal:
        return "continuum"
    if crystal and not continuum:
        return "crystal"
    return "other"


def stage_b(betas: list[float], seeds: list[int]) -> None:
    gates = json.loads(FROZEN.read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    tag = "-".join(f"{b:g}" for b in betas)
    rows = []
    for beta in betas:
        for start in ("R", "X"):
            for seed in seeds:
                if start == "R":
                    rng0 = np.random.default_rng(seed * 7919 + int(round(beta * 10)))
                    pi0 = rng0.permutation(N_ELEMENTS)
                else:
                    pi0 = bipartite_perm(N_ELEMENTS)
                chain_seed = seed * 1000 + int(round(beta * 10)) * 2 + (start == "X")
                samples, acceptance = mcmc_2d_order(
                    pi0, beta=beta, eps=EPS, steps=STEPS, seed=chain_seed,
                    sample_every=SAMPLE_EVERY,
                )
                mean_obs = {
                    key: float(np.mean([s[key] for s in samples]))
                    for key in ("S", "R", "n0", "n1", "n2", "mm_dim", "height")
                }
                row = {"beta": beta, "start": start, "seed": seed, **mean_obs,
                       "acc": float(acceptance),
                       "class": classify(mean_obs, gates)}
                rows.append(row)
                print(
                    f"beta={beta:g} {start} seed={seed}: S={mean_obs['S']:+.1f} "
                    f"n0={mean_obs['n0']:.0f} n12={mean_obs['n1'] + mean_obs['n2']:.0f} "
                    f"h={mean_obs['height']:.1f} acc={acceptance:.3f} -> {row['class']}",
                    flush=True,
                )
    with open(OUT / f"p4_stage_b_{tag}.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote p4_stage_b_{tag}.csv")


def verdicts() -> None:
    """Aggregate all stage-B shard CSVs into the decision registry."""
    gates = json.loads(FROZEN.read_text())
    rows = []
    for path in sorted(OUT.glob("p4_stage_b_*.csv")):
        with open(path) as fh:
            rows.extend(list(csv.DictReader(fh)))
    by_beta: dict[float, list[dict]] = {}
    for row in rows:
        by_beta.setdefault(float(row["beta"]), []).append(row)
    registry = {
        "frozen_constants": str(FROZEN.name),
        "gate_start_gap_max": gates["start_gap_max"],
        "verdict_by_beta": {},
    }
    for beta in sorted(by_beta):
        chains = by_beta[beta]
        classes = [c["class"] for c in chains]
        mean_s = {
            start: float(np.mean([float(c["S"]) for c in chains if c["start"] == start]))
            for start in ("R", "X")
        }
        gap = abs(mean_s["R"] - mean_s["X"])
        if all(c == "continuum" for c in classes) and gap <= gates["start_gap_max"]:
            verdict = "continuum"
        elif all(c == "crystal" for c in classes):
            verdict = "crystal"
        else:
            verdict = "transition/hysteresis"
        registry["verdict_by_beta"][f"{beta:g}"] = {
            "n_chains": len(chains),
            "n_continuum": classes.count("continuum"),
            "n_crystal": classes.count("crystal"),
            "n_other": classes.count("other"),
            "mean_S_R": mean_s["R"],
            "mean_S_X": mean_s["X"],
            "start_gap": gap,
            "verdict": verdict,
        }
        print(f"beta={beta:g}: {registry['verdict_by_beta'][f'{beta:g}']}")
    predictions = gates["predictions"]
    matches = {
        b: registry["verdict_by_beta"][b]["verdict"] == predictions[b]
        for b in predictions
    }
    registry["predictions"] = predictions
    registry["prediction_matches"] = matches
    registry["all_predictions_match"] = all(matches.values())
    with open(OUT / "p4_stage_b_decision_registry.json", "w") as fh:
        json.dump(registry, fh, indent=2)
    print(f"all predictions match: {registry['all_predictions_match']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=["a", "b", "verdict"])
    parser.add_argument("--betas", default="0,1,2,3,4,6")
    parser.add_argument("--seeds", default="100-104")
    args = parser.parse_args()
    if args.stage == "a":
        stage_a()
    elif args.stage == "b":
        betas = [float(b) for b in args.betas.split(",")]
        lo, hi = args.seeds.split("-")
        stage_b(betas, list(range(int(lo), int(hi) + 1)))
    else:
        verdicts()
