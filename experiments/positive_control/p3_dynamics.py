"""P3: does geometry emerge from a geometry-free dynamics?

See docs/prereg/p3_dynamics_emergence.md. P3-A calibrates the order-intrinsic
discriminator on sprinkled 1+1D causal sets (pass) and their column-shuffle
(fail) and proposes gates; P3-B applies the frozen gate to transitive-
percolation (classical sequential growth) orders across a p-sweep and reports,
for each p, whether geometry emerges (pass) or not (block). Reuses the frozen
PC-V1 dissimilarity/fit pipeline; chains are selected from the order.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from pc_common import DEFAULT_OUTPUT_DIR, git_describe, parse_seed_spec, write_rows_csv

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    fit_ordinal_embedding_gradient_descent,
    quadruplet_violation_rate,
)
from causal_spacetime_lab.positive_control.dissimilarity import (
    build_constraint_split,
    margin_from_probe_quantile,
    profile_dissimilarity_matrix,
)
from causal_spacetime_lab.positive_control.dynamics import (
    relation_density,
    transitive_percolation,
)
from causal_spacetime_lab.positive_control.echo_profiles import EchoProfileMatrix
from causal_spacetime_lab.positive_control.order_intrinsic import (
    measure_order_intrinsic_profiles,
    select_bracketed_targets,
    select_disjoint_chains,
)
from causal_spacetime_lab.sprinkling import sprinkle_1p1_causal_diamond

N_SPRINKLE = 1500
N_DYNAMICS = 1500
CHAIN_COUNT = 6
MIN_CHAIN_LEN = 25
MAX_TARGETS = 44
MIN_TARGETS = 20
STEPS = 1500
RESTARTS = 5
TRAIN_C = 3000
HELDOUT_C = 800
P_SWEEP = (0.006, 0.008, 0.010, 0.014, 0.020)
FROZEN_CONSTANTS_PATH = Path("docs/prereg/frozen/p3_test_constants.json")


def _fit_heldout(profiles: EchoProfileMatrix, seed: int):
    diss = profile_dissimilarity_matrix(profiles, 4)
    margin = margin_from_probe_quantile(diss, seed=seed + 3)
    split = build_constraint_split(diss, TRAIN_C, HELDOUT_C, margin, seed=seed + 5)
    coords, _ = fit_ordinal_embedding_gradient_descent(
        profiles.target_count, 1, split.train, steps=STEPS,
        learning_rate=0.05, seed=seed + 100, restarts=RESTARTS,
    )
    return coords, quadruplet_violation_rate(coords, split.heldout)


def _column_shuffle(profiles: EchoProfileMatrix, seed: int) -> EchoProfileMatrix:
    rng = np.random.default_rng(seed)
    delays = profiles.delay_ranks.copy()
    reach = profiles.reachable.copy()
    for col in range(profiles.reference_count):
        perm = rng.permutation(profiles.target_count)
        delays[:, col] = delays[perm, col]
        reach[:, col] = reach[perm, col]
    return EchoProfileMatrix(delays, reach, profiles.target_indices.copy())


def analyze_order(causal, times, coords, seed, want_truth, extra_truth_coords=None):
    """Return a result dict for one causal order (structural blocks recorded)."""

    chains = select_disjoint_chains(causal, times, CHAIN_COUNT, MIN_CHAIN_LEN)
    if len(chains) < CHAIN_COUNT:
        return {"status": f"structural_block: only {len(chains)} chains"}
    targets = select_bracketed_targets(causal, chains, MAX_TARGETS, seed)
    if targets.size < MIN_TARGETS:
        return {"status": f"structural_block: {targets.size} targets"}
    profiles = measure_order_intrinsic_profiles(causal, chains, targets)
    try:
        coords_fit, heldout = _fit_heldout(profiles, seed)
        shuffled = _column_shuffle(profiles, seed + 61)
        _, null_heldout = _fit_heldout(shuffled, seed)
    except ValueError as error:
        return {"status": f"structural_block: {str(error)[:40]}"}
    row = {
        "status": "ok",
        "n_targets": float(targets.size),
        "min_chain_len": float(min(c.size for c in chains)),
        "heldout": float(heldout),
        "null_gap": float(null_heldout - heldout),
    }
    if want_truth and coords is not None:
        row["truth"] = float(
            embedding_distance_order_error(
                coords_fit, coords[targets].reshape(-1, 1),
                num_pair_comparisons=8000, seed=seed + 9,
            )
        )
    for name, values in (extra_truth_coords or {}).items():
        row[f"truth_{name}"] = float(
            embedding_distance_order_error(
                coords_fit, np.asarray(values)[targets].reshape(-1, 1),
                num_pair_comparisons=8000, seed=seed + 9,
            )
        )
    return row


def _rnd(v):
    return float(round(v / 0.05) * 0.05)


def _median(xs):
    xs = sorted(x for x in xs if not (isinstance(x, float) and math.isnan(x)))
    return xs[len(xs) // 2] if xs else None


def cohens_d(a, b):
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
    vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
    pooled = math.sqrt(((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2))
    return (ma - mb) / pooled if pooled else float("inf")


def stage_a(seeds, code_version, output_dir):
    rows, structured, shuffled = [], [], []
    for s in seeds:
        ev = sprinkle_1p1_causal_diamond(N_SPRINKLE, T=2.0, seed=s)
        causal = causal_matrix_1p1(ev)
        res = analyze_order(causal, ev[:, 0], ev[:, 1], s, want_truth=True)
        res.update({"stage": "P3-A", "seed": float(s), "code_version": code_version})
        rows.append(res)
        if res["status"] == "ok":
            structured.append(res)
            # shuffle-only fit for the fail cluster held-out
            shuffled.append(res["heldout"] + res["null_gap"])
        print(f"seed {s}: {res['status']}"
              + (f" heldout={res['heldout']:.3f} null_gap={res['null_gap']:.3f}"
                 f" truth={res.get('truth', float('nan')):.3f}"
                 if res["status"] == "ok" else ""))
    write_rows_csv(output_dir / "p3_stage_a.csv", rows)

    sh = [r["heldout"] for r in structured]
    sng = [r["null_gap"] for r in structured]
    st_ = [r["truth"] for r in structured if "truth" in r]
    summary = {
        "stage": "P3-A", "code_version": code_version,
        "valid": len(structured),
        "median_structured_heldout": _median(sh),
        "median_null_gap": _median(sng),
        "median_truth": _median(st_),
        "effect_size_d": cohens_d(shuffled, sh),
    }
    if structured:
        gate_heldout = min(0.10, _rnd(0.5 * (max(sh) + min(shuffled))))
        gate_nullgap = max(0.10, _rnd(0.5 * min(sng)))
        gate_truth = _rnd(0.5 * (max(st_) + 0.5)) if st_ else None
        summary.update({"proposed_gate_heldout": gate_heldout,
                        "proposed_gate_nullgap": gate_nullgap,
                        "proposed_gate_truth": gate_truth,
                        "structured_heldout_max": max(sh),
                        "shuffle_heldout_min": min(shuffled),
                        "structured_nullgap_min": min(sng)})
    (output_dir / "p3_summary_a.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def stage_b(seeds, code_version, output_dir):
    const = json.loads(FROZEN_CONSTANTS_PATH.read_text(encoding="utf-8"))
    gh, gn = float(const["gate_heldout"]), float(const["gate_nullgap"])
    rows = []
    verdict = {}
    for p in P_SWEEP:
        passes, blocks, valid = 0, 0, 0
        for s in seeds:
            causal, idx = transitive_percolation(N_DYNAMICS, p, s)
            dens = relation_density(causal, idx.astype(float))
            res = analyze_order(causal, idx.astype(float), None, s, want_truth=False)
            res.update({"stage": "P3-B", "p": p, "seed": float(s),
                        "density": dens, "code_version": code_version})
            if res["status"] == "ok":
                valid += 1
                is_pass = res["heldout"] <= gh and res["null_gap"] >= gn
                res["gate_pass"] = float(is_pass)
                passes += int(is_pass)
                blocks += int(not is_pass)
            else:
                res["gate_pass"] = 0.0
                blocks += 1
            rows.append(res)
        verdict[str(p)] = {
            "valid": valid, "pass": passes, "block": blocks,
            "emergence_supported": passes >= 16,
            "no_emergence": blocks >= 16,
        }
        tag = ("EMERGENCE" if passes >= 16
               else ("NO-EMERGENCE" if blocks >= 16 else "inconclusive"))
        print(f"p={p}: valid {valid}/{len(seeds)} pass {passes} block {blocks} "
              f"-> {tag}")
    write_rows_csv(output_dir / "p3_stage_b.csv", rows)
    registry = {"stage": "P3-B", "frozen_commit": const.get("frozen_commit"),
                "code_version": code_version, "gate_heldout": gh, "gate_nullgap": gn,
                "p_sweep": list(P_SWEEP), "verdict_by_p": verdict}
    (output_dir / "p3_stage_b_decision_registry.json").write_text(
        json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    print("\nCommit p3_stage_b_decision_registry.json under docs/prereg/frozen/.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["a", "b"], default="a")
    parser.add_argument("--seeds", default=None)
    args = parser.parse_args()
    stage = args.stage
    seeds = parse_seed_spec(args.seeds or ("0-9" if stage == "a" else "100-119"))
    code_version = git_describe()
    if stage == "b" and not FROZEN_CONSTANTS_PATH.exists():
        raise SystemExit(
            f"frozen P3 constants not found at {FROZEN_CONSTANTS_PATH}; P3-B may "
            "only run after the P3 freeze (preregistration Section 8)")
    if stage == "a":
        stage_a(seeds, code_version, DEFAULT_OUTPUT_DIR)
    else:
        stage_b(seeds, code_version, DEFAULT_OUTPUT_DIR)


if __name__ == "__main__":
    main()
