"""P7 FSS re-scope reconnaissance: do the samplable and instrument-operable
windows overlap?

The standing P7 freeze offers one escape route: re-scope the FSS axis to
the demonstrably samplable range (Wang-Landau tunnels at N <= 80,
projected feasible to N ~ 120) under a fresh prereg. This probe asks the
prerequisite question that must be answered BEFORE any prereg is drafted:
does the geometry instrument that defines G operate at those N at all?

Three measurements, all on the uniform ensemble (beta = 0, i.e. random
permutations -- the maximally geometric positive control, where a working
instrument must score G ~ 1):

1. The FROZEN protocol (6 chains x 25 ticks, >= 20 targets) across
   N = 200..500: where does its own positive control begin to pass?
   An envelope argument says it cannot begin early: the expected longest
   increasing subsequence of a uniform permutation is ~ 2 sqrt(N), so for
   N <= (25/2)^2 ~ 156 even ONE 25-tick chain typically does not exist,
   let alone six disjoint ones -- a selector-independent block. Above
   the envelope the blocks measured here are properties of the frozen
   protocol's GREEDY selector (Greene's theorem permits six disjoint
   25-chains at N = 200..400), so the boundary reported is that of the
   instrument AS FROZEN, not a selector-independent bound.
2. Re-scoped small-N candidate specs (fewer/shorter chains, fewer
   targets and constraints) across N = 100..160: does ANY of them pass
   its own positive control reliably?
3. The bipartite crystal at every probed N: must block structurally
   (height 2), as a sanity check that the re-scoped specs have not
   destroyed the discriminator's negative side.

This is reconnaissance, not a preregistered experiment: nothing is
frozen, no G(beta) is measured, no MCMC is run. Its output is a
feasibility verdict for the re-scope route itself. Deterministic seeds;
the full table is written to the tracked path
docs/p7_fss_rescope_probe_results.json (the audit-trail convention
established in PR #8).

Usage:
    python p7_fss_rescope_probe.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

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
from causal_spacetime_lab.positive_control.echo_profiles import EchoProfileMatrix
from causal_spacetime_lab.positive_control.geometry_score import (
    geometry_order_parameter,
)
from causal_spacetime_lab.positive_control.order_intrinsic import (
    measure_order_intrinsic_profiles,
    select_bracketed_targets,
    select_disjoint_chains,
)
from causal_spacetime_lab.positive_control.two_orders import (
    bipartite_perm,
    perm_to_causal_matrix,
)

ROOT = Path(__file__).resolve().parents[2]
RESULTS_PATH = ROOT / "docs" / "p7_fss_rescope_probe_results.json"

FROZEN_SPEC = {
    "chain_count": 6, "min_length": 25, "max_targets": 44,
    "min_targets": 20, "train_c": 3000, "held_c": 800,
    "steps": 1500, "restarts": 5, "min_cols": 4,
    # the frozen P3 pipeline evaluates the truth gate with 8000 pair
    # comparisons (p3_dynamics.py); the probe must match it exactly to
    # call its result a measurement of the frozen instrument
    "truth_comparisons": 8000,
}
SOFT_SPECS = {
    "3x8_t30": {
        "chain_count": 3, "min_length": 8, "max_targets": 30,
        "min_targets": 10, "train_c": 600, "held_c": 200,
        "steps": 1000, "restarts": 3, "min_cols": 3,
        "truth_comparisons": 4000,
    },
    "4x8_t30": {
        "chain_count": 4, "min_length": 8, "max_targets": 30,
        "min_targets": 10, "train_c": 600, "held_c": 200,
        "steps": 1000, "restarts": 3, "min_cols": 3,
        "truth_comparisons": 4000,
    },
    "4x10_t20": {
        "chain_count": 4, "min_length": 10, "max_targets": 20,
        "min_targets": 8, "train_c": 800, "held_c": 250,
        "steps": 1000, "restarts": 3, "min_cols": 3,
        "truth_comparisons": 4000,
    },
}


def order_inputs(pi: np.ndarray):
    idx = np.arange(pi.size, dtype=float)
    return perm_to_causal_matrix(pi), idx + pi, idx - pi


def _fit_heldout(profiles, seed, spec):
    diss = profile_dissimilarity_matrix(profiles, spec["min_cols"])
    margin = margin_from_probe_quantile(diss, seed=seed + 3)
    split = build_constraint_split(
        diss, spec["train_c"], spec["held_c"], margin, seed=seed + 5
    )
    coords, _ = fit_ordinal_embedding_gradient_descent(
        profiles.target_count, 1, split.train, steps=spec["steps"],
        learning_rate=0.05, seed=seed + 100, restarts=spec["restarts"],
    )
    return coords, quadruplet_violation_rate(coords, split.heldout)


def _column_shuffle(profiles, seed):
    rng = np.random.default_rng(seed)
    delays = profiles.delay_ranks.copy()
    reach = profiles.reachable.copy()
    for col in range(profiles.reference_count):
        perm = rng.permutation(profiles.target_count)
        delays[:, col] = delays[perm, col]
        reach[:, col] = reach[perm, col]
    return EchoProfileMatrix(delays, reach, profiles.target_indices.copy())


def evaluate(pi: np.ndarray, spec: dict, seed: int) -> dict:
    """One instrument evaluation, mirroring p3_dynamics.analyze_order with
    a parameterized spec (the frozen constants are one instance)."""

    causal, times, xcoords = order_inputs(pi)
    chains = select_disjoint_chains(
        causal, times, spec["chain_count"], spec["min_length"]
    )
    if len(chains) < spec["chain_count"]:
        return {"status": "block_chains"}
    targets = select_bracketed_targets(causal, chains, spec["max_targets"], seed)
    if targets.size < spec["min_targets"]:
        return {"status": "block_targets"}
    profiles = measure_order_intrinsic_profiles(causal, chains, targets)
    try:
        coords_fit, heldout = _fit_heldout(profiles, seed, spec)
        _, null_heldout = _fit_heldout(
            _column_shuffle(profiles, seed + 61), seed, spec
        )
    except ValueError as error:
        return {"status": f"block_fit: {str(error)[:40]}"}
    truth = embedding_distance_order_error(
        coords_fit, xcoords[targets].reshape(-1, 1),
        num_pair_comparisons=spec["truth_comparisons"], seed=seed + 9,
    )
    return {
        "status": "ok",
        "n_targets": int(targets.size),
        "heldout": float(heldout),
        "null_gap": float(null_heldout - heldout),
        "truth": float(truth),
        "G": geometry_order_parameter(
            status="ok", heldout=heldout,
            null_gap=null_heldout - heldout, truth_error=truth,
        ),
    }


def scan(spec: dict, n_elements: int, seeds: int) -> dict:
    results = []
    for seed in range(seeds):
        rng = np.random.default_rng(9000 + seed)
        results.append(evaluate(rng.permutation(n_elements), spec, seed))
    ok = [r for r in results if r["status"] == "ok"]
    blocks: dict[str, int] = {}
    for r in results:
        if r["status"] != "ok":
            blocks[r["status"]] = blocks.get(r["status"], 0) + 1
    summary = {
        "n_elements": n_elements,
        "seeds": seeds,
        "ok": len(ok),
        "blocks": blocks,
    }
    if ok:
        summary.update({
            "positive_pass": sum(1 for r in ok if r["G"] >= 0.5),
            "g_median": float(np.median([r["G"] for r in ok])),
            "heldout_median": float(np.median([r["heldout"] for r in ok])),
            "null_gap_median": float(np.median([r["null_gap"] for r in ok])),
            "truth_median": float(np.median([r["truth"] for r in ok])),
        })
    return summary


def main() -> None:
    started = time.perf_counter()
    results: dict = {
        "frozen_spec": FROZEN_SPEC,
        "soft_specs": SOFT_SPECS,
        "frozen_scan": [],
        "soft_scan": {},
        "crystal_control": [],
    }

    print("== frozen protocol on its positive control ==")
    for n in (200, 300, 400, 500):
        row = scan(FROZEN_SPEC, n, seeds=6)
        results["frozen_scan"].append(row)
        print(f"  N={n}: {json.dumps(row)}")

    print("== re-scoped candidate specs on their positive control ==")
    for name, spec in SOFT_SPECS.items():
        results["soft_scan"][name] = []
        for n in (100, 120, 140, 160):
            row = scan(spec, n, seeds=8)
            results["soft_scan"][name].append(row)
            print(f"  N={n} {name}: {json.dumps(row)}")

    print("== crystal negative control (must block structurally) ==")
    crystal_ok = True
    for n in (100, 120, 160, 500):
        causal, times, _ = order_inputs(bipartite_perm(n))
        extracted = len(select_disjoint_chains(causal, times, 3, 8))
        results["crystal_control"].append(
            {"n_elements": n, "chains_extracted": extracted}
        )
        crystal_ok = crystal_ok and extracted == 0
        print(f"  N={n}: chains extracted = {extracted} (want 0)")

    # The verdict this probe exists to compute, under TWO operability
    # criteria reported side by side (so the conclusion is not an
    # artifact of one threshold):
    #  - "generous": at least half the positive-control runs pass.
    #  - "calibration-grade": median G >= 0.5 AND at least 3/4 pass --
    #    still far weaker than the P7 design's frozen calibration
    #    requirement, "G ~ 1 on the uniform ensemble"
    #    (docs/next_experiments_plan.md), which the N = 600 instrument
    #    meets (G median 0.852 here, all ok runs passing).
    # The samplable boundary (N ~ 120) comes from the measured tunneling
    # exponent (docs/p7_enhanced_sampling.md), not from this probe.
    def _generous(row: dict) -> bool:
        return row["ok"] > 0 and row.get("positive_pass", 0) * 2 >= row["ok"]

    def _calibration_grade(row: dict) -> bool:
        return (
            row["ok"] > 0
            and row.get("g_median", 0.0) >= 0.5
            and row.get("positive_pass", 0) * 4 >= row["ok"] * 3
        )

    def _first_operable(rows: list[dict], criterion) -> int | None:
        hits = [row["n_elements"] for row in rows if criterion(row)]
        return min(hits) if hits else None

    all_soft_rows = [
        row for rows in results["soft_scan"].values() for row in rows
    ]
    samplable_rows = [r for r in all_soft_rows if r["n_elements"] <= 120]
    results["verdict"] = {
        "crystal_blocks_everywhere": bool(crystal_ok),
        "frozen_operable_from_generous": _first_operable(
            results["frozen_scan"], _generous
        ),
        "frozen_operable_from_calibration_grade": _first_operable(
            results["frozen_scan"], _calibration_grade
        ),
        "soft_operable_from_generous": _first_operable(
            all_soft_rows, _generous
        ),
        "soft_operable_from_calibration_grade": _first_operable(
            all_soft_rows, _calibration_grade
        ),
        "windows_overlap_generous": bool(
            any(_generous(r) for r in samplable_rows)
        ),
        "windows_overlap_calibration_grade": bool(
            any(_calibration_grade(r) for r in samplable_rows)
        ),
    }
    results["status"] = (
        "reconnaissance; nothing frozen; no G(beta) measured; no MCMC run"
    )
    results["elapsed_seconds"] = round(time.perf_counter() - started, 1)

    RESULTS_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nverdict: {json.dumps(results['verdict'], indent=2)}")
    print(f"wrote {RESULTS_PATH}")


if __name__ == "__main__":
    main()
