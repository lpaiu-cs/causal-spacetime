"""P10 Stage B': the scale-fixed observable, on the untouched instrument.

See docs/prereg/p10_continuum_limit.md Section 9. Stage B failed its own
gate: under both the frozen constants (Stage A, flat) and the
continuum-scaled family (B0, rising), the sprinkled arm's truth error
never fell, so the Section 2 outcome table stayed undecidable. Both
attempts rescaled the INSTRUMENT while keeping the OBSERVABLE — sign
discordance over uniformly sampled target pairs, whose difficulty mix is
set by the fixed continuum distribution of comparison margins and
therefore never eases as resolution improves.

B' inverts that: **the instrument is the frozen P3 discriminator,
verbatim and unscaled** (which also removes B0's unscaled-fit-budget
confound — same constants, same convergence regime at every size), and
only the truth SCORING changes. The scale-fixed observable is

    margin-restricted discordance: the same pair-pair comparisons,
    restricted to quadruples whose TRUE comparison margin
    | |x_a - x_b| - |x_c - x_d| |  (continuum units, x = (i - pi_i)/N)
    is at least DELTA_MARGIN

so every size is asked the same fixed set of questions while the
instrument's per-question precision improves with density. At
``delta = 0`` the scorer reproduces the frozen scorer EXACTLY (same
seed, same pair stream) — the anchor regression that pins the new
observable to the old one at their meeting point.

DELTA_MARGIN is a priori: the 25th percentile of the comparison-margin
distribution of the EXACT continuum model (uniform iid points in the
unit (u, v) square), computed once from pure geometry with a recorded
seed — no experiment data touched it. Eligible fraction is 75% by
construction under the continuum model; the realized fraction over
selected targets is reported per sample and floored.

Stages:
    B'0  the gate, again: arm S through the frozen instrument, BOTH
         scorings per sample. Gate: the restricted discordance's
         top-minus-bottom CI entirely below zero, with full completeness
         and the eligible-comparison floor met everywhere. The
         unrestricted scoring rides along as the replication of Stage
         A's flat curve — same fits, two scorings.
    B'1  frozen in the prereg now, runnable only if B'0 passes; reuses
         the Stage B machinery (block bootstrap, chain validation,
         completeness, dual-start requirement) via import.

Usage:
    python p10_stage_bprime.py --stage b0
    python p10_stage_bprime.py --stage b1 --n 900 --start bipartite
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from p3_dynamics import analyze_order
from p5_two_orders_emergence import order_inputs
from p10_continuum_ladder import LADDER, _difference_ci, _stable_seed
from pc_common import DEFAULT_OUTPUT_DIR, git_describe, write_rows_csv

from causal_spacetime_lab.ordinal_embedding import squared_distance_matrix

#: A priori, from pure continuum geometry (uniform iid points in the
#: unit (u, v) square; margin = ||x_a-x_b| - |x_c-x_d||, x = u - v):
#: the 25th percentile, 4M draws, seed 2718281828, reproducible to
#: 3.5e-6 across seeds. No experiment data was consulted.
DELTA_MARGIN = 0.1399

#: Comparison budget for the restricted scorer. Larger than the frozen
#: 8000 because ~25% of comparisons are rejected by the margin filter;
#: the floor below guards the effective count.
BPRIME_COMPARISONS = 32_000
MIN_ELIGIBLE = 8_000

#: B'0 seeds — fresh: clear of 0-9, 100-119, 400-419, 500-519, 820-959,
#: 1000-1019 (P6-B), 1100-1159 (B0), 9001+, 30000-30379 (B1 envelope).
B0PRIME_SEED_BASE = {600: 40000, 900: 40020, 1200: 40040}
B0PRIME_SAMPLES = 20


def margin_restricted_order_error(
    estimated, reference_x, delta: float,
    num_pair_comparisons: int, seed: int,
) -> tuple[float, int]:
    """Sign-discordance over pair-pair comparisons with true margin
    >= ``delta``; returns (error, n_eligible).

    Mirrors the frozen ``embedding_distance_order_error`` sampling
    stream exactly — same RNG construction, same validity mask — and
    then applies the margin filter. At ``delta = 0`` with the same seed
    and count it must reproduce the frozen scorer bit for bit; a
    regression enforces that anchor, which is what makes this a
    restriction of the frozen observable rather than a second one.
    """

    estimate = np.asarray(estimated, dtype=float)
    if estimate.ndim == 1:
        estimate = estimate.reshape(-1, 1)
    reference = np.asarray(reference_x, dtype=float).reshape(-1, 1)
    if estimate.shape[0] != reference.shape[0]:
        raise ValueError("estimated and reference must have the same count")
    rng = np.random.default_rng(seed)
    n = reference.shape[0]
    first = rng.integers(0, n, size=(num_pair_comparisons, 2))
    second = rng.integers(0, n, size=(num_pair_comparisons, 2))
    valid = (first[:, 0] != first[:, 1]) & (second[:, 0] != second[:, 1])
    first = first[valid]
    second = second[valid]

    ref_d = squared_distance_matrix(reference)
    est_d = squared_distance_matrix(estimate)
    ref_first = ref_d[first[:, 0], first[:, 1]]
    ref_second = ref_d[second[:, 0], second[:, 1]]
    est_first = est_d[first[:, 0], first[:, 1]]
    est_second = est_d[second[:, 0], second[:, 1]]

    margin = np.abs(np.sqrt(ref_first) - np.sqrt(ref_second))
    keep = margin >= delta
    ref_sign = np.sign(ref_first[keep] - ref_second[keep])
    est_sign = np.sign(est_first[keep] - est_second[keep])
    # the frozen scorer excludes exact reference ties from the
    # denominator; the same mask keeps the delta = 0 anchor exact, and
    # for delta > 0 the margin filter already implies ref_sign != 0
    comparable = ref_sign != 0.0
    n_eligible = int(comparable.sum())
    if n_eligible == 0:
        return float("nan"), 0
    return (
        float(np.mean(ref_sign[comparable] != est_sign[comparable])),
        n_eligible,
    )


def evaluate_perm_dual(pi: np.ndarray, seed: int) -> dict:
    """The FROZEN instrument once, two truth scorings on the same fit.

    ``truth`` is the frozen scorer's output, unchanged (the Stage A
    replication arm). ``truth_restricted`` is the scale-fixed scorer on
    the identical fitted coordinates — continuum-normalized true x, the
    a-priori margin, and its own eligible count.
    """

    n = pi.size
    causal, times, coords = order_inputs(pi)
    row, coords_fit, targets = analyze_order(
        causal, times, coords, seed=seed, want_truth=True, return_fit=True,
    )
    if row.get("status") == "ok":
        true_x_continuum = coords[targets] / float(n)
        restricted, n_eligible = margin_restricted_order_error(
            coords_fit, true_x_continuum, DELTA_MARGIN,
            num_pair_comparisons=BPRIME_COMPARISONS, seed=seed + 9,
        )
        row["truth_restricted"] = restricted
        row["n_eligible"] = float(n_eligible)
        row["eligible_floor_met"] = bool(n_eligible >= MIN_ELIGIBLE)
    return row


def bprime_gate_verdict(n_ok_by_rung: dict, floors_met: bool,
                        drop, drop_ci) -> bool:
    """The Section 9 gate: full completeness at every rung, the
    eligible-comparison floor met on every sample, and the restricted
    top-minus-bottom CI entirely below zero. Completeness-conditioned
    verdicts are not verdicts (the B0-gate review lesson, inherited)."""

    complete = all(
        n_ok_by_rung.get(n, 0) == B0PRIME_SAMPLES for n in LADDER
    )
    falls = bool(
        drop is not None and drop_ci[1] is not None and drop_ci[1] < 0.0
    )
    return complete and floors_met and falls


def run_b0prime(output_dir: Path) -> None:
    rows = []
    for n in LADDER:
        base = B0PRIME_SEED_BASE[n]
        for k in range(B0PRIME_SAMPLES):
            pi = np.random.default_rng(base + k).permutation(n)
            row = {
                "arm": "S-frozen-dual", "n": float(n), "start": "uniform",
                "chain_seed": float(base + k), "sample_index": float(k),
                "code_version": git_describe(),
                **evaluate_perm_dual(pi, seed=base + k),
            }
            rows.append(row)
        done = [r for r in rows if r["n"] == n and r.get("status") == "ok"]
        med_u = float(np.median([r["truth"] for r in done]))
        med_r = float(np.median([r["truth_restricted"] for r in done]))
        print(f"B'0 n={n}: {len(done)}/{B0PRIME_SAMPLES} ok | "
              f"unrestricted {med_u:.4f} | restricted {med_r:.4f}",
              flush=True)
    write_rows_csv(output_dir / "p10_bprime0_yardstick.csv", rows)

    ok = [r for r in rows if r.get("status") == "ok"]
    summary: dict = {
        "code_version": git_describe(),
        "delta_margin": DELTA_MARGIN,
        "per_n": {},
    }
    for n in LADDER:
        sub = [r for r in ok if r["n"] == n]
        entry = {"n_ok": len(sub)}
        for key in ("truth", "truth_restricted"):
            vals = [r[key] for r in sub]
            rng = np.random.default_rng(_stable_seed("bprime0", n, key))
            boots = np.median(
                np.asarray(vals)[
                    rng.integers(0, len(vals), size=(4000, len(vals)))
                ], axis=1,
            ) if vals else np.array([])
            entry[f"median_{key}"] = float(np.median(vals)) if vals else None
            entry[f"ci_{key}"] = [
                float(np.percentile(boots, 2.5)),
                float(np.percentile(boots, 97.5)),
            ] if vals else [None, None]
        entry["min_eligible"] = (
            int(min(r["n_eligible"] for r in sub)) if sub else None
        )
        entry["all_floors_met"] = bool(
            sub and all(r["eligible_floor_met"] for r in sub)
        )
        summary["per_n"][str(n)] = entry

    drops = {}
    for key in ("truth", "truth_restricted"):
        drop, ci = _difference_ci(
            [r[key] for r in ok if r["n"] == LADDER[-1]],
            [r[key] for r in ok if r["n"] == LADDER[0]],
            seed=_stable_seed("bprime0-drop", key),
        )
        drops[key] = {"diff": drop, "ci": list(ci)}
    summary["top_minus_bottom"] = drops
    floors = all(
        summary["per_n"][str(n)]["all_floors_met"] for n in LADDER
    )
    summary["yardstick_falls_restricted"] = bprime_gate_verdict(
        {n: summary["per_n"][str(n)]["n_ok"] for n in LADDER},
        floors,
        drops["truth_restricted"]["diff"],
        tuple(drops["truth_restricted"]["ci"]),
    )
    summary["note"] = (
        "The Section 9 gate. The unrestricted scoring rides along as "
        "the replication of Stage A's flat curve on the same fits; only "
        "the restricted scoring gates. If this fails too, the "
        "scale-referenced behaviour is not a property of the question "
        "mix either, and that is the finding."
    )
    (output_dir / "p10_bprime0_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["b0", "b1"], default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.stage == "b0":
        run_b0prime(args.output_dir)
    elif args.stage == "b1":
        raise SystemExit(
            "B'1 is gated on the frozen B'0 record (prereg Section 9) and "
            "its runner is added only after that gate passes -- adding it "
            "now would be a 17-chain-hour path with nothing to authorize "
            "it."
        )
    else:
        raise SystemExit("choose --stage b0/--stage b1")


if __name__ == "__main__":
    main()
