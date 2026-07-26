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

from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    squared_distance_matrix,
)


def _diag_code_version() -> str:
    """git_describe plus a dirty marker. Two diagnostic artifacts were
    produced from uncommitted implementations and stamped a clean HEAD
    that could not reproduce them (review finding); this makes that
    state visible instead of silent."""

    import subprocess
    dirty = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        capture_output=True, text=True,
    ).stdout.strip()
    return git_describe() + ("-dirty" if dirty else "")

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

#: Review round six: the DERIVED scorer seeds base+k+9 sit inside the
#: consecutive chain-seed spans -- 51 of 60 scorer streams coincide
#: with some sample's chain stream (values 40009-40059; e.g. 40020
#: scores N=600 sample 11 AND drives N=900 sample 0's chain). The
#: frozen record stands as run; the seedfix stage re-scores the same
#: deterministic fits under this namespace, disjoint from every seed
#: ever used in the programme (chain seeds 40000-40059, old scorer
#: seeds 40009-40068, and all documented earlier ranges).
SEEDFIX_SCORE_SEED_BASE = {600: 41000, 900: 41020, 1200: 41040}

#: Round eight: the frozen instrument ITSELF derives seed+3 (margin),
#: seed+5 (split), seed+9 (truth scorer), seed+61 (null shuffle) and
#: seed+100 (fit learning) from every chain seed, so consecutive chain
#: seeds share derived streams across rows no matter what the external
#: scorers do -- e.g. chain seed 40000's margin stream (40003) IS chain
#: seed 40003's primary stream. The internal offsets are frozen with
#: the instrument (one pipeline definition, never re-implemented), so
#: the allocation fix is SPACING: stride 200 gives every row the
#: private window [s, s+199], every derived stream stays inside its own
#: row's window, and the spaced scorer seed s+150 does too. New chain
#: seeds mean new samples: the spaced arm is a REPLICATION of the B'0
#: design under a corrected allocation, not a replay.
SPACED_CHAIN_SEED_BASE = {600: 43000, 900: 47000, 1200: 51000}
SPACED_STRIDE = 200
SPACED_SCORE_OFFSET = 150
INSTRUMENT_DERIVED_OFFSETS = (0, 3, 5, 9, 61, 100)


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


#: A-priori bin edges for the POST-HOC mix-matched diagnostic: the
#: conditional quartiles of the continuum margin GIVEN margin >= delta,
#: same pure-geometry stream and seed as DELTA_MARGIN (reproducible to
#: 2.1e-4 or better across seeds). Equal-weighting the four bins asks
#: the same margin MIX at every size, which a common cutoff alone does
#: not (review finding: the surviving distribution above the cutoff can
#: drift with the targets' spatial configuration, and the realized
#: eligible rates 62.6% -> 66.2% show it does).
MARGIN_BIN_EDGES = (0.1399, 0.2633, 0.4160, 0.6367, float("inf"))


def margin_binned_order_error(
    estimated, reference_x, num_pair_comparisons: int, seed: int,
) -> tuple[list, list]:
    """Per-bin discordance over the a-priori margin bins; returns
    (errors_by_bin, counts_by_bin). Same sampling stream as the
    restricted scorer; the equal-weight mean of the four bin errors is
    the mix-matched error."""

    estimate = np.asarray(estimated, dtype=float)
    if estimate.ndim == 1:
        estimate = estimate.reshape(-1, 1)
    reference = np.asarray(reference_x, dtype=float).reshape(-1, 1)
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
    est_sign = np.sign(
        est_d[first[:, 0], first[:, 1]] - est_d[second[:, 0], second[:, 1]]
    )
    ref_sign = np.sign(ref_first - ref_second)
    margin = np.abs(np.sqrt(ref_first) - np.sqrt(ref_second))
    errors, counts = [], []
    for lo, hi in zip(MARGIN_BIN_EDGES[:-1], MARGIN_BIN_EDGES[1:],
                      strict=True):
        keep = (margin >= lo) & (margin < hi) & (ref_sign != 0.0)
        counts.append(int(keep.sum()))
        errors.append(
            float(np.mean(ref_sign[keep] != est_sign[keep]))
            if counts[-1] else float("nan")
        )
    return errors, counts


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


#: 16 a-priori fine bins over the FULL margin range (continuum
#: 16-quantiles, same geometry stream as DELTA_MARGIN, cross-seed delta
#: <= 5e-4). Review escalation on 9.6: four broad bins leave room for
#: WITHIN-bin drift, and error varies sharply with margin, so an
#: aggregate -- however reweighted -- can always be questioned. The
#: pointwise error-versus-margin CURVE cannot: if E_1200(m) sits above
#: E_600(m) bin by fine bin, no mix-matching objection applies, because
#: nothing is being mixed.
FINE_BIN_EDGES = (
    0.0, 0.0330, 0.0673, 0.1028, 0.1399, 0.1788, 0.2198, 0.2633,
    0.3100, 0.3606, 0.4160, 0.4784, 0.5503, 0.6367, 0.7473, 0.9144,
    float("inf"),
)

#: The high-margin pool's threshold IS a bin edge, derived, not a
#: rounded literal: a 0.478 literal once labelled a pool that actually
#: begins at 0.4784 (review finding -- the [0.4780, 0.4784) sliver was
#: excluded from the counts but included in the label).
HIGH_EDGE = FINE_BIN_EDGES[11]


def fine_binned_order_error(
    estimated, reference_x, num_pair_comparisons: int, seed: int,
) -> tuple[list, list]:
    """Per-bin discordance over the 16 fine a-priori bins (full margin
    range); same sampling stream as the other scorers."""

    estimate = np.asarray(estimated, dtype=float)
    if estimate.ndim == 1:
        estimate = estimate.reshape(-1, 1)
    reference = np.asarray(reference_x, dtype=float).reshape(-1, 1)
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
    est_sign = np.sign(
        est_d[first[:, 0], first[:, 1]] - est_d[second[:, 0], second[:, 1]]
    )
    ref_sign = np.sign(ref_first - ref_second)
    margin = np.abs(np.sqrt(ref_first) - np.sqrt(ref_second))
    errors, counts, mean_margins, discordant = [], [], [], []
    disc_margins: list = []
    for lo, hi in zip(FINE_BIN_EDGES[:-1], FINE_BIN_EDGES[1:], strict=True):
        keep = (margin >= lo) & (margin < hi) & (ref_sign != 0.0)
        wrong = keep & (ref_sign != est_sign)
        counts.append(int(keep.sum()))
        discordant.append(int(wrong.sum()))
        errors.append(
            float(discordant[-1] / counts[-1]) if counts[-1]
            else float("nan")
        )
        mean_margins.append(
            float(np.mean(margin[keep])) if counts[-1] else float("nan")
        )
        disc_margins.append([float(m) for m in margin[wrong]])
    return errors, counts, mean_margins, discordant, disc_margins


def run_margin_curve_diagnostic(output_dir: Path) -> None:
    """POST-HOC, second escalation: the pointwise error-vs-margin curve.

    Same 60 samples, deterministic. Per rung and per fine bin, the
    median error over the 20 samples with a bootstrap interval; the
    deliverable is the pointwise comparison of the N = 1200 curve
    against the N = 600 curve. Pointwise domination is immune to mix
    drift at ANY granularity by construction. Like 9.6, this can reword
    the closure's mechanism claim and cannot reopen the frozen gate.
    """

    per_sample: dict = {n: [] for n in LADDER}
    for n in LADDER:
        base = B0PRIME_SEED_BASE[n]
        for k in range(B0PRIME_SAMPLES):
            pi = np.random.default_rng(base + k).permutation(n)
            causal, times, coords = order_inputs(pi)
            row, coords_fit, targets = analyze_order(
                causal, times, coords, seed=base + k, want_truth=True,
                return_fit=True,
            )
            if row.get("status") == "ok":
                (errors, counts, mean_margins, discordant,
                 disc_margins) = fine_binned_order_error(
                    coords_fit, coords[targets] / float(n),
                    num_pair_comparisons=BPRIME_COMPARISONS,
                    seed=base + k + 9,
                )
                per_sample[n].append(
                    (errors, counts, mean_margins, discordant, disc_margins)
                )
        print(f"curve n={n}: {len(per_sample[n])}/{B0PRIME_SAMPLES} ok",
              flush=True)

    n_bins = len(FINE_BIN_EDGES) - 1
    summary: dict = {
        "code_version": _diag_code_version(),
        "fine_bin_edges": list(FINE_BIN_EDGES[:-1]) + ["inf"],
        "curves": {}, "pointwise_1200_vs_600": {},
    }
    curves: dict = {}
    for n in LADDER:
        medians, cis, min_counts = [], [], []
        for b in range(n_bins):
            vals = np.array([
                e[b] for e, *_ in per_sample[n] if not np.isnan(e[b])
            ])
            rng = np.random.default_rng(_stable_seed("curve", n, b))
            boots = np.median(vals[
                rng.integers(0, vals.size, size=(4000, vals.size))
            ], axis=1) if vals.size else np.array([])
            medians.append(float(np.median(vals)) if vals.size else None)
            cis.append([
                float(np.percentile(boots, 2.5)),
                float(np.percentile(boots, 97.5)),
            ] if vals.size else [None, None])
            min_counts.append(min(
                (c[b] for _, c, *_ in per_sample[n]), default=0
            ))
        mean_margin_by_bin = [
            float(np.nanmean([m[b] for _, _, m, *_ in per_sample[n]]))
            for b in range(n_bins)
        ]
        pooled_total_by_bin = [
            int(sum(c[b] for _, c, *_ in per_sample[n]))
            for b in range(n_bins)
        ]
        pooled_discordant_by_bin = [
            int(sum(d[b] for _, _, _, d, _ in per_sample[n]))
            for b in range(n_bins)
        ]
        curves[n] = (medians, cis)
        summary["curves"][str(n)] = {
            "median_by_bin": medians, "ci_by_bin": cis,
            "min_count_by_bin": min_counts,
            "mean_margin_by_bin": mean_margin_by_bin,
            "pooled_total_by_bin": pooled_total_by_bin,
            "pooled_discordant_by_bin": pooled_discordant_by_bin,
        }
    above, separated = [], []
    for b in range(n_bins):
        m600, ci600 = curves[600][0][b], curves[600][1][b]
        m1200, ci1200 = curves[1200][0][b], curves[1200][1][b]
        if m600 is None or m1200 is None:
            continue
        above.append(bool(m1200 > m600))
        separated.append(bool(ci1200[0] > ci600[1]))
    summary["pointwise_1200_vs_600"] = {
        "bins_compared": len(above),
        "bins_where_1200_above_600": int(sum(above)),
        "bins_ci_separated": int(sum(separated)),
    }

    # Heuristic indicator, NOT a bound (review corrected the earlier
    # claim: a mean shift times a neighbour slope bounds nothing --
    # equal means can hide mass reallocation, and the neighbour slope
    # is no Lipschitz constant). Kept as a descriptive indicator only.
    m600 = summary["curves"]["600"]["mean_margin_by_bin"]
    m1200 = summary["curves"]["1200"]["mean_margin_by_bin"]
    e600 = curves[600][0]
    e1200 = curves[1200][0]
    indicator_rows = []
    for b in range(n_bins):
        if e600[b] is None or e1200[b] is None:
            continue
        indicator_rows.append({
            "bin": b, "rise": round(e1200[b] - e600[b], 5),
            "within_bin_mean_margin_shift": round(abs(m1200[b] - m600[b]), 5),
        })
    summary["within_bin_shift_indicator"] = indicator_rows

    # The EXACT pooled analysis, which is what settles the high-margin
    # region without any binning or distributional assumption. Facts
    # about every sampled comparison, not medians of per-sample rates:
    # pooled discordant counts per rung, and the m* construction --
    # m* = the smallest true margin at which the N = 1200 pool contains
    # a discordant comparison above the threshold below; the N = 600
    # pool's comparisons at margins >= m* are then counted exactly.
    pooled = {}
    full_wrong: dict = {}
    for n in LADDER:
        total = sum(
            c[b] for _, c, *_ in per_sample[n] for b in range(n_bins)
            if FINE_BIN_EDGES[b] >= HIGH_EDGE
        )
        # full precision for computation; rounding is display-only in
        # the serialized list (review finding: a 4-decimal threshold
        # leaked into the m* count, making 'exact at >= m*' false off
        # the representable grid)
        full_wrong[n] = sorted(
            m
            for _, _, _, _, dm in per_sample[n]
            for b in range(n_bins) if FINE_BIN_EDGES[b] >= HIGH_EDGE
            for m in dm[b]
        )
        pooled[n] = {"total_at_or_above_high": int(total),
                     "discordant_at_or_above_high": len(full_wrong[n]),
                     "discordant_margins": [round(m, 4)
                                            for m in full_wrong[n][:50]]}
    summary["pooled_high_margin"] = {
        "high_edge": HIGH_EDGE, **{str(n): pooled[n] for n in LADDER}
    }
    star = full_wrong[1200][0] if full_wrong[1200] else None
    if star is not None:
        totals_600_above_star = 0
        wrong_600_above_star = 0
        for _, c, _, _d, dm in per_sample[600]:
            for b in range(n_bins):
                if FINE_BIN_EDGES[b + 1] <= star:
                    continue
                totals_600_above_star += c[b]
                wrong_600_above_star += sum(
                    1 for m in dm[b] if m >= star
                )
        summary["m_star_analysis"] = {
            "m_star": star,
            "note": "totals for N=600 use whole bins overlapping "
                    "[m*, inf); discordant counts are exact at >= m*",
            "n600_comparisons_in_overlapping_bins": int(totals_600_above_star),
            "n600_discordant_at_or_above_m_star": int(wrong_600_above_star),
        }
    (output_dir / "p10_bprime_margincurve_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary["pointwise_1200_vs_600"], indent=2))


def run_mix_matched_diagnostic(output_dir: Path) -> None:
    """POST-HOC robustness for the 9.5 closure, labelled as such.

    Re-runs the exact B'0 samples (same seeds, deterministic) and scores
    each fit per a-priori margin bin; the equal-weight bin mean is the
    mix-matched error, immune to drift in the surviving-margin
    distribution that a common cutoff admits. If the rise survives
    mix-matching, the closure claim strengthens (per-question accuracy
    degrades bin by bin); if it vanishes, the B'0 rise was partly a mix
    artifact and the record must weaken accordingly. Either way it is a
    diagnostic run AFTER the frozen gate failed -- it can reword the
    conclusion, never reopen the gate.
    """

    rows = []
    for n in LADDER:
        base = B0PRIME_SEED_BASE[n]
        for k in range(B0PRIME_SAMPLES):
            pi = np.random.default_rng(base + k).permutation(n)
            causal, times, coords = order_inputs(pi)
            row, coords_fit, targets = analyze_order(
                causal, times, coords, seed=base + k, want_truth=True,
                return_fit=True,
            )
            out = {"n": float(n), "sample_index": float(k),
                   "chain_seed": float(base + k),
                   "code_version": _diag_code_version(),
                   "status": row.get("status")}
            if row.get("status") == "ok":
                errors, counts = margin_binned_order_error(
                    coords_fit, coords[targets] / float(n),
                    num_pair_comparisons=BPRIME_COMPARISONS,
                    seed=base + k + 9,
                )
                for b, (err, cnt) in enumerate(
                        zip(errors, counts, strict=True)):
                    out[f"bin{b}_error"] = err
                    out[f"bin{b}_count"] = float(cnt)
                finite = [e for e in errors if not np.isnan(e)]
                out["mix_matched_error"] = (
                    float(np.mean(errors)) if len(finite) == len(errors)
                    else float("nan")
                )
            rows.append(out)
        done = [r for r in rows if r["n"] == n and r["status"] == "ok"]
        print(f"diag n={n}: {len(done)}/{B0PRIME_SAMPLES} ok, mix-matched "
              f"median {np.median([r['mix_matched_error'] for r in done]):.4f}",
              flush=True)
    write_rows_csv(output_dir / "p10_bprime_mixmatched.csv", rows)

    ok = [r for r in rows if r["status"] == "ok"
          and not np.isnan(r["mix_matched_error"])]
    summary: dict = {"code_version": _diag_code_version(),
                     "bin_edges": list(MARGIN_BIN_EDGES[:-1]) + ["inf"],
                     "per_n": {}}
    for n in LADDER:
        sub = [r for r in ok if r["n"] == n]
        entry = {"n_ok": len(sub)}
        vals = [r["mix_matched_error"] for r in sub]
        rng = np.random.default_rng(_stable_seed("mixmatch", n))
        boots = np.median(np.asarray(vals)[
            rng.integers(0, len(vals), size=(4000, len(vals)))], axis=1)
        entry["median_mix_matched"] = float(np.median(vals))
        entry["ci"] = [float(np.percentile(boots, 2.5)),
                       float(np.percentile(boots, 97.5))]
        for b in range(4):
            bvals = [r[f"bin{b}_error"] for r in sub]
            entry[f"median_bin{b}"] = float(np.median(bvals))
        summary["per_n"][str(n)] = entry
    drop, ci = _difference_ci(
        [r["mix_matched_error"] for r in ok if r["n"] == LADDER[-1]],
        [r["mix_matched_error"] for r in ok if r["n"] == LADDER[0]],
        seed=_stable_seed("mixmatch-drop"),
    )
    summary["top_minus_bottom_mix_matched"] = {"diff": drop, "ci": list(ci)}
    summary["rise_survives_mix_matching"] = bool(
        drop is not None and ci[0] is not None and ci[0] > 0.0
    )
    (output_dir / "p10_bprime_mixmatched_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def _seed_collision_map() -> dict:
    """The round-six diagnosis, computed rather than asserted: which
    old derived scorer seeds (base+k+9) coincide with chain seeds."""

    chain = {b + k for b in B0PRIME_SEED_BASE.values()
             for k in range(B0PRIME_SAMPLES)}
    old_scorer = {b + k + 9 for b in B0PRIME_SEED_BASE.values()
                  for k in range(B0PRIME_SAMPLES)}
    clean = {b + k for b in SEEDFIX_SCORE_SEED_BASE.values()
             for k in range(B0PRIME_SAMPLES)}
    return {
        "chain_seed_span": [min(chain), max(chain)],
        "old_scorer_seed_span": [min(old_scorer), max(old_scorer)],
        "n_old_scorer_seeds_colliding_with_chains":
            len(old_scorer & chain),
        "n_scorer_seeds_total": len(old_scorer),
        "clean_namespace_span": [min(clean), max(clean)],
        "clean_disjoint_from_all": not (clean & (chain | old_scorer)),
    }


def run_seedfix_rescore(output_dir: Path) -> None:
    """POST-HOC, review round six: derived-seed namespace robustness.

    The fits are untouched -- chain seeds do not collide with each
    other, and the fits ARE the frozen objects. This replays the same
    60 fits deterministically and scores every published B' quantity
    under BOTH scorer namespaces: the original base+k+9 (which doubles
    as an exact replay check against the frozen artifacts) and the
    disjoint SEEDFIX namespace. Labelled post hoc: whatever it reads,
    the frozen record stands and the gate stays failed-as-frozen. The
    scorer's own randomness is pair SELECTION only, so per-row shifts
    are bounded in probability by the pair-sampling noise; this run
    measures them instead of arguing that.
    """

    n_bins = len(FINE_BIN_EDGES) - 1
    rows = []
    fine: dict = {"old": {n: [] for n in LADDER},
                  "clean": {n: [] for n in LADDER}}
    for n in LADDER:
        base = B0PRIME_SEED_BASE[n]
        for k in range(B0PRIME_SAMPLES):
            pi = np.random.default_rng(base + k).permutation(n)
            causal, times, coords = order_inputs(pi)
            row, coords_fit, targets = analyze_order(
                causal, times, coords, seed=base + k, want_truth=True,
                return_fit=True,
            )
            out = {"n": float(n), "sample_index": float(k),
                   "chain_seed": float(base + k),
                   "old_score_seed": float(base + k + 9),
                   "clean_score_seed": float(SEEDFIX_SCORE_SEED_BASE[n] + k),
                   "code_version": _diag_code_version(),
                   "status": row.get("status")}
            if row.get("status") == "ok":
                true_x = coords[targets] / float(n)
                # The UNRESTRICTED truth arm (round-seven finding): the
                # frozen instrument scores it internally at seed+9 --
                # the same colliding stream -- so it must be re-scored
                # here too, not copied from the row. The old-namespace
                # recomputation must equal the instrument's own output
                # bit for bit, or the replay is broken.
                truth_ref = coords[targets].reshape(-1, 1)
                truth_old = float(embedding_distance_order_error(
                    coords_fit, truth_ref,
                    num_pair_comparisons=8000, seed=base + k + 9,
                ))
                if truth_old != row["truth"]:
                    raise RuntimeError(
                        f"truth replay mismatch at n={n} k={k}: "
                        f"{truth_old} != {row['truth']}"
                    )
                out["truth_old"] = truth_old
                out["truth_clean"] = float(embedding_distance_order_error(
                    coords_fit, truth_ref,
                    num_pair_comparisons=8000,
                    seed=SEEDFIX_SCORE_SEED_BASE[n] + k,
                ))
                out["delta_truth"] = out["truth_clean"] - out["truth_old"]
                for tag, seed in (
                    ("old", base + k + 9),
                    ("clean", SEEDFIX_SCORE_SEED_BASE[n] + k),
                ):
                    restricted, n_eligible = margin_restricted_order_error(
                        coords_fit, true_x, DELTA_MARGIN,
                        num_pair_comparisons=BPRIME_COMPARISONS, seed=seed,
                    )
                    out[f"restricted_{tag}"] = restricted
                    out[f"eligible_{tag}"] = float(n_eligible)
                    errors, _counts = margin_binned_order_error(
                        coords_fit, true_x,
                        num_pair_comparisons=BPRIME_COMPARISONS, seed=seed,
                    )
                    finite = [e for e in errors if not np.isnan(e)]
                    out[f"mix_matched_{tag}"] = (
                        float(np.mean(errors))
                        if len(finite) == len(errors) else float("nan")
                    )
                    fine[tag][n].append(fine_binned_order_error(
                        coords_fit, true_x,
                        num_pair_comparisons=BPRIME_COMPARISONS, seed=seed,
                    ))
                out["delta_restricted"] = (
                    out["restricted_clean"] - out["restricted_old"]
                )
                out["delta_mix_matched"] = (
                    out["mix_matched_clean"] - out["mix_matched_old"]
                )
            rows.append(out)
        done = [r for r in rows if r["n"] == n and r["status"] == "ok"]
        print(f"seedfix n={n}: {len(done)}/{B0PRIME_SAMPLES} ok | "
              f"max |delta restricted| "
              f"{max(abs(r['delta_restricted']) for r in done):.5f}",
              flush=True)
    write_rows_csv(output_dir / "p10_bprime_seedfix.csv", rows)

    ok = [r for r in rows if r["status"] == "ok"]
    summary: dict = {
        "code_version": _diag_code_version(),
        "collision_map": _seed_collision_map(),
        "delta_margin": DELTA_MARGIN,
        "per_n": {}, "gate": {}, "mix_matched": {}, "curve": {},
    }
    for n in LADDER:
        sub = [r for r in ok if r["n"] == n]
        entry: dict = {"n_ok": len(sub)}
        for tag in ("old", "clean"):
            for quantity in ("restricted", "truth"):
                vals = [r[f"{quantity}_{tag}"] for r in sub]
                rng = np.random.default_rng(
                    _stable_seed("seedfix", quantity, n, tag)
                )
                boots = np.median(np.asarray(vals)[
                    rng.integers(0, len(vals), size=(4000, len(vals)))
                ], axis=1) if vals else np.array([])
                entry[f"median_{quantity}_{tag}"] = (
                    float(np.median(vals)) if vals else None
                )
                entry[f"ci_{quantity}_{tag}"] = [
                    float(np.percentile(boots, 2.5)),
                    float(np.percentile(boots, 97.5)),
                ] if vals else [None, None]
            entry[f"min_eligible_{tag}"] = (
                int(min(r[f"eligible_{tag}"] for r in sub)) if sub else None
            )
        for quantity in ("restricted", "mix_matched", "truth"):
            entry[f"max_abs_delta_{quantity}"] = (
                float(max(abs(r[f"delta_{quantity}"]) for r in sub))
                if sub else None
            )
        summary["per_n"][str(n)] = entry

    for tag in ("old", "clean"):
        drop, ci = _difference_ci(
            [r[f"restricted_{tag}"] for r in ok if r["n"] == LADDER[-1]],
            [r[f"restricted_{tag}"] for r in ok if r["n"] == LADDER[0]],
            seed=_stable_seed("seedfix-drop", tag),
        )
        floors = all(
            r[f"eligible_{tag}"] >= MIN_ELIGIBLE for r in ok
        )
        tdrop, tci = _difference_ci(
            [r[f"truth_{tag}"] for r in ok if r["n"] == LADDER[-1]],
            [r[f"truth_{tag}"] for r in ok if r["n"] == LADDER[0]],
            seed=_stable_seed("seedfix-truthdrop", tag),
        )
        summary["gate"][tag] = {
            "top_minus_bottom_restricted": {"diff": drop, "ci": list(ci)},
            "top_minus_bottom_truth_unrestricted": {
                "diff": tdrop, "ci": list(tci),
            },
            "all_floors_met": bool(floors),
            "yardstick_falls_restricted": bprime_gate_verdict(
                {n: sum(1 for r in ok if r["n"] == n) for n in LADDER},
                floors, drop, tuple(ci),
            ),
        }
        mdrop, mci = _difference_ci(
            [r[f"mix_matched_{tag}"] for r in ok if r["n"] == LADDER[-1]],
            [r[f"mix_matched_{tag}"] for r in ok if r["n"] == LADDER[0]],
            seed=_stable_seed("seedfix-mixdrop", tag),
        )
        summary["mix_matched"][tag] = {
            "top_minus_bottom": {"diff": mdrop, "ci": list(mci)},
            "rise_survives": bool(
                mdrop is not None and mci[0] is not None and mci[0] > 0.0
            ),
        }

    for tag in ("old", "clean"):
        curve_entry: dict = {}
        curves: dict = {}
        for n in LADDER:
            medians, cis = [], []
            for b in range(n_bins):
                vals = np.array([
                    e[b] for e, *_ in fine[tag][n] if not np.isnan(e[b])
                ])
                rng = np.random.default_rng(
                    _stable_seed("seedfix-curve", n, b, tag)
                )
                boots = np.median(vals[
                    rng.integers(0, vals.size, size=(4000, vals.size))
                ], axis=1) if vals.size else np.array([])
                medians.append(
                    float(np.median(vals)) if vals.size else None
                )
                cis.append([
                    float(np.percentile(boots, 2.5)),
                    float(np.percentile(boots, 97.5)),
                ] if vals.size else [None, None])
            curves[n] = (medians, cis)
            curve_entry[str(n)] = {
                "pooled_total_by_bin": [
                    int(sum(c[b] for _, c, *_ in fine[tag][n]))
                    for b in range(n_bins)
                ],
                "pooled_discordant_by_bin": [
                    int(sum(d[b] for _, _, _, d, _ in fine[tag][n]))
                    for b in range(n_bins)
                ],
            }
        above, separated = [], []
        for b in range(n_bins):
            m600, ci600 = curves[600][0][b], curves[600][1][b]
            m1200, ci1200 = curves[1200][0][b], curves[1200][1][b]
            if m600 is None or m1200 is None:
                continue
            above.append(bool(m1200 > m600))
            separated.append(bool(ci1200[0] > ci600[1]))
        curve_entry["pointwise_1200_vs_600"] = {
            "bins_compared": len(above),
            "bins_where_1200_above_600": int(sum(above)),
            "bins_ci_separated": int(sum(separated)),
        }
        curve_entry["high_edge"] = HIGH_EDGE
        curve_entry["last_bin"] = {
            str(n): {
                "total": int(sum(c[-1] for _, c, *_ in fine[tag][n])),
                "discordant": int(
                    sum(d[-1] for _, _, _, d, _ in fine[tag][n])
                ),
            } for n in LADDER
        }
        summary["curve"][tag] = curve_entry

    frozen_dir = (Path(__file__).resolve().parents[2]
                  / "docs" / "prereg" / "frozen" / "p10_stage_bprime")
    replay: dict = {}
    try:
        frozen = json.loads(
            (frozen_dir / "p10_bprime0_summary.json").read_text(
                encoding="utf-8")
        )
        replay["b0prime_median_restricted_max_abs_diff"] = max(
            abs(summary["per_n"][str(n)]["median_restricted_old"]
                - frozen["per_n"][str(n)]["median_truth_restricted"])
            for n in LADDER
        )
        replay["b0prime_median_truth_max_abs_diff"] = max(
            abs(summary["per_n"][str(n)]["median_truth_old"]
                - frozen["per_n"][str(n)]["median_truth"])
            for n in LADDER
        )
    except OSError:
        replay["b0prime_median_restricted_max_abs_diff"] = None
        replay["b0prime_median_truth_max_abs_diff"] = None
    try:
        frozen_curve = json.loads(
            (frozen_dir / "p10_bprime_margincurve_summary.json").read_text(
                encoding="utf-8")
        )
        replay["curve_pooled_counts_match_frozen"] = all(
            summary["curve"]["old"][str(n)][key]
            == frozen_curve["curves"][str(n)][key]
            for n in LADDER
            for key in ("pooled_total_by_bin", "pooled_discordant_by_bin")
        )
    except OSError:
        replay["curve_pooled_counts_match_frozen"] = None
    summary["replay_check"] = replay

    def _de_nan(obj):
        if isinstance(obj, float) and np.isnan(obj):
            return None
        if isinstance(obj, dict):
            return {key: _de_nan(val) for key, val in obj.items()}
        if isinstance(obj, list):
            return [_de_nan(val) for val in obj]
        return obj

    summary = _de_nan(summary)
    (output_dir / "p10_bprime_seedfix_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"collision_map": summary["collision_map"],
                      "gate": summary["gate"],
                      "replay_check": summary["replay_check"]}, indent=2))


def _spaced_seed(n: int, k: int) -> int:
    return SPACED_CHAIN_SEED_BASE[n] + SPACED_STRIDE * k


def _spaced_window_map() -> dict:
    """Every row's private seed window, with the disjointness computed
    rather than asserted: all instrument-derived offsets and the scorer
    offset lie inside the stride, and the 60 windows are pairwise
    disjoint and clear of every earlier namespace."""

    windows = sorted(
        _spaced_seed(n, k)
        for n in LADDER for k in range(B0PRIME_SAMPLES)
    )
    pairwise_disjoint = all(
        b - a >= SPACED_STRIDE
        for a, b in zip(windows, windows[1:], strict=False)
    )
    return {
        "stride": SPACED_STRIDE,
        "instrument_derived_offsets": list(INSTRUMENT_DERIVED_OFFSETS),
        "score_offset": SPACED_SCORE_OFFSET,
        "offsets_inside_stride": bool(
            max(*INSTRUMENT_DERIVED_OFFSETS, SPACED_SCORE_OFFSET)
            < SPACED_STRIDE
        ),
        "window_span": [min(windows), max(windows) + SPACED_STRIDE - 1],
        "pairwise_disjoint": bool(pairwise_disjoint),
    }


def run_spaced_replication(output_dir: Path) -> None:
    """POST-HOC, review round eight: full derived-seed separation.

    Fresh fits at stride-200 chain seeds -- every stream the instrument
    derives (targets, margin, split, truth scorer, null shuffle, fit
    learning) and the external scorer live inside their own row's
    private window, so the bootstrap's independent-replicate basis
    holds by construction, with no exception to argue away. New seeds
    mean new samples: this is a REPLICATION of the B'0 design under a
    corrected allocation, and agreement means directions and verdicts,
    not digits. Labelled post hoc; the frozen record stands as run.
    """

    n_bins = len(FINE_BIN_EDGES) - 1
    rows = []
    fine: dict = {n: [] for n in LADDER}
    for n in LADDER:
        for k in range(B0PRIME_SAMPLES):
            s = _spaced_seed(n, k)
            pi = np.random.default_rng(s).permutation(n)
            causal, times, coords = order_inputs(pi)
            row, coords_fit, targets = analyze_order(
                causal, times, coords, seed=s, want_truth=True,
                return_fit=True,
            )
            out = {"n": float(n), "sample_index": float(k),
                   "chain_seed": float(s),
                   "score_seed": float(s + SPACED_SCORE_OFFSET),
                   "code_version": _diag_code_version(),
                   "status": row.get("status")}
            if row.get("status") == "ok":
                true_x = coords[targets] / float(n)
                # the instrument's internal truth scoring at s+9 is
                # inside this row's private window: collision-free as
                # produced, no re-scoring needed
                out["truth"] = row["truth"]
                restricted, n_eligible = margin_restricted_order_error(
                    coords_fit, true_x, DELTA_MARGIN,
                    num_pair_comparisons=BPRIME_COMPARISONS,
                    seed=s + SPACED_SCORE_OFFSET,
                )
                out["restricted"] = restricted
                out["eligible"] = float(n_eligible)
                out["eligible_floor_met"] = bool(n_eligible >= MIN_ELIGIBLE)
                errors, _counts = margin_binned_order_error(
                    coords_fit, true_x,
                    num_pair_comparisons=BPRIME_COMPARISONS,
                    seed=s + SPACED_SCORE_OFFSET,
                )
                finite = [e for e in errors if not np.isnan(e)]
                out["mix_matched"] = (
                    float(np.mean(errors))
                    if len(finite) == len(errors) else float("nan")
                )
                fine[n].append(fine_binned_order_error(
                    coords_fit, true_x,
                    num_pair_comparisons=BPRIME_COMPARISONS,
                    seed=s + SPACED_SCORE_OFFSET,
                ))
            rows.append(out)
        done = [r for r in rows if r["n"] == n and r["status"] == "ok"]
        print(f"spaced n={n}: {len(done)}/{B0PRIME_SAMPLES} ok | "
              f"median restricted "
              f"{np.median([r['restricted'] for r in done]):.4f}",
              flush=True)
    write_rows_csv(output_dir / "p10_bprime_spaced.csv", rows)

    ok = [r for r in rows if r["status"] == "ok"]
    summary: dict = {
        "code_version": _diag_code_version(),
        "window_map": _spaced_window_map(),
        "delta_margin": DELTA_MARGIN,
        "per_n": {}, "gate": {}, "mix_matched": {}, "curve": {},
    }
    for n in LADDER:
        sub = [r for r in ok if r["n"] == n]
        entry: dict = {"n_ok": len(sub)}
        for quantity in ("restricted", "truth"):
            vals = [r[quantity] for r in sub]
            rng = np.random.default_rng(
                _stable_seed("spaced", quantity, n)
            )
            boots = np.median(np.asarray(vals)[
                rng.integers(0, len(vals), size=(4000, len(vals)))
            ], axis=1) if vals else np.array([])
            entry[f"median_{quantity}"] = (
                float(np.median(vals)) if vals else None
            )
            entry[f"ci_{quantity}"] = [
                float(np.percentile(boots, 2.5)),
                float(np.percentile(boots, 97.5)),
            ] if vals else [None, None]
        entry["min_eligible"] = (
            int(min(r["eligible"] for r in sub)) if sub else None
        )
        summary["per_n"][str(n)] = entry

    floors = all(r["eligible_floor_met"] for r in ok)
    drop, ci = _difference_ci(
        [r["restricted"] for r in ok if r["n"] == LADDER[-1]],
        [r["restricted"] for r in ok if r["n"] == LADDER[0]],
        seed=_stable_seed("spaced-drop"),
    )
    tdrop, tci = _difference_ci(
        [r["truth"] for r in ok if r["n"] == LADDER[-1]],
        [r["truth"] for r in ok if r["n"] == LADDER[0]],
        seed=_stable_seed("spaced-truthdrop"),
    )
    summary["gate"] = {
        "top_minus_bottom_restricted": {"diff": drop, "ci": list(ci)},
        "top_minus_bottom_truth_unrestricted": {
            "diff": tdrop, "ci": list(tci),
        },
        "all_floors_met": bool(floors),
        "yardstick_falls_restricted": bprime_gate_verdict(
            {n: sum(1 for r in ok if r["n"] == n) for n in LADDER},
            floors, drop, tuple(ci),
        ),
    }
    mdrop, mci = _difference_ci(
        [r["mix_matched"] for r in ok if r["n"] == LADDER[-1]],
        [r["mix_matched"] for r in ok if r["n"] == LADDER[0]],
        seed=_stable_seed("spaced-mixdrop"),
    )
    summary["mix_matched"] = {
        "top_minus_bottom": {"diff": mdrop, "ci": list(mci)},
        "rise_survives": bool(
            mdrop is not None and mci[0] is not None and mci[0] > 0.0
        ),
    }

    curves: dict = {}
    curve_entry: dict = {}
    for n in LADDER:
        medians, cis = [], []
        for b in range(n_bins):
            vals = np.array([
                e[b] for e, *_ in fine[n] if not np.isnan(e[b])
            ])
            rng = np.random.default_rng(_stable_seed("spaced-curve", n, b))
            boots = np.median(vals[
                rng.integers(0, vals.size, size=(4000, vals.size))
            ], axis=1) if vals.size else np.array([])
            medians.append(float(np.median(vals)) if vals.size else None)
            cis.append([
                float(np.percentile(boots, 2.5)),
                float(np.percentile(boots, 97.5)),
            ] if vals.size else [None, None])
        curves[n] = (medians, cis)
        curve_entry[str(n)] = {
            "pooled_total_by_bin": [
                int(sum(c[b] for _, c, *_ in fine[n]))
                for b in range(n_bins)
            ],
            "pooled_discordant_by_bin": [
                int(sum(d[b] for _, _, _, d, _ in fine[n]))
                for b in range(n_bins)
            ],
        }
    above, separated = [], []
    for b in range(n_bins):
        m600, ci600 = curves[600][0][b], curves[600][1][b]
        m1200, ci1200 = curves[1200][0][b], curves[1200][1][b]
        if m600 is None or m1200 is None:
            continue
        above.append(bool(m1200 > m600))
        separated.append(bool(ci1200[0] > ci600[1]))
    curve_entry["pointwise_1200_vs_600"] = {
        "bins_compared": len(above),
        "bins_where_1200_above_600": int(sum(above)),
        "bins_ci_separated": int(sum(separated)),
    }
    curve_entry["high_edge"] = HIGH_EDGE
    curve_entry["last_bin"] = {
        str(n): {
            "total": int(sum(c[-1] for _, c, *_ in fine[n])),
            "discordant": int(sum(d[-1] for _, _, _, d, _ in fine[n])),
        } for n in LADDER
    }
    summary["curve"] = curve_entry

    frozen_dir = (Path(__file__).resolve().parents[2]
                  / "docs" / "prereg" / "frozen" / "p10_stage_bprime")
    try:
        frozen = json.loads(
            (frozen_dir / "p10_bprime0_summary.json").read_text(
                encoding="utf-8")
        )
        fr = frozen["top_minus_bottom"]["truth_restricted"]
        ft = frozen["top_minus_bottom"]["truth"]
        summary["comparison_to_frozen"] = {
            "frozen_restricted": fr, "frozen_truth": ft,
            "restricted_sign_replicates": bool(
                drop is not None and fr["diff"] is not None
                and (drop > 0) == (fr["diff"] > 0)
            ),
            "truth_sign_replicates": bool(
                tdrop is not None and ft["diff"] is not None
                and (tdrop > 0) == (ft["diff"] > 0)
            ),
        }
    except OSError:
        summary["comparison_to_frozen"] = None

    def _de_nan(obj):
        if isinstance(obj, float) and np.isnan(obj):
            return None
        if isinstance(obj, dict):
            return {key: _de_nan(val) for key, val in obj.items()}
        if isinstance(obj, list):
            return [_de_nan(val) for val in obj]
        return obj

    summary = _de_nan(summary)
    (output_dir / "p10_bprime_spaced_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"window_map": summary["window_map"],
                      "gate": summary["gate"],
                      "mix_matched": summary["mix_matched"],
                      "comparison_to_frozen":
                          summary["comparison_to_frozen"]}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage",
                        choices=["b0", "b1", "diagnostic", "curve",
                                 "seedfix", "spaced"],
                        default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.stage == "b0":
        run_b0prime(args.output_dir)
    elif args.stage == "diagnostic":
        run_mix_matched_diagnostic(args.output_dir)
    elif args.stage == "curve":
        run_margin_curve_diagnostic(args.output_dir)
    elif args.stage == "seedfix":
        run_seedfix_rescore(args.output_dir)
    elif args.stage == "spaced":
        run_spaced_replication(args.output_dir)
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
