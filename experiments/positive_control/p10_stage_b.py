"""P10 Stage B: the resolution ladder with a scaled instrument.

See docs/prereg/p10_continuum_limit.md Section 8. Stage A found the
frozen-constant instrument scale-covariant: it reads the same relative
window at every size, so the Section 2 outcome table (whose signatures
all presuppose a falling sprinkled-arm curve) was undecidable. Stage B
scales the instrument constants in continuum units, anchored so that
**N = 600 is exactly the frozen instrument** (verified by regression):

    min_chain_len(N) = round(25 * sqrt(N / 600))     chains: fixed length
                                                     in continuum units,
                                                     since the longest
                                                     chain grows ~ 2 sqrt(N)
    max_targets(N)   = round(44 * N / 600)           fixed target density
    min_targets(N)   = round(20 * N / 600)
    train_c(N)       = round(3000 * N / 600)         constraint budget per
    heldout_c(N)     = round(800 * N / 600)          target held fixed --
                                                     the P8 lesson that a
                                                     fixed budget over more
                                                     targets DILUTES the fit
    chain_count      = 6 (fixed: observer count is a choice, not a
                          resolution property)

Stages, in the order the P10-A lessons dictate:

    B0  yardstick pilot -- arm S only, direct sampling, cheap. The
        prerequisite Stage A assumed and then failed: under the scaled
        instrument, does the sprinkled arm's error actually FALL with N?
        If it does not, Stage B stops and that is the finding.
    B1  (only after B0 passes, hypotheses frozen in the prereg first)
        arm E chains, 48 retained samples per chain so the inherited
        ESS >= 20 screen (P7's actual bar at P7's actual m) means
        something, plus TOST equivalence against the preregistered
        margin.

Usage:
    python p10_stage_b.py --stage b0
    python p10_stage_b.py --stage b1 --n 900 --start bipartite
    python p10_stage_b.py --aggregate
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from p3_dynamics import analyze_order
from p5_two_orders_emergence import order_inputs
from p10_continuum_ladder import LADDER, STARTS, _difference_ci, _stable_seed
from pc_common import DEFAULT_OUTPUT_DIR, git_describe, write_rows_csv

from causal_spacetime_lab.positive_control.mcmc_diagnostics import (
    integrated_autocorrelation,
)
from causal_spacetime_lab.positive_control.two_orders import (
    bipartite_perm,
    mcmc_2d_order_fast,
)

BETA = 2.0
EPS_TIMES_N = 12.0
ANCHOR_N = 600

#: Stage B1 chain schedule: 48 retained samples so the inherited screen
#: (ESS >= 20 at m = 48, P7's actual bar at P7's actual m) has room to
#: mean something; 25k spacing was measured decorrelated in Stage 0.
BURN_STEPS = 100_000
SAMPLE_SPACING = 25_000
SAMPLES_PER_CHAIN = 48
TOTAL_STEPS = BURN_STEPS + SAMPLES_PER_CHAIN * SAMPLE_SPACING  # 1.3M
MIN_ESS = 20.0

#: Frozen seed table (prereg Section 8, corrected). The first table used
#: 1000-1050 and was NOT fresh: P6-B's uniform reference orders are
#: exactly ``default_rng(1000..1019).permutation(600)``, so the
#: (600, random) chain's ``seed + 1`` start would have literally reused
#: a P6-B reference order (review finding). 30000-30050 is verified
#: clean against every seed range in the programme, including the
#: derived streams (``seed + 1`` starts, ``seed + 7k`` discriminator
#: seeds up to 30379).
B1_CHAIN_SEEDS = {
    (600, "random"): 30000, (600, "bipartite"): 30010,
    (900, "random"): 30020, (900, "bipartite"): 30030,
    (1200, "random"): 30040, (1200, "bipartite"): 30050,
}
B0_SEED_BASE = {600: 1100, 900: 1120, 1200: 1140}
B0_SAMPLES = 20


def scaled_constants(n: int) -> dict:
    """The instrument family; at ANCHOR_N these ARE the frozen values."""

    ratio = n / ANCHOR_N
    return {
        "chain_count": 6,
        "min_chain_len": round(25 * math.sqrt(ratio)),
        "max_targets": round(44 * ratio),
        "min_targets": round(20 * ratio),
        "train_c": round(3000 * ratio),
        "heldout_c": round(800 * ratio),
    }


def evaluate_perm_scaled(pi: np.ndarray, seed: int) -> dict:
    """One permutation through the SAME pipeline definition, at the
    scaled operating point for its size. P7's G is computed with its
    frozen margins unchanged -- the score definition does not scale."""

    n = pi.size
    causal, times, coords = order_inputs(pi)
    row = analyze_order(
        causal, times, coords, seed=seed, want_truth=True,
        **scaled_constants(n),
    )
    if row.get("status") == "ok":
        m_h = (0.10 - row["heldout"]) / 0.10
        m_n = (row["null_gap"] - 0.10) / 0.10
        m_t = (0.40 - row["truth"]) / 0.40
        row["G"] = float(np.clip(0.5 + min(m_h, m_n, m_t), 0.0, 1.0))
    else:
        row["G"] = 0.0
    return row


def b0_gate_verdict(n_ok_by_rung: dict, drop, drop_ci) -> bool:
    """The 8.2 gate, as a pure function.

    The gate may pass only when BOTH hold: every rung retained its full
    `B0_SAMPLES` successful measurements, and the top-minus-bottom CI
    sits entirely below zero. Without the completeness half, a run in
    which the instrument structurally blocked most samples could still
    authorize B1 from whatever it happened to measure — even one
    surviving sample per endpoint yields a degenerate 'passing' CI
    (review finding). A gate conditioned on the cases the instrument
    could measure is not the preregistered gate.
    """

    complete = all(
        n_ok_by_rung.get(n, 0) == B0_SAMPLES for n in LADDER
    )
    falls = bool(
        drop is not None and drop_ci[1] is not None and drop_ci[1] < 0.0
    )
    return complete and falls


def run_b0(output_dir: Path) -> None:
    """Arm S under the scaled instrument: the yardstick pilot."""

    rows = []
    for n in LADDER:
        base = B0_SEED_BASE[n]
        for k in range(B0_SAMPLES):
            pi = np.random.default_rng(base + k).permutation(n)
            row = {
                "arm": "S-scaled", "n": float(n), "start": "uniform",
                "chain_seed": float(base + k), "sample_index": float(k),
                "code_version": git_describe(),
                **{f"const_{key}": float(value)
                   for key, value in scaled_constants(n).items()},
                **evaluate_perm_scaled(pi, seed=base + k),
            }
            rows.append(row)
        done = [r for r in rows if r["n"] == n and r.get("status") == "ok"]
        med = float(np.median([r["truth"] for r in done])) if done else None
        print(f"B0 n={n}: {len(done)}/{B0_SAMPLES} ok, median truth "
              f"{med if med is None else round(med, 4)}", flush=True)
    write_rows_csv(output_dir / "p10_b0_yardstick.csv", rows)

    ok = [r for r in rows if r.get("status") == "ok"]
    summary: dict = {"code_version": git_describe(), "per_n": {}}
    for n in LADDER:
        vals = [r["truth"] for r in ok if r["n"] == n]
        rng = np.random.default_rng(_stable_seed("b0", n))
        boots = np.median(
            np.asarray(vals)[rng.integers(0, len(vals), size=(4000, len(vals)))],
            axis=1,
        ) if vals else np.array([])
        summary["per_n"][str(n)] = {
            "n_ok": len(vals),
            "median_truth": float(np.median(vals)) if vals else None,
            "ci_truth": [
                float(np.percentile(boots, 2.5)),
                float(np.percentile(boots, 97.5)),
            ] if vals else [None, None],
            "constants": scaled_constants(n),
        }
    top = summary["per_n"][str(LADDER[-1])]
    bottom = summary["per_n"][str(LADDER[0])]
    drop, drop_ci = _difference_ci(
        [r["truth"] for r in ok if r["n"] == LADDER[-1]],
        [r["truth"] for r in ok if r["n"] == LADDER[0]],
        seed=_stable_seed("b0-drop"),
    )
    summary["top_minus_bottom_truth"] = {"diff": drop, "ci": drop_ci}
    summary["all_rungs_complete"] = all(
        summary["per_n"][str(n)]["n_ok"] == B0_SAMPLES for n in LADDER
    )
    summary["yardstick_falls"] = b0_gate_verdict(
        {n: summary["per_n"][str(n)]["n_ok"] for n in LADDER}, drop, drop_ci
    )
    summary["note"] = (
        "B0 is the prerequisite Stage A assumed and then failed: the "
        "scaled instrument must make the sprinkled arm's error FALL "
        "with N (top-minus-bottom CI entirely below zero) before any "
        "B1 hypothesis is frozen or any E chain runs. top/bottom "
        f"medians: {top['median_truth']} vs {bottom['median_truth']}."
    )
    (output_dir / "p10_b0_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


#: The authoritative B0 record. B1 is gated on it (prereg 8.2), and the
#: gate is enforced here in code, the same way P8-B refuses to run
#: without frozen constants -- a documented rule a CLI can bypass is not
#: a rule.
B0_FROZEN_DIR = Path("docs/prereg/frozen/p10_stage_b")
B0_FROZEN_SUMMARY = B0_FROZEN_DIR / "p10_b0_summary.json"
B0_FROZEN_CSV = B0_FROZEN_DIR / "p10_b0_yardstick.csv"


def require_frozen_chain(n: int, start: str, chain_rows: list) -> None:
    """Every aggregated chain must BE one of the six frozen chains.

    Grouping by ``(n, start)`` alone would let a stale or extra CSV --
    a correctly named start produced with a different seed, or a
    seventh chain -- slip into the screen and the hypotheses (review
    finding). Membership in the frozen table and a uniform, matching
    ``chain_seed`` on every row are required; anything else exits.
    """

    if (n, start) not in B1_CHAIN_SEEDS:
        raise SystemExit(
            f"chain ({n}, {start!r}) is not one of the six frozen B1 "
            "chains; remove the extra CSV before aggregating"
        )
    expected = float(B1_CHAIN_SEEDS[(n, start)])
    seeds = {float(r["chain_seed"]) for r in chain_rows}
    if seeds != {expected}:
        raise SystemExit(
            f"chain ({n}, {start!r}) carries chain_seed(s) "
            f"{sorted(seeds)} but the frozen table says {expected}; "
            "this CSV is not the preregistered chain"
        )


def require_all_chains(present: set) -> None:
    """All six preregistered chains, before any hypothesis is evaluated.

    A run aggregated from random-start shards alone would auto-satisfy
    the melt criterion (which needs the bipartite chain to test) and
    report the frozen hypotheses from half the design, unable even to
    detect start disagreement (review finding). Presence is required for
    every (rung, start); whether a present chain *survives* is then the
    screen's job -- the dual-start design demands both be run and
    checked, not that both pass.
    """

    expected = {(n, start) for n in LADDER for start in STARTS}
    missing = sorted(expected - present)
    if missing:
        raise SystemExit(
            "B1 aggregation requires all six preregistered chains before "
            f"any hypothesis is evaluated; missing: {missing}"
        )


def _require_b0_gate() -> None:
    if not B0_FROZEN_SUMMARY.exists():
        raise SystemExit(
            f"frozen B0 record not found at {B0_FROZEN_SUMMARY}; B1 may only "
            "run after B0 has been run, frozen, and PASSED (prereg 8.2)"
        )
    record = json.loads(B0_FROZEN_SUMMARY.read_text(encoding="utf-8"))
    if record.get("yardstick_falls") is not True:
        raise SystemExit(
            "the frozen B0 record says yardstick_falls = "
            f"{record.get('yardstick_falls')!r}: the 8.2 gate FAILED, no B1 "
            "chain may run, and Stage B is closed at 8.5. This refusal is "
            "the preregistration operating, not an error."
        )


def chain_is_complete(chain_rows: list) -> bool:
    """True iff a B1 chain retained every scheduled sample.

    The frozen sampler draws ``(i, j)`` BEFORE its sampling check and
    ``continue``s on ``i == j``, so a scheduled retention point silently
    vanishes with probability ``1/N`` — about a 4-8% event per 48-sample
    chain at these sizes (review finding). The sampler cannot be touched
    (P5's exact-replay validation pins its trajectory), so completeness
    is enforced here: a shortened chain no longer carries the ``m = 48``
    diagnostic basis the screen claims, and its nominal step labels are
    shifted from the true iterations, so it must not enter the frozen
    hypotheses.
    """

    indices = sorted(
        int(float(r.get("sample_index", -1))) for r in chain_rows
    )
    return (
        len(chain_rows) == SAMPLES_PER_CHAIN
        # exactly 0..47, each once: a row count alone would accept a
        # duplicated row papering over a missing scheduled state, which
        # biases the medians while wearing the right length (review
        # finding)
        and indices == list(range(SAMPLES_PER_CHAIN))
        and all(str(r.get("chain_complete")) in ("True", "1.0", "1")
                for r in chain_rows)
    )


def run_b1_chain(n: int, start: str, output_dir: Path) -> None:
    seed = B1_CHAIN_SEEDS[(n, start)]
    eps = EPS_TIMES_N / n
    pi0 = (bipartite_perm(n) if start == "bipartite"
           else np.random.default_rng(seed + 1).permutation(n))
    samples, acceptance, perms = mcmc_2d_order_fast(
        pi0, beta=BETA, eps=eps, steps=TOTAL_STEPS, seed=seed,
        sample_every=SAMPLE_SPACING, burn_frac=BURN_STEPS / TOTAL_STEPS,
        collect_perms=True,
    )
    complete = len(samples) == SAMPLES_PER_CHAIN
    if not complete:
        print(f"B1 n={n} {start}: SHORT CHAIN -- {len(samples)} of "
              f"{SAMPLES_PER_CHAIN} scheduled samples retained (an i == j "
              "draw landed on a retention point); nominal steps are "
              "unreliable and the screen will exclude this chain",
              flush=True)
    rows = []
    for k, (obs, pi) in enumerate(zip(samples, perms, strict=True)):
        row = {
            "arm": "E-scaled", "n": float(n), "start": start,
            "chain_seed": float(seed), "sample_index": float(k),
            # the nominal schedule; only trustworthy when the chain is
            # complete, which chain_complete records per row
            "step": (float(BURN_STEPS + k * SAMPLE_SPACING)
                     if complete else float("nan")),
            "chain_complete": bool(complete),
            "n0": float(obs["n0"]), "action": float(obs["S"]),
            "acceptance": float(acceptance),
            "code_version": git_describe(),
            **evaluate_perm_scaled(pi, seed=seed + 7 * k),
        }
        rows.append(row)
        print(f"B1 n={n} {start} sample {k}: {row.get('status')} "
              f"G={row['G']:.3f}", flush=True)
    write_rows_csv(output_dir / f"p10_b1_n{n}_{start}.csv", rows)


#: The frozen 8.3 equivalence margin: P2/P9's standing saturation
#: tolerance, chosen a priori (see the prereg for the provenance).
TOST_MARGIN = 0.05


def _truth_block_length(series: list) -> int:
    """Block length for one chain's truth series: ceil(2 tau), tau from
    the shared Geyer diagnostic ON THE TRUTH SERIES itself.

    The mixing screen is frozen on n0 and stays frozen; but n0 mixing
    does not certify truth mixing, and an iid bootstrap over correlated
    truth values narrows the CI by roughly sqrt(tau) -- enough to
    falsely support an equivalence hypothesis (review finding). The
    dependence is handled where the uncertainty is computed instead.
    """

    tau = integrated_autocorrelation(list(series))[0]
    return max(1, math.ceil(2.0 * tau))


def _block_resample_median(chains: list, rng, draws: int) -> np.ndarray:
    """Bootstrap medians preserving within-chain dependence.

    Circular moving-block bootstrap per chain -- blocks never cross a
    chain boundary -- with each chain's block length set from its own
    truth autocorrelation. Chains of iid values reduce to the ordinary
    bootstrap (block length 1).
    """

    lengths = [_truth_block_length(c) for c in chains]
    out = np.empty(draws)
    arrays = [np.asarray(c, dtype=float) for c in chains]
    for d in range(draws):
        pieces = []
        for arr, b in zip(arrays, lengths, strict=True):
            m = arr.size
            n_blocks = math.ceil(m / b)
            starts = rng.integers(0, m, size=n_blocks)
            idx = (starts[:, None] + np.arange(b)[None, :]) % m
            pieces.append(arr[idx.ravel()][:m])
        out[d] = float(np.median(np.concatenate(pieces)))
    return out


def _dependent_difference_ci(e_chains: list, s_values: list, seed,
                             draws: int = 4000):
    """median(E) - median(S) with arm-E dependence preserved.

    Arm E enters as per-chain series and is block-resampled; arm S is
    direct sampling (genuinely iid) and is resampled ordinarily.
    """

    flat_e = [v for c in e_chains for v in c]
    s = np.asarray(s_values, dtype=float)
    if not flat_e or s.size == 0:
        return None, (None, None)
    rng = np.random.default_rng(seed)
    e_meds = _block_resample_median(e_chains, rng, draws)
    s_meds = np.median(
        s[rng.integers(0, s.size, size=(draws, s.size))], axis=1
    )
    diffs = e_meds - s_meds
    return float(np.median(flat_e) - np.median(s)), (
        float(np.percentile(diffs, 2.5)),
        float(np.percentile(diffs, 97.5)),
    )


def _dependent_e_vs_e_ci(top_chains: list, bottom_chains: list, seed,
                         draws: int = 4000):
    """H-DEEPEN's top-minus-bottom, block-resampled on both sides."""

    flat_top = [v for c in top_chains for v in c]
    flat_bottom = [v for c in bottom_chains for v in c]
    if not flat_top or not flat_bottom:
        return None, (None, None)
    rng = np.random.default_rng(seed)
    diffs = (
        _block_resample_median(top_chains, rng, draws)
        - _block_resample_median(bottom_chains, rng, draws)
    )
    return float(np.median(flat_top) - np.median(flat_bottom)), (
        float(np.percentile(diffs, 2.5)),
        float(np.percentile(diffs, 97.5)),
    )


def evaluate_frozen_hypotheses(
    e_chains_by_n: dict, s_truth_by_n: dict, seed_fn=_stable_seed
) -> dict:
    """The preregistered 8.3 decision block, as a pure function.

    ``e_chains_by_n`` maps each rung to the surviving arm-E chains,
    each a list of truth values IN SAMPLE ORDER -- the chain structure
    is required, not decoration, because arm-E uncertainty is computed
    by a within-chain circular block bootstrap (block length from each
    chain's own truth autocorrelation). Flat pooling with an iid
    bootstrap narrows the CI by roughly sqrt(tau) under the dependence
    the ESS >= 20 screen still admits, which could falsely support an
    equivalence hypothesis (review finding). ``s_truth_by_n`` holds the
    B0 arm-S values, genuinely iid by direct sampling.

    * H-TRACK: at EVERY rung, the E - S median-difference 95% CI must
      lie inside [-TOST_MARGIN, +TOST_MARGIN]. A rung with no surviving
      E samples cannot be evaluated and H-TRACK is NOT supported.
    * H-DEEPEN: arm E's own top-rung minus bottom-rung median difference
      has its 95% CI entirely below zero.
    """

    rungs = sorted(e_chains_by_n)
    per_rung = {}
    track_passes = []
    for n in rungs:
        diff, ci = _dependent_difference_ci(
            e_chains_by_n[n], s_truth_by_n.get(n, []),
            seed=seed_fn("b1-tost", n),
        )
        evaluable = diff is not None
        tost_pass = bool(
            evaluable
            and -TOST_MARGIN <= ci[0]
            and ci[1] <= TOST_MARGIN
        )
        per_rung[str(n)] = {
            "E_minus_S_diff": diff, "ci": list(ci),
            "evaluable": evaluable, "tost_pass": tost_pass,
        }
        track_passes.append(tost_pass)
    h_track = bool(rungs) and all(track_passes)

    top, bottom = (rungs[-1], rungs[0]) if rungs else (None, None)
    deepen_diff, deepen_ci = _dependent_e_vs_e_ci(
        e_chains_by_n.get(top, []), e_chains_by_n.get(bottom, []),
        seed=seed_fn("b1-deepen"),
    ) if rungs else (None, (None, None))
    h_deepen = bool(
        deepen_diff is not None and deepen_ci[1] is not None
        and deepen_ci[1] < 0.0
    )
    return {
        "tost_margin": TOST_MARGIN,
        "per_rung": per_rung,
        "h_track_supported": h_track,
        "h_deepen_diff": deepen_diff,
        "h_deepen_ci": list(deepen_ci),
        "h_deepen_supported": h_deepen,
        "consistent_with_continuum_limit": bool(h_track and h_deepen),
    }


def aggregate(output_dir: Path) -> None:
    """B1 aggregation: inherited diagnostics, strict band, both frozen
    criteria combined -- every Stage A review lesson pre-applied -- and
    the frozen 8.3 hypotheses evaluated via the pure decision block
    above. Gated on the B0 record exactly as running B1 is: evaluating
    hypotheses on post-gate-invalid data is as bad as generating it."""

    import csv

    _require_b0_gate()
    rows = []
    for path in sorted(output_dir.glob("p10_b1_n*.csv")):
        rows.extend(csv.DictReader(path.open(encoding="utf-8")))
    if not rows:
        raise SystemExit("no B1 chains found; run --stage b1 first")

    summary: dict = {"code_version": git_describe(), "chains": [],
                     "per_n": {}}
    screen_pass: dict = {}
    for n in LADDER:
        chains: dict = {}
        for r in rows:
            if int(float(r["n"])) == n:
                chains.setdefault(r["start"], []).append(r)
        if not chains:
            continue
        random_tail = [
            float(r["n0"])
            for r in sorted(chains.get("random", []),
                            key=lambda r: float(r["sample_index"]))
            [SAMPLES_PER_CHAIN // 2:]
        ]
        band = (min(random_tail), max(random_tail)) if random_tail else None
        for start, chain_rows in chains.items():
            require_frozen_chain(n, start, chain_rows)
            chain_rows.sort(key=lambda r: float(r["sample_index"]))
            n0s = [float(r["n0"]) for r in chain_rows]
            ess = integrated_autocorrelation(n0s)[1]
            melted = (
                band is not None and band[0] <= n0s[0] <= band[1]
            ) if start == "bipartite" else True
            complete = chain_is_complete(chain_rows)
            # Retention completeness is not measurement completeness: a
            # chain can keep all 48 MCMC states while the discriminator
            # structurally blocks on some of them, silently shrinking
            # the basis the hypotheses are computed from (review
            # finding). Every retained state must also have evaluated.
            n_evals_ok = sum(
                1 for r in chain_rows if r.get("status") == "ok"
            )
            evals_ok = n_evals_ok == SAMPLES_PER_CHAIN
            screen_pass[(n, start)] = (
                bool(ess >= MIN_ESS) and bool(melted)
                and complete and evals_ok
            )
            summary["chains"].append({
                "n": n, "start": start, "ess_n0": round(float(ess), 1),
                "ess_pass": bool(ess >= MIN_ESS),
                "first_sample_in_random_band": bool(melted),
                "chain_complete": complete,
                "n_retained": len(chain_rows),
                "n_evals_ok": n_evals_ok,
                "evals_complete": evals_ok,
                "screen_pass": screen_pass[(n, start)],
                "acceptance": float(chain_rows[0]["acceptance"]),
            })

    import csv as _csv

    # The FROZEN B0 artifacts, not whatever sits in output_dir: the gate
    # authorizes against the frozen record, so the hypotheses must be
    # evaluated against the same pilot that passed it -- a stale or
    # rerun output-directory copy would silently decouple the two
    # (review finding).
    require_all_chains(set(screen_pass))
    b0 = json.loads(B0_FROZEN_SUMMARY.read_text(encoding="utf-8"))
    b0_rows = list(_csv.DictReader(
        B0_FROZEN_CSV.open(encoding="utf-8")
    ))
    e_chains_by_n: dict = {}
    s_truth_by_n: dict = {}
    for n in LADDER:
        # chain structure preserved, in sample order: the block
        # bootstrap must never resample across a chain boundary
        e_chains_by_n[n] = []
        for start in STARTS:
            if not screen_pass.get((n, start), False):
                continue
            chain_rows = sorted(
                (r for r in rows
                 if int(float(r["n"])) == n and r["start"] == start
                 and r.get("status") == "ok"),
                key=lambda r: float(r["sample_index"]),
            )
            if chain_rows:
                e_chains_by_n[n].append(
                    [float(r["truth"]) for r in chain_rows]
                )
        s_truth_by_n[n] = [
            float(r["truth"]) for r in b0_rows
            if int(float(r["n"])) == n and r.get("status") == "ok"
        ]
        flat = [v for c in e_chains_by_n[n] for v in c]
        summary["per_n"][str(n)] = {
            "E_n_ok": len(flat),
            "E_median_truth": float(np.median(flat)) if flat else None,
            "S_median_truth": b0["per_n"][str(n)]["median_truth"],
            "E_truth_block_lengths": [
                _truth_block_length(c) for c in e_chains_by_n[n]
            ],
        }
    summary["frozen_hypotheses"] = evaluate_frozen_hypotheses(
        e_chains_by_n, s_truth_by_n
    )

    (output_dir / "p10_b1_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["b0", "b1"], default=None)
    parser.add_argument("--n", type=int, choices=LADDER, default=None)
    parser.add_argument("--start", choices=STARTS, default=None)
    parser.add_argument("--aggregate", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.aggregate:
        aggregate(args.output_dir)
    elif args.stage == "b0":
        run_b0(args.output_dir)
    elif args.stage == "b1":
        if args.n is None or args.start is None:
            raise SystemExit("--stage b1 requires --n and --start")
        _require_b0_gate()
        run_b1_chain(args.n, args.start, args.output_dir)
    else:
        raise SystemExit("choose --stage b0/--stage b1/--aggregate")


if __name__ == "__main__":
    main()
