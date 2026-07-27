"""P11: the metric instrument for the continuum limit — stage runners.

Implements docs/prereg/p11_continuum_metric.md (v1.9 at this
writing) and nothing else.
Every estimator goes through the repository's shared definitions
(longest_chain_length, estimate_tau_from_longest_chain_1p1,
estimate_tau_from_interval_count) — the shared-definition rule; the
frozen coordinate convention is the unit-diamond rank grid
(u, v) = (i, pi(i)) / N with density rho = N/2 in (t, x).

Stage order is gated: verify (completeness pin + wall times) must pass
before pilot; pilot must be feasible before Stage A. Stage B and C
refuse by design until their addenda are frozen.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from p5_two_orders_emergence import order_inputs
from p10_continuum_ladder import _stable_seed
from p10_stage_bprime import _diag_code_version
from pc_common import DEFAULT_OUTPUT_DIR, write_rows_csv

from causal_spacetime_lab.chains import longest_chain_length
from causal_spacetime_lab.estimators import estimate_tau_from_longest_chain_1p1
from causal_spacetime_lab.metrics import estimate_tau_from_interval_count

#: Ladder and pair protocol (prereg Sections 1.1 and 4).
P11_LADDER = (600, 1200, 2400)
TAU_BAND = (0.35, 0.45)
K_PAIRS = 6
DRAW_REJECTION_CAP = 200
SKIP_CAP = 20

#: Power design (prereg Section 1.2; the literals ARE the frozen spec).
DELTA_STAR = -0.2007
DELTA_EQ = 0.067
N_SUP_COEFF = 260.9
N_EQ_COEFF = 2895.2
#: v1.8: the v1.7 chi-square factor was a Gaussian pivot — invalid
#: for the non-Gaussian y this experiment scores. Each rung's
#: variance now takes its kurtosis-adaptive Bonett upper bound at
#: one-sided 95% (z below), and the two bounds are summed so the
#: pair-level guarantee is >= 90% by Bonferroni, no normal model
#: anywhere. The pilot grows to 200 per endpoint rung to resolve the
#: kurtosis (~20 seconds by the frozen verification record).
BONETT_Z = 1.6449
N_FLOOR = 12
N_CAP = 60
PROJECTION_LIMIT_HOURS = 12.0

#: Seed windows (prereg Section 5, v1.6): stride-200 private windows,
#: pilots fill 12 of 16 slots, stage rungs fill up to 60 of 80.
STRIDE = 200
SCORE_OFFSET = 150
PILOT_BLOCKS = {600: (200000, 220), 2400: (244000, 220)}
STAGE_A_BLOCKS = {600: (288000, 80), 1200: (304000, 80),
                  2400: (320000, 80)}
PILOT_SAMPLES = 200

#: Verification block (prereg Sections 4 and 5): consecutive seeds,
#: single generator per sample, discarded after the pin.
VERIFY_BASE = {600: 190000, 1200: 192000, 2400: 194000}
VERIFY_COUNT = 2000
VERIFY_PIN = 1998

VERIFY_ARTIFACT = "p11_verification_summary.json"
PILOT_ARTIFACT = "p11_pilot_summary.json"

#: Stage B (the Section 10 addendum): dual-box chain estimator on
#: spacelike pairs, same band, same protocol, derived
#: Delta*_B == DELTA_STAR (same LIS law on a rank box of expected
#: count N d^2/4). Frozen windows from Section 5; the B verification
#: block from Section 10.
PILOT_B_BLOCKS = {600: (336000, 220), 2400: (380000, 220)}
STAGE_B_BLOCKS = {600: (424000, 80), 1200: (440000, 80),
                  2400: (456000, 80)}
VERIFY_B_BASE = {600: 600000, 1200: 602000, 2400: 604000}
VERIFY_B_ARTIFACT = "p11_verification_b_summary.json"
PILOT_B_ARTIFACT = "p11_pilot_b_summary.json"


def continuum_uv(pi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The frozen convention: (u, v) = (i, pi(i)) / N on the unit square.

    Derived from the shared dictionary's (t, x) = (i + pi, i - pi) so
    the construction cannot drift from order_inputs: u = (t + x) / 2N,
    v = (t - x) / 2N.
    """

    _causal, t, x = order_inputs(pi)
    n = float(pi.size)
    return (t + x) / (2.0 * n), (t - x) / (2.0 * n)


def tau_true_pairs(u: np.ndarray, v: np.ndarray, first, second) -> np.ndarray:
    du = u[second] - u[first]
    dv = v[second] - v[first]
    return 2.0 * np.sqrt(du * dv)


def eligible_pool(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """All ordered related pairs whose true proper time lies in the
    frozen band — the precomputed candidate pool of Section 4 (v1.4).
    Returns an array of (first, second) index pairs."""

    du = u[None, :] - u[:, None]
    dv = v[None, :] - v[:, None]
    related = (du > 0) & (dv > 0)
    tau = np.where(related, 2.0 * np.sqrt(np.abs(du * dv)), np.nan)
    keep = related & (tau >= TAU_BAND[0]) & (tau <= TAU_BAND[1])
    first, second = np.nonzero(keep)
    return np.column_stack([first, second])


def supports_disjoint(box_a, box_b) -> bool:
    """Stage A supports are the pairs' (u, v) bounding boxes (== the
    closed causal intervals); two boxes are disjoint iff they are
    disjoint in u OR in v."""

    au0, au1, av0, av1 = box_a
    bu0, bu1, bv0, bv1 = box_b
    u_disjoint = (au1 < bu0) or (bu1 < au0)
    v_disjoint = (av1 < bv0) or (bv1 < av0)
    return u_disjoint or v_disjoint


def draw_disjoint_pairs(pool: np.ndarray, u: np.ndarray, v: np.ndarray,
                        rng: np.random.Generator):
    """Greedy accept from a uniform without-replacement ordering of the
    eligible pool; the frozen cap counts support-overlap REJECTIONS
    only (Section 4, v1.4). Returns (pairs, complete, rejections)."""

    accepted: list = []
    boxes: list = []
    rejections = 0
    if pool.shape[0] == 0:
        return accepted, False, 0
    for idx in rng.permutation(pool.shape[0]):
        i, j = int(pool[idx, 0]), int(pool[idx, 1])
        # v-normalized so the same box works for timelike pairs
        # (v[i] < v[j]) and Stage B's spacelike duals (v[i] > v[j]);
        # for timelike pairs this is the identity
        box = (u[i], u[j], min(v[i], v[j]), max(v[i], v[j]))
        if all(supports_disjoint(box, other) for other in boxes):
            accepted.append((i, j))
            boxes.append(box)
            if len(accepted) == K_PAIRS:
                return accepted, True, rejections
        else:
            rejections += 1
            if rejections >= DRAW_REJECTION_CAP:
                return accepted, False, rejections
    return accepted, False, rejections


def score_pair(u: np.ndarray, v: np.ndarray, i: int, j: int) -> dict:
    """Both estimators on one pair, through the shared definitions."""

    n = u.size
    inside = ((u > u[i]) & (u < u[j]) & (v > v[i]) & (v < v[j]))
    interior = np.flatnonzero(inside)
    members = np.concatenate(([i], interior, [j]))
    su, sv = u[members], v[members]
    causal_sub = (su[:, None] < su[None, :]) & (sv[:, None] < sv[None, :])
    chain = longest_chain_length(
        causal_sub, start=0, end=members.size - 1, event_times=su,
    )
    rho = n / 2.0
    tau_chain = estimate_tau_from_longest_chain_1p1(chain, rho=rho)
    tau_vol = estimate_tau_from_interval_count(int(interior.size), rho=rho)
    du, dv = u[j] - u[i], v[j] - v[i]
    tau_ref = 2.0 * float(np.sqrt(du * dv))
    # endpoint-conditioned expected interior count (Section 3, v1.5):
    # rank gaps a = N du, b = N dv
    a, b = round(n * du), round(n * dv)
    m_cond = (a - 1) * (b - 1) / (n - 2)
    return {
        "tau_true": tau_ref,
        "tau_chain": float(tau_chain),
        "tau_vol": float(tau_vol),
        "relerr_chain": abs(float(tau_chain) - tau_ref) / tau_ref,
        "relerr_vol": abs(float(tau_vol) - tau_ref) / tau_ref,
        "m_open": int(interior.size),
        "m_conditioned": float(m_cond),
        "chain_length": int(chain),
    }


def run_sample(n: int, seed: int, single_stream: bool = False):
    """One sample: permutation, pool, draw, score. Returns (record,
    complete). ``single_stream`` is the verification mode of Section 5:
    one generator for permutation and draws, no derived offsets."""

    rng = np.random.default_rng(seed)
    pi = rng.permutation(n)
    u, v = continuum_uv(pi)
    pool = eligible_pool(u, v)
    draw_rng = rng if single_stream else np.random.default_rng(
        seed + SCORE_OFFSET
    )
    pairs, complete, rejections = draw_disjoint_pairs(pool, u, v, draw_rng)
    record = {
        "n": float(n), "seed": float(seed),
        "pool_size": float(pool.shape[0]),
        "rejections": float(rejections),
        "complete": bool(complete),
    }
    if not complete:
        # y is never computed for an incomplete sample -- the frozen
        # order of operations that keeps completeness decisions ahead
        # of any estimator evaluation (Section 4, v1.6)
        return record, False
    scored = [score_pair(u, v, i, j) for i, j in pairs]
    record["y"] = float(np.log10(
        np.median([s["relerr_chain"] for s in scored])
    ))
    record["y_vol"] = float(np.log10(
        np.median([s["relerr_vol"] for s in scored])
    ))
    record["median_relerr_chain"] = float(
        np.median([s["relerr_chain"] for s in scored])
    )
    record["mean_m_conditioned"] = float(
        np.mean([s["m_conditioned"] for s in scored])
    )
    record["min_m_open"] = float(min(s["m_open"] for s in scored))
    return record, True


def fill_block(n: int, base: int, slots: int, needed: int,
               sample_fn=None):
    """Seeds in window order, first ``needed`` complete (Section 4,
    v1.6). Returns (records, skipped_seeds, filled)."""

    runner = sample_fn if sample_fn is not None else run_sample
    records, skipped = [], []
    for k in range(slots):
        seed = base + STRIDE * k
        record, complete = runner(n, seed)
        if complete:
            records.append(record)
            if len(records) == needed:
                return records, skipped, True
        else:
            skipped.append(seed)
            if len(skipped) > SKIP_CAP:
                return records, skipped, False
    return records, skipped, False


def verdict(lo: float, hi: float, flat_available: bool) -> str:
    """The frozen table of Section 1.4: rows in order, first match."""

    if hi < 0.0:
        return "IMPROVES"
    if lo > 0.0:
        return "DEGRADES"
    if flat_available and (-DELTA_EQ < lo) and (hi < DELTA_EQ):
        return "FLAT-WITHIN-MARGIN"
    return "INCONCLUSIVE" if flat_available else "UNRESOLVED"


def bonett_variance_bound(y: np.ndarray, z: float = BONETT_Z) -> dict:
    """Kurtosis-adaptive upper bound on the variance (Section 1.2) —
    Bonett's log-variance interval at nominal one-sided 95%. This is
    an ASYMPTOTIC approximation, not a finite-sample guarantee (v1.9,
    review); the calibrated wrapper below validates it against the
    realized distribution."""

    y = np.asarray(y, dtype=float)
    n = y.size
    centered = y - y.mean()
    s2 = float(np.sum(centered ** 2) / (n - 1))
    g4 = float(n * np.sum(centered ** 4) / np.sum(centered ** 2) ** 2)
    se = float(np.sqrt((g4 - (n - 3) / n) / (n - 1)))
    return {"s2": s2, "g4": g4, "se": se,
            "bound": float(s2 * np.exp(z * se))}


CALIBRATION_RESAMPLES = 2000
CALIBRATION_TARGET = 0.95


def calibrated_variance_bound(y: np.ndarray, label: str) -> dict:
    """Section 1.2 (v1.9): the Bonett bound, VALIDATED on the realized
    distribution instead of trusted by construction. Bootstrap
    coverage check: resample the pilot, ask how often the resample's
    bound covers the full-pilot variance; if the nominal z falls short
    of the target, raise z until coverage is met (monotone). The
    residual limitation is stated in the prereg: tail mass absent from
    the 200 pilot observations is beyond any construction's reach."""

    y = np.asarray(y, dtype=float)
    n = y.size
    base = bonett_variance_bound(y)
    rng = np.random.default_rng(_stable_seed("p11-bonett-cal", label))
    resamples = [
        y[rng.integers(0, n, size=n)] for _ in range(CALIBRATION_RESAMPLES)
    ]

    def coverage(z: float) -> float:
        hits = sum(
            1 for r in resamples
            if bonett_variance_bound(r, z=z)["bound"] >= base["s2"]
        )
        return hits / CALIBRATION_RESAMPLES

    coverage_nominal = coverage(BONETT_Z)
    z_used = BONETT_Z
    if coverage_nominal < CALIBRATION_TARGET:
        lo, hi = BONETT_Z, 8.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if coverage(mid) >= CALIBRATION_TARGET:
                hi = mid
            else:
                lo = mid
        z_used = hi
    return {
        **bonett_variance_bound(y, z=z_used),
        "coverage_at_nominal": coverage_nominal,
        "z_used": float(z_used),
        "calibrated": bool(z_used > BONETT_Z),
    }


def power_requirements(y_b: np.ndarray, y_t: np.ndarray) -> dict:
    """Section 1.2 (v1.9): per-rung calibrated Bonett bounds summed
    (Bonferroni at the nominal levels — labelled calibrated-
    approximate, not exact), frozen literals, the selection rule, and
    the ex-ante flat-availability declaration."""

    bound_b = calibrated_variance_bound(y_b, "bottom")
    bound_t = calibrated_variance_bound(y_t, "top")
    s2 = bound_b["s2"] + bound_t["s2"]
    s2_90 = bound_b["bound"] + bound_t["bound"]
    calibration = {
        "bottom": {k: bound_b[k] for k in
                   ("coverage_at_nominal", "z_used", "calibrated")},
        "top": {k: bound_t[k] for k in
                ("coverage_at_nominal", "z_used", "calibrated")},
    }
    n_sup = int(np.ceil(N_SUP_COEFF * s2_90))
    n_eq = int(np.ceil(N_EQ_COEFF * s2_90))
    clamp = lambda k: int(min(max(k, N_FLOOR), N_CAP))  # noqa: E731
    if n_sup > N_CAP:
        return {"s2": s2, "s2_90": s2_90, "calibration": calibration,
                "n_sup": n_sup, "n_eq": n_eq,
                "n_per_rung": None, "flat_available": False,
                "infeasible": True}
    if n_eq <= N_CAP:
        return {"s2": s2, "s2_90": s2_90, "calibration": calibration,
                "n_sup": n_sup, "n_eq": n_eq,
                "n_per_rung": clamp(max(n_sup, n_eq)),
                "flat_available": True, "infeasible": False}
    return {"s2": s2, "s2_90": s2_90, "calibration": calibration,
            "n_sup": n_sup, "n_eq": n_eq,
            "n_per_rung": clamp(n_sup), "flat_available": False,
            "infeasible": False}


def run_verify(output_dir: Path) -> None:
    """The completeness pin and the measured per-rung wall times
    (Sections 4 and 1.3). Outcomes are computed to price the pipeline
    and then discarded; only completeness and time are recorded."""

    stamp = _preflight_clean()
    summary: dict = {"code_version": stamp,
                     "pin_required": VERIFY_PIN, "per_n": {}}
    for n in P11_LADDER:
        base = VERIFY_BASE[n]
        complete_count = 0
        start = time.perf_counter()
        for k in range(VERIFY_COUNT):
            _record, complete = run_sample(n, base + k, single_stream=True)
            complete_count += int(complete)
        elapsed = time.perf_counter() - start
        summary["per_n"][str(n)] = {
            "complete": complete_count, "total": VERIFY_COUNT,
            "mean_seconds_per_sample": elapsed / VERIFY_COUNT,
        }
        print(f"verify n={n}: {complete_count}/{VERIFY_COUNT} complete | "
              f"{elapsed / VERIFY_COUNT:.3f} s/sample", flush=True)
    summary["pin_passed"] = bool(all(
        summary["per_n"][str(n)]["complete"] >= VERIFY_PIN
        for n in P11_LADDER
    ))
    (output_dir / VERIFY_ARTIFACT).write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pin_passed": summary["pin_passed"]}, indent=2))


def _preflight_clean() -> str:
    """Stages run only from clean commits (prereg Section 8); a dirty
    stamp aborts BEFORE any quarantined seed window is consumed
    (review: running dirty burned windows and merely marked the
    artifact, which is a record of the damage, not a guard)."""

    stamp = _diag_code_version()
    if stamp.endswith("-dirty"):
        raise SystemExit(
            f"working tree is dirty ({stamp}) -- stages run only from "
            "clean commits; commit first, then run."
        )
    if not stamp or stamp.startswith("unknown"):
        # git failed (source archive, copied directory): "unknown"
        # carries no provenance, and two unknown stamps would even
        # pass the prerequisite equality check (review) -- refuse.
        raise SystemExit(
            "no git provenance available (stamp 'unknown') -- stages "
            "run only from a clean git checkout."
        )
    return stamp


def _load_gate_artifact(output_dir: Path, name: str, stamp: str) -> dict:
    """Prerequisites must exist AND have been produced by the current
    implementation: the artifact's stamp must equal this clean HEAD's.
    This also catches a crashed rerun leaving a stale predecessor --
    its old stamp cannot equal the current one."""

    path = output_dir / name
    if not path.exists():
        raise SystemExit(
            f"{name} not found -- the prerequisite stage has not run; "
            "the stage order is frozen (verify -> pilot -> a)."
        )
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("code_version") != stamp:
        raise SystemExit(
            f"{name} was produced at {artifact.get('code_version')} but "
            f"HEAD is {stamp} -- the implementation changed after the "
            "prerequisite ran; regenerate it (verify and pilot are "
            "cheap by the frozen wall-time record)."
        )
    return artifact


def run_pilot(output_dir: Path) -> None:
    """Stage P (Section 1.3, v1.6): BOTH endpoint rungs, variance and
    wall time only. Cross-rung statistics are forbidden: nothing here
    computes a mean difference, and the artifact holds per-rung SDs
    and times only."""

    stamp = _preflight_clean()
    verification = _load_gate_artifact(output_dir, VERIFY_ARTIFACT, stamp)
    if not verification.get("pin_passed"):
        raise SystemExit("verification pin failed -- pilot refuses to run")

    summary: dict = {"code_version": stamp, "per_rung": {}}
    ys: dict = {}
    for n, (base, slots) in PILOT_BLOCKS.items():
        start = time.perf_counter()
        records, skipped, filled = fill_block(n, base, slots, PILOT_SAMPLES)
        elapsed = time.perf_counter() - start
        if not filled:
            raise SystemExit(
                f"pilot block n={n} could not fill {PILOT_SAMPLES} "
                f"complete samples (skips: {len(skipped)}) -- "
                "INFEASIBLE-INCOMPLETE"
            )
        ys[n] = np.array([r["y"] for r in records])
        rung_bound = bonett_variance_bound(ys[n])
        summary["per_rung"][str(n)] = {
            "n_samples": len(records),
            "variance": rung_bound["s2"],
            "kurtosis_g4": rung_bound["g4"],
            "variance_bound_95": rung_bound["bound"],
            "skipped_seeds": skipped,
            "mean_seconds_per_sample": elapsed / len(records),
            # the raw per-sample statistics behind the variance and
            # calibration (review: without them the artifact could
            # not be audited or recalibrated)
            "y": [float(val) for val in ys[n]],
        }
        print(f"pilot n={n}: var={rung_bound['s2']:.5f} "
              f"g4={rung_bound['g4']:.2f} "
              f"bound={rung_bound['bound']:.5f} | skips={len(skipped)}",
              flush=True)

    power = power_requirements(ys[600], ys[2400])
    times = verification["per_n"]
    projected_hours = None
    if power["n_per_rung"] is not None:
        projected_hours = power["n_per_rung"] * sum(
            times[str(n)]["mean_seconds_per_sample"] for n in P11_LADDER
        ) / 3600.0
    summary["power"] = power
    summary["projected_stage_a_hours"] = projected_hours
    summary["feasible"] = bool(
        not power["infeasible"]
        and projected_hours is not None
        and projected_hours <= PROJECTION_LIMIT_HOURS
    )
    (output_dir / PILOT_ARTIFACT).write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"power": power,
                      "projected_stage_a_hours": projected_hours,
                      "feasible": summary["feasible"]}, indent=2))


def run_stage_a(output_dir: Path) -> None:
    """Stage A (Sections 1.4 and 6): the primary gate, at the pilot's
    n, with the frozen verdict table; slope, constant-level, and
    middle-rung checks ride along labelled, never gating."""

    stamp = _preflight_clean()
    pilot = _load_gate_artifact(output_dir, PILOT_ARTIFACT, stamp)
    if not pilot.get("feasible"):
        raise SystemExit("pilot declared the design infeasible -- "
                         "Stage A refuses to run")
    n_per_rung = int(pilot["power"]["n_per_rung"])
    flat_available = bool(pilot["power"]["flat_available"])

    rows, per_rung_y, per_rung_y_vol = [], {}, {}
    skip_counts, skipped_seeds = {}, {}
    for n, (base, slots) in STAGE_A_BLOCKS.items():
        records, skipped, filled = fill_block(n, base, slots, n_per_rung)
        if not filled:
            raise SystemExit(
                f"stage A block n={n} could not fill {n_per_rung} "
                f"complete samples (skips: {len(skipped)}) -- "
                "INFEASIBLE-INCOMPLETE"
            )
        for r in records:
            r["stage"] = "A"
            r["code_version"] = stamp
        rows.extend(records)
        per_rung_y[n] = np.array([r["y"] for r in records])
        per_rung_y_vol[n] = np.array([r["y_vol"] for r in records])
        skip_counts[n] = len(skipped)
        skipped_seeds[n] = [int(s) for s in skipped]
        print(f"stage A n={n}: {len(records)} complete | "
              f"mean y {per_rung_y[n].mean():.4f} | skips {len(skipped)}",
              flush=True)
    write_rows_csv(output_dir / "p11_stage_a.csv", rows)

    delta = float(per_rung_y[2400].mean() - per_rung_y[600].mean())
    rng = np.random.default_rng(_stable_seed("p11-a-delta"))
    boots = []
    for _ in range(4000):
        top = rng.choice(per_rung_y[2400], size=n_per_rung, replace=True)
        bot = rng.choice(per_rung_y[600], size=n_per_rung, replace=True)
        boots.append(top.mean() - bot.mean())
    lo, hi = (float(np.percentile(boots, 2.5)),
              float(np.percentile(boots, 97.5)))

    log_n = np.log10(np.array(P11_LADDER, dtype=float))
    means = np.array([per_rung_y[n].mean() for n in P11_LADDER])

    def _slope_with_ci(per_rung: dict, label: str):
        point = float(np.polyfit(
            log_n, [per_rung[n].mean() for n in P11_LADDER], 1
        )[0])
        srng = np.random.default_rng(_stable_seed(label))
        boots = []
        for _ in range(4000):
            resampled = [
                srng.choice(per_rung[n], size=n_per_rung,
                            replace=True).mean()
                for n in P11_LADDER
            ]
            boots.append(float(np.polyfit(log_n, resampled, 1)[0]))
        return point, [float(np.percentile(boots, 2.5)),
                       float(np.percentile(boots, 97.5))]

    slope, slope_ci = _slope_with_ci(per_rung_y, "p11-a-slope")
    slope_vol, slope_vol_ci = _slope_with_ci(per_rung_y_vol,
                                             "p11-a-slope-vol")

    summary = {
        "code_version": stamp,
        "n_per_rung": n_per_rung,
        "flat_available": flat_available,
        "skip_counts": {str(n): skip_counts[n] for n in P11_LADDER},
        # the identities, not just the counts (review): which seeds
        # were excluded is what a selection audit needs
        "skipped_seeds": {str(n): skipped_seeds[n] for n in P11_LADDER},
        # any realized skip conditions the estimand on pair packing
        # (prereg Section 4, v1.7) -- carried beside the verdict
        "selection_caveat": bool(any(
            skip_counts[n] > 0 for n in P11_LADDER
        )),
        "delta": delta, "delta_ci": [lo, hi],
        "verdict": verdict(lo, hi, flat_available),
        "mean_y_by_rung": {str(n): float(per_rung_y[n].mean())
                           for n in P11_LADDER},
        "middle_rung_between_endpoints": bool(
            min(means[0], means[2]) <= means[1] <= max(means[0], means[2])
        ),
        "labelled_checks": {
            "slope_mean_y_vs_log10N": slope,
            "slope_ci": slope_ci,
            "predicted_slope_chain": -1.0 / 3.0,
            "slope_mean_y_vol_vs_log10N": slope_vol,
            "slope_vol_ci": slope_vol_ci,
            "predicted_slope_vol": -1.0 / 2.0,
            "constant_level_by_rung": {
                str(n): {
                    "median_relerr_chain": float(np.median(
                        [r["median_relerr_chain"] for r in rows
                         if r["n"] == n]
                    )),
                    "predicted_0.89_m_cond^-1/3": float(0.89 * np.mean(
                        [r["mean_m_conditioned"] for r in rows
                         if r["n"] == n]
                    ) ** (-1.0 / 3.0)),
                } for n in P11_LADDER
            },
        },
    }
    (output_dir / "p11_stage_a_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"delta": delta, "delta_ci": [lo, hi],
                      "verdict": summary["verdict"]}, indent=2))


def eligible_pool_spacelike(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """All ordered spacelike pairs (u_i < u_j, v_i > v_j) whose true
    separation lies in the frozen band (Section 10: same band as
    Stage A)."""

    du = u[None, :] - u[:, None]
    dv = v[None, :] - v[:, None]
    space = (du > 0) & (dv < 0)
    dist = np.where(space, 2.0 * np.sqrt(np.abs(du * dv)), np.nan)
    keep = space & (dist >= TAU_BAND[0]) & (dist <= TAU_BAND[1])
    first, second = np.nonzero(keep)
    return np.column_stack([first, second])


def score_pair_spacelike(u: np.ndarray, v: np.ndarray,
                         i: int, j: int) -> dict:
    """The Section 10 dual-box chain estimator: the unanchored
    interior LIS of the meet-join rank box, no endpoint correction
    (the box contains no endpoint events), through the shared
    definitions with chain_counts_endpoints=False."""

    n = u.size
    u_lo, u_hi = u[i], u[j]
    v_lo, v_hi = v[j], v[i]
    inside = (u > u_lo) & (u < u_hi) & (v > v_lo) & (v < v_hi)
    interior = np.flatnonzero(inside)
    if interior.size:
        su, sv = u[interior], v[interior]
        causal_sub = (su[:, None] < su[None, :]) & (sv[:, None] < sv[None, :])
        chain = longest_chain_length(causal_sub, event_times=su)
    else:
        chain = 0
    rho = n / 2.0
    d_chain = estimate_tau_from_longest_chain_1p1(
        chain, rho=rho, chain_counts_endpoints=False,
    )
    d_vol = estimate_tau_from_interval_count(int(interior.size), rho=rho)
    d_ref = 2.0 * float(np.sqrt((u_hi - u_lo) * (v_hi - v_lo)))
    a, b = round(n * (u_hi - u_lo)), round(n * (v_hi - v_lo))
    m_cond = (a - 1) * (b - 1) / (n - 2)
    return {
        "tau_true": d_ref,
        "tau_chain": float(d_chain),
        "tau_vol": float(d_vol),
        "relerr_chain": abs(float(d_chain) - d_ref) / d_ref,
        "relerr_vol": abs(float(d_vol) - d_ref) / d_ref,
        "m_open": int(interior.size),
        "m_conditioned": float(m_cond),
        "chain_length": int(chain),
    }


def dual_box_order_certificate(u: np.ndarray, v: np.ndarray,
                               i: int, j: int,
                               comp: np.ndarray | None = None) -> bool:
    """Order-only certificate that the causal order FORCES every
    interior member of the dual box (Section 10 amendment; review:
    the forcing regression sampled pairs, production must certify
    every scored pair).

    Anchored percolation to the fixed point: under the two Gamma
    rules every derived pair keeps one endpoint of its parent, so the
    full forcing's x-anchored closure equals iterating
    S <- S ∪ {c : c || x, c comparable to some b in S} from S = {y}
    (and dually T from {x} for the y-anchor). This is EXACT relative
    to Gamma-forcing for the anchored pairs, not a k-step truncation.
    An event certified on both sides is a box member by soundness, so
    certifying all members verifies exactly what the scorer uses."""

    n = u.size
    if comp is None:
        causal = (u[:, None] < u[None, :]) & (v[:, None] < v[None, :])
        comp = causal | causal.T
    members = np.flatnonzero(
        (u > u[i]) & (u < u[j]) & (v > v[j]) & (v < v[i])
    )
    if members.size == 0:
        return True
    comp_u8 = comp.astype(np.uint8)
    idx = np.arange(n)
    incomp_x = ~comp[i] & (idx != i)
    incomp_y = ~comp[j] & (idx != j)
    reach_x = np.zeros(n, dtype=bool)
    reach_x[j] = True
    reach_y = np.zeros(n, dtype=bool)
    reach_y[i] = True
    while True:
        grow_x = reach_x | (incomp_x
                            & (comp_u8 @ reach_x.astype(np.uint8) > 0))
        grow_y = reach_y | (incomp_y
                            & (comp_u8 @ reach_y.astype(np.uint8) > 0))
        if (grow_x == reach_x).all() and (grow_y == reach_y).all():
            break
        reach_x, reach_y = grow_x, grow_y
    return bool(np.all(reach_x[members]) and np.all(reach_y[members]))


def run_sample_spacelike(n: int, seed: int, single_stream: bool = False):
    """One Stage B sample; mirrors run_sample with the spacelike pool
    and the dual-box scorer, same record schema, plus the per-pair
    order certificate."""

    rng = np.random.default_rng(seed)
    pi = rng.permutation(n)
    u, v = continuum_uv(pi)
    pool = eligible_pool_spacelike(u, v)
    draw_rng = rng if single_stream else np.random.default_rng(
        seed + SCORE_OFFSET
    )
    pairs, complete, rejections = draw_disjoint_pairs(pool, u, v, draw_rng)
    record = {
        "n": float(n), "seed": float(seed),
        "pool_size": float(pool.shape[0]),
        "rejections": float(rejections),
        "complete": bool(complete),
    }
    if not complete:
        return record, False
    scored = [score_pair_spacelike(u, v, i, j) for i, j in pairs]
    causal_full = (u[:, None] < u[None, :]) & (v[:, None] < v[None, :])
    comp_full = causal_full | causal_full.T
    record["box_order_certified"] = bool(all(
        dual_box_order_certificate(u, v, i, j, comp=comp_full)
        for i, j in pairs
    ))
    record["y"] = float(np.log10(
        np.median([s["relerr_chain"] for s in scored])
    ))
    record["y_vol"] = float(np.log10(
        np.median([s["relerr_vol"] for s in scored])
    ))
    record["median_relerr_chain"] = float(
        np.median([s["relerr_chain"] for s in scored])
    )
    record["mean_m_conditioned"] = float(
        np.mean([s["m_conditioned"] for s in scored])
    )
    record["min_m_open"] = float(min(s["m_open"] for s in scored))
    return record, True


FROZEN_P11_DIR = (Path(__file__).resolve().parents[2]
                  / "docs" / "prereg" / "frozen" / "p11")


def _require_stage_pass(name: str, stage_label: str) -> dict:
    """Cross-stage gate (Section 6): a later stage runs only after the
    earlier stage's FROZEN record reads IMPROVES. The frozen copy is
    the repo-canonical one, and its stamp must be REACHABLE from HEAD
    (ancestry, not equality — the earlier stage's code legitimately
    predates later commits; a squash-merge that orphaned the stamp
    fails here loudly, review)."""

    import subprocess

    path = FROZEN_P11_DIR / name
    if not path.exists():
        raise SystemExit(
            f"{name} not found in the frozen record -- {stage_label} "
            "has not run and been frozen; the stage order of Section 6 "
            "is frozen."
        )
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("verdict") != "IMPROVES":
        raise SystemExit(
            f"{stage_label} verdict is {artifact.get('verdict')!r}, not "
            "IMPROVES -- the campaign stops at that stage (Section 6)."
        )
    stamp = str(artifact.get("code_version", ""))
    reachable = subprocess.run(
        ["git", "merge-base", "--is-ancestor", stamp, "HEAD"],
        capture_output=True,
    ).returncode == 0
    if not reachable:
        raise SystemExit(
            f"{stage_label}'s frozen stamp {stamp} is not an ancestor "
            "of HEAD -- its producing commit is unreachable from this "
            "history, so its provenance cannot be audited here."
        )
    return artifact


def run_verify_b(output_dir: Path) -> None:
    """Stage B completeness pin and wall times (Section 10)."""

    stamp = _preflight_clean()
    summary: dict = {"code_version": stamp,
                     "pin_required": VERIFY_PIN, "per_n": {}}
    for n in P11_LADDER:
        base = VERIFY_B_BASE[n]
        complete_count = 0
        start = time.perf_counter()
        for k in range(VERIFY_COUNT):
            _record, complete = run_sample_spacelike(
                n, base + k, single_stream=True,
            )
            complete_count += int(complete)
        elapsed = time.perf_counter() - start
        summary["per_n"][str(n)] = {
            "complete": complete_count, "total": VERIFY_COUNT,
            "mean_seconds_per_sample": elapsed / VERIFY_COUNT,
        }
        print(f"verify-b n={n}: {complete_count}/{VERIFY_COUNT} complete | "
              f"{elapsed / VERIFY_COUNT:.3f} s/sample", flush=True)
    summary["pin_passed"] = bool(all(
        summary["per_n"][str(n)]["complete"] >= VERIFY_PIN
        for n in P11_LADDER
    ))
    (output_dir / VERIFY_B_ARTIFACT).write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pin_passed": summary["pin_passed"]}, indent=2))


def run_pilot_b(output_dir: Path) -> None:
    """Stage P-B (Sections 1.3 and 10): both endpoint rungs, 200
    spacelike samples each, calibrated bounds, Delta*_B = Delta*
    (derived in Section 10). Cross-rung statistics forbidden."""

    stamp = _preflight_clean()
    # Section 6: Stage P-B runs only after Stage A PASSES -- checking
    # the verification pin alone let the sequence be bypassed (review)
    _require_stage_pass("p11_stage_a_summary.json", "Stage A")
    verification = _load_gate_artifact(output_dir, VERIFY_B_ARTIFACT, stamp)
    if not verification.get("pin_passed"):
        raise SystemExit("B verification pin failed -- pilot-b refuses")

    summary: dict = {"code_version": stamp, "per_rung": {}}
    ys: dict = {}
    for n, (base, slots) in PILOT_B_BLOCKS.items():
        start = time.perf_counter()
        records, skipped, filled = fill_block(
            n, base, slots, PILOT_SAMPLES, sample_fn=run_sample_spacelike,
        )
        elapsed = time.perf_counter() - start
        if not filled:
            raise SystemExit(
                f"pilot-b block n={n} could not fill {PILOT_SAMPLES} "
                f"complete samples (skips: {len(skipped)}) -- "
                "INFEASIBLE-INCOMPLETE"
            )
        ys[n] = np.array([r["y"] for r in records])
        rung_bound = bonett_variance_bound(ys[n])
        summary["per_rung"][str(n)] = {
            "n_samples": len(records),
            "variance": rung_bound["s2"],
            "kurtosis_g4": rung_bound["g4"],
            "variance_bound_95": rung_bound["bound"],
            "skipped_seeds": skipped,
            "mean_seconds_per_sample": elapsed / len(records),
            "y": [float(val) for val in ys[n]],
        }
        print(f"pilot-b n={n}: var={rung_bound['s2']:.5f} "
              f"g4={rung_bound['g4']:.2f} | skips={len(skipped)}",
              flush=True)

    power = power_requirements(ys[600], ys[2400])
    times = verification["per_n"]
    projected_hours = None
    if power["n_per_rung"] is not None:
        projected_hours = power["n_per_rung"] * sum(
            times[str(n)]["mean_seconds_per_sample"] for n in P11_LADDER
        ) / 3600.0
    summary["power"] = power
    summary["projected_stage_b_hours"] = projected_hours
    summary["feasible"] = bool(
        not power["infeasible"]
        and projected_hours is not None
        and projected_hours <= PROJECTION_LIMIT_HOURS
    )
    (output_dir / PILOT_B_ARTIFACT).write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"power": {k: power[k] for k in
                                ("s2", "s2_90", "n_sup", "n_eq",
                                 "n_per_rung", "flat_available",
                                 "infeasible")},
                      "projected_stage_b_hours": projected_hours,
                      "feasible": summary["feasible"]}, indent=2))


def run_stage_b(output_dir: Path) -> None:
    """Stage B (Sections 1.4, 6, 10): the spacelike gate at pilot-b's
    n, frozen verdict table, labelled checks riding along."""

    stamp = _preflight_clean()
    pilot = _load_gate_artifact(output_dir, PILOT_B_ARTIFACT, stamp)
    if not pilot.get("feasible"):
        raise SystemExit("pilot-b declared the design infeasible -- "
                         "Stage B refuses to run")
    n_per_rung = int(pilot["power"]["n_per_rung"])
    flat_available = bool(pilot["power"]["flat_available"])

    rows, per_rung_y, per_rung_y_vol = [], {}, {}
    skip_counts, skipped_seeds = {}, {}
    for n, (base, slots) in STAGE_B_BLOCKS.items():
        records, skipped, filled = fill_block(
            n, base, slots, n_per_rung, sample_fn=run_sample_spacelike,
        )
        if not filled:
            raise SystemExit(
                f"stage B block n={n} could not fill {n_per_rung} "
                f"complete samples (skips: {len(skipped)}) -- "
                "INFEASIBLE-INCOMPLETE"
            )
        for r in records:
            r["stage"] = "B"
            r["code_version"] = stamp
        rows.extend(records)
        per_rung_y[n] = np.array([r["y"] for r in records])
        per_rung_y_vol[n] = np.array([r["y_vol"] for r in records])
        skip_counts[n] = len(skipped)
        skipped_seeds[n] = [int(s) for s in skipped]
        print(f"stage B n={n}: {len(records)} complete | "
              f"mean y {per_rung_y[n].mean():.4f} | skips {len(skipped)}",
              flush=True)
    write_rows_csv(output_dir / "p11_stage_b.csv", rows)

    delta = float(per_rung_y[2400].mean() - per_rung_y[600].mean())
    rng = np.random.default_rng(_stable_seed("p11-b-delta"))
    boots = []
    for _ in range(4000):
        top = rng.choice(per_rung_y[2400], size=n_per_rung, replace=True)
        bot = rng.choice(per_rung_y[600], size=n_per_rung, replace=True)
        boots.append(top.mean() - bot.mean())
    lo, hi = (float(np.percentile(boots, 2.5)),
              float(np.percentile(boots, 97.5)))

    log_n = np.log10(np.array(P11_LADDER, dtype=float))
    means = np.array([per_rung_y[n].mean() for n in P11_LADDER])

    def _slope_with_ci(per_rung: dict, label: str):
        point = float(np.polyfit(
            log_n, [per_rung[n].mean() for n in P11_LADDER], 1
        )[0])
        srng = np.random.default_rng(_stable_seed(label))
        sboots = []
        for _ in range(4000):
            resampled = [
                srng.choice(per_rung[n], size=n_per_rung,
                            replace=True).mean()
                for n in P11_LADDER
            ]
            sboots.append(float(np.polyfit(log_n, resampled, 1)[0]))
        return point, [float(np.percentile(sboots, 2.5)),
                       float(np.percentile(sboots, 97.5))]

    slope, slope_ci = _slope_with_ci(per_rung_y, "p11-b-slope")
    slope_vol, slope_vol_ci = _slope_with_ci(per_rung_y_vol,
                                             "p11-b-slope-vol")

    summary = {
        "code_version": stamp,
        "n_per_rung": n_per_rung,
        "flat_available": flat_available,
        "skip_counts": {str(n): skip_counts[n] for n in P11_LADDER},
        "skipped_seeds": {str(n): skipped_seeds[n] for n in P11_LADDER},
        "selection_caveat": bool(any(
            skip_counts[n] > 0 for n in P11_LADDER
        )),
        # every scored pair carries the order-only membership
        # certificate (Section 10 amendment); a count above zero
        # scopes the verdict as coordinate-assisted for those samples
        "box_uncertified_samples": int(sum(
            1 for r in rows if not r.get("box_order_certified", False)
        )),
        "delta": delta, "delta_ci": [lo, hi],
        "verdict": verdict(lo, hi, flat_available),
        "mean_y_by_rung": {str(n): float(per_rung_y[n].mean())
                           for n in P11_LADDER},
        "middle_rung_between_endpoints": bool(
            min(means[0], means[2]) <= means[1] <= max(means[0], means[2])
        ),
        "labelled_checks": {
            "slope_mean_y_vs_log10N": slope,
            "slope_ci": slope_ci,
            "predicted_slope_chain": -1.0 / 3.0,
            "slope_mean_y_vol_vs_log10N": slope_vol,
            "slope_vol_ci": slope_vol_ci,
            "predicted_slope_vol": -1.0 / 2.0,
            "constant_level_by_rung": {
                str(n): {
                    "median_relerr_chain": float(np.median(
                        [r["median_relerr_chain"] for r in rows
                         if r["n"] == n]
                    )),
                    "predicted_0.89_m_cond^-1/3": float(0.89 * np.mean(
                        [r["mean_m_conditioned"] for r in rows
                         if r["n"] == n]
                    ) ** (-1.0 / 3.0)),
                } for n in P11_LADDER
            },
        },
    }
    (output_dir / "p11_stage_b_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"delta": delta, "delta_ci": [lo, hi],
                      "verdict": summary["verdict"]}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage",
                        choices=["verify", "pilot", "a",
                                 "verify-b", "pilot-b", "b", "c"],
                        default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.stage == "verify":
        run_verify(args.output_dir)
    elif args.stage == "pilot":
        run_pilot(args.output_dir)
    elif args.stage == "a":
        run_stage_a(args.output_dir)
    elif args.stage == "verify-b":
        run_verify_b(args.output_dir)
    elif args.stage == "pilot-b":
        run_pilot_b(args.output_dir)
    elif args.stage == "b":
        run_stage_b(args.output_dir)
    elif args.stage == "c":
        raise SystemExit(
            "stage c is gated on its frozen addendum (prereg Section 6) "
            "and has no runner until that addendum lands -- this "
            "refusal is the preregistration operating."
        )
    else:
        raise SystemExit("choose --stage verify/pilot/a/verify-b/pilot-b/b")


if __name__ == "__main__":
    main()
