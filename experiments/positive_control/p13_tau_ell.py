"""P13: how far does curvature let the flat-normalized reading go?

Implements docs/prereg/p13_tau_ell_ladder.md v1.1 and nothing else.
P12 swept density at fixed tau/ell; P13 sweeps tau/ell at FIXED
discreteness, which is why every rung carries its own patch and its
own intensity and why the m-matching is gated rather than assumed.

Everything the estimator does is P11's and P12's, imported: the
sprinkler's thinning, the exact dS_2 truth, the eligibility pair
(band AND box-inside-patch), the greedy disjoint packing, the
calibrated Bonett power machinery, and the shared chain estimator.
Only the per-rung parameterization and the inverted verdict table are
new.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from p10_continuum_ladder import _stable_seed
from p11_metric import (
    DRAW_REJECTION_CAP,
    K_PAIRS,
    PILOT_SAMPLES,
    PROJECTION_LIMIT_HOURS,
    SKIP_CAP,
    STRIDE,
    VERIFY_COUNT,
    VERIFY_PIN,
    _de_nan,
    _load_gate_artifact,
    _preflight_clean,
    _require_stage_pass,
    bonett_variance_bound,
    calibrated_variance_bound,
    score_pair,
    supports_disjoint,
)
from pc_common import DEFAULT_OUTPUT_DIR, write_rows_csv

from causal_spacetime_lab.estimators import (
    estimate_tau_from_longest_chain_1p1,
)

ELL = 1.0

#: Frozen per-rung constants (prereg Section 2): (eta_lo, X, rho).
#: The X values grow steeply because a pair at tau/ell = 1.5 spans an
#: eta interval comparable to the patch, so six disjoint boxes can
#: only stack along x; the rho values are CALIBRATED so the realized
#: interval count matches at every rung, since at fixed rho*tau^2 the
#: realized m drifts (that drift is the curvature correction to the
#: diamond volume).
RUNGS = (0.30, 0.60, 1.00, 1.50)
PATCH = {
    0.30: (-1.8, 1.8, 1707.0),
    0.60: (-2.6, 4.0, 434.8),
    1.00: (-4.0, 8.0, 163.0),
    1.50: (-7.0, 17.0, 76.1),
}
BAND_REL = 0.10
M_TARGET = 76.0
M_TOLERANCE = 0.05          # Section 5 gate: +/- 5% of the grand mean

#: Power design (Section 1). The cap is 200 here, not P11's 60: the
#: design-check variance puts n_eq around 170, and a design that
#: cannot buy its own expected verdict is not a design.
DELTA_DETECT = 0.15
DELTA_EQ = 0.05
N_FLOOR = 12
N_CAP = 200
N_SUP_COEFF = (1.960 + 1.282) ** 2 / DELTA_DETECT ** 2
N_EQ_COEFF = (1.960 + 1.645) ** 2 / DELTA_EQ ** 2

#: Seed windows (Section 6). Design-check space is 7000000+ and is
#: never touched by a runner.
VERIFY_BASE = {0.30: 6000000, 0.60: 6002000, 1.00: 6004000,
               1.50: 6006000}
PILOT_BLOCKS = {0.30: (6020000, 220), 1.50: (6064000, 220)}
STAGE_A_BLOCKS = {0.30: (6120000, 220), 0.60: (6164000, 220),
                  1.00: (6208000, 220), 1.50: (6252000, 220)}
TWIN_BLOCKS = {0.30: (6300000, 220), 0.60: (6344000, 220),
               1.00: (6388000, 220), 1.50: (6432000, 220)}
SCORE_OFFSET = 150

#: P13's prerequisite record lives in P12's frozen directory.
FROZEN_P12_DIR = (Path(__file__).resolve().parents[2]
                  / "docs" / "prereg" / "frozen" / "p12")

VERIFY_ARTIFACT = "p13_verification_summary.json"
PILOT_ARTIFACT = "p13_pilot_summary.json"


def patch_proper_volume(eta_lo: float, xhalf: float) -> float:
    """ell^2 (1 - 1/|eta_lo|) * 2X -- the closed form of the frozen
    quadrature, and the denominator of the Section 3 rho convention."""

    return ELL ** 2 * (1.0 - 1.0 / abs(eta_lo)) * 2.0 * xhalf


def tau_curved(eta_p, x_p, eta_q, x_q):
    """The exact dS_2 geodesic proper time (P12 Section 3), pinned
    there against an independent geodesic integration at 1e-9."""

    z = ((np.asarray(eta_p) ** 2 + np.asarray(eta_q) ** 2
          - (np.asarray(x_p) - np.asarray(x_q)) ** 2)
         / (2.0 * np.asarray(eta_p) * np.asarray(eta_q)))
    return ELL * np.arccosh(np.where(z >= 1.0, z, np.nan))


def sprinkle(eta_lo: float, xhalf: float, rho: float,
             rng: np.random.Generator, flat: bool = False):
    """Inhomogeneous Poisson by thinning (P12 Section 2). ``flat``
    makes the intensity uniform -- the twin control of Section 5."""

    if flat:
        n_prop = int(rng.poisson(
            rho * (abs(eta_lo) - 1.0) * 2.0 * xhalf
        ))
        return (rng.uniform(eta_lo, -1.0, n_prop),
                rng.uniform(-xhalf, xhalf, n_prop))
    omega2_max = 1.0
    n_prop = int(rng.poisson(
        rho * (abs(eta_lo) - 1.0) * 2.0 * xhalf * omega2_max
    ))
    eta = rng.uniform(eta_lo, -1.0, n_prop)
    x = rng.uniform(-xhalf, xhalf, n_prop)
    keep = rng.random(n_prop) < (ELL ** 2 / eta ** 2) / omega2_max
    return eta[keep], x[keep]


def eligible_pairs(eta, x, eta_lo, xhalf, band, flat: bool = False):
    """Both frozen conditions: separation in band, causal box inside
    the patch. Coordinates only, decided before any measurement.

    For the curved rungs the band is on tau_curved; for the flat twin
    it is on tau_flat = sqrt(dU dV), centred on the curved rung's
    design-check mean box area so the packing pressure matches.
    """

    u, v = eta + x, eta - x
    du = u[None, :] - u[:, None]
    dv = v[None, :] - v[:, None]
    related = (du > 0) & (dv > 0)
    if flat:
        sep = np.sqrt(np.where(related, du * dv, np.nan))
    else:
        sep = tau_curved(eta[:, None], x[:, None],
                         eta[None, :], x[None, :])
    keep = related & (sep >= band[0]) & (sep <= band[1])
    a, b = np.nonzero(keep)
    if a.size == 0:
        return np.empty((0, 2), dtype=int)
    e1, x1 = (u[b] + v[a]) / 2.0, (u[b] - v[a]) / 2.0
    e2, x2 = (u[a] + v[b]) / 2.0, (u[a] - v[b]) / 2.0
    inside = (
        (e1 >= eta_lo) & (e1 <= -1.0) & (np.abs(x1) <= xhalf)
        & (e2 >= eta_lo) & (e2 <= -1.0) & (np.abs(x2) <= xhalf)
    )
    return np.column_stack([a[inside], b[inside]])


def draw_disjoint(pool, u, v, rng):
    accepted, boxes, rejections = [], [], 0
    if pool.shape[0] == 0:
        return accepted, False, 0
    for idx in rng.permutation(pool.shape[0]):
        i, j = int(pool[idx, 0]), int(pool[idx, 1])
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


def run_sample(tau_c: float, seed: int, flat: bool = False,
               single_stream: bool = False):
    """One P13 sample at a rung. ``flat`` runs the twin control."""

    eta_lo, xhalf, rho = PATCH[tau_c]
    band = (tau_c * (1.0 - BAND_REL), tau_c * (1.0 + BAND_REL))
    if flat:
        centre = TWIN_BAND_CENTRE[tau_c]
        band = (centre * (1.0 - BAND_REL), centre * (1.0 + BAND_REL))
        rho = TWIN_RHO[tau_c]
    rng = np.random.default_rng(seed)
    eta, x = sprinkle(eta_lo, xhalf, rho, rng, flat=flat)
    # Section 3 convention: realized count over the patch's PROPER
    # volume (not the sprinkler's nominal intensity). For the flat
    # twin that proper volume IS the coordinate area -- using the
    # curved patch's volume there mis-normalizes tau_hat, which the
    # first implementation did (twin relative error read 39-70%
    # instead of ~18%).
    volume = (patch_proper_volume(eta_lo, xhalf) if not flat
              else (abs(eta_lo) - 1.0) * 2.0 * xhalf)
    rho_hat = eta.size / volume
    pool = eligible_pairs(eta, x, eta_lo, xhalf, band, flat=flat)
    draw_rng = rng if single_stream else np.random.default_rng(
        seed + SCORE_OFFSET
    )
    u, v = eta + x, eta - x
    pairs, complete, rejections = draw_disjoint(pool, u, v, draw_rng)
    record = {
        "tau_ell": float(tau_c), "seed": float(seed),
        "flat_twin": bool(flat),
        "n_realized": float(eta.size), "rho_hat": float(rho_hat),
        "pool_size": float(pool.shape[0]),
        "rejections": float(rejections), "complete": bool(complete),
    }
    if not complete:
        return record, False
    errors, counts = [], []
    for i, j in pairs:
        u_half, v_half = u / 2.0, v / 2.0
        scored = score_pair(u_half, v_half, i, j)
        tau_hat = float(estimate_tau_from_longest_chain_1p1(
            scored["chain_length"], rho=rho_hat
        ))
        if flat:
            truth = float(np.sqrt((u[j] - u[i]) * (v[j] - v[i])))
        else:
            truth = float(tau_curved(eta[i], x[i], eta[j], x[j]))
        errors.append(abs(tau_hat - truth) / truth)
        counts.append(int(scored["m_open"]))
    record["y"] = float(np.log10(np.median(errors)))
    record["median_relerr"] = float(np.median(errors))
    record["mean_m"] = float(np.mean(counts))
    record["sd_m"] = float(np.std(counts, ddof=1))
    return record, True


#: The flat twin's two calibrated constants per rung (prereg Section
#: 5): the band CENTRE, measured as the curved rung's realized mean
#: box side sqrt(dU dV) in design-check space, so the box-to-patch
#: ratio -- the packing pressure -- matches; and the uniform intensity
#: rho_twin = 2 * M_TARGET / centre^2, which puts the twin's realized
#: mean m on the same target. Both are MEASURED with this runner
#: rather than guessed: a first implementation guessed a scale factor
#: and the twin's m came out 172 and 872 against a target of 76.
TWIN_BAND_CENTRE = {0.30: 0.4016, 0.60: 0.9725,
                    1.00: 2.0515, 1.50: 4.2262}
TWIN_RHO = {tau_c: 2.0 * M_TARGET / centre ** 2
            for tau_c, centre in TWIN_BAND_CENTRE.items()}


def verdict(lo: float, hi: float, robust_available: bool) -> str:
    """The frozen table of Section 1.2, by precedence. The polarity is
    inverted from P11/P12: here the null is no effect and the
    interesting outcome is POSITIVE."""

    if lo > 0.0:
        return "CURVATURE-DEGRADES"
    if hi < 0.0:
        return "CURVATURE-HELPS"
    if robust_available and (-DELTA_EQ < lo) and (hi < DELTA_EQ):
        return "CURVATURE-ROBUST"
    return "INCONCLUSIVE" if robust_available else "UNRESOLVED"


def power_requirements(y_bottom: np.ndarray, y_top: np.ndarray) -> dict:
    """Section 1.3: calibrated Bonett bounds summed, then the two
    requirements, then the frozen selection rule with cap 200."""

    bottom = calibrated_variance_bound(y_bottom, "bottom")
    top = calibrated_variance_bound(y_top, "top")
    s2 = bottom["s2"] + top["s2"]
    s2_90 = bottom["bound"] + top["bound"]
    n_sup = int(np.ceil(N_SUP_COEFF * s2_90))
    n_eq = int(np.ceil(N_EQ_COEFF * s2_90))
    calibration = {
        side: {k: data[k] for k in
               ("coverage_at_nominal", "z_used", "calibrated")}
        for side, data in (("bottom", bottom), ("top", top))
    }

    def clamp(k):
        return int(min(max(k, N_FLOOR), N_CAP))

    if n_sup > N_CAP:
        return {"s2": s2, "s2_90": s2_90, "calibration": calibration,
                "n_sup": n_sup, "n_eq": n_eq, "n_per_rung": None,
                "robust_available": False, "infeasible": True}
    if n_eq <= N_CAP:
        return {"s2": s2, "s2_90": s2_90, "calibration": calibration,
                "n_sup": n_sup, "n_eq": n_eq,
                "n_per_rung": clamp(max(n_sup, n_eq)),
                "robust_available": True, "infeasible": False}
    return {"s2": s2, "s2_90": s2_90, "calibration": calibration,
            "n_sup": n_sup, "n_eq": n_eq, "n_per_rung": clamp(n_sup),
            "robust_available": False, "infeasible": False}


def _fill_block(tau_c, base, slots, needed, flat=False):
    records, skipped = [], []
    for k in range(slots):
        seed = base + STRIDE * k
        record, complete = run_sample(tau_c, seed, flat=flat)
        if complete:
            records.append(record)
            if len(records) == needed:
                return records, skipped, True
        else:
            skipped.append(seed)
            if len(skipped) > SKIP_CAP:
                return records, skipped, False
    return records, skipped, False


def run_verify(output_dir: Path) -> None:
    stamp = _preflight_clean()
    summary: dict = {"code_version": stamp, "pin_required": VERIFY_PIN,
                     "per_rung": {}}
    for tau_c in RUNGS:
        base = VERIFY_BASE[tau_c]
        done, start = 0, time.perf_counter()
        for k in range(VERIFY_COUNT):
            _r, ok = run_sample(tau_c, base + k, single_stream=True)
            done += int(ok)
        elapsed = time.perf_counter() - start
        summary["per_rung"][str(tau_c)] = {
            "complete": done, "total": VERIFY_COUNT,
            "mean_seconds_per_sample": elapsed / VERIFY_COUNT,
        }
        print(f"verify-13 tau/ell={tau_c}: {done}/{VERIFY_COUNT} | "
              f"{elapsed / VERIFY_COUNT:.3f} s/sample", flush=True)
    summary["pin_passed"] = bool(all(
        summary["per_rung"][str(t)]["complete"] >= VERIFY_PIN
        for t in RUNGS
    ))
    (output_dir / VERIFY_ARTIFACT).write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pin_passed": summary["pin_passed"]}, indent=2))


def run_pilot(output_dir: Path) -> None:
    stamp = _preflight_clean()
    _require_stage_pass("p12_stage_a_summary.json", "P12 Stage A",
                        directory=FROZEN_P12_DIR)
    verification = _load_gate_artifact(output_dir, VERIFY_ARTIFACT, stamp)
    if not verification.get("pin_passed"):
        raise SystemExit("verification-13 pin failed -- pilot refuses")

    summary: dict = {"code_version": stamp, "per_rung": {}}
    ys: dict = {}
    for tau_c, (base, slots) in PILOT_BLOCKS.items():
        start = time.perf_counter()
        records, skipped, filled = _fill_block(
            tau_c, base, slots, PILOT_SAMPLES)
        elapsed = time.perf_counter() - start
        if not filled:
            raise SystemExit(
                f"pilot-13 rung {tau_c} could not fill {PILOT_SAMPLES} "
                f"(skips: {len(skipped)}) -- INFEASIBLE-INCOMPLETE")
        ys[tau_c] = np.array([r["y"] for r in records])
        nominal = bonett_variance_bound(ys[tau_c])
        calibrated = calibrated_variance_bound(
            ys[tau_c], "bottom" if tau_c == RUNGS[0] else "top")
        summary["per_rung"][str(tau_c)] = {
            "n_samples": len(records),
            "variance": nominal["s2"],
            "kurtosis_g4": nominal["g4"],
            "variance_bound_95_nominal": nominal["bound"],
            "variance_bound_95_calibrated": calibrated["bound"],
            "calibration_z_used": calibrated["z_used"],
            "calibration_coverage_at_nominal":
                calibrated["coverage_at_nominal"],
            "mean_m": float(np.mean([r["mean_m"] for r in records])),
            "skipped_seeds": skipped,
            "mean_seconds_per_sample": elapsed / len(records),
            "y": [float(val) for val in ys[tau_c]],
        }
        print(f"pilot-13 tau/ell={tau_c}: var={nominal['s2']:.5f} "
              f"g4={nominal['g4']:.2f} | skips={len(skipped)}", flush=True)

    power = power_requirements(ys[RUNGS[0]], ys[RUNGS[-1]])
    times = verification["per_rung"]
    projected = None
    if power["n_per_rung"] is not None:
        # four curved rungs plus four twin rungs
        projected = 2.0 * power["n_per_rung"] * sum(
            times[str(t)]["mean_seconds_per_sample"] for t in RUNGS
        ) / 3600.0
    summary["power"] = power
    summary["projected_stage_a_hours"] = projected
    summary["feasible"] = bool(
        not power["infeasible"] and projected is not None
        and projected <= PROJECTION_LIMIT_HOURS
    )
    (output_dir / PILOT_ARTIFACT).write_text(
        json.dumps(_de_nan(summary), indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"power": power, "feasible": summary["feasible"],
                      "projected_stage_a_hours": projected}, indent=2))


def run_stage_a(output_dir: Path) -> None:
    stamp = _preflight_clean()
    pilot = _load_gate_artifact(output_dir, PILOT_ARTIFACT, stamp)
    if not pilot.get("feasible"):
        raise SystemExit("pilot-13 declared the design infeasible")
    n_per_rung = int(pilot["power"]["n_per_rung"])
    robust_available = bool(pilot["power"]["robust_available"])

    rows, per_rung, twin_rows, twin_y = [], {}, [], {}
    skip_counts, skipped_seeds = {}, {}
    for tau_c in RUNGS:
        base, slots = STAGE_A_BLOCKS[tau_c]
        records, skipped, filled = _fill_block(
            tau_c, base, slots, n_per_rung)
        if not filled:
            raise SystemExit(
                f"stage A-13 rung {tau_c} could not fill {n_per_rung} "
                f"(skips: {len(skipped)}) -- INFEASIBLE-INCOMPLETE")
        for r in records:
            r["stage"] = "A13"
            r["code_version"] = stamp
        rows.extend(records)
        per_rung[tau_c] = np.array([r["y"] for r in records])
        skip_counts[tau_c] = len(skipped)
        skipped_seeds[tau_c] = [int(s) for s in skipped]
        print(f"stage A-13 tau/ell={tau_c}: {len(records)} complete | "
              f"mean y {per_rung[tau_c].mean():.4f} | "
              f"skips {len(skipped)}", flush=True)

        tbase, tslots = TWIN_BLOCKS[tau_c]
        trecords, tskipped, tfilled = _fill_block(
            tau_c, tbase, tslots, n_per_rung, flat=True)
        if not tfilled:
            raise SystemExit(
                f"flat twin rung {tau_c} could not fill {n_per_rung} "
                f"(skips: {len(tskipped)}) -- INFEASIBLE-INCOMPLETE")
        for r in trecords:
            r["stage"] = "A13-twin"
            r["code_version"] = stamp
        twin_rows.extend(trecords)
        twin_y[tau_c] = np.array([r["y"] for r in trecords])
        print(f"    twin tau/ell={tau_c}: {len(trecords)} complete | "
              f"mean y {twin_y[tau_c].mean():.4f}", flush=True)
    write_rows_csv(output_dir / "p13_stage_a.csv", rows + twin_rows)

    def contrast(per, label):
        delta = float(per[RUNGS[-1]].mean() - per[RUNGS[0]].mean())
        rng = np.random.default_rng(_stable_seed(label))
        boots = [
            float(rng.choice(per[RUNGS[-1]], size=n_per_rung,
                             replace=True).mean()
                  - rng.choice(per[RUNGS[0]], size=n_per_rung,
                               replace=True).mean())
            for _ in range(4000)
        ]
        return delta, [float(np.percentile(boots, 2.5)),
                       float(np.percentile(boots, 97.5))]

    delta, (lo, hi) = contrast(per_rung, "p13-delta")
    twin_delta, twin_ci = contrast(twin_y, "p13-twin-delta")

    # Section 5 gate 1: m matching
    mean_m = {t: float(np.mean([r["mean_m"] for r in rows
                                if r["tau_ell"] == t])) for t in RUNGS}
    grand = float(np.mean(list(mean_m.values())))
    m_ok = all(abs(mean_m[t] - grand) / grand <= M_TOLERANCE
               for t in RUNGS)
    # Section 5 gate 2: the twin's contrast must be within delta_eq
    twin_ok = abs(twin_delta) <= DELTA_EQ

    # labelled: the (tau/ell)^2 trend across all four rungs
    x2 = np.array(RUNGS, dtype=float) ** 2
    means = np.array([per_rung[t].mean() for t in RUNGS])
    trend = float(np.polyfit(x2, means, 1)[0])
    trng = np.random.default_rng(_stable_seed("p13-trend"))
    trend_boots = [
        float(np.polyfit(x2, [
            trng.choice(per_rung[t], size=n_per_rung, replace=True).mean()
            for t in RUNGS
        ], 1)[0])
        for _ in range(4000)
    ]

    summary = {
        "code_version": stamp,
        "n_per_rung": n_per_rung,
        "robust_available": robust_available,
        "delta": delta, "delta_ci": [lo, hi],
        "delta_eq": DELTA_EQ,
        "verdict": (verdict(lo, hi, robust_available)
                    if (m_ok and twin_ok) else "CONFOUNDED"),
        "confound_gates": {
            "m_matched": m_ok, "mean_m_by_rung": mean_m,
            "m_grand_mean": grand, "m_tolerance": M_TOLERANCE,
            "flat_twin_delta": twin_delta, "flat_twin_ci": twin_ci,
            "flat_twin_within_delta_eq": twin_ok,
        },
        "mean_y_by_rung": {str(t): float(per_rung[t].mean())
                           for t in RUNGS},
        "skip_counts": {str(t): skip_counts[t] for t in RUNGS},
        "skipped_seeds": {str(t): skipped_seeds[t] for t in RUNGS},
        "selection_caveat": bool(any(skip_counts[t] > 0 for t in RUNGS)),
        "labelled_checks": {
            "trend_vs_tau_ell_squared": trend,
            "trend_ci": [float(np.percentile(trend_boots, 2.5)),
                         float(np.percentile(trend_boots, 97.5))],
            "median_relerr_by_rung": {
                str(t): float(np.median(
                    [r["median_relerr"] for r in rows
                     if r["tau_ell"] == t]
                )) for t in RUNGS
            },
            "sd_m_by_rung": {
                str(t): float(np.mean([r["sd_m"] for r in rows
                                       if r["tau_ell"] == t]))
                for t in RUNGS
            },
        },
    }
    (output_dir / "p13_stage_a_summary.json").write_text(
        json.dumps(_de_nan(summary), indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"delta": delta, "delta_ci": [lo, hi],
                      "verdict": summary["verdict"],
                      "m_matched": m_ok,
                      "flat_twin_delta": twin_delta}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["verify", "pilot", "a"],
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
    else:
        raise SystemExit("choose --stage verify/pilot/a")


if __name__ == "__main__":
    main()
