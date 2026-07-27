"""P12: does the instrument read CURVED geometry? — stage runners.

Implements docs/prereg/p12_curved_sprinkling.md v1.1a and nothing
else. The estimator is P11's, unchanged and imported; what changes is
the ambient geometry, which in 1+1D enters only through the sprinkling
density (the conformal factor cannot alter causal structure).

Convention, pinned by the prereg and asserted by a test: this module
works in FULL-null coordinates U = eta + x, V = eta - x, where the
flat proper time is sqrt(dU dV); every call into a P11 routine passes
U/2, V/2, the half-null convention in which the same number reads
2 sqrt(du dv).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from p10_continuum_ladder import _stable_seed
from p11_metric import (
    BONETT_Z,
    DRAW_REJECTION_CAP,
    K_PAIRS,
    N_CAP,
    N_EQ_COEFF,
    N_FLOOR,
    N_SUP_COEFF,
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
    power_requirements,
    score_pair,
    supports_disjoint,
    verdict,
)
from pc_common import DEFAULT_OUTPUT_DIR, write_rows_csv

#: Frozen patch constants (prereg Section 2, v1.1).
ELL = 1.0
ETA_RANGE = (-2.0, -1.0)
X_HALF = 1.5
RICCI_SCALAR = 2.0 / ELL ** 2

#: Ladder, band, protocol (Sections 1.1 and 4).
P12_LADDER = (600, 1200, 2400)
TAU_BAND = (0.28, 0.34)

#: Seed windows (Section 6). Design-check space is 2000000-5999999 and
#: is never touched by a runner.
VERIFY_BASE = {600: 1000000, 1200: 1002000, 2400: 1004000}
PILOT_BLOCKS = {600: (1020000, 220), 2400: (1064000, 220)}
STAGE_A_BLOCKS = {600: (1120000, 80), 1200: (1136000, 80),
                  2400: (1152000, 80)}
SCORE_OFFSET = 150

VERIFY_ARTIFACT = "p12_verification_summary.json"
PILOT_ARTIFACT = "p12_pilot_summary.json"


def patch_proper_volume() -> float:
    """Exact proper 2-volume of the frozen patch, by quadrature:
    integral of Omega^2 = ell^2/eta^2 over the coordinate rectangle."""

    eta = np.linspace(ETA_RANGE[0], ETA_RANGE[1], 200001)
    return float(np.trapezoid(ELL ** 2 / eta ** 2, eta) * (2.0 * X_HALF))


PATCH_VOLUME = patch_proper_volume()


def tau_curved(eta_p, x_p, eta_q, x_q):
    """Exact dS_2 geodesic proper time via the de Sitter invariant
    (prereg Section 3). Vectorized; NaN where not timelike."""

    z = ((np.asarray(eta_p) ** 2 + np.asarray(eta_q) ** 2
          - (np.asarray(x_p) - np.asarray(x_q)) ** 2)
         / (2.0 * np.asarray(eta_p) * np.asarray(eta_q)))
    return ELL * np.arccosh(np.where(z >= 1.0, z, np.nan))


def sprinkle(n_expected: int, rng: np.random.Generator):
    """A genuine inhomogeneous Poisson sprinkling of intensity
    proportional to Omega^2 = ell^2/eta^2, by thinning a uniform
    proposal (prereg Section 2). Returns (eta, x)."""

    omega2_max = ELL ** 2 / ETA_RANGE[1] ** 2
    accept = (PATCH_VOLUME
              / ((ETA_RANGE[1] - ETA_RANGE[0]) * 2.0 * X_HALF * omega2_max))
    n_prop = int(rng.poisson(n_expected / accept))
    eta = rng.uniform(ETA_RANGE[0], ETA_RANGE[1], n_prop)
    x = rng.uniform(-X_HALF, X_HALF, n_prop)
    keep = rng.random(n_prop) < (ELL ** 2 / eta ** 2) / omega2_max
    return eta[keep], x[keep]


def eligible_pairs(eta: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Both frozen eligibility conditions (Section 4, v1.1): curved
    proper time in the band, AND the causal box inside the patch.
    Both are conditions on coordinates, decided before any
    measurement, so completeness never conditions on a measured
    value."""

    u, v = eta + x, eta - x
    du = u[None, :] - u[:, None]
    dv = v[None, :] - v[:, None]
    related = (du > 0) & (dv > 0)
    tau = tau_curved(eta[:, None], x[:, None], eta[None, :], x[None, :])
    keep = related & (tau >= TAU_BAND[0]) & (tau <= TAU_BAND[1])
    a, b = np.nonzero(keep)
    if a.size == 0:
        return np.empty((0, 2), dtype=int)
    # the box's two null corners, in (eta, x)
    e1, x1 = (u[b] + v[a]) / 2.0, (u[b] - v[a]) / 2.0
    e2, x2 = (u[a] + v[b]) / 2.0, (u[a] - v[b]) / 2.0
    inside = (
        (e1 >= ETA_RANGE[0]) & (e1 <= ETA_RANGE[1]) & (np.abs(x1) <= X_HALF)
        & (e2 >= ETA_RANGE[0]) & (e2 <= ETA_RANGE[1])
        & (np.abs(x2) <= X_HALF)
    )
    return np.column_stack([a[inside], b[inside]])


def draw_disjoint(pool: np.ndarray, u: np.ndarray, v: np.ndarray,
                  rng: np.random.Generator):
    """P11's greedy disjoint-support rule, supports being the pairs'
    (U, V) boxes; the cap counts overlap rejections only."""

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


def score_curved_pair(eta, x, i, j, rho):
    """The P11 Stage A estimator, unchanged, against curved truth.

    The half-null mapping u = U/2, v = V/2 is what makes score_pair's
    convention agree with this module's; score_pair normalizes by its
    own n, so it is called on a scaled copy and the chain length is
    renormalized here with the frozen rho.
    """

    u_half, v_half = (eta + x) / 2.0, (eta - x) / 2.0
    scored = score_pair(u_half, v_half, i, j)
    chain = scored["chain_length"]
    tau_hat = max(chain - 2, 0) / np.sqrt(2.0 * rho)
    tau_true = float(tau_curved(eta[i], x[i], eta[j], x[j]))
    tau_flat = float(np.sqrt(
        ((eta[j] + x[j]) - (eta[i] + x[i]))
        * ((eta[j] - x[j]) - (eta[i] - x[i]))
    ))
    return {
        "tau_hat": float(tau_hat),
        "tau_curved": tau_true,
        "tau_flat": tau_flat,
        "relerr_curved": abs(tau_hat - tau_true) / tau_true,
        "relerr_flat_arm": abs(tau_flat - tau_true) / tau_true,
        "m_open": int(scored["m_open"]),
        "chain_length": int(chain),
    }


def run_sample(n_expected: int, seed: int, single_stream: bool = False):
    rng = np.random.default_rng(seed)
    eta, x = sprinkle(n_expected, rng)
    rho = eta.size / PATCH_VOLUME
    pool = eligible_pairs(eta, x)
    draw_rng = rng if single_stream else np.random.default_rng(
        seed + SCORE_OFFSET
    )
    u, v = eta + x, eta - x
    pairs, complete, rejections = draw_disjoint(pool, u, v, draw_rng)
    record = {
        "n_expected": float(n_expected), "seed": float(seed),
        "n_realized": float(eta.size), "rho": float(rho),
        "pool_size": float(pool.shape[0]),
        "rejections": float(rejections), "complete": bool(complete),
    }
    if not complete:
        return record, False
    scored = [score_curved_pair(eta, x, i, j, rho) for i, j in pairs]
    record["y"] = float(np.log10(
        np.median([s["relerr_curved"] for s in scored])
    ))
    record["median_relerr_curved"] = float(
        np.median([s["relerr_curved"] for s in scored])
    )
    record["median_relerr_flat_arm"] = float(
        np.median([s["relerr_flat_arm"] for s in scored])
    )
    record["min_m_open"] = float(min(s["m_open"] for s in scored))
    record["mean_m_open"] = float(np.mean([s["m_open"] for s in scored]))
    return record, True


def _fill_block(n: int, base: int, slots: int, needed: int):
    records, skipped = [], []
    for k in range(slots):
        seed = base + STRIDE * k
        record, complete = run_sample(n, seed)
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
                     "patch_volume": PATCH_VOLUME, "per_n": {}}
    for n in P12_LADDER:
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
        print(f"verify-12 n={n}: {complete_count}/{VERIFY_COUNT} complete "
              f"| {elapsed / VERIFY_COUNT:.3f} s/sample", flush=True)
    summary["pin_passed"] = bool(all(
        summary["per_n"][str(n)]["complete"] >= VERIFY_PIN
        for n in P12_LADDER
    ))
    (output_dir / VERIFY_ARTIFACT).write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pin_passed": summary["pin_passed"]}, indent=2))


def run_pilot(output_dir: Path) -> None:
    stamp = _preflight_clean()
    # P12's instrument IS P11's Stage A estimator, so its flat
    # certification is the prerequisite (prereg Section 6)
    _require_stage_pass("p11_stage_a_summary.json", "P11 Stage A")
    verification = _load_gate_artifact(output_dir, VERIFY_ARTIFACT, stamp)
    if not verification.get("pin_passed"):
        raise SystemExit("verification-12 pin failed -- pilot refuses")

    summary: dict = {"code_version": stamp, "per_rung": {}}
    ys: dict = {}
    for n, (base, slots) in PILOT_BLOCKS.items():
        start = time.perf_counter()
        records, skipped, filled = _fill_block(n, base, slots, PILOT_SAMPLES)
        elapsed = time.perf_counter() - start
        if not filled:
            raise SystemExit(
                f"pilot-12 block n={n} could not fill {PILOT_SAMPLES} "
                f"(skips: {len(skipped)}) -- INFEASIBLE-INCOMPLETE")
        ys[n] = np.array([r["y"] for r in records])
        nominal = bonett_variance_bound(ys[n])
        calibrated = calibrated_variance_bound(
            ys[n], "bottom" if n == P12_LADDER[0] else "top")
        summary["per_rung"][str(n)] = {
            "n_samples": len(records),
            "variance": nominal["s2"],
            "kurtosis_g4": nominal["g4"],
            "variance_bound_95_nominal": nominal["bound"],
            "variance_bound_95_calibrated": calibrated["bound"],
            "calibration_z_used": calibrated["z_used"],
            "calibration_coverage_at_nominal":
                calibrated["coverage_at_nominal"],
            "skipped_seeds": skipped,
            "mean_seconds_per_sample": elapsed / len(records),
            "y": [float(val) for val in ys[n]],
        }
        print(f"pilot-12 n={n}: var={nominal['s2']:.5f} "
              f"g4={nominal['g4']:.2f} | skips={len(skipped)}", flush=True)

    power = power_requirements(ys[600], ys[2400])
    times = verification["per_n"]
    projected = None
    if power["n_per_rung"] is not None:
        projected = power["n_per_rung"] * sum(
            times[str(n)]["mean_seconds_per_sample"] for n in P12_LADDER
        ) / 3600.0
    summary["power"] = power
    summary["projected_stage_a_hours"] = projected
    summary["feasible"] = bool(
        not power["infeasible"] and projected is not None
        and projected <= PROJECTION_LIMIT_HOURS
    )
    (output_dir / PILOT_ARTIFACT).write_text(
        json.dumps(_de_nan(summary), indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"power": power,
                      "projected_stage_a_hours": projected,
                      "feasible": summary["feasible"]}, indent=2))


def run_stage_a(output_dir: Path) -> None:
    stamp = _preflight_clean()
    pilot = _load_gate_artifact(output_dir, PILOT_ARTIFACT, stamp)
    if not pilot.get("feasible"):
        raise SystemExit("pilot-12 declared the design infeasible -- "
                         "Stage A refuses to run")
    n_per_rung = int(pilot["power"]["n_per_rung"])
    flat_available = bool(pilot["power"]["flat_available"])

    rows, per_rung_y, skip_counts, skipped_seeds = [], {}, {}, {}
    for n, (base, slots) in STAGE_A_BLOCKS.items():
        records, skipped, filled = _fill_block(n, base, slots, n_per_rung)
        if not filled:
            raise SystemExit(
                f"stage A-12 block n={n} could not fill {n_per_rung} "
                f"(skips: {len(skipped)}) -- INFEASIBLE-INCOMPLETE")
        for r in records:
            r["stage"] = "A12"
            r["code_version"] = stamp
        rows.extend(records)
        per_rung_y[n] = np.array([r["y"] for r in records])
        skip_counts[n] = len(skipped)
        skipped_seeds[n] = [int(s) for s in skipped]
        print(f"stage A-12 n={n}: {len(records)} complete | mean y "
              f"{per_rung_y[n].mean():.4f} | skips {len(skipped)}",
              flush=True)
    write_rows_csv(output_dir / "p12_stage_a.csv", rows)

    delta = float(per_rung_y[2400].mean() - per_rung_y[600].mean())
    rng = np.random.default_rng(_stable_seed("p12-a-delta"))
    boots = [
        float(rng.choice(per_rung_y[2400], size=n_per_rung,
                         replace=True).mean()
               - rng.choice(per_rung_y[600], size=n_per_rung,
                            replace=True).mean())
        for _ in range(4000)
    ]
    lo, hi = (float(np.percentile(boots, 2.5)),
              float(np.percentile(boots, 97.5)))

    log_n = np.log10(np.array(P12_LADDER, dtype=float))
    means = np.array([per_rung_y[n].mean() for n in P12_LADDER])
    slope = float(np.polyfit(log_n, means, 1)[0])
    srng = np.random.default_rng(_stable_seed("p12-a-slope"))
    slope_boots = [
        float(np.polyfit(log_n, [
            srng.choice(per_rung_y[n], size=n_per_rung,
                        replace=True).mean() for n in P12_LADDER
        ], 1)[0])
        for _ in range(4000)
    ]

    summary = {
        "code_version": stamp,
        "n_per_rung": n_per_rung,
        "flat_available": flat_available,
        "ricci_scalar": RICCI_SCALAR,
        "skip_counts": {str(n): skip_counts[n] for n in P12_LADDER},
        "skipped_seeds": {str(n): skipped_seeds[n] for n in P12_LADDER},
        "selection_caveat": bool(any(skip_counts[n] > 0
                                     for n in P12_LADDER)),
        "delta": delta, "delta_ci": [lo, hi],
        "verdict": verdict(lo, hi, flat_available),
        "mean_y_by_rung": {str(n): float(per_rung_y[n].mean())
                           for n in P12_LADDER},
        "middle_rung_between_endpoints": bool(
            min(means[0], means[2]) <= means[1] <= max(means[0], means[2])
        ),
        "labelled_checks": {
            "slope_mean_y_vs_log10N": slope,
            "slope_ci": [float(np.percentile(slope_boots, 2.5)),
                         float(np.percentile(slope_boots, 97.5))],
            "predicted_slope": -1.0 / 3.0,
            "median_relerr_curved_by_rung": {
                str(n): float(np.median(
                    [r["median_relerr_curved"] for r in rows
                     if r["n_expected"] == n]
                )) for n in P12_LADDER
            },
            "flat_arm_median_relerr_by_rung": {
                str(n): float(np.median(
                    [r["median_relerr_flat_arm"] for r in rows
                     if r["n_expected"] == n]
                )) for n in P12_LADDER
            },
            "flat_arm_frozen_prediction": [0.42, 0.49],
            "mean_m_open_by_rung": {
                str(n): float(np.mean(
                    [r["mean_m_open"] for r in rows if r["n_expected"] == n]
                )) for n in P12_LADDER
            },
        },
    }
    (output_dir / "p12_stage_a_summary.json").write_text(
        json.dumps(_de_nan(summary), indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"delta": delta, "delta_ci": [lo, hi],
                      "verdict": summary["verdict"]}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["verify", "pilot", "a", "b"],
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
    elif args.stage == "b":
        raise SystemExit(
            "stage b is gated on its frozen addendum (prereg Section 5), "
            "which must budget |R_hat - R|/R per rung and use flat-twin "
            "differencing -- this refusal is the preregistration "
            "operating."
        )
    else:
        raise SystemExit("choose --stage verify/pilot/a")


if __name__ == "__main__":
    main()

# Keep the unused-import guard honest about constants the prereg pins
# but the runners reach only through P11's helpers.
_PINNED = (BONETT_Z, N_SUP_COEFF, N_EQ_COEFF, N_FLOOR, N_CAP)
