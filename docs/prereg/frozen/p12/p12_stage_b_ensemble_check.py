"""Stage B design check on the ENSEMBLE: what the budget actually is.

Section 5 item 4 demands the expected `|R_hat - R| / R` be stated per
rung BEFORE running, precisely so "a gate that passes while recovering
nothing" is a visible outcome rather than a surprise. This script
measures it, together with the two other quantities the addendum has to
freeze: the per-pair dispersion of the observable, and whether six
disjoint pairs are realizable at the BOTTOM rung (item 3).

THE OBSERVABLE. Per pair, order-plus-count only:

    g = m_open / (rho_hat * tau_hat^2)

with `m_open` the open interval cardinality (order data), `rho_hat` the
frozen density calibration, and `tau_hat` the unchanged P11 chain
estimator. In the continuum `g -> Vol/tau^2 = G(tau/2ell)`, which is
1/2 in flat space and less than 1/2 under positive curvature.

THE CONSTRUCTION (Section 5 item 5, the flat-twin difference, sharpened
into a RATIO). Both `m` and `tau_hat` carry discreteness bias, and the
bias is a function of `m` alone. Matching `m` between a curved arm and
a flat twin therefore cancels it to first order IN THE RATIO

    Q = g_curved / g_flat  ->  G(s) / (1/2),     s = tau / (2 ell)

so that `1 - Q = 1 - 2 G(s)`, and inverting G gives

    R_hat = 8 s^2 / tau_hat^2,    s = G^-1( Q / 2 ).

A difference would leave the bias in the units; the ratio removes it and
leaves a dimensionless quantity, which is why the ratio is what this
design freezes.

WHY tau/ell = 1.5 AND NOT 0.3. The volume check measured the
inversion's amplification -- relative error in `g` per relative error in
`R_hat` -- at 268.8x for tau/ell = 0.30 and 12.9x for 1.50. P12 Stage A
ran at 0.3 because that was where the chain law was known to hold; P13
has since certified the same estimator out to tau/ell = 1.5
(CURVATURE-ROBUST, its Section 13), so Stage B can take the rung where
the signal is 22x larger and the amplification 21x smaller. The patch
constants are P13's 1.50 rung verbatim.

Seeds: P12's design-check space `2000000-5999999` (Section 6), disjoint
from every experimental window. Run from the repository root.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "experiments/positive_control")
import p13_tau_ell as P  # noqa: E402

OUT = Path("docs/prereg/frozen/p12/p12_stage_b_ensemble_check.json")

TAU_ELL = 1.50                     # P13's top rung, certified there
ETA_LO, XHALF, RHO_TOP = P.PATCH[TAU_ELL]
TWIN_CENTRE = P.TWIN_BAND_CENTRE[TAU_ELL]

#: Density ladder at FIXED geometry, 4x end to end as P12 Stage A's is.
#: RHO_TOP is P13's calibrated intensity for mean m ~ 76.
LADDER = {"600": RHO_TOP / 4.0, "1200": RHO_TOP / 2.0, "2400": RHO_TOP}
TWIN_RHO_TOP = P.TWIN_RHO[TAU_ELL]
TWIN_LADDER = {"600": TWIN_RHO_TOP / 4.0, "1200": TWIN_RHO_TOP / 2.0,
               "2400": TWIN_RHO_TOP}

N_SAMPLES = 400
BASE = {"600": 2_000_000, "1200": 2_100_000, "2400": 2_200_000}
TWIN_BASE = {"600": 2_300_000, "1200": 2_400_000, "2400": 2_500_000}
DESIGN_FLOOR, DESIGN_CEIL = 2_000_000, 5_999_999


def g_of_s(s):
    if s < 1e-8:
        return 0.5 - s ** 2 / 12.0
    return math.log(math.cosh(s)) / s ** 2


def invert_g(g):
    """s = G^-1(g); None when g >= 1/2, i.e. when the measurement is at
    or beyond the flat value and no positive R is representable. Section
    5 item 4's failure mode, counted rather than hidden."""

    if not (0.0 < g < 0.5):
        return None
    lo, hi = 1e-8, 1.0
    while g_of_s(hi) > g:
        hi *= 2.0
        if hi > 1e6:
            return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g_of_s(mid) > g:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def sample_pairs(rho, seed, flat):
    """One sample through P13's own code path, returning the per-pair
    (g, tau_hat, m) rather than P13's aggregated y."""

    band_centre = TWIN_CENTRE if flat else TAU_ELL
    band = (band_centre * 0.9, band_centre * 1.1)
    rng = np.random.default_rng(seed)
    eta, x = P.sprinkle(ETA_LO, XHALF, rho, rng, flat=flat)
    volume = (P.patch_proper_volume(ETA_LO, XHALF) if not flat
              else (abs(ETA_LO) - 1.0) * 2.0 * XHALF)
    rho_hat = eta.size / volume
    pool = P.eligible_pairs(eta, x, ETA_LO, XHALF, band, flat=flat)
    u, v = eta + x, eta - x
    draw = np.random.default_rng(seed + P.SCORE_OFFSET)
    pairs, complete, _rej = P.draw_disjoint(pool, u, v, draw)
    if not complete:
        return None
    out = []
    for i, j in pairs:
        scored = P.score_pair(u / 2.0, v / 2.0, i, j)
        tau_hat = float(P.estimate_tau_from_longest_chain_1p1(
            scored["chain_length"], rho=rho_hat))
        m_open = int(scored["m_open"])
        if tau_hat <= 0.0:
            continue
        out.append({"g": m_open / (rho_hat * tau_hat ** 2),
                    "tau_hat": tau_hat, "m": m_open})
    return out


def run_arm(ladder, bases, flat):
    arm = {}
    for rung, rho in ladder.items():
        gs, taus, ms, complete = [], [], [], 0
        for k in range(N_SAMPLES):
            got = sample_pairs(rho, bases[rung] + P.STRIDE * k, flat)
            if got is None:
                continue
            complete += 1
            gs.append([p["g"] for p in got])
            taus.extend(p["tau_hat"] for p in got)
            ms.extend(p["m"] for p in got)
        flatg = np.array([g for row in gs for g in row])
        arm[rung] = {
            "rho": float(rho),
            "n_attempted": N_SAMPLES, "n_complete": complete,
            "completion": complete / N_SAMPLES,
            "pairs": int(flatg.size),
            "mean_g": float(flatg.mean()),
            "sd_g": float(flatg.std(ddof=1)),
            "rel_sd_g": float(flatg.std(ddof=1) / flatg.mean()),
            "se_mean_g": float(flatg.std(ddof=1) / math.sqrt(flatg.size)),
            "mean_m": float(np.mean(ms)),
            "mean_tau_hat": float(np.mean(taus)),
            "per_sample_mean_g": [float(np.mean(row)) for row in gs],
        }
        a = arm[rung]
        print(f"  {'twin' if flat else 'curved'} rung {rung}:"
              f" completion {a['completion']:.3f} | pairs {a['pairs']}"
              f" | mean g {a['mean_g']:.5f} (rel sd {a['rel_sd_g']:.3f},"
              f" se {a['se_mean_g']:.5f}) | mean m {a['mean_m']:.2f}"
              f" | mean tau_hat {a['mean_tau_hat']:.4f}", flush=True)
    return arm


def main() -> None:
    for base in list(BASE.values()) + list(TWIN_BASE.values()):
        assert DESIGN_FLOOR <= base <= DESIGN_CEIL, base
        assert base + P.STRIDE * N_SAMPLES <= DESIGN_CEIL

    results: dict = {
        "tau_over_ell": TAU_ELL,
        "patch": {"eta_lo": ETA_LO, "xhalf": XHALF},
        "twin_band_centre": TWIN_CENTRE,
        "ladder_rho": {k: float(v) for k, v in LADDER.items()},
        "twin_ladder_rho": {k: float(v) for k, v in TWIN_LADDER.items()},
        "n_samples_per_rung": N_SAMPLES,
        "seed_space": [DESIGN_FLOOR, DESIGN_CEIL],
        "g_true_continuum": g_of_s(TAU_ELL / 2.0),
        "one_minus_Q_true": 1.0 - g_of_s(TAU_ELL / 2.0) / 0.5,
    }
    print("=== continuum targets ===")
    print(f"  G(s) = {results['g_true_continuum']:.6f}"
          f" | 1 - Q = {results['one_minus_Q_true']:.6f}", flush=True)

    print("=== curved arm ===")
    results["curved"] = run_arm(LADDER, BASE, flat=False)
    print("=== flat twin ===")
    results["twin"] = run_arm(TWIN_LADDER, TWIN_BASE, flat=True)

    print("=== the budget Section 5 item 4 asks for ===")
    budget = {}
    for rung in LADDER:
        c, t = results["curved"][rung], results["twin"][rung]
        q = c["mean_g"] / t["mean_g"]
        s = invert_g(q / 2.0)
        pooled = None
        if s is not None:
            r_hat = 8.0 * s ** 2 / c["mean_tau_hat"] ** 2
            pooled = abs(r_hat - 2.0) / 2.0
        # per-SAMPLE recovery: each sample's own six pairs against the
        # twin's rung mean, which is how the campaign statistic would be
        # built if it followed the P11/P12 per-sample shape
        per_sample, undefined = [], 0
        for gm in c["per_sample_mean_g"]:
            s_i = invert_g((gm / t["mean_g"]) / 2.0)
            if s_i is None:
                undefined += 1
                continue
            r_i = 8.0 * s_i ** 2 / c["mean_tau_hat"] ** 2
            per_sample.append(abs(r_i - 2.0) / 2.0)
        # the se on Q, and what it implies for the pooled recovery
        rel_se_q = math.sqrt((c["se_mean_g"] / c["mean_g"]) ** 2
                             + (t["se_mean_g"] / t["mean_g"]) ** 2)
        budget[rung] = {
            "Q": q, "one_minus_Q": 1.0 - q,
            "R_hat_pooled": (8.0 * s ** 2 / c["mean_tau_hat"] ** 2
                             if s is not None else None),
            "pooled_rel_error": pooled,
            "rel_se_Q": rel_se_q,
            "se_one_minus_Q": rel_se_q * q,
            "signal_to_noise_one_minus_Q":
                (1.0 - q) / (rel_se_q * q) if q > 0 else None,
            "per_sample_median_rel_error":
                float(np.median(per_sample)) if per_sample else None,
            "per_sample_undefined_count": undefined,
            "per_sample_undefined_fraction":
                undefined / len(c["per_sample_mean_g"]),
        }
        b = budget[rung]
        pooled_txt = ("undefined" if b["pooled_rel_error"] is None
                      else f"{b['pooled_rel_error']:.3f}")
        per_txt = ("n/a" if b["per_sample_median_rel_error"] is None
                   else f"{b['per_sample_median_rel_error']:.2f}")
        print(f"  rung {rung}: 1-Q {b['one_minus_Q']:+.5f}"
              f" (true {results['one_minus_Q_true']:.5f}, se"
              f" {b['se_one_minus_Q']:.5f}, S/N"
              f" {b['signal_to_noise_one_minus_Q']:.1f})"
              f" | pooled |R_hat-R|/R {pooled_txt}"
              f" | per-sample median {per_txt}"
              f" | undefined {b['per_sample_undefined_fraction']:.1%}",
              flush=True)
    results["budget"] = budget

    # --- Section 5 item 2: Delta*_B derived, not borrowed ------------
    # The FORM is theory: g = m / (rho tau_hat^2) inherits a Poisson
    # 1/m relative variance from the count and a BDJ m^(-2/3) one from
    # the chain (sd(L) ~ m^(1/6) against E L ~ 2 sqrt(m), doubled by
    # the square). The two COEFFICIENTS are calibrated here, because m
    # and L are positively correlated -- a fuller box has a longer
    # chain -- so an analytic-only coefficient would overstate the
    # dispersion of their ratio. The correlation is measured alongside
    # so the discrepancy is documented rather than absorbed.
    ms = np.array([results["curved"][r]["mean_m"] for r in LADDER])
    sds = np.array([results["curved"][r]["rel_sd_g"] for r in LADDER])
    design = np.column_stack([1.0 / ms, ms ** (-2.0 / 3.0)])
    coef, *_ = np.linalg.lstsq(design, sds ** 2, rcond=None)
    a_coef, b_coef = float(coef[0]), float(coef[1])

    def rel_sd_model(m):
        return math.sqrt(a_coef / m + b_coef * m ** (-2.0 / 3.0))

    m_bot, m_top = float(ms[0]), float(ms[-1])
    delta_star_b = math.log10(rel_sd_model(m_top) / rel_sd_model(m_bot))
    results["delta_star_b"] = {
        "poisson_coefficient_a": a_coef,
        "chain_coefficient_b": b_coef,
        "analytic_b_if_uncorrelated": (2.0 * 0.4509) ** 2,
        "model": "rel_sd(g)^2 = a/m + b m^(-2/3)",
        "fitted_vs_measured": [
            {"m": float(m), "measured": float(s),
             "model": rel_sd_model(float(m))}
            for m, s in zip(ms, sds, strict=True)],
        "m_bottom": m_bot, "m_top": m_top,
        "delta_star_b": delta_star_b,
        "chain_only_would_be": -math.log10(4.0) / 3.0,
    }
    print("=== Delta*_B, derived from the propagated input rates ===")
    print(f"  rel_sd(g)^2 = {a_coef:.4f}/m + {b_coef:.4f} m^(-2/3)"
          f"  (analytic b if uncorrelated: "
          f"{(2.0 * 0.4509) ** 2:.4f})")
    for row in results["delta_star_b"]["fitted_vs_measured"]:
        print(f"    m {row['m']:6.2f}: measured {row['measured']:.4f}"
              f" | model {row['model']:.4f}")
    print(f"  Delta*_B = {delta_star_b:+.4f} dex"
          f"  (the chain-only rate would be"
          f" {-math.log10(4.0) / 3.0:+.4f})")

    # --- power: the per-sample statistic and its variance ------------
    # y_B = log10(|R_hat - R| / R) with R_hat the BOUNDARY-CONSTRAINED
    # estimate: the positive-curvature domain is g_ratio < 1/2, and a
    # sample outside it has recovered no curvature, so R_hat = 0 and
    # the relative error is exactly 1. That is the boundary of the
    # parameter space, not a clamp over noise, and the boundary-hit
    # fraction is published per rung.
    ys = {}
    for rung in LADDER:
        c, t = results["curved"][rung], results["twin"][rung]
        vals, hits = [], 0
        for gm in c["per_sample_mean_g"]:
            s_i = invert_g((gm / t["mean_g"]) / 2.0)
            if s_i is None:
                hits += 1
                rel = 1.0
            else:
                rel = abs(8.0 * s_i ** 2 / c["mean_tau_hat"] ** 2
                          - 2.0) / 2.0
            vals.append(math.log10(rel) if rel > 0 else -12.0)
        arr = np.array(vals)
        ys[rung] = {
            "n": int(arr.size), "mean_y": float(arr.mean()),
            "variance": float(arr.var(ddof=1)),
            "boundary_hits": hits,
            "boundary_fraction": hits / arr.size,
            "median_rel_error": float(10 ** np.median(arr)),
        }
        print(f"  rung {rung}: mean y_B {ys[rung]['mean_y']:+.4f}"
              f" | var {ys[rung]['variance']:.5f}"
              f" | boundary hits {ys[rung]['boundary_fraction']:.1%}")
    results["per_sample_y_B"] = ys

    # --- the DIMENSIONLESS recovery, which is what Q alone gives -----
    # R_hat needs a length and the only one available is tau_hat, whose
    # chain bias then enters once. The combination R tau^2 = 8 s^2 needs
    # NO length at all, so it separates what the ratio recovers from
    # what the normalization spoils. Reported per rung so the record can
    # quote the split instead of asserting it.
    dimensionless = {}
    r_tau2_true = 2.0 * TAU_ELL ** 2
    for rung in LADDER:
        c, t = results["curved"][rung], results["twin"][rung]
        s = invert_g((c["mean_g"] / t["mean_g"]) / 2.0)
        r_tau2 = 8.0 * s ** 2 if s is not None else 0.0
        dimensionless[rung] = {
            "R_tau2_hat": r_tau2, "R_tau2_true": r_tau2_true,
            "rel_error": abs(r_tau2 - r_tau2_true) / r_tau2_true,
            "tau_hat_over_tau_true": c["mean_tau_hat"] / TAU_ELL,
        }
        d = dimensionless[rung]
        print(f"  rung {rung}: R tau^2 hat {r_tau2:.4f}"
              f" (true {r_tau2_true:.4f}) -> rel err"
              f" {d['rel_error']:.3f} | tau_hat/tau"
              f" {d['tau_hat_over_tau_true']:.3f}")
    results["dimensionless_recovery"] = dimensionless

    # --- power, on P11's conventions ---------------------------------
    # P11 sizes n_sup to detect the DERIVED target itself and sets
    # delta_eq = |Delta*|/3. Both are carried over with Delta*_B in
    # place of Delta*_A. The 1.645 in the equivalence slot is
    # Phi^-1(0.95), i.e. ALREADY the two-sided quantile for beta = 0.10
    # -- noted because P13's review C2 caught exactly this slot being
    # "upgraded" with a one-sided quantile.
    s2 = ys["600"]["variance"] + ys["2400"]["variance"]
    delta_eq_b = abs(delta_star_b) / 3.0
    n_sup = math.ceil(s2 * (1.960 + 1.282) ** 2 / delta_star_b ** 2)
    n_eq = math.ceil(s2 * (1.960 + 1.645) ** 2 / delta_eq_b ** 2)
    results["power_projection"] = {
        "s2_uncalibrated_sum": s2,
        "delta_star_b": delta_star_b,
        "delta_eq_b": delta_eq_b,
        "n_sup_uncalibrated": n_sup,
        "n_eq_uncalibrated": n_eq,
        "note": "UNCALIBRATED and indicative only. Stage P-B computes "
                "the calibrated Bonett bounds on fresh pilot seeds and "
                "the frozen formulas decide n; these numbers exist so "
                "the addendum can state affordability before freezing.",
    }
    print(f"  indicative (UNCALIBRATED) s2 {s2:.4f}"
          f" | detect {delta_star_b:+.4f} -> n_sup {n_sup}"
          f" | delta_eq {delta_eq_b:.4f} -> n_eq {n_eq}")

    # --- why b is below its uncorrelated value ----------------------
    results["m_chain_correlation_note"] = (
        "measured b sits below the uncorrelated analytic value because "
        "m and L rise together within a box, so their ratio disperses "
        "less than the independent sum would")

    OUT.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
