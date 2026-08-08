"""P14 §8 P3-E: discriminability, exploration stage.

EXPLORATORY. No gate, no threshold, no verdict, no confirmation seed
window appears in this file, and none may be added to it (design
§8.2). P3-E explores and CHOOSES -- the primary operating point, the
primary discriminator, the fit ranges -- and everything it chooses is
frozen as input to P3-C, which runs on a fresh seed block and is the
ONLY stage whose composite equivalence gate may close anything. The
frozen termination sentence (P3-C's, never P3-E's):

    "동결된 운영점에서, 3a의 (i) 표준화 평균이동,
    (ii) 순위 판별(AUC), (iii) 중점-역치 분류 --
    이 세 규칙으로는 분리가 해상되지 않았다."

It closes the chosen operating point under those three rules only --
not 3a's full distribution, not any other pure-order statistic (the
N(0,1)-vs-N(0,4) counterexample: mean shift, rank discrimination, and
a midpoint threshold are all blind to a pure variance separation).

Protocol frozen in the P3 design review (five rounds, in-session):

- **Ladder-E** (exponent): ONE box (slice-a1.0) and ONE density
  (E[N] = 300); only `w` varies, `A = w^2`. The P1 seven-point table
  varies box and density together and is NEVER used for exponent
  fits. Regression axis is `log A` (slope in `w` would read 2 where
  the heuristic says O(|A|)). Shared point sets across rungs, so the
  ladder is one sprinkling cluster: every CI is a sprinkling-level
  cluster bootstrap. `D = gained + lost` and `dr = gained - lost` are
  regressed SEPARATELY: the review's structure is D = O(A) while
  dr = O(A^2), and the pilot's gained/lost ratio falls toward 1 as
  A -> 0 (72:1 at A = 4 is finite-amplitude imbalance, not a
  leading-order signature). Fits: D over the full ladder and over
  A <= 0.5; dr asymptotic over A <= 0.0625 only (the pilot shows
  (gained-lost)/A^2 still drifting through A ~ 0.5) plus a
  finite-range A <= 0.5 fit reported as an effective exponent. The
  three lowest rungs carry TOTAL n = 800 sprinklings (400 shared +
  400 low-A-only), sized for a 95% slope-CI half-width <= 0.15.
  Bootstrap replicates where a rung mean is <= 0 are recorded as fit
  failures, never discarded; the drop-highest-rung refit is a REFIT
  DIAGNOSTIC (a two-point slope), not an independent validation.
  `w = 0` is a pipeline null only, excluded from every log fit.
  The -A contrast (polarization rotated 90 deg) is the x<->y swap on
  the square transverse box -- the geometry does not admit negative
  `w`; D and gained/lost must match the +A DISTRIBUTIONS (the axis
  swap itself is statistic #4's business, not gained<->lost).

- **Ladder-G** (operating-point map): the P1 seven points, n = 250
  each (CV <= 0.31 -> RSE(mean D) <= 2%). D as the §5.1 interval
  [D, D + ambiguous], gained/lost, per-arm 3a, and the paired
  difference. AUC and s are reported as REALIZED precision, never as
  sizing bases (SE(s) ~ sqrt(2/n) holds only near s = 0; at the
  pilot's s ~ 5-12 the realized SE is ~0.3-0.6).

- **Class C characterization** (descriptive ONLY, this PR): at
  aniso-a1.0, densities E[N] in {1200, 2400}, a FIXED n = 72
  sprinklings each, both arms. Reported: causal-interval counts,
  cardinality profiles, opportunity counts sum C(m_i, k) for
  k = 2, 3, 4, raw chain counts, opportunity-weighted
  C_k = sum(chains_k) / sum(C(m_i, k)) -- the frozen PRIMARY
  normalization (no cardinality cutoff; m_i < k contributes zero to
  both sides); C2, C3 are primary-candidate, C4 exploratory.
  Per-arm AND union zero-denominator rates with exact CP intervals,
  and the FULL-pipeline wall-clock (census + interval extraction +
  chain enumeration + aggregation) median and p95 with the CPU
  recorded. NO pass/fail, NO density freeze here: the 0/72 rule was
  rejected in review (a true 1% zero-rate passes 0/72 with
  probability ~48.5%, and "expand n on a failure" is optional
  stopping) -- the freeze gate is a later, power-sized decision.

- **P3-C preliminary calibration** (code only, no campaign): margins
  FROZEN at eps_s = 0.0806, eps_AUC = 0.0233, eps_BA = 0.0285 (the
  accuracy standard never moves to meet a power result); the joint
  90% power is met by SAMPLE SIZE. Preliminary candidate
  n = 4800/arm under the normal design model: joint pass 18273/20000
  = 0.9136, exact CP 95% [0.9097, 0.9175], lower bound >= 0.90. The
  artifact stores the RAW integers, reps, seed path, full-precision
  margins, and CP endpoints. After P3-E, the power model is refit on
  the measured 3a distribution AT THE SAME margins, the final n is
  recomputed and committed before P3-C runs, and P3-C is separately
  approved.

Seeds: one stream per component (the consumption unit; substreams
via `default_rng([seed, tag])`). Burned and never reused for any
campaign: 20260808 (P1), 777 (P2 design), 778/781 (P3 pilots),
779/780 (independent review), 20260811-14 (P2 campaign).

Run:  python experiments/positive_control/p14_probe_p3e.py
"""

from __future__ import annotations

import json
import math
import platform
import time
from pathlib import Path

import numpy as np
from p14_plane_wave import Slab, arms, sprinkle
from p14_probe_p1 import clopper_pearson, element_eligible, relation_census

# ---------------------------------------------------------------------
# Frozen protocol constants (P3 design review, 2026-08-09)
# ---------------------------------------------------------------------

LADDER_SLAB = dict(du=1.0, dv=0.2, dx=6.0, dy=6.0)  # slice-a1.0, square
LADDER_EN = 300
W_LADDER = (0.125, 0.177, 0.25, 0.354, 0.5, 0.707, 1.0, 1.414, 2.0)
LOW_A_RUNGS = W_LADDER[:3]          # A <= 0.0625: the asymptotic-fit set
N_LADDER = 400                      # sprinklings computing ALL rungs
N_LOW_EXTRA = 400                   # + low-A-only -> TOTAL 800 per low rung
W_SWAP = 1.0                        # the -A contrast rung
FIT_DR_ASYMPTOTIC_MAX_A = 0.0625
FIT_FINITE_MAX_A = 0.5

#: (label, w, du, dv, dx, dy) -- P1's seven points, operating-point
#: MAP only; never an exponent ladder (box and density co-vary).
LADDER_G_POINTS = (
    ("slice-a0.3", 1.0, 0.3, 0.2, 6.0, 6.0),
    ("slice-a0.6", 1.0, 0.6, 0.2, 6.0, 6.0),
    ("slice-a1.0", 1.0, 1.0, 0.2, 6.0, 6.0),
    ("aniso-a1.0", 1.0, 1.0, 1.0, 2.0, 6.0),
    ("roomy-a0.2", 1.0, 0.2, 16.0, 6.0, 6.0),
    ("high-a2.0", 1.0, 2.0, 0.5, 2.0, 2.0),
    ("edge-a2.4", 1.0, 2.4, 0.2, 1.2, 0.8),
)
LADDER_G_EN = 300
N_G = 250

CLASSC_POINT = ("aniso-a1.0", 1.0, 1.0, 1.0, 2.0, 6.0)
CLASSC_DENSITIES = (1200, 2400)
N_CLASSC = 72                       # FIXED; descriptive, no pass/fail
CHAIN_KS = (2, 3, 4)                # C2, C3 primary-candidate; C4 expl.

B_BOOT = 4_000

#: P3-C margins, FROZEN (never rescaled to meet power); the
#: preliminary candidate n and its certification record.
P3C_MARGINS = {
    "eps_s": 3.605 * math.sqrt(2.0 / 4000.0),
    "eps_auc": 3.605 * math.sqrt((2 * 4000 + 1) / (12.0 * 4000 * 4000)),
    "eps_ba": 3.605 / (2.0 * math.sqrt(4000.0)),
}
P3C_N_CANDIDATE = 4_800             # PRELIMINARY (normal design model);
P3C_POWER_REPS = 20_000             # final n recomputed from P3-E dist.
P3C_POWER_SEED = (781, 4820)

SEEDS = {
    "ladder_e": 20260821,
    "low_a": 20260822,
    "ladder_g": 20260823,
    "class_c": 20260824,
}
BURNED_SEEDS = (20260808, 777, 778, 779, 780, 781,
                20260811, 20260812, 20260813, 20260814)

_ARTIFACT = (Path(__file__).resolve().parents[2]
             / "docs" / "prereg" / "p14_probe_p3e_results.json")


def assert_seed_layout() -> None:
    seeds = list(SEEDS.values())
    assert len(set(seeds)) == len(seeds), "seed collision"
    clash = set(seeds) & set(BURNED_SEEDS)
    assert not clash, f"burned seeds reused: {sorted(clash)}"


# ---------------------------------------------------------------------
# Pair statistics from two censuses of the SAME points
# ---------------------------------------------------------------------


def pair_flips(c_flat, c_curved) -> tuple[int, int, int]:
    """(gained, lost, ambiguous_union) between the two relations."""

    gained = int((c_curved.related & ~c_flat.related).sum())
    lost = int((~c_curved.related & c_flat.related).sum())
    amb = len(set(c_flat.ambiguous_pairs) | set(c_curved.ambiguous_pairs))
    return gained, lost, amb


# ---------------------------------------------------------------------
# Ladder-E
# ---------------------------------------------------------------------


def run_ladder_e() -> dict:
    """Raw per-sprinkling rung records; two strata (full / low-A)."""

    slab = Slab(**LADDER_SLAB)
    rho = LADDER_EN / slab.coordinate_volume
    geoms = {w: arms(slab, w)[0] for w in W_LADDER}
    flat = arms(slab, 1.0)[1]

    def one_stratum(seed: int, n: int, rungs) -> dict:
        rng = np.random.default_rng(seed)
        rec = {w: [] for w in rungs}
        swap = []
        for _ in range(n):
            pts = sprinkle(flat, rho, rng)
            npts = len(pts)
            pairs = npts * (npts - 1) // 2
            c0 = relation_census(flat, pts)
            for w in rungs:
                cA = relation_census(geoms[w], pts)
                g, lo, amb = pair_flips(c0, cA)
                rec[w].append((g / pairs, lo / pairs, amb / pairs))
            if W_SWAP in rungs:
                cS = relation_census(geoms[W_SWAP], pts[:, [0, 1, 3, 2]])
                g, lo, _a = pair_flips(c0, cS)
                swap.append((g / pairs, lo / pairs))
        return {"rungs": {str(w): rec[w] for w in rungs}, "swap": swap,
                "n": n}

    full = one_stratum(SEEDS["ladder_e"], N_LADDER, W_LADDER)
    low = one_stratum(SEEDS["low_a"], N_LOW_EXTRA, LOW_A_RUNGS)
    return {"full": full, "low": low}


def _rung_means(raw: dict, rng: np.random.Generator | None) -> dict:
    """Per-rung mean gained/lost fractions; optionally a bootstrap
    resample WITHIN each stratum (the two strata are independent
    clusters and are resampled independently)."""

    out = {}
    for w in W_LADDER:
        key = str(w)
        rows = list(raw["full"]["rungs"][key])
        if key in raw["low"]["rungs"]:
            rows += list(raw["low"]["rungs"][key])
        arr = np.asarray(rows)
        if rng is not None:
            n_full = len(raw["full"]["rungs"][key])
            i_full = rng.integers(0, n_full, n_full)
            parts = [np.asarray(raw["full"]["rungs"][key])[i_full]]
            if key in raw["low"]["rungs"]:
                n_low = len(raw["low"]["rungs"][key])
                i_low = rng.integers(0, n_low, n_low)
                parts.append(np.asarray(raw["low"]["rungs"][key])[i_low])
            arr = np.concatenate(parts)
        g, lo = arr[:, 0].mean(), arr[:, 1].mean()
        out[w] = {"d": g + lo, "dr": g - lo, "gained": g, "lost": lo,
                  "n": len(arr)}
    return out


def _fit_slope(means: dict, quantity: str, max_a: float,
               min_a: float = 0.0) -> float | None:
    """Log-log LS slope of `quantity` vs A over rungs in (min_a,
    max_a]; None (a fit failure) if any rung mean is <= 0."""

    xs, ys = [], []
    for w in W_LADDER:
        a = w * w
        if min_a < a <= max_a:
            v = means[w][quantity]
            if v <= 0.0:
                return None
            xs.append(math.log(a))
            ys.append(math.log(v))
    return float(np.polyfit(xs, ys, 1)[0])


#: (name, quantity, max_A, drop_highest_for_diagnostic)
_FITS = (
    ("d_full", "d", 4.0),
    ("d_low", "d", FIT_FINITE_MAX_A),
    ("dr_asymptotic", "dr", FIT_DR_ASYMPTOTIC_MAX_A),
    ("dr_finite_range", "dr", FIT_FINITE_MAX_A),
)


def ladder_e_fits(raw: dict) -> dict:
    """Point fits + sprinkling-cluster bootstrap CIs + failure counts
    + the drop-highest-rung REFIT DIAGNOSTIC (not a validation)."""

    means = _rung_means(raw, None)
    rng = np.random.default_rng([SEEDS["ladder_e"], 999])
    boots: dict = {name: [] for name, *_ in _FITS}
    failures = {name: 0 for name, *_ in _FITS}
    for _ in range(B_BOOT):
        bm = _rung_means(raw, rng)
        for name, qty, max_a in _FITS:
            s = _fit_slope(bm, qty, max_a)
            if s is None:
                failures[name] += 1
            else:
                boots[name].append(s)
    out = {"rungs": {str(w): means[w] for w in W_LADDER}}
    for name, qty, max_a in _FITS:
        bs = np.asarray(boots[name])
        point = _fit_slope(means, qty, max_a)
        # refit diagnostic: drop the highest rung of THIS fit range
        in_range = [w for w in W_LADDER if w * w <= max_a]
        drop = _fit_slope(means, qty, in_range[-2] ** 2) \
            if len(in_range) >= 3 else None
        out[name] = {
            "slope": point,
            "ci95": ([float(np.percentile(bs, 2.5)),
                      float(np.percentile(bs, 97.5))]
                     if len(bs) else None),
            "boot_reps": B_BOOT,
            "fit_failures": failures[name],
            "refit_diagnostic_drop_highest": drop,
        }
    # the -A contrast, descriptive: distribution match of D
    sw = np.asarray(raw["full"]["swap"])
    plus = np.asarray(raw["full"]["rungs"][str(W_SWAP)])[:len(sw)]
    d_p = plus[:, 0] + plus[:, 1]
    d_m = sw[:, 0] + sw[:, 1]
    out["swap_contrast"] = {
        "w": W_SWAP,
        "mean_d_plus": float(d_p.mean()), "mean_d_minus": float(d_m.mean()),
        "paired_diff": float((d_p - d_m).mean()),
        "paired_diff_se": float((d_p - d_m).std(ddof=1)
                                / math.sqrt(len(d_p))),
        "mean_gained_minus": float(sw[:, 0].mean()),
        "mean_lost_minus": float(sw[:, 1].mean()),
    }
    return out


# ---------------------------------------------------------------------
# Ladder-G
# ---------------------------------------------------------------------


def run_ladder_g() -> list[dict]:
    rows = []
    for k, (label, w, du, dv, dx, dy) in enumerate(LADDER_G_POINTS):
        slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
        rho = LADDER_G_EN / slab.coordinate_volume
        curved, flat = arms(slab, w)
        rng = np.random.default_rng([SEEDS["ladder_g"], k])
        rec = []
        for _ in range(N_G):
            pts = sprinkle(flat, rho, rng)
            npts = len(pts)
            pairs = npts * (npts - 1) // 2
            c0 = relation_census(flat, pts)
            cA = relation_census(curved, pts)
            g, lo, amb = pair_flips(c0, cA)
            rec.append((g / pairs, lo / pairs, amb / pairs,
                        cA.related.sum() / pairs,
                        c0.related.sum() / pairs))
        arr = np.asarray(rec)
        d = arr[:, 0] + arr[:, 1]
        fa, f0 = arr[:, 3], arr[:, 4]
        df = fa - f0
        sd_arm = math.sqrt(0.5 * (fa.var(ddof=1) + f0.var(ddof=1)))
        s_real = float(df.mean() / sd_arm) if sd_arm > 0 else 0.0
        # realized precision for s: the full-variance formula, never
        # the near-zero approximation (design review round 4)
        se_s = math.sqrt(2.0 / N_G + s_real ** 2 / (4.0 * (N_G - 1)))
        # descriptive AUC over per-sprinkling 3a values, stratified
        # bootstrap (per-arm independent resampling of sprinklings)
        auc = _mann_whitney_auc(fa, f0)
        brng = np.random.default_rng([SEEDS["ladder_g"], 500 + k])
        bs = np.empty(B_BOOT)
        for b in range(B_BOOT):
            bs[b] = _mann_whitney_auc(
                fa[brng.integers(0, N_G, N_G)],
                f0[brng.integers(0, N_G, N_G)])
        rows.append({
            "label": label, "w": w, "du": du, "dv": dv, "dx": dx,
            "dy": dy, "n": N_G,
            "raw": {"sum_g": float(arr[:, 0].sum()),
                    "sum_l": float(arr[:, 1].sum()),
                    "sum_amb": float(arr[:, 2].sum()),
                    # per-sprinkling 3a samples: the DISTRIBUTION the
                    # P3-C power model is refit on (fixed margins,
                    # recomputed n, committed before P3-C runs)
                    "f_curved": [float(v) for v in fa],
                    "f_flat": [float(v) for v in f0]},
            "mean_d": float(d.mean()),
            "rse_d": float(d.std(ddof=1) / math.sqrt(N_G) / d.mean())
            if d.mean() > 0 else None,
            "d_interval": [float(d.mean()),
                           float((d + arr[:, 2]).mean())],
            "mean_gained": float(arr[:, 0].mean()),
            "mean_lost": float(arr[:, 1].mean()),
            "mean_f_curved": float(fa.mean()), "mean_f_flat": float(f0.mean()),
            "delta_f": float(df.mean()),
            "se_delta_f": float(df.std(ddof=1) / math.sqrt(N_G)),
            "s_realized": s_real, "se_s_realized": float(se_s),
            "auc": float(auc),
            "auc_ci95_boot": [float(np.percentile(bs, 2.5)),
                              float(np.percentile(bs, 97.5))],
        })
        print(f"  ladder-G {label} done", flush=True)
    return rows


def _mann_whitney_auc(a: np.ndarray, b: np.ndarray) -> float:
    """AUC = P(a > b) + 0.5 P(a = b), midrank implementation."""

    both = np.concatenate([a, b])
    order = both.argsort(kind="mergesort")
    sorted_v = both[order]
    ranks_sorted = np.arange(1, len(both) + 1, dtype=float)
    i = 0
    while i < len(both):
        j = i
        while j + 1 < len(both) and sorted_v[j + 1] == sorted_v[i]:
            j += 1
        if j > i:
            ranks_sorted[i:j + 1] = 0.5 * (i + 1 + j + 1)
        i = j + 1
    ranks = np.empty_like(ranks_sorted)
    ranks[order] = ranks_sorted
    ra = ranks[:len(a)].sum()
    return (ra - len(a) * (len(a) + 1) / 2.0) / (len(a) * len(b))


# ---------------------------------------------------------------------
# Class C characterization (descriptive only)
# ---------------------------------------------------------------------


def _interval_chain_stats(census, elig_sorted: np.ndarray) -> dict:
    """One arm's causal intervals (eligible related endpoint pairs)
    with interior cardinalities, opportunity counts, and chain counts
    by matrix powers of the interior relation."""

    rel = census.related
    idx = np.where(elig_sorted)[0]
    m_list = []
    chains = {k: 0 for k in CHAIN_KS}
    opp = {k: 0 for k in CHAIN_KS}
    n_intervals = 0
    for ii in range(len(idx)):
        i = idx[ii]
        for jj in range(ii + 1, len(idx)):
            j = idx[jj]
            if not rel[i, j]:
                continue
            n_intervals += 1
            interior = np.where(rel[i, :] & rel[:, j])[0]
            m = len(interior)
            m_list.append(m)
            for k in CHAIN_KS:
                opp[k] += math.comb(m, k)
            if m >= 2:
                mm = rel[np.ix_(interior, interior)].astype(np.int64)
                p = mm
                chains[2] += int(p.sum())
                for k in (3, 4):
                    p = p @ mm
                    chains[k] += int(p.sum())
    m_arr = np.asarray(m_list) if m_list else np.zeros(0, dtype=int)
    return {"n_intervals": n_intervals,
            "m_profile": {
                "mean": float(m_arr.mean()) if len(m_arr) else 0.0,
                "max": int(m_arr.max()) if len(m_arr) else 0,
                "n_ge2": int((m_arr >= 2).sum()),
                "n_ge4": int((m_arr >= 4).sum()),
            },
            "chains": chains, "opportunities": opp}


def run_class_c() -> list[dict]:
    label, w, du, dv, dx, dy = CLASSC_POINT
    slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
    curved, flat = arms(slab, w)
    out = []
    for target_n in CLASSC_DENSITIES:
        rho = target_n / slab.coordinate_volume
        rng = np.random.default_rng([SEEDS["class_c"], target_n])
        per_arm = {"curved": [], "flat": []}
        seconds = []
        for _ in range(N_CLASSC):
            pts = sprinkle(curved, rho, rng)
            t0 = time.perf_counter()
            order = np.argsort(pts[:, 0], kind="stable")
            spts = pts[order]
            # ONE eligibility rule, the curved guard's, for both arms
            # (Class C eligibility is a property of the frozen
            # geometry, not of the arm being read -- P1 §4.6)
            elig = np.array([element_eligible(curved, p) for p in spts])
            for arm_name, geom in (("curved", curved), ("flat", flat)):
                cen = relation_census(geom, pts)
                per_arm[arm_name].append(
                    _interval_chain_stats(cen, elig))
            seconds.append(time.perf_counter() - t0)
        rec = {"label": label, "e_n": target_n, "n_sprinklings": N_CLASSC,
               "pipeline_seconds_median": float(np.median(seconds)),
               "pipeline_seconds_p95": float(np.percentile(seconds, 95)),
               "cpu": platform.processor()}
        for arm_name in ("curved", "flat"):
            rows = per_arm[arm_name]
            arm = {"mean_intervals": float(np.mean(
                [r["n_intervals"] for r in rows])),
                "m_profile": {
                    "mean_of_means": float(np.mean(
                        [r["m_profile"]["mean"] for r in rows])),
                    "max": int(max(r["m_profile"]["max"] for r in rows)),
                    "mean_n_ge2": float(np.mean(
                        [r["m_profile"]["n_ge2"] for r in rows])),
                    "mean_n_ge4": float(np.mean(
                        [r["m_profile"]["n_ge4"] for r in rows])),
                }}
            for k in CHAIN_KS:
                opp = np.array([r["opportunities"][k] for r in rows])
                ch = np.array([r["chains"][k] for r in rows])
                zero = int((opp == 0).sum())
                lo, hi = clopper_pearson(zero, N_CLASSC)
                arm[f"k{k}"] = {
                    "sum_opportunities": int(opp.sum()),
                    "sum_chains": int(ch.sum()),
                    "c_hat_pooled": (float(ch.sum() / opp.sum())
                                     if opp.sum() else None),
                    "zero_denominator_sprinklings": zero,
                    "zero_rate_ci95": [lo, hi],
                }
            rec[arm_name] = arm
        for k in CHAIN_KS:
            zc = np.array([r["opportunities"][k]
                           for r in per_arm["curved"]])
            zf = np.array([r["opportunities"][k] for r in per_arm["flat"]])
            union_zero = int(((zc == 0) | (zf == 0)).sum())
            lo, hi = clopper_pearson(union_zero, N_CLASSC)
            rec[f"union_zero_k{k}"] = {
                "sprinklings": union_zero, "ci95": [lo, hi]}
        out.append(rec)
        print(f"  class-C E[N]={target_n} done "
              f"(median {rec['pipeline_seconds_median']:.1f}s)", flush=True)
    return out


# ---------------------------------------------------------------------
# P3-C preliminary joint-power record (code + artifact; NO campaign)
# ---------------------------------------------------------------------


def p3c_joint_power(n_arm: int = P3C_N_CANDIDATE,
                    reps: int = P3C_POWER_REPS) -> dict:
    """Null-effect joint power of the composite gate at the FROZEN
    margins, normal design model, full pipeline per replicate
    (threshold and direction re-learned). Returns RAW integers with
    the exact CP interval -- the certification is the CP LOWER bound
    >= 0.90, never the point estimate."""

    z = 1.959964
    se_s = math.sqrt(2.0 / n_arm)
    se_auc = math.sqrt((2 * n_arm + 1) / (12.0 * n_arm * n_arm))
    se_ba = 1.0 / (2.0 * math.sqrt(n_arm))
    rng = np.random.default_rng(list(P3C_POWER_SEED))
    passes = 0
    done = 0
    while done < reps:
        k = min(200, reps - done)
        a = rng.normal(size=(k, n_arm))
        b = rng.normal(size=(k, n_arm))
        s = ((a.mean(1) - b.mean(1))
             / np.sqrt(0.5 * (a.var(1, ddof=1) + b.var(1, ddof=1))))
        both = np.concatenate([a, b], axis=1)
        ranks = both.argsort(axis=1).argsort(axis=1) + 1
        auc = ((ranks[:, :n_arm].sum(1) - n_arm * (n_arm + 1) / 2.0)
               / (n_arm * n_arm))
        h = n_arm // 2
        ta, tb = a[:, :h].mean(1), b[:, :h].mean(1)
        thr = 0.5 * (ta + tb)
        d = np.sign(ta - tb)
        acc_a = np.where(d > 0, (a[:, h:] > thr[:, None]).mean(1),
                         np.where(d < 0, (a[:, h:] < thr[:, None]).mean(1),
                                  0.0))
        acc_b = np.where(d > 0, (b[:, h:] <= thr[:, None]).mean(1),
                         np.where(d < 0, (b[:, h:] >= thr[:, None]).mean(1),
                                  1.0))
        ba = 0.5 * (acc_a + acc_b)
        ok = ((np.abs(s) + z * se_s < P3C_MARGINS["eps_s"])
              & (np.abs(auc - 0.5) + z * se_auc < P3C_MARGINS["eps_auc"])
              & (np.abs(ba - 0.5) + z * se_ba < P3C_MARGINS["eps_ba"]))
        passes += int(ok.sum())
        done += k
    lo, hi = clopper_pearson(passes, reps)
    return {
        "status": "PRELIMINARY candidate under the normal design "
                  "model; final n recomputed from the P3-E "
                  "distribution at the SAME margins and separately "
                  "approved before P3-C runs",
        "n_arm": n_arm, "reps": reps,
        "passes_raw": passes,
        "joint_power": passes / reps,
        "joint_power_ci95_exact": [lo, hi],
        "certified_lower_bound_rule": "CP 95% lower bound >= 0.90",
        "margins": dict(P3C_MARGINS),
        "seed_path": list(P3C_POWER_SEED),
    }


# ---------------------------------------------------------------------
# Artifact + rendered tables
# ---------------------------------------------------------------------


def write_artifact(art: dict) -> Path:
    with open(_ARTIFACT, "w", encoding="utf-8", newline="\n") as f:
        json.dump(art, f, indent=2, sort_keys=True)
        f.write("\n")
    return _ARTIFACT


def ladder_e_table(art: dict) -> str:
    fits = art["ladder_e_fits"]
    lines = ["| fit | range | slope | boot 95% CI | failures | "
             "drop-highest refit |",
             "|---|---|---|---|---|---|"]
    ranges = {"d_full": "A ≤ 4", "d_low": "A ≤ 0.5",
              "dr_asymptotic": "A ≤ 0.0625",
              "dr_finite_range": "A ≤ 0.5 (effective)"}
    for name, *_ in _FITS:
        f = fits[name]
        ci = (f"[{f['ci95'][0]:.3f}, {f['ci95'][1]:.3f}]"
              if f["ci95"] else "—")
        drop = (f"{f['refit_diagnostic_drop_highest']:.3f}"
                if f["refit_diagnostic_drop_highest"] is not None else "—")
        lines.append(
            f"| {name} | {ranges[name]} | {f['slope']:.3f} | {ci} "
            f"| {f['fit_failures']}/{f['boot_reps']} | {drop} |")
    return "\n".join(lines)


def ladder_g_table(art: dict) -> str:
    lines = ["| point | mean D | RSE | D interval | g/l | Δf̄ ± SE | "
             "s (realized ± SE) | AUC [boot 95%] |",
             "|---|---|---|---|---|---|---|---|"]
    for r in art["ladder_g"]:
        gl = (f"{r['mean_gained'] / max(r['mean_lost'], 1e-12):.1f}"
              if r["mean_lost"] > 0 else "inf")
        lines.append(
            f"| {r['label']} | {r['mean_d']:.5f} | "
            f"{100 * r['rse_d']:.1f}% | [{r['d_interval'][0]:.5f}, "
            f"{r['d_interval'][1]:.5f}] | {gl} | "
            f"{r['delta_f']:+.5f} ± {r['se_delta_f']:.5f} | "
            f"{r['s_realized']:.2f} ± {r['se_s_realized']:.2f} | "
            f"{r['auc']:.3f} [{r['auc_ci95_boot'][0]:.3f}, "
            f"{r['auc_ci95_boot'][1]:.3f}] |")
    return "\n".join(lines)


def class_c_table(art: dict) -> str:
    lines = ["| E[N] | arm | intervals/spr | k | Σopp | Σchains | Ĉ_k | "
             "zero-denom | union-zero |",
             "|---|---|---|---|---|---|---|---|---|"]
    for rec in art["class_c"]:
        for arm_name in ("curved", "flat"):
            arm = rec[arm_name]
            for k in CHAIN_KS:
                kk = arm[f"k{k}"]
                ch = (f"{kk['c_hat_pooled']:.4f}"
                      if kk["c_hat_pooled"] is not None else "—")
                uz = rec[f"union_zero_k{k}"]
                lines.append(
                    f"| {rec['e_n']} | {arm_name} | "
                    f"{arm['mean_intervals']:.1f} | {k} | "
                    f"{kk['sum_opportunities']} | {kk['sum_chains']} | "
                    f"{ch} | {kk['zero_denominator_sprinklings']}"
                    f"/{rec['n_sprinklings']} | "
                    f"{uz['sprinklings']}/{rec['n_sprinklings']} |")
    return "\n".join(lines)


def main() -> None:
    assert_seed_layout()
    print("P14 probe P3-E -- exploratory; seeds "
          f"{sorted(SEEDS.values())}\n")
    t0 = time.perf_counter()
    art: dict = {"script": "experiments/positive_control/p14_probe_p3e.py",
                 "seeds": dict(SEEDS),
                 "burned_seeds": list(BURNED_SEEDS)}

    print("ladder-E ...", flush=True)
    raw_e = run_ladder_e()
    art["ladder_e_raw"] = raw_e
    art["ladder_e_fits"] = ladder_e_fits(raw_e)
    print(f"  done ({time.perf_counter() - t0:.0f}s)")

    print("ladder-G ...", flush=True)
    art["ladder_g"] = run_ladder_g()

    print("class-C characterization ...", flush=True)
    art["class_c"] = run_class_c()

    print("P3-C preliminary joint power ...", flush=True)
    art["p3c_preliminary"] = p3c_joint_power()

    print("\n== ladder-E fits ==")
    print(ladder_e_table(art))
    print("\n== ladder-G ==")
    print(ladder_g_table(art))
    print("\n== class-C characterization (descriptive) ==")
    print(class_c_table(art))
    p = art["p3c_preliminary"]
    print(f"\np3c preliminary: n={p['n_arm']}/arm, "
          f"{p['passes_raw']}/{p['reps']} = {p['joint_power']:.4f}, "
          f"CP [{p['joint_power_ci95_exact'][0]:.4f}, "
          f"{p['joint_power_ci95_exact'][1]:.4f}]")
    print(f"\nartifact written: {write_artifact(art)} "
          f"({time.perf_counter() - t0:.0f}s total)")


if __name__ == "__main__":
    main()
