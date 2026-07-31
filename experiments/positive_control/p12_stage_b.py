"""P12 Stage B: recover the CURVATURE itself, from order plus count.

Implements the frozen addendum, docs/prereg/p12_curved_sprinkling.md
Section 10, and nothing else. Stage A asked whether the chain reading
follows curved proper time; Stage B asks whether the scalar curvature
`R` can be read back out of the causal order together with the interval
count. It is the first quantity in this programme that is not a length
or a time.

THE OBSERVABLE, order-plus-count only (Section 10.1):

    g_hat = m_open / (rho_hat * tau_hat^2)   ->   Vol/tau^2 = G(s)

with `s = tau / (2 ell)` and, exactly rather than by expansion,

    Vol = 4 ell^2 ln cosh(tau / (2 ell)),    G(s) = ln cosh(s) / s^2,

strictly decreasing from `G(0+) = 1/2`. So a flat twin normalizes the
discreteness bias away IN A RATIO and the inversion returns curvature:

    Q_hat = g_curved / g_flat,   s_hat = G^-1(Q_hat / 2),
    R_hat = 8 s_hat^2 / tau_hat^2.

WHY A RATIO AND NOT SECTION 5's DIFFERENCE. Each arm's `g` is biased by
53-111% because `tau_hat` underestimates `tau` and enters squared. A
difference leaves that bias in the units of the answer; the ratio
cancels it to first order in matched `m`. Section 10.3 has the measured
case.

WHAT IS IMPORTED RATHER THAN REIMPLEMENTED (Section 10.11). P13's
per-rung sprinkler, exact dS_2 truth, eligibility pair and greedy
disjoint packing for the 1.50 rung; P11's `score_pair`, calibrated
Bonett machinery, verdict table and preflight/gate helpers; the shared
chain estimator. The only new code here is `G`, its bisection inverse
with the boundary convention, `g_hat`, and the paired twin-ratio
aggregation.

WHY THIS IS A SEPARATE MODULE. Section 10.11 says to extend
`p12_curved.py` with a Stage B module; `p10_stage_b.py` is the house
precedent for that phrasing, and it keeps Stage A's code — which
produced a frozen IMPROVES record this stage's own gate reads — out of
the diff.

Run through `p12_curved.py --stage b`, or directly:

    python experiments/positive_control/p12_stage_b.py --stage verify-b
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import p13_tau_ell as P13
from p10_continuum_ladder import _stable_seed
from p11_metric import (
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
    verdict,
)
from pc_common import DEFAULT_OUTPUT_DIR, write_rows_csv

from causal_spacetime_lab.estimators import (
    estimate_tau_from_longest_chain_1p1,
)

ELL = 1.0

# --------------------------------------------------------------------
# Frozen constants (Section 10.8). The operating point is P13's 1.50
# rung VERBATIM, taken from P13's own tables so the two cannot drift
# apart silently: Section 10.8 moved Stage B off Stage A's band because
# at tau/ell = 0.3 the signal is 0.37% and the inversion amplifies
# every input error 269-fold, against 8.17% and 12.9x here.
# --------------------------------------------------------------------
TAU_ELL = 1.50
ETA_LO, XHALF, RHO_TOP = P13.PATCH[TAU_ELL]
BAND_REL = P13.BAND_REL
TWIN_CENTRE = P13.TWIN_BAND_CENTRE[TAU_ELL]

#: Truth. R = 2/ell^2 in dS_2; the DIMENSIONLESS combination R tau^2 is
#: what `Q_hat` alone determines, needing no length at all.
R_TRUE = 2.0 / ELL ** 2
R_TAU2_TRUE = R_TRUE * TAU_ELL ** 2

#: Density ladder at FIXED geometry, 4x end to end as Stage A's was.
P12B_LADDER = (600, 1200, 2400)
LADDER_RHO = {600: RHO_TOP / 4.0, 1200: RHO_TOP / 2.0, 2400: RHO_TOP}

#: The twin's intensities are FROZEN LITERALS from the design check's
#: calibration (Section 10.8), not recalibrated at run time. Two
#: reasons, both of which a re-calibrating runner would violate:
#: recalibration would consume experimental seeds to re-derive a
#: preregistered constant, and it would let the campaign tune its own
#: normalizer after seeing its own data. The cross-arm gate below is
#: what checks that the frozen values still land where they should.
TWIN_RHO = {600: 2.1948, 1200: 4.3370, 2400: 8.6990}

#: Per-rung realized `m` from the design check (Section 10.6 table:
#: 18.87 / 37.80 / 75.77), used as the WITHIN-ARM gate's targets.
M_TARGET = {600: 18.8681, 1200: 37.8014, 2400: 75.7675}
#: P13's Section 5 tolerance, reused for the within-arm gate.
M_TOLERANCE = P13.M_TOLERANCE
#: Section 10.6 (2): the twin-vs-curved `m` ratio, gated tightly
#: because `Q_hat` divides one arm by the other.
CROSS_ARM_TOLERANCE = 0.01

#: The offset the design check actually left in each rung, from
#: `p12_stage_b_ensemble_check.json` (`cross_arm_m[rung]`). This is the
#: campaign's EXPECTED cross-arm offset on correctly calibrated arms,
#: so it is what the spurious-trip disclosure has to be centred on --
#: review C28 found one literal `+0.55%`, the bottom rung's, standing in
#: for all three, which overstated the upper rungs' risk by 1.5x and
#: 50x. Literals here rather than a run-time read of a 1 MB artifact,
#: on VERIFY_PIN_B's precedent, and
#: test_frozen_cross_arm_offsets_match_the_design_check recomputes them
#: from that artifact rather than trusting these.
FROZEN_CROSS_ARM_OFFSET = {600: 0.00552373710017573,
                           1200: -0.004819043378003873,
                           2400: -0.0014848054904806895}
assert set(FROZEN_CROSS_ARM_OFFSET) == set(P12B_LADDER), (
    "every rung the campaign runs needs its own frozen offset -- this "
    "disclosure is assembled at the END of a full run, so a missing key "
    "would abort after the compute instead of before it")

# --------------------------------------------------------------------
# Power (Section 10.4 and 10.10). Delta*_B is DERIVED there, not
# borrowed from P11: g_hat inherits a Poisson 1/m from the count and a
# BDJ m^(-2/3) from the chain, and propagating the fitted rates over
# the frozen m ladder gives -0.2508 dex -- steeper than the chain-only
# -0.2007, which is the evidence it was derived.
#
# The 1.645 in the equivalence slot is Phi^-1(0.95), i.e. ALREADY the
# two-sided quantile for beta = 0.10, because a within-margin verdict
# needs BOTH interval bounds inside the margin and so has power
# 2 Phi(z) - 1. P13's review C2 caught exactly this slot being
# "upgraded" with a one-sided quantile; test_equivalence_coefficient_
# is_two_sided derives the delivered power back out of the constant.
# --------------------------------------------------------------------
DELTA_STAR_B = -0.2508
DELTA_EQ_B = 0.0836
N_SUP_COEFF = (1.960 + 1.282) ** 2 / DELTA_STAR_B ** 2
N_EQ_COEFF = (1.960 + 1.645) ** 2 / DELTA_EQ_B ** 2
N_FLOOR = 12
N_CAP = 1300

#: Section 10.7's co-requirement: the top rung's DIMENSIONLESS recovery.
RECOVERY_THRESHOLD = 0.25

#: The design check's fitted dispersion model (Section 10.4), kept so
#: test 5 can recompute Delta*_B from it rather than trusting the
#: literal above.
REL_SD_A = 1.4923
REL_SD_B = 0.4586

# --------------------------------------------------------------------
# Seed windows (Section 10.10). Section 6's 1200000-1999999 reservation
# is SUPERSEDED there: six blocks at these sizes need ~1.6M seeds
# against a reservation of 800k, and 2000000+ is P12's design-check
# space, which the addendum's own checks consumed. Slots are cap + 20
# at STRIDE. Design-check space from here on is 40000000+.
# --------------------------------------------------------------------
WINDOW_FLOOR, WINDOW_CEIL = 30_000_000, 33_263_999
VERIFY_B_BASE = {600: 30_000_000, 1200: 30_002_000, 2400: 30_004_000}
PILOT_B_BLOCKS = {600: (30_100_000, 320), 2400: (30_200_000, 320)}
PILOT_B_TWIN_BLOCKS = {600: (30_300_000, 320), 2400: (30_400_000, 320)}
STAGE_B_BLOCKS = {600: (31_000_000, 1320), 1200: (31_400_000, 1320),
                  2400: (31_800_000, 1320)}
STAGE_B_TWIN_BLOCKS = {600: (32_200_000, 1320), 1200: (32_600_000, 1320),
                       2400: (33_000_000, 1320)}
SCORE_OFFSET = P13.SCORE_OFFSET

#: Every range this programme has already spent, for the freshness test
#: (Section 10.11 item 7 asks it be asserted against the RANGES rather
#: than against a floor, so a future window that lands *below* a spent
#: range is caught too).
SPENT_RANGES = (
    ("P11/P12 experimental", 1_000_000, 1_999_999),
    ("P12 design check", 2_000_000, 5_999_999),
    ("P13 campaign v1", 6_000_000, 6_475_999),
    ("P13 campaign v2", 8_000_000, 8_795_999),
    ("P13 11.3 diagnostic", 14_000_000, 14_549_999),
    ("P13 campaign v3", 15_000_000, 21_503_999),
    ("P13 design-check reservation", 22_000_000, 29_999_999),
)

# --------------------------------------------------------------------
# Verification-B's completeness pin, PER RUNG. Section 10.10 first wrote
# P11's single `VERIFY_PIN = 1998`; that is amended, and the amendment is
# disclosed rather than applied quietly because a preregistered threshold
# is being moved.
#
# WHAT WAS WRONG. 1998 of 2000 demands completion 0.999, and 10.5's own
# frozen artifact had already measured the bottom rung at 0.997 (3988 of
# 4000). Worse, 1998 sits ON the median of Bin(2000, 0.999): at exactly
# the rate it demands it fails half the time, and it needs better than
# 0.9995 to be reliably clearable. No rung was ever measured there. So
# 10.5 and 10.10 contradicted each other inside one section -- 10.5 says
# "12 samples in 4000 fail to pack and the frozen fill rule replaces them
# from reserve slots, well inside the skip cap", and even notes "the true
# rate is not zero", while 10.10 required it to be nearly zero.
#
# The fragility was not new and had already been grazed at THIS operating
# point: P13 campaign v1's tau/ell = 1.5 verification read 1998 of 2000,
# clearing the pin by exactly zero, and v3's 1.0 rung by one.
#
# WHAT THE PIN IS FOR, which is what fixes the calibration. It is a smoke
# test run before any quarantined window is spent. It is NOT the
# feasibility gate: the fill rule is, and the fill rule is exact -- 1320
# slots, SKIP_CAP 20, needing 1300, aborting INFEASIBLE-INCOMPLETE on the
# real campaign. At n = 2000 the healthy rate and the rate at which the
# campaign starts aborting (~0.990) are 2.68 sigma apart, so NO pin on
# 2000 samples can be both clearable and a feasibility gate. 1998 was
# calibrated as though it were one.
#
# THE RULE, and it may only ever loosen:
#
#     pin(rung) = min( 1998, largest k with
#                      P(X >= k | X ~ Bin(2000, p_lo)) >= 0.999 )
#
# with p_lo the one-sided 95% Clopper-Pearson lower bound on that rung's
# frozen CURVED completion (verification runs the curved arm; the twin's
# frozen rate is at least as high at every rung, so the curved arm binds).
# Two levels, both conventional rather than tuned: 95% on the rate so the
# pin survives the frozen estimate sitting at the pessimistic end of its
# own sample, and 99.9% clearance so a healthy pipeline fails at most once
# in a thousand runs.
#
# The inputs are frozen data ONLY. 3988/4000, 2000/2000 and 1000/1000
# predate any Stage B run, so the pins below are not tuned to a measured
# verification -- and test_verification_pins_are_derived_from_frozen_data
# recomputes them from the artifact rather than trusting these literals,
# which is review C25's lesson applied before it can recur.
#
# WHAT IT COSTS, stated because it is not nothing. Against a drift to
# 0.990-0.995 the derived pin has little power where 1998 had nearly all;
# against a collapse to 0.985 or below, which is where the campaign
# actually aborts, it retains 0.95-1.00. The sensitivity given up is to
# drifts that do not threaten the campaign, and the fill rule catches
# those exactly at fill time.
VERIFY_PIN_B = {600: 1979, 1200: 1990, 2400: 1985}
assert all(pin <= VERIFY_PIN for pin in VERIFY_PIN_B.values()), (
    "the pin derivation may only LOOSEN the inherited 1998, never tighten "
    "it -- otherwise it could be used to make a later stage easier to "
    "reach")

VERIFY_ARTIFACT = "p12_verification_b_summary.json"
PILOT_ARTIFACT = "p12_pilot_b_summary.json"
STAGE_ARTIFACT = "p12_stage_b_summary.json"

#: The two frozen records Section 10.10's cross-stage gate reads. P13's
#: CURVATURE-ROBUST is campaign v3 ("Stage A-13C"), which lives in its
#: own directory: `frozen/p13/` holds campaign v1's CONFOUNDED record,
#: so naming the directory matters and is asserted by a test.
FROZEN_P12_DIR = (Path(__file__).resolve().parents[2]
                  / "docs" / "prereg" / "frozen" / "p12")
FROZEN_P13V3_DIR = (Path(__file__).resolve().parents[2]
                    / "docs" / "prereg" / "frozen" / "p13v3")


# ====================================================================
# The new mathematics: G, its inverse, and the boundary convention
# ====================================================================

_LOG_2 = math.log(2.0)

#: Above this the asymptote below is exact to the last bit, and well
#: below the `|x| ~ 710` where `math.cosh` overflows.
_LOG_COSH_ASYMPTOTIC_FROM = 20.0


def log_cosh(x: float) -> float:
    """`ln cosh(x)`, finite wherever `x` is.

    Review C29: `math.cosh` RAISES above `|x| ~ 710`, so the direct
    quotient did not merely lose accuracy for strongly curved readings,
    it aborted on them -- and `invert_g`'s bracket doubles through
    `s = 1024` on any `g` below about `0.002`, so the documented domain
    `0 < g < 1/2` was not covered and the `hi` ceiling that was supposed
    to catch it could never run. Both ends need their own form: the
    direct one keeps the significant figures the `s -> 0` limit needs,
    and `ln((1 + e^-2|x|)/2)` carries the large-`x` tail without ever
    forming `cosh` itself.
    """

    ax = abs(x)
    if ax < _LOG_COSH_ASYMPTOTIC_FROM:
        return math.log(math.cosh(ax))
    return ax + math.log1p(math.exp(-2.0 * ax)) - _LOG_2


def g_of_s(s: float) -> float:
    """`G(s) = ln cosh(s) / s^2`, with the `s -> 0` limit filled in.

    Strictly decreasing on `(0, inf)` from `G(0+) = 1/2`, which is what
    makes the inversion single-valued; the series `1/2 - s^2/12` is used
    near zero because the quotient loses all its significant figures
    there.
    """

    if s < 1e-8:
        return 0.5 - s ** 2 / 12.0
    return log_cosh(s) / s ** 2


#: The bracket may not double past this. `G(s) ~ 1/s` out here, and the
#: doubling stops at the last power of two under the ceiling, so the
#: reachable domain is `g >= G(2^39) = 1.8e-12` -- nine orders below
#: the smallest `Q_hat/2` these ladders can produce. It bounds the
#: loop; it is NOT a second boundary convention, which is why hitting
#: it raises rather than returning 10.2's None (C29).
S_MAX = 1e12


def invert_g(g: float) -> float | None:
    """`s = G^-1(g)` by bisection, or None at the boundary.

    THE BOUNDARY CONVENTION (Section 10.2), frozen. The
    positive-curvature domain is `g < 1/2`. A measurement at or above
    the flat value has recovered no curvature, so it returns None and
    the caller records `R_hat = 0` with relative error exactly 1. That
    is the EDGE OF THE PARAMETER SPACE, not a clamp over noise, and the
    boundary-hit fraction is published per rung rather than absorbed.

    10.2 gives that None exactly ONE meaning, so the other end of the
    domain may not borrow it: a `g` too small for `S_MAX` is the most
    strongly curved reading the instrument can report, and returning
    None would file it as `R_hat = 0` -- no curvature at all. It raises
    instead.
    """

    if not (0.0 < g < 0.5):
        return None
    lo, hi = 1e-8, 1.0
    while g_of_s(hi) > g:
        hi *= 2.0
        if hi > S_MAX:
            raise ValueError(
                f"G^-1({g!r}) lies past s = {S_MAX:g}. That is the "
                "strongly curved end, NOT 10.2's flat boundary, and the "
                "two may not share a return value.")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g_of_s(mid) > g:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def volume_closed_form(tau: float, ell: float = ELL) -> float:
    """`4 ell^2 ln cosh(tau / 2 ell)` -- the exact dS_2 diamond volume.

    Section 5 refused the small-diamond expansion because its
    coefficient is "the kind of unverified constant that killed Stage C
    v1". This is the closed form the volume check pinned against
    quadrature at 1e-14; position drops out entirely, as maximal
    symmetry requires.
    """

    return 4.0 * ell ** 2 * log_cosh(tau / (2.0 * ell))


def g_hat(m_open: int, rho_hat: float, tau_hat: float) -> float:
    """The per-pair observable: order data over count-calibrated area.

    Nothing here knows a length: `m_open` is the open interval
    cardinality, `rho_hat` is the realized count over the patch's
    proper volume, and `tau_hat` is the shared chain estimator's
    reading. That is the whole of "order plus number".
    """

    return m_open / (rho_hat * tau_hat ** 2)


# ====================================================================
# Sampling: P13's code path, returning per-pair g rather than P13's y
# ====================================================================

def run_sample(rung: int, seed: int, flat: bool = False,
               single_stream: bool = False):
    """One Stage B sample at a density rung. Returns (record, complete).

    The geometry is fixed and only `rho` moves down the ladder, which is
    P12's sweep; `flat` runs the twin, whose band centre and intensity
    are its own frozen constants.
    """

    rho = TWIN_RHO[rung] if flat else LADDER_RHO[rung]
    centre = TWIN_CENTRE if flat else TAU_ELL
    band = (centre * (1.0 - BAND_REL), centre * (1.0 + BAND_REL))
    rng = np.random.default_rng(seed)
    eta, x = P13.sprinkle(ETA_LO, XHALF, rho, rng, flat=flat)
    # Section 3's convention: realized count over the patch's PROPER
    # volume. For the flat twin that proper volume IS the coordinate
    # area -- using the curved patch's volume there mis-normalizes
    # tau_hat, which P13's first implementation did.
    volume = (P13.patch_proper_volume(ETA_LO, XHALF) if not flat
              else (abs(ETA_LO) - 1.0) * 2.0 * XHALF)
    rho_hat = eta.size / volume
    pool = P13.eligible_pairs(eta, x, ETA_LO, XHALF, band, flat=flat)
    draw_rng = rng if single_stream else np.random.default_rng(
        seed + SCORE_OFFSET
    )
    u, v = eta + x, eta - x
    pairs, complete, rejections = P13.draw_disjoint(pool, u, v, draw_rng)
    record = {
        "rung": float(rung), "seed": float(seed), "flat_twin": bool(flat),
        "rho": float(rho), "n_realized": float(eta.size),
        "rho_hat": float(rho_hat), "pool_size": float(pool.shape[0]),
        "rejections": float(rejections), "complete": bool(complete),
    }
    if not complete:
        # completeness is decided before any estimator is evaluated --
        # the frozen order of operations (P11 Section 4)
        return record, False
    gs, taus, ms = [], [], []
    for i, j in pairs:
        scored = score_pair(u / 2.0, v / 2.0, i, j)
        tau_hat = float(estimate_tau_from_longest_chain_1p1(
            scored["chain_length"], rho=rho_hat
        ))
        if tau_hat <= 0.0:
            # a zero-length chain carries no reading; g would be
            # infinite. Counted, not silently dropped.
            continue
        gs.append(g_hat(int(scored["m_open"]), rho_hat, tau_hat))
        taus.append(tau_hat)
        ms.append(int(scored["m_open"]))
    if not gs:
        record["complete"] = False
        return record, False
    record["mean_g"] = float(np.mean(gs))
    record["mean_tau_hat"] = float(np.mean(taus))
    record["mean_m"] = float(np.mean(ms))
    record["sd_m"] = float(np.std(ms, ddof=1)) if len(ms) > 1 else 0.0
    record["pairs_scored"] = float(len(gs))
    return record, True


def _fill_block(rung: int, base: int, slots: int, needed: int,
                flat: bool = False):
    """Seeds in window order, first `needed` complete (Section 6)."""

    records, skipped = [], []
    for k in range(slots):
        seed = base + STRIDE * k
        record, complete = run_sample(rung, seed, flat=flat)
        if complete:
            records.append(record)
            if len(records) == needed:
                return records, skipped, True
        else:
            skipped.append(seed)
            if len(skipped) > SKIP_CAP:
                return records, skipped, False
    return records, skipped, False


# ====================================================================
# The paired twin-ratio aggregation
# ====================================================================

#: A relative error of exactly zero would take log10 to -inf; the
#: design check's floor is kept so the two agree. Measure-zero for
#: continuous g, and recorded rather than hidden.
Y_FLOOR = 1e-12


def paired_recovery(curved: list[dict], twin: list[dict]) -> dict:
    """`y_B = log10(|R_hat - R| / R)`, PAIRED sample by sample.

    Review C10 (on the design check) is why this is paired rather than
    divided by the twin's rung mean. With a single estimated
    denominator every `y_B` shares one random quantity, so the stored
    variance contains none of the twin arm's sampling variation while
    the power calculation treats the values as independent -- which
    makes the projected `n` too optimistic. Pairing by index is
    arbitrary between two independent seed streams, and therefore
    harmless, and it makes each `y_B` a genuine iid draw carrying BOTH
    arms' variation. Nothing rung-level leaks in: `tau_hat` is the
    sample's own.

    Review C13 is why the boundary population is here rather than
    excluded: every summary below is over the same boundary-constrained
    population, where the first version quoted medians that excluded
    undefined inversions while the power block included them at error 1
    -- two summaries of one population disagreeing inside one script.
    """

    n = min(len(curved), len(twin))
    rels, boundary = [], 0
    for i in range(n):
        q = curved[i]["mean_g"] / twin[i]["mean_g"]
        s = invert_g(q / 2.0)
        if s is None:
            boundary += 1
            rels.append(1.0)
            continue
        r_hat = 8.0 * s ** 2 / curved[i]["mean_tau_hat"] ** 2
        rels.append(abs(r_hat - R_TRUE) / R_TRUE)
    rel = np.array(rels, dtype=float)
    y = np.log10(np.where(rel > 0.0, rel, Y_FLOOR))
    return {
        "y": y, "rel_error": rel, "n": n,
        "boundary_hits": boundary,
        "boundary_fraction": boundary / n if n else float("nan"),
        "median_rel_error": float(np.median(rel)) if n else float("nan"),
        "count_below_one": int((rel < 1.0).sum()),
        "fraction_below_one": float((rel < 1.0).mean()) if n else
        float("nan"),
    }


def dimensionless_recovery(curved: list[dict], twin: list[dict],
                           label: str) -> dict:
    """`R tau^2 = 8 s^2` from `Q_hat` alone -- Section 10.7's gate (ii).

    This needs NO length, so it separates what the ratio recovers from
    what normalizing by `tau_hat^2` spoils. The interval is a JOINT
    bootstrap over both arms, because resampling only the curved arm
    would understate it exactly as a fixed denominator did (C10).
    """

    cg = np.array([r["mean_g"] for r in curved], dtype=float)
    tg = np.array([r["mean_g"] for r in twin], dtype=float)
    ct = np.array([r["mean_tau_hat"] for r in curved], dtype=float)

    def recover(q: float) -> float:
        s = invert_g(q / 2.0)
        return 0.0 if s is None else 8.0 * s ** 2

    q_hat = float(cg.mean() / tg.mean())
    r_tau2 = recover(q_hat)
    rng = np.random.default_rng(_stable_seed(f"p12-b-recovery-{label}"))
    boots = []
    for _ in range(4000):
        ic = rng.integers(0, cg.size, cg.size)
        it = rng.integers(0, tg.size, tg.size)
        boots.append(recover(float(cg[ic].mean() / tg[it].mean())))
    boots = np.array(boots, dtype=float)
    rel_boots = np.abs(boots - R_TAU2_TRUE) / R_TAU2_TRUE
    return {
        "Q_hat": q_hat, "one_minus_Q": 1.0 - q_hat,
        "R_tau2_hat": r_tau2, "R_tau2_true": R_TAU2_TRUE,
        "rel_error": abs(r_tau2 - R_TAU2_TRUE) / R_TAU2_TRUE,
        "rel_error_ci": [float(np.percentile(rel_boots, 2.5)),
                         float(np.percentile(rel_boots, 97.5))],
        "R_tau2_ci": [float(np.percentile(boots, 2.5)),
                      float(np.percentile(boots, 97.5))],
        "tau_hat_over_tau_true": float(ct.mean() / TAU_ELL),
    }


# ====================================================================
# Gates
# ====================================================================

def m_gate(records_by_rung: dict, label: str) -> dict:
    """The WITHIN-ARM `m` gate, per rung against its frozen target.

    P13's form of this gate does not transfer literally and the reason
    is worth stating rather than quietly adapting. P13 swept `tau/ell`
    at FIXED discreteness, so it could require every rung within 5% of
    the arm's grand mean. Stage B sweeps DENSITY, so `m` runs 19 -> 76
    by design and a grand-mean test would fail on the design itself.
    The transferable content of the gate is what it protects against --
    an arm drifting away from the discreteness the design was priced at
    -- so each rung is checked against ITS OWN frozen target.

    Applied to BOTH arms. That is P13's review C5 ported: Section 5
    gives the twin the same eligibility, K, cap, fill rule and m-gate,
    and P13's first implementation gated only the curved rows, so a
    twin rung could drift in discreteness and still be used as the
    normalizer. Here the twin IS inside the estimator, so an ungated
    twin is worse than decoration -- it is a silent bias.
    """

    per_rung = {}
    for rung, records in records_by_rung.items():
        realized = float(np.mean([r["mean_m"] for r in records]))
        target = M_TARGET[rung]
        offset = realized / target - 1.0
        per_rung[str(rung)] = {
            "mean_m": realized, "target_m": target, "offset": offset,
            "passed": bool(abs(offset) <= M_TOLERANCE),
        }
    return {
        "arm": label, "tolerance": M_TOLERANCE, "per_rung": per_rung,
        "passed": all(r["passed"] for r in per_rung.values()),
    }


def cross_arm_m_gate(curved_by_rung: dict, twin_by_rung: dict,
                     enforce: bool = True) -> dict:
    """The CROSS-ARM `m` gate at +/- 1% (Section 10.6, review C11).

    `Q_hat` divides one arm by the other, so a consistent offset
    between the arms passes both within-arm gates while numerator and
    denominator keep different discreteness biases. Review C11 measured
    what that costs: the twin sat 2-3.4% below the curved arm, giving a
    ~1.5% bias mismatch, and with `Q/(1-Q) ~ 11` that is ~17% injected
    straight into the recovery. This is the case a within-arm gate
    cannot see.

    The bias mismatch is published alongside, read off each arm's own
    `g` against its own continuum value, because it is the systematic
    Section 10.6 prices the whole recovery with.

    WHY `enforce` EXISTS, and it is 10.6's own principle applied to this
    gate rather than only to the calibration probe. 10.6 says of the twin
    calibration: "the probe must be able to SEE the tolerance it
    calibrates to ... calibrating below the probe's own resolution is not
    calibration". The same holds for gating. At Stage P-B's 200 samples
    the offset's own standard error is about 0.94%, i.e. the tolerance
    itself, so a 1% gate there would abort roughly a third of the time on
    perfectly calibrated arms -- it would be noise with a verdict
    attached. So the pilot MEASURES the offset and publishes it with its
    standard error, and Stage B, at about 1300 samples and a standard
    error near 0.37%, is where the gate binds. `offset_se` and
    `resolution_ratio` are recorded either way so the gate's
    meaningfulness is auditable from the artifact instead of argued in a
    docstring.
    """

    g_curved_true = g_of_s(TAU_ELL / 2.0)
    per_rung = {}
    for rung in curved_by_rung:
        c, t = curved_by_rung[rung], twin_by_rung[rung]
        m_c_all = np.array([r["mean_m"] for r in c], dtype=float)
        m_t_all = np.array([r["mean_m"] for r in t], dtype=float)
        m_c, m_t = float(m_c_all.mean()), float(m_t_all.mean())
        se_c = float(m_c_all.std(ddof=1) / math.sqrt(m_c_all.size))
        se_t = float(m_t_all.std(ddof=1) / math.sqrt(m_t_all.size))
        g_c = float(np.mean([r["mean_g"] for r in c]))
        g_t = float(np.mean([r["mean_g"] for r in t]))
        ratio = m_t / m_c
        offset = ratio - 1.0
        # arms run on disjoint seed blocks, so independent propagation
        offset_se = ratio * math.hypot(se_t / m_t, se_c / m_c)
        bias_c, bias_t = g_c / g_curved_true, g_t / 0.5
        per_rung[str(rung)] = {
            "mean_m_curved": m_c, "mean_m_twin": m_t,
            "se_mean_m_curved": se_c, "se_mean_m_twin": se_t,
            "cross_arm_offset": offset,
            "cross_arm_offset_se": offset_se,
            "resolution_ratio": (CROSS_ARM_TOLERANCE / offset_se
                                 if offset_se > 0 else None),
            "bias_factor_curved": bias_c, "bias_factor_twin": bias_t,
            "bias_mismatch": bias_t / bias_c - 1.0,
            "passed": bool(abs(offset) <= CROSS_ARM_TOLERANCE),
        }
    return {
        "name": "cross-arm m",
        "tolerance": CROSS_ARM_TOLERANCE,
        "enforced": bool(enforce),
        "per_rung": per_rung,
        # when the gate is not enforced its verdict must not be read as
        # one, so `passed` is reported and the abort list ignores it
        "passed": (all(r["passed"] for r in per_rung.values())
                   if enforce else True),
        "measured_passed": all(r["passed"] for r in per_rung.values()),
        "note": (
            "ENFORCED: a failure aborts before anything is written."
            if enforce else
            "MEASURED ONLY at this stage: the offset's standard error is "
            "comparable to the tolerance at pilot sample sizes, so a "
            "gate here would be noise with a verdict attached (10.6's "
            "resolution principle). Stage B enforces it."),
    }


def _resolution_ratio_text(ratio: float | None) -> str:
    """Render `resolution_ratio`, which the gate types `float | None`.

    Found by the first test that runs `run_stage_b` end to end (C27):
    both runners printed this field as a bare float while the line
    directly beneath tested it for None, so a zero standard error --
    perfect resolution, the harmless end -- crashed the campaign in its
    own progress log. One renderer rather than the same special case
    written twice.
    """

    return "n/a (se = 0)" if ratio is None else f"{ratio:.2f}"


def cross_arm_spurious_trip_probability(offset: float, se: float) -> float:
    """P(the cross-arm gate trips) on arms calibrated to `offset`.

    Two-tailed, because the gate reads `|realized| > tolerance`. The
    realized offset is centred on what the design check left in THAT
    rung; review C28 found the bottom rung's `+0.55%` centring all
    three, which is the one substitution that can only ever inflate the
    upper rungs -- their frozen offsets are negative and smaller.
    """

    scale = se * math.sqrt(2.0)
    return (0.5 * (1.0 - math.erf((CROSS_ARM_TOLERANCE - offset) / scale))
            + 0.5 * (1.0 + math.erf((-CROSS_ARM_TOLERANCE - offset) / scale)))


def _abort_on_gate_failure(gates: list[dict]) -> None:
    """Gates ABORT; they do not record their own failure and continue.

    Review C21 found the first fix for C11 doing exactly that, which is
    review C12's shape reappearing inside the repair for C11: a guard
    that writes `passed: false` into an artifact and then writes the
    artifact is a record of the damage, not a guard. Section 10.11's
    test 4b requires asserting the ABORT rather than the flag.
    """

    failures = []
    for gate in gates:
        if gate.get("passed"):
            continue
        name = gate.get("arm") or gate.get("name") or "gate"
        for rung, row in gate["per_rung"].items():
            if not row["passed"]:
                key = ("cross_arm_offset" if "cross_arm_offset" in row
                       else "offset")
                failures.append(
                    f"{name} rung {rung}: {row[key]:+.2%} outside "
                    f"{gate['tolerance']:.0%}")
    if failures:
        raise SystemExit(
            "Stage B ABORTS on its frozen m gates (nothing written): "
            + "; ".join(failures))


def power_requirements_b(y_bottom: np.ndarray,
                         y_top: np.ndarray) -> dict:
    """Section 10.10's frozen formulas, with 10.4's `Delta*_B`.

    Calibrated Bonett bounds per rung, summed (Bonferroni at the
    nominal levels -- calibrated-approximate, not exact), then the
    frozen selection rule against the cap of 1300. The
    FLAT-WITHIN-MARGIN row is available only if `n_eq` fits.
    """

    bottom = calibrated_variance_bound(y_bottom, "stage-b-bottom")
    top = calibrated_variance_bound(y_top, "stage-b-top")
    s2 = bottom["s2"] + top["s2"]
    s2_90 = bottom["bound"] + top["bound"]
    n_sup = int(math.ceil(N_SUP_COEFF * s2_90))
    n_eq = int(math.ceil(N_EQ_COEFF * s2_90))
    calibration = {
        side: {k: data[k] for k in
               ("coverage_at_nominal", "z_used", "calibrated")}
        for side, data in (("bottom", bottom), ("top", top))
    }

    def clamp(k: int) -> int:
        return int(min(max(k, N_FLOOR), N_CAP))

    common = {"s2": s2, "s2_90": s2_90, "calibration": calibration,
              "n_sup": n_sup, "n_eq": n_eq,
              "delta_star_b": DELTA_STAR_B, "delta_eq_b": DELTA_EQ_B,
              "n_cap": N_CAP}
    if n_sup > N_CAP:
        return {**common, "n_per_rung": None, "flat_available": False,
                "infeasible": True}
    if n_eq <= N_CAP:
        return {**common, "n_per_rung": clamp(max(n_sup, n_eq)),
                "flat_available": True, "infeasible": False}
    return {**common, "n_per_rung": clamp(n_sup),
            "flat_available": False, "infeasible": False}


def seed_blocks() -> list[tuple[str, int, int]]:
    """Every block this stage consumes, as (label, lo, hi) inclusive.

    Review C17 (on the design check) is why these are enumerated rather
    than argued: the previous layout put one rung's calibration
    iteration on exactly another rung's, and 14 cross-rung collisions
    followed behind a comment claiming separation. The eye is not a
    disjointness checker.
    """

    blocks: list[tuple[str, int, int]] = []
    for rung, base in VERIFY_B_BASE.items():
        # verification is single-stream over CONSECUTIVE seeds
        blocks.append((f"verify-{rung}", base, base + VERIFY_COUNT - 1))
    for name, table in (("pilot-curved", PILOT_B_BLOCKS),
                        ("pilot-twin", PILOT_B_TWIN_BLOCKS),
                        ("stage-curved", STAGE_B_BLOCKS),
                        ("stage-twin", STAGE_B_TWIN_BLOCKS)):
        for rung, (base, slots) in table.items():
            blocks.append((f"{name}-{rung}", base,
                           base + STRIDE * (slots - 1)))
    return blocks


def assert_windows_disjoint_and_fresh() -> list[tuple[str, int, int]]:
    """Pairwise disjointness, containment, and freshness. Aborts."""

    blocks = seed_blocks()
    for label, lo, hi in blocks:
        if not (WINDOW_FLOOR <= lo and hi <= WINDOW_CEIL):
            raise SystemExit(
                f"seed block {label} [{lo}..{hi}] leaves Stage B's "
                f"window [{WINDOW_FLOOR}..{WINDOW_CEIL}]")
        for name, slo, shi in SPENT_RANGES:
            if not (hi < slo or shi < lo):
                raise SystemExit(
                    f"seed block {label} [{lo}..{hi}] overlaps the "
                    f"already-spent {name} [{slo}..{shi}]")
    for i, (la, loa, hia) in enumerate(blocks):
        for lb, lob, hib in blocks[i + 1:]:
            if not (hia < lob or hib < loa):
                raise SystemExit(
                    f"seed blocks overlap: {la} [{loa}..{hia}] vs "
                    f"{lb} [{lob}..{hib}]")
    return blocks


def cross_stage_gate() -> dict:
    """Section 10.10: BOTH prerequisites, with reachable stamps.

    P12 Stage A's IMPROVES because the estimator is that one, and
    P13's CURVATURE-ROBUST because the OPERATING POINT is that one --
    Stage B does not run at Stage A's band, so it cannot inherit Stage
    A's certification of its own rung and inherits P13's instead
    (Section 10.8).

    P13's record is campaign v3's, in `frozen/p13v3/`. `frozen/p13/`
    holds campaign v1's CONFOUNDED record, so the directory is named
    explicitly and a test pins it.
    """

    stage_a = _require_stage_pass(
        "p12_stage_a_summary.json", "P12 Stage A",
        directory=FROZEN_P12_DIR, expected="IMPROVES")
    p13 = _require_stage_pass(
        "p13_stage_a_summary.json", "P13 Stage A-13C",
        directory=FROZEN_P13V3_DIR, expected="CURVATURE-ROBUST")
    return {
        "p12_stage_a": {"verdict": stage_a["verdict"],
                        "code_version": stage_a["code_version"]},
        "p13_stage_a_13c": {"verdict": p13["verdict"],
                            "code_version": p13["code_version"],
                            "control_result": p13.get("control_result"),
                            # P13 v3's own record carries a selection
                            # caveat; an inherited gate inherits its
                            # caveats too, so it rides along here
                            # rather than being dropped at the border.
                            "selection_caveat":
                                bool(p13.get("selection_caveat"))},
    }


# ====================================================================
# Stage runners
# ====================================================================

def run_verify_b(output_dir: Path) -> None:
    """Verification-B (Section 10.10): completeness pin and wall times.

    2000 single-stream samples per rung against the PER-RUNG pins
    derived above. Outcomes are computed to price the pipeline and then
    discarded; only completeness and time are recorded.
    """

    stamp = _preflight_clean()
    assert_windows_disjoint_and_fresh()
    summary: dict = {
        "code_version": stamp, "tau_over_ell": TAU_ELL,
        "pin_required": {str(r): VERIFY_PIN_B[r] for r in P12B_LADDER},
        # the inherited value is recorded beside the derived ones so the
        # amendment is legible in the artifact, not only in the prose
        "pin_inherited_single_value": VERIFY_PIN,
        "pin_note": (
            "PER-RUNG, amended: the inherited 1998 demanded completion "
            "0.999 while 10.5's frozen artifact measured 0.997 at the "
            "bottom rung, and 1998 sits on the median of Bin(2000, "
            "0.999). Derived from frozen completions only; the pin is a "
            "smoke test and the fill rule is the feasibility gate."),
        "per_rung": {},
    }
    for rung in P12B_LADDER:
        base = VERIFY_B_BASE[rung]
        done, start = 0, time.perf_counter()
        for k in range(VERIFY_COUNT):
            _r, ok = run_sample(rung, base + k, single_stream=True)
            done += int(ok)
        elapsed = time.perf_counter() - start
        summary["per_rung"][str(rung)] = {
            "complete": done, "total": VERIFY_COUNT,
            "pin": VERIFY_PIN_B[rung],
            "pin_margin": done - VERIFY_PIN_B[rung],
            "completion": done / VERIFY_COUNT,
            "mean_seconds_per_sample": elapsed / VERIFY_COUNT,
        }
        print(f"verify-12B rung={rung}: {done}/{VERIFY_COUNT} complete "
              f"(pin {VERIFY_PIN_B[rung]}, margin "
              f"{done - VERIFY_PIN_B[rung]:+d}) "
              f"| {elapsed / VERIFY_COUNT:.3f} s/sample", flush=True)
    summary["pin_passed"] = bool(all(
        summary["per_rung"][str(r)]["complete"] >= VERIFY_PIN_B[r]
        for r in P12B_LADDER
    ))
    (output_dir / VERIFY_ARTIFACT).write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pin_passed": summary["pin_passed"]}, indent=2))


def run_pilot_b(output_dir: Path) -> None:
    """Stage P-B (Section 10.10): both endpoint rungs, both arms.

    200 samples each. Cross-rung statistics are FORBIDDEN: nothing here
    computes a mean difference, and the artifact holds per-rung bounds
    and times only. The twin is piloted at the same `n` because
    Section 10.6 measures its per-pair dispersion equal to the curved
    arm's, so equal `n` equalizes the two contributions to `se(Q_hat)`.
    """

    stamp = _preflight_clean()
    assert_windows_disjoint_and_fresh()
    gate = cross_stage_gate()
    verification = _load_gate_artifact(output_dir, VERIFY_ARTIFACT, stamp)
    if not verification.get("pin_passed"):
        raise SystemExit("verification-12B pin failed -- pilot refuses")

    summary: dict = {"code_version": stamp, "cross_stage_gate": gate,
                     "per_rung": {}}
    ys: dict = {}
    curved_by_rung, twin_by_rung = {}, {}
    for rung, (base, slots) in PILOT_B_BLOCKS.items():
        tbase, tslots = PILOT_B_TWIN_BLOCKS[rung]
        start = time.perf_counter()
        curved, cskip, cok = _fill_block(rung, base, slots, PILOT_SAMPLES)
        twin, tskip, tok = _fill_block(rung, tbase, tslots, PILOT_SAMPLES,
                                       flat=True)
        elapsed = time.perf_counter() - start
        for ok, arm, skips in ((cok, "curved", cskip), (tok, "twin", tskip)):
            if not ok:
                raise SystemExit(
                    f"pilot-B {arm} rung {rung} could not fill "
                    f"{PILOT_SAMPLES} (skips: {len(skips)}) -- "
                    "INFEASIBLE-INCOMPLETE")
        curved_by_rung[rung], twin_by_rung[rung] = curved, twin
        paired = paired_recovery(curved, twin)
        ys[rung] = paired["y"]
        nominal = bonett_variance_bound(paired["y"])
        calibrated = calibrated_variance_bound(
            paired["y"], "stage-b-bottom" if rung == P12B_LADDER[0]
            else "stage-b-top")
        summary["per_rung"][str(rung)] = {
            "n_samples": paired["n"],
            "variance": nominal["s2"],
            "kurtosis_g4": nominal["g4"],
            "variance_bound_95_nominal": nominal["bound"],
            "variance_bound_95_calibrated": calibrated["bound"],
            "calibration_z_used": calibrated["z_used"],
            "calibration_coverage_at_nominal":
                calibrated["coverage_at_nominal"],
            "boundary_fraction": paired["boundary_fraction"],
            "median_rel_error": paired["median_rel_error"],
            "mean_m_curved": float(np.mean([r["mean_m"] for r in curved])),
            "mean_m_twin": float(np.mean([r["mean_m"] for r in twin])),
            "skipped_seeds_curved": cskip,
            "skipped_seeds_twin": tskip,
            "mean_seconds_per_sample": elapsed / (2 * PILOT_SAMPLES),
            "y": [float(val) for val in paired["y"]],
        }
        print(f"pilot-12B rung={rung}: var={nominal['s2']:.5f} "
              f"g4={nominal['g4']:.2f} | boundary "
              f"{paired['boundary_fraction']:.1%} | skips "
              f"{len(cskip)}/{len(tskip)}", flush=True)

    # The WITHIN-arm gates run here: their tolerance is 5% against a
    # pilot resolution near 0.7%, so they can see what they gate on, and
    # an arm that has drifted off its design target should stop the
    # sequence where it costs 400 samples instead of 7800.
    #
    # The CROSS-arm gate does NOT bind here. Its tolerance is 1% and the
    # offset's standard error at 200 samples is about 0.94%, so gating on
    # it would abort roughly a third of the time on correctly calibrated
    # arms -- 10.6's own resolution principle, applied to the gate rather
    # than only to the calibration probe. It is measured and published
    # with its standard error, and Stage B enforces it.
    cross = cross_arm_m_gate(curved_by_rung, twin_by_rung, enforce=False)
    gates = [m_gate(curved_by_rung, "curved"), m_gate(twin_by_rung, "twin")]
    summary["m_gate_curved"], summary["m_gate_twin"] = gates[0], gates[1]
    summary["cross_arm_m_gate_measured"] = cross
    for rung, row in cross["per_rung"].items():
        print(f"  cross-arm rung {rung}: offset "
              f"{row['cross_arm_offset']:+.2%} +/- "
              f"{row['cross_arm_offset_se']:.2%} (tolerance "
              f"{CROSS_ARM_TOLERANCE:.0%}, resolution ratio "
              f"{_resolution_ratio_text(row['resolution_ratio'])})"
              " -- measured, not gating", flush=True)
    _abort_on_gate_failure(gates)

    power = power_requirements_b(ys[P12B_LADDER[0]], ys[P12B_LADDER[-1]])
    times = verification["per_rung"]
    projected = None
    if power["n_per_rung"] is not None:
        # three rungs, both arms
        per_sample = sum(
            times[str(r)]["mean_seconds_per_sample"] for r in P12B_LADDER
        )
        projected = 2.0 * power["n_per_rung"] * per_sample / 3600.0
    summary["power"] = power
    summary["n_twin"] = power["n_per_rung"]
    summary["projected_stage_b_hours"] = projected
    summary["feasible"] = bool(
        not power["infeasible"] and projected is not None
        and projected <= PROJECTION_LIMIT_HOURS
    )
    (output_dir / PILOT_ARTIFACT).write_text(
        json.dumps(_de_nan(summary), indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"power": power,
                      "projected_stage_b_hours": projected,
                      "feasible": summary["feasible"]}, indent=2))


def stage_b_record(delta: float, ci: list[float], flat_available: bool,
                   recovery: dict) -> dict:
    """Section 10.7's co-requirement, as a function so a test can drive
    it without a campaign.

    > A pass requires BOTH. (i) the Section 1.4 table returns IMPROVES
    > on Delta_B, and (ii) the top rung's DIMENSIONLESS recovery
    > satisfies |R tau^2 hat - R tau^2| / (R tau^2) <= 0.25. If (i)
    > passes and (ii) fails, the record reads RATE-ONLY and states in
    > the same sentence that no curvature was recovered.

    The reason this exists at all: the budget alone makes the failure
    visible, but it does not stop the VERDICT WORD from being read as
    recovery, which is what P13's Section 11 record cost.
    """

    # 10.4 derives this stage's own margin, `delta_eq_B = |Delta*_B|/3 =
    # 0.0836`, and 10.10 SIZES the equivalence row with it -- so reading
    # the 1.4 table at P11's 0.067, as this did until review C30, would
    # have returned INCONCLUSIVE on an interval the campaign had powered
    # itself to call FLAT-WITHIN-MARGIN. The table is shared; the margin
    # is this stage's.
    rate = verdict(ci[0], ci[1], flat_available, delta_eq=DELTA_EQ_B)
    recovered = bool(recovery["rel_error"] <= RECOVERY_THRESHOLD)
    if rate == "IMPROVES" and not recovered:
        outcome = "RATE-ONLY"
        statement = (
            "RATE-ONLY: the rate gate passed and NO CURVATURE WAS "
            f"RECOVERED -- the top rung's dimensionless recovery is "
            f"{recovery['rel_error']:.3f} against a threshold of "
            f"{RECOVERY_THRESHOLD:.2f}, so the gate measured a decaying "
            "bias and nothing else.")
    elif rate == "IMPROVES" and recovered:
        outcome = "RECOVERS-CURVATURE"
        statement = (
            "RECOVERS-CURVATURE: the rate gate passed and the top "
            f"rung's dimensionless recovery is {recovery['rel_error']:.3f} "
            f"<= {RECOVERY_THRESHOLD:.2f}. Scope, which the record must "
            "carry in the same breath: 1+1D metrics are conformally "
            "flat, so curvature reaches the causal order only through "
            "the volume this estimator reads; d >= 3 is where the "
            "general question lives.")
    else:
        outcome = rate
        statement = (
            f"{rate}: the rate gate did not return IMPROVES, so the "
            "co-requirement's second half is reported but not reached.")
    return {
        "delta_b": delta, "delta_b_ci": ci,
        "rate_verdict": rate,
        "recovery_rel_error": recovery["rel_error"],
        "recovery_rel_error_ci": recovery["rel_error_ci"],
        "recovery_threshold": RECOVERY_THRESHOLD,
        "recovery_passed": recovered,
        "outcome": outcome,
        "statement": statement,
    }


def run_stage_b(output_dir: Path) -> None:
    """Stage B (Sections 10.7 and 10.10): the campaign."""

    stamp = _preflight_clean()
    assert_windows_disjoint_and_fresh()
    gate = cross_stage_gate()
    pilot = _load_gate_artifact(output_dir, PILOT_ARTIFACT, stamp)
    if not pilot.get("feasible"):
        raise SystemExit("pilot-12B declared the design infeasible -- "
                         "Stage B refuses to run")
    n_per_rung = int(pilot["power"]["n_per_rung"])
    n_twin = int(pilot["n_twin"])
    flat_available = bool(pilot["power"]["flat_available"])

    rows: list[dict] = []
    curved_by_rung, twin_by_rung, ys = {}, {}, {}
    skips: dict = {}
    for rung in P12B_LADDER:
        base, slots = STAGE_B_BLOCKS[rung]
        tbase, tslots = STAGE_B_TWIN_BLOCKS[rung]
        curved, cskip, cok = _fill_block(rung, base, slots, n_per_rung)
        twin, tskip, tok = _fill_block(rung, tbase, tslots, n_twin,
                                       flat=True)
        for ok, arm, skipped, need in ((cok, "curved", cskip, n_per_rung),
                                       (tok, "twin", tskip, n_twin)):
            if not ok:
                raise SystemExit(
                    f"stage B {arm} rung {rung} could not fill {need} "
                    f"(skips: {len(skipped)}) -- INFEASIBLE-INCOMPLETE")
        for r in curved:
            r["stage"], r["code_version"] = "B12", stamp
        for r in twin:
            r["stage"], r["code_version"] = "B12-twin", stamp
        rows.extend(curved)
        rows.extend(twin)
        curved_by_rung[rung], twin_by_rung[rung] = curved, twin
        paired = paired_recovery(curved, twin)
        ys[rung] = paired["y"]
        # review C14 ported: BOTH arms' skipped identities are kept, not
        # just counted in an error message the run never reaches
        skips[rung] = {
            "curved": [int(s) for s in cskip],
            "twin": [int(s) for s in tskip],
            "paired": paired,
        }
        print(f"stage B-12 rung={rung}: {paired['n']} paired | mean y_B "
              f"{paired['y'].mean():+.4f} | boundary "
              f"{paired['boundary_fraction']:.1%} | below 1 "
              f"{paired['count_below_one']}/{paired['n']} | skips "
              f"{len(cskip)}/{len(tskip)}", flush=True)

    cross = cross_arm_m_gate(curved_by_rung, twin_by_rung, enforce=True)
    # This is where the cross-arm gate binds, so this is where its
    # resolution has to be adequate. A gate whose tolerance is inside its
    # own standard error is noise with a verdict attached, and it would be
    # worse here than in the pilot because it would abort a full campaign.
    # Refusing loudly beats gating meaninglessly.
    for rung, row in cross["per_rung"].items():
        ratio = row["resolution_ratio"]
        print(f"  cross-arm rung {rung}: offset "
              f"{row['cross_arm_offset']:+.2%} +/- "
              f"{row['cross_arm_offset_se']:.2%}"
              f" | resolution ratio {_resolution_ratio_text(ratio)}",
              flush=True)
        if ratio is not None and ratio < 2.0:
            raise SystemExit(
                f"cross-arm gate at rung {rung} cannot see its own "
                f"tolerance: {CROSS_ARM_TOLERANCE:.0%} against a standard "
                f"error of {row['cross_arm_offset_se']:.2%} (ratio "
                f"{ratio:.2f} < 2). Gating on it would be noise with a "
                "verdict attached -- 10.6's resolution principle.")
    gates = [m_gate(curved_by_rung, "curved"), m_gate(twin_by_rung, "twin"),
             cross]
    _abort_on_gate_failure(gates)

    # ONLY here. Both the gate's note and the abort's own message say
    # "nothing written", and review C27 found the campaign CSV going out
    # above the resolution check and all three m gates -- so an aborted
    # run left a complete, freezable-looking `p12_stage_b.csv` on disk
    # and the two statements were false about their own function. This
    # is C21's shape one layer out: there the guard wrote its failure
    # and continued, here the guard was honest and the write beat it.
    # Everything that can abort now runs above this line.
    write_rows_csv(output_dir / "p12_stage_b.csv", rows)

    bottom, top = P12B_LADDER[0], P12B_LADDER[-1]
    delta = float(ys[top].mean() - ys[bottom].mean())
    rng = np.random.default_rng(_stable_seed("p12-b-delta"))
    boots = [
        float(rng.choice(ys[top], size=ys[top].size, replace=True).mean()
              - rng.choice(ys[bottom], size=ys[bottom].size,
                           replace=True).mean())
        for _ in range(4000)
    ]
    ci = [float(np.percentile(boots, 2.5)),
          float(np.percentile(boots, 97.5))]

    recovery = {str(r): dimensionless_recovery(
        curved_by_rung[r], twin_by_rung[r], str(r)) for r in P12B_LADDER}
    decision = stage_b_record(delta, ci, flat_available, recovery[str(top)])

    summary = {
        "code_version": stamp,
        "cross_stage_gate": gate,
        "n_per_rung": n_per_rung, "n_twin": n_twin,
        "flat_available": flat_available,
        "tau_over_ell": TAU_ELL,
        "ricci_scalar_true": R_TRUE, "r_tau2_true": R_TAU2_TRUE,
        **decision,
        # the verdict word alone is not the record: Section 10.9 forbids
        # reading it as accuracy, so the size sits beside it
        "verdict": decision["outcome"],
        "dimensionless_recovery_by_rung": recovery,
        "mean_y_b_by_rung": {str(r): float(ys[r].mean())
                             for r in P12B_LADDER},
        "per_rung": {
            str(r): {
                "n_paired": skips[r]["paired"]["n"],
                "boundary_hits": skips[r]["paired"]["boundary_hits"],
                "boundary_fraction":
                    skips[r]["paired"]["boundary_fraction"],
                "median_rel_error":
                    skips[r]["paired"]["median_rel_error"],
                "count_below_one":
                    skips[r]["paired"]["count_below_one"],
                "fraction_below_one":
                    skips[r]["paired"]["fraction_below_one"],
                "variance": float(ys[r].var(ddof=1)),
            } for r in P12B_LADDER
        },
        "m_gate_curved": gates[0], "m_gate_twin": gates[1],
        "cross_arm_m_gate": cross,
        # Disclosed pre-run rather than discovered on an abort: each
        # rung's frozen offset (+0.55% / -0.48% / -0.15%) sits inside a
        # 1% tolerance by a margin comparable to this arm's standard
        # error, so the campaign has a real chance of tripping its own
        # cross-arm gate on correctly calibrated arms -- concentrated at
        # the bottom rung, where the frozen offset is both the largest
        # and the only positive one. The tolerance is NOT loosened for
        # it -- 10.6 chose 1% because a 1% offset injects ~17% into the
        # recovery -- and if it trips, the halt is the record.
        "cross_arm_spurious_trip_probability": {
            rung: float(cross_arm_spurious_trip_probability(
                FROZEN_CROSS_ARM_OFFSET[int(rung)],
                row["cross_arm_offset_se"]))
            for rung, row in cross["per_rung"].items()
            if row["cross_arm_offset_se"] > 0
        },
        "skip_counts": {str(r): {"curved": len(skips[r]["curved"]),
                                 "twin": len(skips[r]["twin"])}
                        for r in P12B_LADDER},
        "skipped_seeds": {str(r): {"curved": skips[r]["curved"],
                                   "twin": skips[r]["twin"]}
                          for r in P12B_LADDER},
        # spans BOTH arms (P13 review C14) and inherits P13's own
        "selection_caveat": bool(
            any(skips[r]["curved"] or skips[r]["twin"]
                for r in P12B_LADDER)
            or gate["p13_stage_a_13c"]["selection_caveat"]),
        "scope": (
            "Supports: the scalar curvature of the ambient geometry is "
            "recoverable from causal order plus the interval count, "
            "OVER THIS LADDER, to the accuracy reported. Does NOT "
            "support the general claim: 1+1D metrics are conformally "
            "flat, so curvature reaches the order only through the "
            "volume this estimator reads. d >= 3 is where that "
            "question lives."),
        "construction_limit": (
            "A flat twin cannot reproduce a curved box's intra-box "
            "intensity gradient, so the residual bias mismatch is a "
            "property of the CONSTRUCTION, not of the instrument. "
            "Removing it needs a control matched in the box-area "
            "distribution rather than only in m, which is a different "
            "design."),
        "labelled_checks": {
            "mean_m_by_rung_curved": {
                str(r): float(np.mean([x["mean_m"] for x in
                                       curved_by_rung[r]]))
                for r in P12B_LADDER},
            "mean_m_by_rung_twin": {
                str(r): float(np.mean([x["mean_m"] for x in
                                       twin_by_rung[r]]))
                for r in P12B_LADDER},
            "mean_tau_hat_by_rung": {
                str(r): float(np.mean([x["mean_tau_hat"] for x in
                                       curved_by_rung[r]]))
                for r in P12B_LADDER},
        },
        "seed_blocks": [{"label": la, "lo": lo, "hi": hi}
                        for la, lo, hi in seed_blocks()],
    }
    (output_dir / STAGE_ARTIFACT).write_text(
        json.dumps(_de_nan(summary), indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "delta_b": delta, "delta_b_ci": ci,
        "rate_verdict": decision["rate_verdict"],
        "recovery_rel_error": decision["recovery_rel_error"],
        "recovery_passed": decision["recovery_passed"],
        "outcome": decision["outcome"],
    }, indent=2))
    print("\n" + decision["statement"])


STAGES = {"verify-b": run_verify_b, "pilot-b": run_pilot_b,
          "b": run_stage_b}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=sorted(STAGES), default=None)
    parser.add_argument("--output-dir", type=Path,
                        default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args(argv)
    if args.stage is None:
        raise SystemExit("choose --stage " + "/".join(sorted(STAGES)))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    STAGES[args.stage](args.output_dir)


if __name__ == "__main__":
    main()
