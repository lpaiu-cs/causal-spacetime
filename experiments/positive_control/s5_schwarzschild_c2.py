"""S5: preregistered Schwarzschild C2-unpaired discrimination.

The frozen rule is docs/prereg/p14_s5_schwarzschild_c2.md. Two
INDEPENDENT arms (fresh Schwarzschild sprinkles; fresh flat sprinkles)
on the common measure sampler; the per-reading global relation
fraction f is the classifier score, direction frozen (curved BELOW
flat, the S3/S4 sign).

Margin declaration (order matters -- review ruling): the exploration
results are a DESIGN INPUT, and AUC 0.60 is frozen as the minimum
practically useful single-poset discrimination INDEPENDENTLY of the
budget, of the selected n, and of the future confirmation data
(confirmation-data-independent, not anchor-independent).

Primary metric and gate (AUC, DeLong, strict comparisons, unclipped):
  DETECTED             AUC CI95 lower > 0.60
  EQUIVALENT-AT-MARGIN CI95 entirely inside (0.40, 0.60)
  DIRECTION-REVERSED   CI95 upper < 0.40
  INCONCLUSIVE         otherwise
Secondary gate (BA, out-of-sample, deterministic CI): threshold
learned on the TRAIN half only (midpoint of train arm means,
direction frozen, ties to the flat side); per-arm TEST accuracies get
exact Clopper-Pearson intervals at 97.5% two-sided coverage each
(Bonferroni joint 95%), and the BA interval is the average of the two
bounds; gate: BA lower > 0.60 (eps_BA = 0.10). No inner bootstrap, so
the outer power certification is single-layer: each replicate
retrains the threshold and evaluates the deterministic CI.

Sample size: n/arm is chosen by the frozen rule "the minimum
n in {300, 600, 1000} for which EVERY required branch (AUC effect,
AUC equivalence, BA effect) certifies CP95 lower bound >= 0.90"; if
all three candidates fail, the stage is blocked/review-reopened --
no new n is added silently.

Power sources (frozen; AUC is a rank statistic, so a common positive
rescaling of both arms cannot change it -- no conservativeness claim
is attached to scaling in the null):
  null:   BOTH arms resampled independently from the S3 curved
          marginal empirical distribution, unscaled.
  effect: the S3 flat and curved marginals resampled independently,
          each arm's residuals about its own mean scaled by that
          arm's chi-square one-sided 95% SD upper factor and
          re-centered at the ORIGINAL arm mean -- spreads widen
          (conservative for detection), the mean gap is preserved
          exactly.
Streams: entropy vectors [781, 4860, n, j] (row j at candidate n),
design lineage, never campaign streams.

Campaign seeds: curved arm 40_000_251, flat arm 40_000_261 (fresh,
ledger-asserted at entry). Freeze/result commit split, exact-checkout
execution, entry/exit git state and content-addressed manifest
verification as in S4; the results commit must move both seeds to
OBSERVED, flip this runner to replay, and preserve the executed
freeze manifest snapshot immediately.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
import subprocess
import time
from pathlib import Path

import numpy as np
import s1_schwarzschild_cost as s1
from probe_seed_ledger import (
    FRESH_PROBE_SCALARS,
    S3_SMOKE_SEED,
    replay_scalar,
)
from s3_schwarzschild_probe import flat_relation, reading

_REPO = Path(__file__).resolve().parents[2]
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_s5_results.json"
_REPLAY_ARTIFACT = _REPO / "docs" / "prereg" / "p14_s5_replay_results.json"
_FREEZE_MANIFEST = _REPO / "docs" / "prereg" / "p14_s5_freeze_manifest.json"
_S3_ARTIFACT = _REPO / "docs" / "prereg" / "p14_s3_probe_results.json"

# ---------------------------------------------------------------------
# FROZEN margins and rules
# ---------------------------------------------------------------------

EPS_AUC = 0.10               # minimum useful discrimination: AUC 0.60
EPS_BA = 0.10                # BA secondary threshold 0.60
Z975 = 1.959964              # two-sided 95% normal critical value
CANDIDATE_N = (300, 600, 1000)
POWER_TARGET = 0.90
POWER_B = 20_000
POWER_SEED_BASE = (781, 4860)

#: SELECTED by the frozen minimum-n rule (certification pinned in the
#: prereg doc and recomputed by tests): n = 300 certifies every
#: required branch (20000/20000, 19564/20000, 20000/20000).
N_ARM = 300
E_N = 300
SEED_CURVED = 40_000_251     # OBSERVED (outcome DETECTED) -- reruns
SEED_FLAT = 40_000_261       # of the campaign path are replays
SMOKE_SEED = S3_SMOKE_SEED   # 40_000_221, observed smoke stream

SENTENCES = {
    "DETECTED": ("동결된 Schwarzschild 도메인·밀도에서, 단일 causal set의 "
                 "global relation fraction은 flat/Schwarzschild 앙상블을 "
                 "우연 수준보다 판별하는 정보를 운반한다 (AUC CI95 하한 > "
                 "0.60, 프로그램 내부 진술)."),
    "EQUIVALENT-AT-MARGIN": ("AUC가 chance의 동결 ±0.10 무시가능성 밴드 "
                             "안으로 해상됐다 — 판별력 부재의 단언이 "
                             "아니다."),
    "DIRECTION-REVERSED": ("판별이 동결 방향과 반대로 해상됐다 — 동결 "
                           "방향 규칙 위반, stage 실패로 기록하고 원인을 "
                           "재검토한다."),
    "INCONCLUSIVE": "어느 분기 조건도 충족되지 않아 판정을 유보한다.",
    "BA_PASS": ("out-of-sample balanced accuracy의 결합 95% 하한이 "
                "0.60을 넘어, 학습-외 판별이 확인됐다 (secondary)."),
    "BA_FAIL": ("out-of-sample balanced accuracy는 0.60 문턱을 넘지 "
                "못했다 (secondary; stage 판정에 불참)."),
}


# ---------------------------------------------------------------------
# DeLong AUC with midrank ties, unclipped two-sided 95% CI
# ---------------------------------------------------------------------


def auc_brute(curved: np.ndarray, flat: np.ndarray) -> float:
    """Exact pairwise AUC = P(curved < flat) + half ties -- the
    validation reference for the DeLong implementation."""

    less = (curved[:, None] < flat[None, :]).sum()
    ties = (curved[:, None] == flat[None, :]).sum()
    return float((less + 0.5 * ties) / (len(curved) * len(flat)))


def delong_auc_ci(curved: np.ndarray, flat: np.ndarray
                  ) -> tuple[float, float, float]:
    """(auc, lower, upper): midrank DeLong variance, CI unclipped.

    Placement structure: V10_i = mean_j psi(c_i, f_j),
    V01_j = mean_i psi(c_i, f_j), psi = 1(c<f) + 0.5 * 1(c==f);
    var = S10/m + S01/n with sample variances (ddof=1)."""

    m, n = len(curved), len(flat)
    psi = ((curved[:, None] < flat[None, :]).astype(float)
           + 0.5 * (curved[:, None] == flat[None, :]))
    auc = float(psi.mean())
    v10 = psi.mean(axis=1)
    v01 = psi.mean(axis=0)
    var = v10.var(ddof=1) / m + v01.var(ddof=1) / n
    half = Z975 * math.sqrt(var)
    return auc, auc - half, auc + half


def auc_outcome(ci: tuple[float, float]) -> str:
    """The frozen four-way partition, strict comparisons."""

    lo, hi = ci
    if lo > 0.5 + EPS_AUC:
        return "DETECTED"
    if 0.5 - EPS_AUC < lo and hi < 0.5 + EPS_AUC:
        return "EQUIVALENT-AT-MARGIN"
    if hi < 0.5 - EPS_AUC:
        return "DIRECTION-REVERSED"
    return "INCONCLUSIVE"


# ---------------------------------------------------------------------
# Exact Clopper-Pearson bounds (cached per test-arm size)
# ---------------------------------------------------------------------


def _binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Bin(n, p), exact log-sum."""

    if k <= 0:
        return 1.0
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    lp, lq = math.log(p), math.log(1.0 - p)
    terms = [math.lgamma(n + 1) - math.lgamma(i + 1)
             - math.lgamma(n - i + 1) + i * lp + (n - i) * lq
             for i in range(k, n + 1)]
    mx = max(terms)
    return math.exp(mx) * sum(math.exp(t - mx) for t in terms)


@functools.lru_cache(maxsize=64)
def _cp_table(n: int, alpha_each: float) -> tuple:
    """Clopper-Pearson [lower, upper] for every k in 0..n at
    two-sided coverage 1 - 2*alpha_each (alpha_each per side)."""

    rows = []
    for k in range(n + 1):
        if k == 0:
            lo = 0.0
        else:
            a, b = 0.0, 1.0
            for _ in range(80):
                mid = 0.5 * (a + b)
                if _binom_sf(k, n, mid) < alpha_each:
                    a = mid
                else:
                    b = mid
            lo = 0.5 * (a + b)
        if k == n:
            hi = 1.0
        else:
            a, b = 0.0, 1.0
            for _ in range(80):
                mid = 0.5 * (a + b)
                if 1.0 - _binom_sf(k + 1, n, mid) < alpha_each:
                    b = mid
                else:
                    a = mid
            hi = 0.5 * (a + b)
        rows.append((lo, hi))
    return tuple(rows)


def ba_cp_bonferroni(curved_test: np.ndarray, flat_test: np.ndarray,
                     threshold: float) -> tuple[float, float, float]:
    """(ba, lower, upper): per-arm test accuracies with exact CP
    intervals at 97.5% two-sided coverage each (alpha 0.0125 per
    side; Bonferroni joint 95%), BA interval = average of bounds.
    Ties at the threshold go to the flat side (frozen direction rule:
    curved means f < threshold)."""

    kc = int((curved_test < threshold).sum())
    kf = int((flat_test >= threshold).sum())
    tc = _cp_table(len(curved_test), 0.0125)
    tf = _cp_table(len(flat_test), 0.0125)
    ba = 0.5 * (kc / len(curved_test) + kf / len(flat_test))
    return (ba, 0.5 * (tc[kc][0] + tf[kf][0]),
            0.5 * (tc[kc][1] + tf[kf][1]))


def ba_gate(ba_ci: tuple[float, float]) -> bool:
    return ba_ci[0] > 0.5 + EPS_BA


def train_test_split(arm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Frozen split: first half by reading index is TRAIN."""

    half = len(arm) // 2
    return arm[:half], arm[half:]


def learn_threshold(curved_train: np.ndarray,
                    flat_train: np.ndarray) -> float:
    """Frozen rule: midpoint of the two TRAIN arm means."""

    return 0.5 * (float(curved_train.mean()) + float(flat_train.mean()))


# ---------------------------------------------------------------------
# Power certification (single-layer: deterministic CIs per replicate)
# ---------------------------------------------------------------------


def _gammp(a: float, x: float) -> float:
    """Regularized lower incomplete gamma P(a, x), exact to double
    precision: series for x < a + 1, continued fraction otherwise
    (the standard gammp/gammq pair). Verified independently by the
    identities P(1/2, x) = erf(sqrt(x)) and P(1, x) = 1 - e^-x in the
    contract tests."""

    if x < 0.0 or a <= 0.0:
        raise ValueError("gammp domain")
    if x == 0.0:
        return 0.0
    lg = math.lgamma(a)
    if x < a + 1.0:
        ap, term, total = a, 1.0 / a, 1.0 / a
        for _ in range(500):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-16:
                break
        return total * math.exp(-x + a * math.log(x) - lg)
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return 1.0 - math.exp(-x + a * math.log(x) - lg) * h


@functools.lru_cache(maxsize=32)
def chi2_quantile(p: float, df: int) -> float:
    """EXACT chi-square quantile by bisection on the CDF
    P(df/2, q/2) -- the frozen replacement for the Wilson-Hilferty
    approximation the review caught (relative error -5.65e-7 at
    df = 299, optimistic direction)."""

    lo, hi = 0.0, df * 10.0 + 100.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _gammp(df / 2.0, mid / 2.0) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _sd_upper_factor(df: int) -> float:
    """Exact chi-square one-sided 95% upper-bound factor for a sample
    SD: sqrt(df / chi2_quantile(0.05, df)). At df = 299 this is
    1.0724931824796697 (pinned by test against the review's
    independent computation)."""

    factor = math.sqrt(df / chi2_quantile(0.05, df))
    assert factor > 1.0
    return factor


def power_sources(s3_curved: np.ndarray, s3_flat: np.ndarray) -> dict:
    """The frozen resampling sources (module docstring)."""

    mc, mf = float(s3_curved.mean()), float(s3_flat.mean())
    fc = _sd_upper_factor(len(s3_curved) - 1)
    ff = _sd_upper_factor(len(s3_flat) - 1)
    return {
        "null_pool": s3_curved,
        "effect_curved": mc + (s3_curved - mc) * fc,
        "effect_flat": mf + (s3_flat - mf) * ff,
    }


def power_row(row: str, n_arm: int, s3_curved: np.ndarray,
              s3_flat: np.ndarray, b: int = POWER_B) -> int:
    """Successes over b replicates for one certification row at one
    candidate n. Rows: 'auc_effect' (DETECTED), 'auc_null'
    (EQUIVALENT-AT-MARGIN), 'ba_effect' (BA gate). Stream
    [781, 4860, n_arm, j] with j = frozen row index."""

    src = power_sources(s3_curved, s3_flat)
    j = {"auc_effect": 0, "auc_null": 1, "ba_effect": 2}[row]
    rng = np.random.default_rng([*POWER_SEED_BASE, n_arm, j])
    ok = 0
    for _ in range(b):
        if row == "auc_null":
            pool = src["null_pool"]
            c = pool[rng.integers(0, len(pool), n_arm)]
            f = pool[rng.integers(0, len(pool), n_arm)]
        else:
            ec, ef = src["effect_curved"], src["effect_flat"]
            c = ec[rng.integers(0, len(ec), n_arm)]
            f = ef[rng.integers(0, len(ef), n_arm)]
        if row == "auc_effect":
            _, lo, hi = delong_auc_ci(c, f)
            ok += auc_outcome((lo, hi)) == "DETECTED"
        elif row == "auc_null":
            _, lo, hi = delong_auc_ci(c, f)
            ok += auc_outcome((lo, hi)) == "EQUIVALENT-AT-MARGIN"
        else:
            ct, cs = train_test_split(c)
            ft, fs = train_test_split(f)
            thr = learn_threshold(ct, ft)
            _, lo, _ = ba_cp_bonferroni(cs, fs, thr)
            ok += lo > 0.5 + EPS_BA
    return ok


def select_n(s3_curved: np.ndarray, s3_flat: np.ndarray,
             b: int = POWER_B) -> tuple[int | None, dict]:
    """The frozen selection rule: minimum candidate n whose three
    required rows all certify CP95 lower >= POWER_TARGET; None means
    blocked/review-reopened (no silent extra candidates)."""

    table = {}
    chosen = None
    for n_arm in CANDIDATE_N:
        rows = {}
        for row in ("auc_effect", "auc_null", "ba_effect"):
            k = power_row(row, n_arm, s3_curved, s3_flat, b)
            lo = cp_lower(k, b)
            rows[row] = {"successes": k, "b": b, "cp95_lower": lo}
        table[n_arm] = rows
        if chosen is None and all(r["cp95_lower"] >= POWER_TARGET
                                  for r in rows.values()):
            chosen = n_arm
            break
    return chosen, table


def cp_lower(k: int, n: int, conf: float = 0.95) -> float:
    """Exact Clopper-Pearson lower bound for k successes of n."""

    if k == n:
        return ((1.0 - conf) / 2.0) ** (1.0 / n)
    a, bb = 0.0, 1.0
    for _ in range(80):
        mid = 0.5 * (a + bb)
        if _binom_sf(k, n, mid) < (1.0 - conf) / 2.0:
            a = mid
        else:
            bb = mid
    return 0.5 * (a + bb)


def negative_controls() -> dict:
    """The frozen falsifiability battery: constructed arm pairs on
    stream [781, 4860, 0, 99] for which each gate must fail or the
    DIRECTION-REVERSED branch must be reached (never hidden in
    INCONCLUSIVE)."""

    r = json.loads(_S3_ARTIFACT.read_text(encoding="utf-8"))
    s3c = np.array(r["f_schwarzschild_lower"]["per_reading"])
    s3f = np.array(r["f_flat"]["per_reading"])
    rng = np.random.default_rng([*POWER_SEED_BASE, 0, 99])

    def draw(pool, n=N_ARM):
        return pool[rng.integers(0, len(pool), n)]

    # (i) no effect: both arms one pool -> DETECTED must not fire
    c0, f0 = draw(s3c), draw(s3c)
    out_null = auc_outcome(delong_auc_ci(c0, f0)[1:])
    # (ii) real effect -> EQUIVALENT-AT-MARGIN must not fire
    ce, fe = draw(s3c), draw(s3f)
    out_eff = auc_outcome(delong_auc_ci(ce, fe)[1:])
    # (iii) reversed strong effect -> DIRECTION-REVERSED reached
    out_rev = auc_outcome(delong_auc_ci(fe, ce)[1:])
    # (iv) BA gate fails with no effect (threshold from its own train)
    ct, cs = train_test_split(c0)
    ft, fs = train_test_split(f0)
    ba_null = ba_gate(ba_cp_bonferroni(
        cs, fs, learn_threshold(ct, ft))[1:])
    return {"detected_on_null": out_null == "DETECTED",
            "equivalent_on_effect": out_eff == "EQUIVALENT-AT-MARGIN",
            "effect_outcome": out_eff,
            "reversed_outcome": out_rev,
            "ba_gate_on_null": ba_null}


# ---------------------------------------------------------------------
# The campaign: two independent arms
# ---------------------------------------------------------------------


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_freeze(stage: str) -> None:
    """Content-addressed freeze identity, verified at entry and
    re-verified at exit (the S4 mechanism, unchanged in kind)."""

    manifest = json.loads(_FREEZE_MANIFEST.read_text(encoding="utf-8"))
    for rel, want in manifest["files"].items():
        got = _sha256(_REPO / rel)
        if got != want:
            raise SystemExit(
                f"S5 {stage}: frozen protocol surface drifted: {rel}\n"
                f"  frozen {want}\n  found  {got}")


def _git_state() -> dict:
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=_REPO,
                         capture_output=True, text=True,
                         check=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=_REPO,
                           capture_output=True, text=True,
                           check=True).stdout.strip()
    return {"rev": rev, "dirty": bool(dirty)}


def flat_reading(rng: np.random.Generator, n_events: int) -> float:
    """Flat-arm census: exact chord predicate only (no ambiguity)."""

    pts = s1.sample_events(n_events, rng)
    pts = pts[np.argsort(pts[:, 0], kind="stable")]
    rel = 0
    for i in range(n_events - 1):
        for j in range(i + 1, n_events):
            rel += flat_relation(pts[i], pts[j])
    return rel / (n_events * (n_events - 1) // 2)


def curved_arm(rng: np.random.Generator, n_readings: int,
               e_n: int) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Curved-arm censuses via the tested S3 reading (its flat side is
    computed and DISCARDED -- negligible cost, and the arm's sprinkles
    are its own, so nothing is paired). Returns (f_lower, f_upper,
    ambiguous, escalated)."""

    lo = np.empty(n_readings)
    hi = np.empty(n_readings)
    amb = esc = 0
    start = time.perf_counter()
    for k in range(n_readings):
        n = int(rng.poisson(e_n))
        if n < 2:
            raise SystemExit(f"S5: degenerate sprinkle n={n}")
        fm_lo, fm_hi, _f0_discarded, a, e = reading(rng, n)
        lo[k], hi[k] = fm_lo, fm_hi
        amb += a
        esc += e
        if (k + 1) % 10 == 0 or k == 0 or (k + 1) == n_readings:
            dt = time.perf_counter() - start
            eta = dt / (k + 1) * (n_readings - k - 1)
            print(f"  s5 curved {k + 1}/{n_readings}  "
                  f"elapsed {dt / 60:.1f} min, eta {eta / 60:.1f} min",
                  flush=True)
    return lo, hi, amb, esc


def flat_arm(rng: np.random.Generator, n_readings: int,
             e_n: int) -> np.ndarray:
    out = np.empty(n_readings)
    for k in range(n_readings):
        n = int(rng.poisson(e_n))
        if n < 2:
            raise SystemExit(f"S5: degenerate sprinkle n={n}")
        out[k] = flat_reading(rng, n)
    return out


def stage_outcome(c_lower: np.ndarray, c_upper: np.ndarray,
                  f: np.ndarray) -> tuple[str, dict]:
    """Identified-agreement rule (frozen): every AUC quantity is
    evaluated on BOTH curved bound series; a branch fires only if
    BOTH series give it, else INCONCLUSIVE. With zero ambiguity the
    series coincide and this is the plain rule."""

    ci_lo = delong_auc_ci(c_lower, f)
    ci_hi = delong_auc_ci(c_upper, f)
    o1, o2 = auc_outcome(ci_lo[1:]), auc_outcome(ci_hi[1:])
    outcome = o1 if o1 == o2 else "INCONCLUSIVE"
    detail = {"auc_lower_series": ci_lo, "auc_upper_series": ci_hi,
              "outcomes": [o1, o2]}
    return outcome, detail


def main() -> None:
    # Fail-closed CLI (PR #60 review): the ONLY accepted argument is
    # --smoke. argparse rejects unknown arguments and handles --help
    # by exiting, so a typo or a help invocation can never fall
    # through into the fresh-seed campaign.
    parser = argparse.ArgumentParser(
        description="S5 preregistered campaign runner "
                    "(docs/prereg/p14_s5_schwarzschild_c2.md); the "
                    "no-argument form runs the ONE official campaign "
                    "and consumes the fresh seeds.",
        allow_abbrev=False)
    parser.add_argument("--smoke", action="store_true",
                        help="validation run on the observed smoke "
                             "stream; writes no artifact")
    smoke = parser.parse_args().smoke
    n_readings = 6 if smoke else N_ARM
    e_n = 40 if smoke else E_N

    verify_freeze("entry")
    if smoke:
        seed_c = seed_f = SMOKE_SEED
    else:
        seed_c = replay_scalar("s5_curved")
        seed_f = replay_scalar("s5_flat")
        if (SEED_CURVED, SEED_FLAT) != (seed_c, seed_f):
            raise SystemExit("S5: seeds drifted from the observed "
                             "streams; a rerun must replay them.")
    state_start = _git_state()
    if state_start["dirty"] and not smoke:
        raise SystemExit(
            "S5: refusing a campaign run from a dirty tree "
            f"(rev {state_start['rev']}).")

    c_lo, c_hi, amb, esc = curved_arm(
        np.random.default_rng(seed_c), n_readings, e_n)
    f = flat_arm(np.random.default_rng(seed_f), n_readings, e_n)

    state_end = _git_state()
    if not smoke and state_end != state_start:
        raise SystemExit("S5: refusing to write -- tree changed "
                         f"({state_start} -> {state_end}).")
    if not smoke:
        verify_freeze("exit")

    outcome, detail = stage_outcome(c_lo, c_hi, f)
    # BA secondary on the unfavorable (upper) curved series
    ct, cs = train_test_split(c_hi)
    ft, fs = train_test_split(f)
    thr = learn_threshold(ct, ft)
    ba, ba_lo, ba_hi = ba_cp_bonferroni(cs, fs, thr)
    ba_pass = ba_gate((ba_lo, ba_hi))

    sentences = [SENTENCES[outcome],
                 SENTENCES["BA_PASS" if ba_pass else "BA_FAIL"]]
    if smoke:
        run_kind = "smoke"
    elif ("s5_curved" in FRESH_PROBE_SCALARS
          and "s5_flat" in FRESH_PROBE_SCALARS):
        run_kind = "fresh_observation"
    else:
        run_kind = "replay"
    result = {
        "stage": "S5 Schwarzschild C2-unpaired discrimination",
        "rule": "docs/prereg/p14_s5_schwarzschild_c2.md",
        "run_kind": run_kind,
        "params": {"n_arm": n_readings, "e_n": e_n,
                   "seed_curved": seed_c, "seed_flat": seed_f,
                   "smoke": smoke},
        "margins": {"eps_auc": EPS_AUC, "eps_ba": EPS_BA},
        "outcome": outcome,
        "auc": detail,
        "ba": {"threshold_train_midpoint": thr, "ba": ba,
               "ci95_cp_bonferroni": [ba_lo, ba_hi], "pass": ba_pass},
        "sentences": sentences,
        "curved_f_lower": {"mean": float(c_lo.mean()),
                           "per_reading": [float(v) for v in c_lo]},
        "curved_f_upper": {"mean": float(c_hi.mean()),
                           "per_reading": [float(v) for v in c_hi]},
        "flat_f": {"mean": float(f.mean()),
                   "per_reading": [float(v) for v in f]},
        "ambiguity": {"ambiguous": amb, "escalated": esc},
        "code": {"start": state_start, "end": state_end},
    }

    print(f"\nS5: outcome {outcome}  "
          f"AUC(lower-series) {detail['auc_lower_series'][0]:.4f} "
          f"CI [{detail['auc_lower_series'][1]:.4f}, "
          f"{detail['auc_lower_series'][2]:.4f}]  "
          f"BA {ba:.4f} [{ba_lo:.4f}, {ba_hi:.4f}] pass={ba_pass}  "
          f"ambiguous {amb}, escalated {esc}")

    if smoke:
        print("smoke run -- artifact NOT written")
        return
    # Write-ownership boundary (the S4/PR #58 rule): replay output is
    # owned by its own path; a 'fresh' run while the fresh artifact
    # exists refuses instead of replacing the executed lineage.
    if run_kind == "replay":
        target = _REPLAY_ARTIFACT
        result["replay_of"] = ("observed streams 40000251/40000261; "
                               "the fresh observation is "
                               "docs/prereg/p14_s5_results.json "
                               "(lineage 86e3674) and this replay may "
                               "not replace it")
    elif _ARTIFACT.exists():
        raise SystemExit(
            "S5: the fresh-observation artifact already exists; any "
            "rerun is a replay and may not replace it.")
    else:
        target = _ARTIFACT
    target.write_text(json.dumps(result, indent=2, ensure_ascii=False)
                      + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
