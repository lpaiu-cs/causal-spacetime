"""S4: preregistered Schwarzschild C1-paired confirmation.

The frozen rule is docs/prereg/p14_s4_schwarzschild_c1.md; this runner
implements it and nothing else. Protocol (domain, measure, Poisson-N
sampling, predicates, tol ladder, undecided contract) is S3's frozen
protocol, imported from the S3 module rather than duplicated.

Gates (all identified quantities, all STRICT comparisons -- an exact
boundary value does not pass):
  Gate A (primary, C1 detection):  S4 identified CI95 top < -EPS_DET
  Gate B (secondary, replication): identified Welch CI95 of
    theta_S4 - theta_S3 -> REPLICATED / DISCORDANT / B-INCONCLUSIVE
  Stage verdict: CONFIRMED = A and REPLICATED;
    DETECTED-NOT-REPLICATED = A and not REPLICATED;
    NO-DETECTION = S4 identified CI95 inside (-EPS_DET, +EPS_DET);
    else INCONCLUSIVE.
B=REPLICATED forces A at these margins (doc section 4), so CONFIRMED
is numerically equivalent to REPLICATED; both verdicts are still
reported by their own frozen sentences.

Power certification (doc section 7) is reproduced by `power_row` and
`negative_controls`, contract-tested against the pinned table: both
blocks independently resampled from the centered SD_cons-scaled S3
empirical distribution, row j on entropy-vector stream [781, 4840, j]
(S4 block drawn before S3 block in each replicate).

Seed discipline: the campaign draws the fresh allocation 40_000_241
exactly once; --smoke draws only the observed smoke stream 40_000_221.
The results commit must move the seed to OBSERVED and flip this
runner's entry to the replay path.
"""

from __future__ import annotations

import functools
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from p14_probe_p2 import student_t_crit
from probe_seed_ledger import S3_SMOKE_SEED, S4_SEED, assert_fresh_scalar
from s3_schwarzschild_probe import reading

_REPO = Path(__file__).resolve().parents[2]
_S3_ARTIFACT = _REPO / "docs" / "prereg" / "p14_s3_probe_results.json"
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_s4_results.json"
_FREEZE_MANIFEST = _REPO / "docs" / "prereg" / "p14_s4_freeze_manifest.json"


def _sha256(path: Path) -> str:
    """Raw-byte digest; every manifest path is .gitattributes-pinned
    to LF, so the working-tree bytes equal the committed blob bytes
    on every platform (the PR #55 storage-boundary rule)."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_freeze(stage: str) -> None:
    """Refuse unless every protocol-surface file matches the frozen
    digest table -- the content-addressed freeze identity, which
    survives merge commits (PR #57 review: a clean-but-LATER tree
    must not consume the once-only fresh seed on a drifted protocol;
    a clean-tree check alone cannot see that). Verified at entry and
    re-verified at exit before the artifact is written."""

    manifest = json.loads(_FREEZE_MANIFEST.read_text(encoding="utf-8"))
    for rel, want in manifest["files"].items():
        got = _sha256(_REPO / rel)
        if got != want:
            raise SystemExit(
                f"S4 {stage}: frozen protocol surface drifted: {rel}\n"
                f"  frozen {want}\n  found  {got}\n"
                "The campaign must run from the exact freeze content; "
                "any change requires review and a new freeze commit.")

# ---------------------------------------------------------------------
# FROZEN constants (docs/prereg/p14_s4_schwarzschild_c1.md)
# ---------------------------------------------------------------------

N_READINGS = 300
E_N = 300
EPS_DET = 0.0036             # ~9.98% of |anchor| ("about 10%")
EPS_REP = 0.0012             # rounded from the S3 SD 0.0012108375
KAPPA = 3.0e-4               # assurance region half-width (eps_rep/4)
SD_CONS = 0.00129861497      # chi2 one-sided 95% upper bound, df=299
POWER_SEED = (781, 4840)     # entropy vector, design lineage
POWER_B = 20_000
SEED = S4_SEED               # 40_000_241, fresh allocation
SMOKE_SEED = S3_SMOKE_SEED   # 40_000_221, observed smoke stream

SENTENCES = {
    "GATE_A": ("동결된 Schwarzschild 좌표·도메인(M=1, r∈[10,20], 극관각 캡 "
               "1.0, T=40)의 공통 측도 위에서, paired 앙상블 평균 이동의 "
               "identified CI95가 동결 방향으로 검출문턱 ε_det = 0.0036을 "
               "초과했다 — 유한밀도 인과 census가 type D 진공 곡률의 "
               "빛원뿔 변형을 C1-급으로 검출했다(프로그램 내부 진술)."),
    "REPLICATED": ("S4 블록과 S3 탐색 블록의 독립 두-표본 차이의 identified "
                   "Welch CI95가 ±ε_rep = ±0.0012 안에 들어, 탐색 효과가 "
                   "정량적으로 재현됐다."),
    "DISCORDANT": ("차이의 identified Welch CI95가 ±ε_rep 밖에 전부 놓여, "
                   "S4 효과 크기가 탐색 효과와 불일치한다."),
    "B-INCONCLUSIVE": ("차이의 identified Welch CI95가 ±ε_rep 경계에 걸쳐 "
                       "재현 판정을 유보한다."),
    "DETECTED-NOT-REPLICATED": ("검출문턱은 초과했으나 탐색 효과의 정량 "
                                "재현에는 이르지 못했다 — 검출 판정은 "
                                "유효하되 효과 크기는 재심의 대상이다."),
    "NO-DETECTION": ("동결 조건과 검출문턱 ε_det에서 평균 이동이 0과 "
                     "동등한 것으로 판정됐다."),
    "INCONCLUSIVE": "어느 분기 조건도 충족되지 않아 판정을 유보한다.",
}


# ---------------------------------------------------------------------
# Frozen statistics: identified intervals and gates
# ---------------------------------------------------------------------


@functools.lru_cache(maxsize=1024)
def _t_crit_int(df: int) -> float:
    return student_t_crit(df)


def _t_crit(df: float) -> float:
    """Two-sided 95% Student-t critical value at floor(df) -- the
    frozen conservative rule for non-integer Welch df."""

    return _t_crit_int(max(1, int(math.floor(df))))


def s3_block() -> np.ndarray:
    r = json.loads(_S3_ARTIFACT.read_text(encoding="utf-8"))
    return np.array(r["delta_lower"]["per_reading"])


def identified_ci(lower: np.ndarray, upper: np.ndarray) -> tuple[float, float]:
    """S4-side identified CI95: the union bound over the two 97.5%
    one-sided endpoints (doc section 6). With zero ambiguity
    lower == upper and this is the ordinary Student-t CI."""

    n = len(lower)
    t = _t_crit(n - 1)
    lo = float(lower.mean()) - t * float(lower.std(ddof=1)) / math.sqrt(n)
    hi = float(upper.mean()) + t * float(upper.std(ddof=1)) / math.sqrt(n)
    return lo, hi


def welch_identified_ci(lower: np.ndarray, upper: np.ndarray,
                        s3d: np.ndarray) -> tuple[float, float]:
    """Identified Welch CI95 of theta_S4 - theta_S3 (doc section 6):
    per-bound Welch SE and Welch-Satterthwaite df, floor(df) rule."""

    n4, n3 = len(lower), len(s3d)
    m3 = float(s3d.mean())
    v3 = float(s3d.var(ddof=1))
    out = []
    for series, is_lo in ((lower, True), (upper, False)):
        v4 = float(series.var(ddof=1))
        a, c = v4 / n4, v3 / n3
        se = math.sqrt(a + c)
        df = (a + c) ** 2 / (a * a / (n4 - 1) + c * c / (n3 - 1))
        t = _t_crit(df)
        d = float(series.mean()) - m3
        out.append(d - t * se if is_lo else d + t * se)
    return out[0], out[1]


def gate_a(ci: tuple[float, float]) -> bool:
    return ci[1] < -EPS_DET


def gate_b(wci: tuple[float, float]) -> str:
    lo, hi = wci
    if -EPS_REP < lo and hi < EPS_REP:
        return "REPLICATED"
    if lo > EPS_REP or hi < -EPS_REP:
        return "DISCORDANT"
    return "B-INCONCLUSIVE"


def gate_null(ci: tuple[float, float]) -> bool:
    return -EPS_DET < ci[0] and ci[1] < EPS_DET


def stage_verdict(ci: tuple[float, float],
                  wci: tuple[float, float]) -> tuple[str, str]:
    a = gate_a(ci)
    b = gate_b(wci)
    if a and b == "REPLICATED":
        return "CONFIRMED", b
    if a:
        return "DETECTED-NOT-REPLICATED", b
    if gate_null(ci):
        return "NO-DETECTION", b
    return "INCONCLUSIVE", b


# ---------------------------------------------------------------------
# Frozen power certification (doc section 7)
# ---------------------------------------------------------------------

#: (kind, theta_S4, theta_S3) per table row; anchor filled at call.
_POWER_ROWS = ("A@anchor", "B@diff0", "B@diff-kappa", "B@diff+kappa",
               "null@0")


def _power_source() -> tuple[np.ndarray, float]:
    d3 = s3_block()
    m3 = float(d3.mean())
    sd3 = float(d3.std(ddof=1))
    return (d3 - m3) * (SD_CONS / sd3), m3


def power_row(j: int, b: int = POWER_B) -> int:
    """Successes for row j of the frozen table: both blocks resampled
    independently on stream [781, 4840, j], S4 drawn before S3."""

    src, m3 = _power_source()
    n3 = len(src)
    kind, th4, th3 = [
        ("A", m3, None), ("B", m3, m3), ("B", m3 - KAPPA, m3),
        ("B", m3 + KAPPA, m3), ("null", 0.0, None)][j]
    rng = np.random.default_rng([*POWER_SEED, j])
    ok = 0
    for _ in range(b):
        x4 = src[rng.integers(0, n3, N_READINGS)] + th4
        if kind == "A":
            ok += gate_a(identified_ci(x4, x4))
        elif kind == "null":
            ok += gate_null(identified_ci(x4, x4))
        else:
            x3 = src[rng.integers(0, n3, n3)] + th3
            ok += gate_b(welch_identified_ci(x4, x4, x3)) == "REPLICATED"
    return ok


def cp_lower_all_success(n: int, conf: float = 0.95) -> float:
    """Exact Clopper-Pearson lower bound for n successes out of n."""

    return ((1.0 - conf) / 2.0) ** (1.0 / n)


def negative_controls() -> dict:
    """The frozen falsifiability battery (doc section 7): constructed
    blocks on stream [781, 4840, 99] for which each gate must FAIL or
    the A-only stage path must be reached."""

    src, m3 = _power_source()
    n3 = len(src)
    rng = np.random.default_rng([*POWER_SEED, 99])

    def block(theta, n=N_READINGS):
        return src[rng.integers(0, n3, n)] + theta

    x_at = block(-EPS_DET)
    x_in = block(-0.5 * EPS_DET)
    x_off = block(m3 + 3 * EPS_REP)
    x3_ref = block(m3, n3)
    x_large = block(2 * EPS_DET)
    x_aonly = block(m3 + 3 * EPS_REP)
    aonly_ci = identified_ci(x_aonly, x_aonly)
    aonly_wci = welch_identified_ci(x_aonly, x_aonly, x3_ref)
    verdict, b_label = stage_verdict(aonly_ci, aonly_wci)
    return {
        "gate_a_at_threshold": gate_a(identified_ci(x_at, x_at)),
        "gate_a_inside": gate_a(identified_ci(x_in, x_in)),
        "gate_b_offband": gate_b(welch_identified_ci(x_off, x_off, x3_ref)),
        "null_at_2eps": gate_null(identified_ci(x_large, x_large)),
        "a_only_stage": verdict,
        "a_only_b_label": b_label,
    }


# ---------------------------------------------------------------------
# The campaign
# ---------------------------------------------------------------------


def _git_state() -> dict:
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=_REPO,
                         capture_output=True, text=True,
                         check=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=_REPO,
                           capture_output=True, text=True,
                           check=True).stdout.strip()
    return {"rev": rev, "dirty": bool(dirty)}


def main() -> None:
    smoke = "--smoke" in sys.argv
    n_readings = 3 if smoke else N_READINGS
    e_n = 40 if smoke else E_N

    verify_freeze("entry")
    if smoke:
        seed = SMOKE_SEED
    else:
        seed = assert_fresh_scalar("s4_campaign")
    state_start = _git_state()
    if state_start["dirty"] and not smoke:
        raise SystemExit(
            "S4: refusing a campaign run from a dirty tree "
            f"(rev {state_start['rev']}); the rule requires the exact "
            "freeze content.")

    rng = np.random.default_rng(seed)
    fm_lo = np.empty(n_readings)
    fm_hi = np.empty(n_readings)
    f0 = np.empty(n_readings)
    counts = np.empty(n_readings, dtype=int)
    ambiguous = escalated = 0
    start = time.perf_counter()
    for k in range(n_readings):
        n = int(rng.poisson(e_n))
        if n < 2:
            raise SystemExit(f"S4: degenerate sprinkle n={n}")
        counts[k] = n
        fm_lo[k], fm_hi[k], f0[k], amb, esc = reading(rng, n)
        ambiguous += amb
        escalated += esc
        if (k + 1) % 10 == 0 or k == 0 or (k + 1) == n_readings:
            dt = time.perf_counter() - start
            eta = dt / (k + 1) * (n_readings - k - 1)
            print(f"  s4 {k + 1}/{n_readings}  "
                  f"delta so far {np.mean(fm_lo[:k + 1] - f0[:k + 1]):+.6f}"
                  f"  elapsed {dt / 60:.1f} min, eta {eta / 60:.1f} min",
                  flush=True)
    runtime = time.perf_counter() - start

    state_end = _git_state()
    if not smoke and state_end != state_start:
        raise SystemExit(
            "S4: refusing to write -- the tree changed during the run "
            f"({state_start} -> {state_end}).")
    if not smoke:
        verify_freeze("exit")

    delta_lower = fm_lo - f0
    delta_upper = fm_hi - f0
    s3d = s3_block()
    ci = identified_ci(delta_lower, delta_upper)
    wci = welch_identified_ci(delta_lower, delta_upper, s3d)
    verdict, b_label = stage_verdict(ci, wci)
    sentences = []
    if gate_a(ci):
        sentences.append(SENTENCES["GATE_A"])
    sentences.append(SENTENCES[b_label])
    if verdict in ("DETECTED-NOT-REPLICATED", "NO-DETECTION",
                   "INCONCLUSIVE"):
        sentences.append(SENTENCES[verdict])

    result = {
        "stage": "S4 Schwarzschild C1-paired confirmation",
        "rule": "docs/prereg/p14_s4_schwarzschild_c1.md",
        "run_kind": "fresh_observation" if not smoke else "smoke",
        "params": {"e_n": e_n, "n_readings": n_readings, "seed": seed,
                   "smoke": smoke},
        "margins": {"eps_det": EPS_DET, "eps_rep": EPS_REP,
                    "kappa": KAPPA, "sd_cons": SD_CONS},
        "event_counts": {"mean": float(np.mean(counts)),
                         "min": int(counts.min()),
                         "max": int(counts.max()),
                         "per_reading": [int(c) for c in counts]},
        "identified_ci95": list(ci),
        "welch_identified_ci95": list(wci),
        "gate_a": gate_a(ci),
        "gate_b": b_label,
        "verdict": verdict,
        "sentences": sentences,
        "delta_lower": {"mean": float(delta_lower.mean()),
                        "sd": float(delta_lower.std(ddof=1)),
                        "per_reading": [float(v) for v in delta_lower]},
        "delta_upper": {"mean": float(delta_upper.mean()),
                        "sd": float(delta_upper.std(ddof=1)),
                        "per_reading": [float(v) for v in delta_upper]},
        "f_flat": {"mean": float(f0.mean()),
                   "per_reading": [float(v) for v in f0]},
        "f_schwarzschild_lower": {"mean": float(fm_lo.mean()),
                                  "per_reading": [float(v) for v in fm_lo]},
        "f_schwarzschild_upper": {"mean": float(fm_hi.mean()),
                                  "per_reading": [float(v) for v in fm_hi]},
        "ambiguity": {"ambiguous": ambiguous, "escalated": escalated},
        "runtime_seconds": runtime,
        "code": {"start": state_start, "end": state_end},
    }

    print(f"\nS4: identified CI95 [{ci[0]:+.6f}, {ci[1]:+.6f}]  "
          f"welch [{wci[0]:+.6f}, {wci[1]:+.6f}]  "
          f"gate_a {gate_a(ci)}  gate_b {b_label}  -> {verdict}  "
          f"ambiguous {ambiguous}, escalated {escalated}  "
          f"({runtime / 3600:.2f} h)")

    if smoke:
        print("smoke run -- artifact NOT written")
        return
    # newline="\n": the artifact must be LF in the WORKING TREE too,
    # not only in the blob after git normalizes -- a CRLF working copy
    # is what broke the S3-artifact digest across platforms (PR #57
    # R2: the storage-boundary rule applies at write time).
    _ARTIFACT.write_text(json.dumps(result, indent=2, ensure_ascii=False)
                         + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {_ARTIFACT}")


if __name__ == "__main__":
    main()
