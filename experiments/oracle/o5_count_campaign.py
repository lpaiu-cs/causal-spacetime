"""O5 Poisson-count campaign -- the frozen runner (freeze, NO results).

EXECUTION IS NOT APPROVED BY THIS FREEZE. Per the delegation ruling,
the arc STOPS after this freeze merges: the campaign executes only
after ONE integrated PI approval, once, from the clean exact checkout
of the freeze branch head (the merge commit's second parent) -- and
NEVER twice. No rerun, no seed reuse, no sample extension, no cap
raise.

WHAT THE CAMPAIGN MEASURES -- the prediction-anchored count. A
genuine 4D Poisson process of frozen intensity A per unit certified
4-volume is sprinkled into the L4 box (the SAME box measure the
ambiguity pilot scanned): N ~ Poisson(A * SCALE) proposal points, each
asked the tri-state membership question "is x causally between the
frozen anchors?" (leg1 P -> x, then leg2 x -> Q, early exit). The
estimand is V_op - V_true, where V_op = A^-1 E[K_op] is the measure of
the wrapper-accepted region; K ~ Poisson(A * V_true) holds ONLY under
the defect-free alternative, which is what the power sentence
conditions on -- the verdict itself assumes nothing.

THE FROZEN VERDICT RULE -- general in (K_certain, U_amb), 96-bit
end to end:

    L        = exact Garwood lower bound at alpha/2 from K_certain;
    U        = exact Garwood upper bound at alpha/2 from
               K_certain + U_amb          (ONE outer interval: every
               ambiguous point may or may not be a member);
    C        = [L / A, U / A]             (volume units);
    D        = [C_lo - V'_hi, C_hi - V'_lo]   (identified discrepancy
               against the O3' certification, interval arithmetic);
    B        = tau * V_ref'               (the equivalence band);
    verdict  = CONCORDANT     if D  is contained in [-B, B]
               DISCORDANT     if D  is disjoint from [-B, B]
               INCONCLUSIVE   otherwise (published as-is).

TRI-STATE AMBIGUITY AND THE PILOT. The certified U ceiling comes from
the ambiguity pilot (docs/prereg/p14_o5_amb_pilot.json, FEASIBLE at
k = 0): P(U > 30) <= 1e-3 at 99% confidence under the PILOT's
ambiguity rule. The campaign's rule is A FORTIORI NARROWER: the pilot
counted an event ambiguous when ANY asked leg returned None with leg2
always asked after a decided leg1, while the campaign early-exits on
leg1 False (the point is out regardless of leg2), so
{campaign-ambiguous} is a subset of {pilot-ambiguous} pointwise on the
same measure and the ceiling transfers.

THE POWER SENTENCE (design, frozen here; conditions on the defect-free
alternative): acceptance [52,752, 53,980] on K_true; worst endpoint
power 0.930794 at V'_hi; at U = 30 the guaranteed window keeps
0.911855; joint completion-and-concordance
    >= 0.911855 - 0.001 (pilot tail budget) - 0.010 (pilot alpha)
       - exp(-179,855)  (exact Chernoff on P(2N > 60,000,000))
    =  0.900855 >= 0.90.
Wall-cap risk is NOT included in that bound (the wall cap is a start
rule, not a probability statement). Three independent engines
reproduced the acceptance integers and the endpoint powers to 6
decimals before this freeze; the contract tests re-derive them from
THIS module's own 96-bit rule.

BUDGET. max_calls = 60,000,000 + 6,537: the scan's hard cap (P(2N >
cap) <= exp(-179,855) under the alternative) plus the deterministic
full-wrapper G3a preflight, which runs BEFORE the reservation and is
charged to the same budget. E[calls] <= 2 * A * SCALE = 53,667,552
(~11.5 h at the projected price; the pilot ran ~19% slower, ~13.7 h
-- both far under the 24 h wall). The caps are frozen: RAISING THEM
DURING OR AFTER THE RUN IS FORBIDDEN; a cap termination is an
incident published as-is.

ORDER. static preflight (digests, environment lock, clean tree at the
exact approved SHA, artifact absence, fresh-seed assertion, retained
refs x3, namespace probe) -> metered G3a wrapper preflight ->
reservation claim on `refs/o5/reservation` (the claim RETURNS its
object; any uncertain push outcome is SEED POSSIBLY SPENT,
fail-closed) -> RNG construction -> N ~ Poisson(A * SCALE) drawn once
-> chunked scan with an atomic non-verdict checkpoint per chunk ->
reservation re-read -> write-once publication. Every failure path
files a write-once incident that names its `failure_point` and
preserves the partial tallies.

Run (ONLY after the separate integrated execution approval, from a
clean exact checkout of the freeze branch head):

    python experiments/oracle/o5_count_campaign.py --preflight \
        --freeze-rev <full 40-hex freeze branch head>
    python experiments/oracle/o5_count_campaign.py \
        --freeze-rev <full 40-hex freeze branch head>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parents[0] / "positive_control"))

import gmpy2  # noqa: E402
import numpy as np  # noqa: E402
import o4_sizing as sz  # noqa: E402
import o4b_budget  # noqa: E402
import o4b_g3a as g3a  # noqa: E402
import o4b_meter as meter  # noqa: E402
import o5_reservation as reservation  # noqa: E402
import s1_schwarzschild_cost as s1  # noqa: E402
from probe_seed_ledger import assert_fresh_scalar  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
_MANIFEST = (_REPO / "docs" / "prereg"
             / "p14_o5_count_freeze_manifest.json")
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_o5_count.json"
_INCIDENT = _REPO / "docs" / "prereg" / "p14_o5_count_incident.json"
_CHECKPOINT = (_REPO / "docs" / "prereg"
               / "p14_o5_count_checkpoint.json")

#: The deterministic full-wrapper preflight's exact metered call
#: count, frozen by the O4b sizing artifact and re-asserted at run
#: time after the preflight completes.
G3A_PREFLIGHT_CALLS = 6_537

#: The frozen configuration. The oracle endpoints are the ACTUAL O3'
#: artifact values (docs/prereg/p14_o3p_volume.json, target-met at
#: ratio 0.0049995), adopted verbatim; V_ref' is that artifact's
#: recommendation (the standalone midpoint), whose adoption the O3'
#: results PR explicitly deferred to THIS freeze. Caps are part of
#: the freeze: no auto-raise, ever.
FROZEN = {
    "a_intensity": 940.0,
    "tau": 0.025,
    "alpha": 0.05,
    "v_lo": 56.49295916680885,
    "v_hi": 57.06066656647874,
    "v_ref": 56.776812866643795,
    "u_power_ceiling": 30,
    "seed_name": "o5_campaign",
    "tol": s1.DEFAULT_TOL,
    "chunk": 65_536,
    "scan_call_cap": 60_000_000,
    "max_calls": 60_006_537,
    "max_wall_s": 86_400.0,
}
assert FROZEN["max_calls"] == (FROZEN["scan_call_cap"]
                               + G3A_PREFLIGHT_CALLS)
assert FROZEN["v_lo"] < FROZEN["v_ref"] < FROZEN["v_hi"]

#: The design acceptance window on K_true (documentation of the power
#: sentence, not the verdict rule -- the verdict is the interval
#: containment in `decide`). The contract tests prove the two agree:
#: at U_amb = 0, decide(k, 0) is CONCORDANT exactly for k in this
#: window.
ACCEPTANCE = (52_752, 53_980)

#: The anchors, exactly the pilot's and the O4b campaign's.
_P = np.array([0.0, sz.R_IN, 0.0, 0.0])
_Q = np.array([sz.DT, sz.R_OUT, 0.0, 0.0])


def _import_gate() -> None:
    """The two artifacts this freeze quotes must SAY what the frozen
    literals claim -- checked at import, refused before anything runs.
    The manifest pins their bytes; this pins their meaning."""

    o3p = json.loads((_REPO / "docs" / "prereg"
                      / "p14_o3p_volume.json").read_text(
                          encoding="utf-8"))
    r = o3p["result"]
    if (r["v_lo"] != FROZEN["v_lo"] or r["v_hi"] != FROZEN["v_hi"]
            or r["status"] != "target-met"):
        raise SystemExit(
            f"o5_count_campaign: the frozen oracle endpoints "
            f"({FROZEN['v_lo']}, {FROZEN['v_hi']}) do not match the "
            f"O3' artifact ({r['v_lo']}, {r['v_hi']}, {r['status']})")
    v_ref = o3p["v_ref_prime_recommendation"]["value"]
    if v_ref != FROZEN["v_ref"]:
        raise SystemExit(
            f"o5_count_campaign: frozen v_ref {FROZEN['v_ref']} is "
            f"not the O3' recommendation {v_ref}")
    pilot = json.loads((_REPO / "docs" / "prereg"
                        / "p14_o5_amb_pilot.json").read_text(
                            encoding="utf-8"))
    d = pilot["decision"]
    if (d["verdict"] != "FEASIBLE"
            or d["u_max"] != FROZEN["u_power_ceiling"]):
        raise SystemExit(
            f"o5_count_campaign: the ambiguity pilot's decision "
            f"({d['verdict']} at u_max={d['u_max']}) does not certify "
            f"the frozen U ceiling {FROZEN['u_power_ceiling']} -- the "
            f"power sentence has no basis")


_import_gate()


# ------------------------------------------------ the frozen rule
#
# THE DECIDING ARITHMETIC RUNS AT 96-BIT MPFR PRECISION, END TO END
# (the pilot's PR #83 R1 lesson, applied without the pilot's one
# double step): the Garwood bounds, the division by A, the band and
# every comparison in `decide` stay mpfr until serialization. gmpy2
# at a fixed precision is part of the locked instrument; the
# double-precision engines survive only in the contract tests, as
# independent cross-checks.

_PREC = 96


def _ctx():
    ctx = gmpy2.get_context()
    ctx.precision = _PREC
    return ctx


def pois_cdf(k: int, lam) -> gmpy2.mpfr:
    """P(Poisson(lam) <= k) at 96-bit, via the regularized upper
    incomplete gamma: Q(k+1, lam)."""

    _ctx()
    if k < 0:
        return gmpy2.mpfr(0)
    lam = gmpy2.mpfr(lam)
    if lam <= 0:
        return gmpy2.mpfr(1)
    a = gmpy2.mpfr(k + 1)
    return gmpy2.gamma_inc(a, lam) / gmpy2.gamma(a)


def garwood_low(k: int, alpha: float) -> gmpy2.mpfr:
    """The exact lower Garwood bound: the lam with
    P(Poisson(lam) >= k) = alpha/2; 0 at k = 0."""

    _ctx()
    if k < 0:
        raise ValueError(f"garwood_low({k})")
    if k == 0:
        return gmpy2.mpfr(0)
    a2 = gmpy2.mpfr(alpha) / 2
    lo = gmpy2.mpfr(0)
    hi = gmpy2.mpfr(k + 10.0 * math.sqrt(k) + 20.0)
    for _ in range(300):
        mid = (lo + hi) / 2
        if 1 - pois_cdf(k - 1, mid) < a2:
            lo = mid
        else:
            hi = mid
        if hi - lo <= gmpy2.mpfr("1e-25") * max(gmpy2.mpfr(1), hi):
            break
    return (lo + hi) / 2


def garwood_up(k: int, alpha: float) -> gmpy2.mpfr:
    """The exact upper Garwood bound: the lam with
    P(Poisson(lam) <= k) = alpha/2."""

    _ctx()
    if k < 0:
        raise ValueError(f"garwood_up({k})")
    a2 = gmpy2.mpfr(alpha) / 2
    lo = gmpy2.mpfr(max(0.0, k - 1.0))
    hi = gmpy2.mpfr(k + 10.0 * math.sqrt(k + 1.0) + 30.0)
    for _ in range(300):
        mid = (lo + hi) / 2
        if pois_cdf(k, mid) > a2:
            lo = mid
        else:
            hi = mid
        if hi - lo <= gmpy2.mpfr("1e-25") * max(gmpy2.mpfr(1), hi):
            break
    return (lo + hi) / 2


def decide(k_certain: int, u_amb: int, cfg: dict = FROZEN) -> dict:
    """The frozen verdict, computed by the GENERAL rule for whatever
    (K_certain, U_amb) the scan produced -- one outer Garwood
    interval, identified-discrepancy containment, three ways out."""

    if k_certain < 0 or u_amb < 0:
        raise ValueError(f"decide({k_certain}, {u_amb})")
    _ctx()
    a = gmpy2.mpfr(cfg["a_intensity"])
    lam_lo = garwood_low(k_certain, cfg["alpha"])
    lam_hi = garwood_up(k_certain + u_amb, cfg["alpha"])
    c_lo, c_hi = lam_lo / a, lam_hi / a
    v_lo, v_hi = gmpy2.mpfr(cfg["v_lo"]), gmpy2.mpfr(cfg["v_hi"])
    band = gmpy2.mpfr(cfg["tau"]) * gmpy2.mpfr(cfg["v_ref"])
    d_lo, d_hi = c_lo - v_hi, c_hi - v_lo
    contained = d_lo >= -band and d_hi <= band
    disjoint = d_hi < -band or d_lo > band
    verdict = ("CONCORDANT" if contained
               else "DISCORDANT" if disjoint
               else "INCONCLUSIVE")
    return {
        "k_certain": k_certain,
        "u_ambiguous": u_amb,
        "alpha": cfg["alpha"],
        "garwood_low_counts": float(lam_lo),
        "garwood_up_counts": float(lam_hi),
        "c_lo_volume": float(c_lo),
        "c_hi_volume": float(c_hi),
        "d_lo": float(d_lo),
        "d_hi": float(d_hi),
        "band": float(band),
        "oracle": {"v_lo": cfg["v_lo"], "v_hi": cfg["v_hi"],
                   "v_ref": cfg["v_ref"], "tau": cfg["tau"]},
        "verdict": verdict,
        "rule": ("ONE outer Garwood interval at alpha/2 per side: "
                 "L(K_certain), U(K_certain + U_amb); C = counts / A; "
                 "D = [C_lo - V'_hi, C_hi - V'_lo] against "
                 "[-tau*V_ref', +tau*V_ref']; contained -> CONCORDANT,"
                 " disjoint -> DISCORDANT, else INCONCLUSIVE -- "
                 "general in (K_certain, U_amb), 96-bit end to end, "
                 "published as computed"),
    }


# ------------------------------------------------ the freeze gates

def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _environment() -> dict:
    return {
        "python": ".".join(platform.python_version_tuple()[:2]),
        "gmpy2": gmpy2.version(),
        "mpfr": gmpy2.mpfr_version(),
        "gmp": gmpy2.mp_version(),
        "numpy": np.__version__,
    }


def verify_freeze(stage: str) -> dict:
    manifest = json.loads(_MANIFEST.read_text(encoding="utf-8"))
    bad = [rel for rel, want in manifest["files"].items()
           if _sha256(_REPO / rel) != want]
    if bad:
        raise SystemExit(
            f"freeze verification ({stage}): digest mismatch on "
            f"{bad} -- refusing to run on a drifted protocol surface")
    env, locked = _environment(), manifest["environment"]
    drift = {kk: (locked[kk], env.get(kk)) for kk in locked
             if env.get(kk) != locked[kk]}
    if drift:
        raise SystemExit(
            f"freeze verification ({stage}): environment drift "
            f"{drift} -- the solver and rounding instrument are part "
            f"of the frozen apparatus")
    return manifest


def _git_state() -> dict:
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=_REPO,
                         check=True, capture_output=True,
                         text=True).stdout.strip()
    porcelain = subprocess.run(["git", "status", "--porcelain"],
                               cwd=_REPO, check=True,
                               capture_output=True,
                               text=True).stdout.splitlines()
    mine = {p.relative_to(_REPO).as_posix()
            for p in (_ARTIFACT, _INCIDENT, _CHECKPOINT)}
    dirt = [ln for ln in porcelain
            if ln.strip() and ln[3:].strip().strip('"') not in mine]
    return {"rev": rev, "dirty": bool(dirt), "dirt": dirt}


def require_full_sha(freeze_rev: str) -> str:
    """The manifest cannot certify itself (the O4b/O3' rule): the run
    demands the full 40-hex SHA of the freeze branch head."""

    want = (freeze_rev or "").strip().lower()
    if len(want) != 40 or any(c not in "0123456789abcdef"
                              for c in want):
        raise SystemExit(
            f"--freeze-rev must be the full 40-hex SHA of the freeze "
            f"branch head, got {freeze_rev!r} -- a short or malformed "
            f"rev is refused before anything runs")
    return want


def preflight(expected_rev: str) -> dict:
    """Everything before a single seed is touched, observing nothing.
    `--preflight` and the real execution pass through this same
    check."""

    want = require_full_sha(expected_rev)
    manifest = verify_freeze("preflight")
    state = _git_state()
    if state["dirty"]:
        raise SystemExit(
            f"preflight: working tree is dirty {state['dirt']} -- the "
            f"campaign runs from a clean exact checkout only")
    if state["rev"] != want:
        raise SystemExit(
            f"preflight: HEAD is {state['rev']}, not the approved "
            f"freeze SHA {want} -- a clean tree with matching digests "
            f"is NOT enough; the campaign runs from the exact "
            f"approved commit only")
    for path in (_ARTIFACT, _INCIDENT):
        if path.exists():
            raise SystemExit(
                f"preflight: {path.name} already exists -- the "
                f"campaign is write-once and the first attempt stands")
    if _CHECKPOINT.exists():
        raise SystemExit(
            f"preflight: {_CHECKPOINT.name} already exists -- it is a "
            f"partial from an earlier attempt; a new run must not "
            f"start next to one")
    reservation.verify_retained()
    if (already := reservation.held()) is not None:
        raise SystemExit(
            f"preflight: {reservation.REF} is already held by "
            f"{already} -- the stream was opened by some checkout and "
            f"is spent regardless of whether a result was published")
    reservation.probe_namespace()
    seed = assert_fresh_scalar(FROZEN["seed_name"])
    return {"manifest": manifest, "git": state, "seed": seed}


# ------------------------------------------------ the sampler

def draw_points(rng, k: int) -> np.ndarray:
    """`k` proposal points from the 4D box measure, exactly the
    pilot's sampler: r by the inverse CDF of r^2 dr on [R_LO, R_HI];
    cos(psi) uniform on [cos(PSI_MAX), 1]; phi uniform on [0, 2pi) --
    drawn to consume the stream by the literal measure even though the
    predicate is axisymmetric in phi; t uniform on [0, DT]. Returns
    (k, 4) rows (t, r, psi, phi)."""

    u = rng.random((k, 4))
    r3_lo, r3_hi = sz.R_LO ** 3, sz.R_HI ** 3
    r = np.cbrt(r3_lo + u[:, 0] * (r3_hi - r3_lo))
    cmin = math.cos(sz.PSI_MAX)
    psi = np.arccos(cmin + u[:, 1] * (1.0 - cmin))
    phi = u[:, 2] * (2.0 * math.pi)
    t = u[:, 3] * sz.DT
    return np.column_stack([t, r, psi, phi])


def point_membership(t: float, r: float, psi: float,
                     tol: float) -> tuple[str, int]:
    """The tri-state membership with early exit:
    ('in' | 'out' | 'ambiguous', calls). leg1 None -> ambiguous;
    leg1 False -> out (leg2 never asked); leg1 True -> ask leg2.
    NARROWER than the pilot's ambiguity rule (which always asked leg2
    after a decided leg1), so the pilot's certified U ceiling
    transfers a fortiori."""

    x = np.array([t, r, psi, 0.0])
    leg1 = s1.causal_relation(_P, x, s1.M, tol)
    if leg1 is None:
        return "ambiguous", 1
    if not leg1:
        return "out", 1
    leg2 = s1.causal_relation(x, _Q, s1.M, tol)
    if leg2 is None:
        return "ambiguous", 2
    return ("in" if leg2 else "out"), 2


# ------------------------------------------------ checkpoints

_CK_REQUIRED = ("freeze_sha", "manifest_digest", "seed", "rng_stream",
                "rng_position", "n_total", "points_done", "k_certain",
                "u_ambiguous", "calls", "chunk_provenance", "budget")


def write_checkpoint(payload: dict) -> None:
    """Atomic (same-directory temp + fsync + os.replace), stamped
    partial/non_verdict by THIS writer so no caller can claim
    otherwise -- the O4b checkpoint discipline, with the stream named
    (rng_stream) and every required key enforced at the write."""

    missing = [kk for kk in _CK_REQUIRED if kk not in payload]
    if missing:
        raise ValueError(f"checkpoint missing {missing}")
    record = {
        "kind": "checkpoint", "stage": "chunk", **payload,
        "partial": True, "non_verdict": True,
        "why": ("progress record, not a result: the sprinkled scan "
                "has not finished, so no count here is a verdict"),
    }
    _CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(_CHECKPOINT.parent),
                               prefix=_CHECKPOINT.name + ".",
                               suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as f:
            json.dump(record, f, ensure_ascii=False, indent=1)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, _CHECKPOINT)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def publish_write_once(path: Path, payload: dict,
                       receipt: dict | None = None) -> None:
    """Atomic no-clobber publication (fsynced temp + os.link).

    `receipt["committed"]` is stamped the instant os.link returns
    (the pilot's PR #83 R1, the o4b R26-R29 pattern): whether THIS
    call published is a fact the incident boundary must be able to
    read, because an interrupt one line later must not file an
    incident beside a published result. Cleanup cannot raise."""

    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as f:
            json.dump(payload, f, ensure_ascii=False, indent=1)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        try:
            os.link(tmp, path)
        except FileExistsError:
            raise SystemExit(
                f"publish: {path.name} already exists -- write-once; "
                f"the first observation stands") from None
        if receipt is not None:
            receipt["committed"] = True
    finally:
        try:
            tmp.unlink(missing_ok=True)
        except BaseException:                   # noqa: BLE001
            pass


def artifact_names_this_run(claim_object: str | None) -> bool:
    """Was the published artifact written by THIS run? Read from the
    file, not from a flag (the o4b R29 rule)."""

    if claim_object is None or not _ARTIFACT.exists():
        return False
    try:
        rec = json.loads(_ARTIFACT.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return rec.get("reservation", {}).get("object") == claim_object


# ------------------------------------------------ the run

def run_scan(rng, budget, cfg: dict, on_chunk=None,
             live: dict | None = None) -> dict:
    """The sprinkled scan: N ~ Poisson(A * SCALE) drawn ONCE from the
    stream, then N points in chunks. Every chunk commits an atomic
    checkpoint before the next begins; `on_chunk` receives the running
    tallies (tests use it).

    `live`, when given, is updated IN PLACE from the moment N is
    drawn and after EVERY completed observation (review PR #87 R1):
    a failed run spends the seed permanently with no rerun, so the
    incident must preserve the actual last observation point -- not
    the last successful chunk, and never `n_total: null` after N was
    already drawn under a claimed reservation. Checkpoints stay
    chunk-granular; the incident does not."""

    if live is None:
        live = {}
    n_total = int(rng.poisson(cfg["a_intensity"] * sz.SCALE))
    k_certain, u_amb, done, calls = 0, 0, 0, 0
    live.update(n_total=n_total, points_done=0, k_certain=0,
                u_ambiguous=0, calls=0)
    with meter.metered(budget):
        while done < n_total:
            size = min(cfg["chunk"], n_total - done)
            state_before = rng.bit_generator.state
            points = draw_points(rng, size)
            consumed = 0
            for t, r, psi, _phi in points:
                status, c = point_membership(float(t), float(r),
                                             float(psi), cfg["tol"])
                calls += c
                if status == "in":
                    k_certain += 1
                elif status == "ambiguous":
                    u_amb += 1
                done += 1
                consumed += 1
                live["points_done"] = done
                live["k_certain"] = k_certain
                live["u_ambiguous"] = u_amb
                live["calls"] = calls
            yield_state = {
                "freeze_sha": cfg["_freeze_sha"],
                "manifest_digest": cfg["_digest"],
                "seed": cfg["_seed"],
                "rng_stream": cfg["seed_name"],
                "rng_position": rng.bit_generator.state,
                "n_total": n_total,
                "points_done": done,
                "k_certain": k_certain,
                "u_ambiguous": u_amb,
                "calls": calls,
                "chunk_provenance": {
                    "state_before_draw": state_before,
                    "size": size, "consumed": consumed,
                },
                "budget": budget.state(),
            }
            write_checkpoint(yield_state)
            if on_chunk is not None:
                on_chunk(yield_state)
    return {"n_total": n_total, "points_done": done,
            "k_certain": k_certain, "u_ambiguous": u_amb,
            "calls": calls}


def file_incident(failure_point: str, reason: str, preserved: dict,
                  context: dict) -> None:
    record = {
        "kind": "incident",
        "run_kind": "o5_count",
        "failure_point": failure_point,
        "termination_reason": reason,
        "verdict": None,
        "why_no_verdict": ("fail-closed: the sprinkled scan or its "
                           "publication did not complete, so the "
                           "frozen rule was never applied"),
        "preserved": preserved,
        **context,
    }
    publish_write_once(_INCIDENT, record)


def main() -> None:
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("--preflight", action="store_true",
                        help="run every pre-execution check and "
                             "observe nothing")
    parser.add_argument("--freeze-rev", default="",
                        help="REQUIRED for preflight and execution: "
                             "the full 40-hex SHA of the freeze "
                             "branch head")
    parser.add_argument("--write-manifest", action="store_true",
                        help="regenerate the freeze manifest; a "
                             "freeze action, never part of a run")
    args = parser.parse_args()

    if args.write_manifest:
        print(f"wrote {write_manifest()}")
        return

    checks = preflight(args.freeze_rev)
    if args.preflight:
        print("preflight PASS: freeze digests, environment lock, "
              f"clean tree at the approved {checks['git']['rev'][:9]},"
              f" seed {checks['seed']} fresh, retained refs x3 "
              "verified, namespace writable, no artifacts. Nothing "
              "was observed.")
        return

    start = checks["git"]                  # the VERIFIED baseline
    approved = start["rev"]
    digest = _sha256(_MANIFEST)
    budget = o4b_budget.Budget(FROZEN["max_calls"],
                               FROZEN["max_wall_s"])
    cfg = dict(FROZEN, _freeze_sha=approved, _digest=digest,
               _seed=checks["seed"])

    claimed: dict = {"object": None}
    receipt: dict = {}
    t0 = time.perf_counter()

    def context() -> dict:
        return {
            "freeze_sha": approved, "manifest_digest": digest,
            "seed": {FROZEN["seed_name"]: checks["seed"]},
            "reservation_claimed": (
                "uncertain" if claimed.get("uncertain")
                else claimed["object"] is not None),
            "reservation_object": claimed["object"],
            "reservation_uncertainty": claimed.get("uncertain"),
            # distinct from claim-time uncertainty (the pilot's PR #85
            # R1): a confirmed claim whose pre-publication ownership
            # re-check could not complete stays a confirmed claim
            "publish_reverify_uncertainty": claimed.get(
                "reverify_uncertain"),
            "seed_spent": (claimed["object"] is not None
                           or bool(claimed.get("uncertain"))),
            "environment": _environment(),
        }

    tallies = {"n_total": None, "points_done": 0, "k_certain": 0,
               "u_ambiguous": 0, "calls": 0}
    try:
        # metered deterministic full-wrapper preflight, BEFORE the
        # reservation: a failure here has spent nothing
        budget.enter("g3a_preflight")
        with meter.metered(budget):
            wrapper = g3a.run_preflight(FROZEN["tol"])
        if not wrapper["passed"]:
            raise RuntimeError(
                f"full-wrapper preflight failed "
                f"({wrapper.get('failures')})")
        if budget.reserved != G3A_PREFLIGHT_CALLS:
            raise RuntimeError(
                f"wrapper preflight spent {budget.reserved} metered "
                f"calls, not the frozen {G3A_PREFLIGHT_CALLS}")

        # the point of no return -- the claim RETURNS its object.
        # claim()'s exception contract (the pilot's PR #83 R2 + PR
        # #85): every SystemExit it raises is PRE-PUSH by construction
        # -- this run attempted no write and spent nothing, so it may
        # not file the campaign's GLOBAL incident. Everything AFTER
        # the push is ClaimUncertain, which does file.
        budget.enter("reservation")
        try:
            claimed["object"] = reservation.claim({
                "campaign": "o5_count",
                "freeze_rev": approved,
                "manifest_sha256": digest,
                "seed": {FROZEN["seed_name"]: checks["seed"]},
            })
        except reservation.ClaimUncertain as uncertain:
            claimed["uncertain"] = uncertain.as_record()
            raise
        except SystemExit:
            claimed["clean_refusal"] = True
            raise

        budget.enter("scan")
        rng = np.random.default_rng(checks["seed"])

        # `tallies` IS the live run state (review PR #87 R1): run_scan
        # updates it per observation, so the incident boundary reads
        # the actual last observation point on any mid-chunk failure
        scan = run_scan(rng, budget, cfg, live=tallies)
        tallies.update(scan)

        verdict = decide(scan["k_certain"], scan["u_ambiguous"],
                         FROZEN)

        budget.enter("publish")
        try:
            reservation.verify_still_held(claimed["object"])
        except reservation.ClaimUncertain as uncertain:
            # the publish-time re-verify can be uncertain too; its
            # record lands under its OWN key (the pilot's PR #85 R1)
            claimed["reverify_uncertain"] = uncertain.as_record()
            raise
        end = _git_state()
        if end["rev"] != approved or end["dirty"]:
            raise RuntimeError(
                f"exit lineage: HEAD is {end['rev']} "
                f"(dirty={end['dirty']}), not the approved "
                f"{approved}")
        result = {
            "kind": "results",
            "run_kind": "o5_count",
            "stage": ("O5 Poisson-count campaign: sprinkled tri-state "
                      "membership count against the O3' certified "
                      "volume"),
            "frozen_config": {kk: v for kk, v in FROZEN.items()},
            "g3a_preflight": {
                "passed": True,
                "metered_calls": G3A_PREFLIGHT_CALLS,
            },
            "decision": verdict,
            "scan": scan,
            "reservation": {
                "ref": reservation.REF,
                "authority": reservation.CANONICAL_AUTHORITY,
                "object": claimed["object"],
                "seed_spent": True,
                "verified_at_exit": True,
            },
            "seed": {FROZEN["seed_name"]: checks["seed"]},
            "freeze_sha": approved,
            "manifest_digest": digest,
            "code": {"start": start, "end": end},
            "budget": budget.state(),
            "environment": _environment(),
            "total_wall_s": time.perf_counter() - t0,
        }
        publish_write_once(_ARTIFACT, result, receipt=receipt)
    except o4b_budget.CapReached as cap:
        file_incident(
            "scan", f"completion budget exhausted ({cap.reason})",
            dict(tallies, budget=budget.state()), context())
        raise SystemExit(
            f"cap fired ({cap.reason}); incident written to "
            f"{_INCIDENT.name}, no verdict published") from None
    except BaseException as exc:
        # THIS RUN's publication commit is decided two ways (the
        # pilot's PR #83 R1, the o4b R26-R29 pattern): the receipt
        # stamped the instant os.link returned, and the artifact asked
        # whether it names this run's nonce-unique claim object.
        if (receipt.get("committed")
                or artifact_names_this_run(claimed["object"])):
            raise SystemExit(
                f"{_ARTIFACT.name} was published; the failure came "
                f"after the commit ({type(exc).__name__}: {exc}) and "
                f"no incident is filed over a published result"
            ) from exc
        if claimed.get("clean_refusal"):
            # a pre-push claim refusal: nothing was written, nothing
            # was spent by THIS run, and the winning checkout's record
            # is the authoritative one -- propagate without an incident
            raise
        # KeyboardInterrupt and SystemExit included: once the
        # reservation is claimed the seed is spent, and a stop with
        # neither a result nor a record is the one observability
        # failure this program has already outlawed.
        # the STAGE names the failure point (the pilot's post-merge
        # review): a ClaimUncertain out of the publish-time re-verify
        # comes after the whole scan
        point = ("publish" if budget.stage == "publish"
                 else "reservation_claim"
                 if isinstance(exc, reservation.ClaimUncertain)
                 else budget.stage)
        file_incident(
            point, f"unhandled {type(exc).__name__}: {exc}",
            dict(tallies, budget=budget.state()), context())
        raise SystemExit(
            f"{point}: {type(exc).__name__}; incident written to "
            f"{_INCIDENT.name}, no verdict published") from exc
    # success logging OUTSIDE the incident boundary: a BrokenPipeError
    # here must not be able to file an incident over the published
    # result
    scan = tallies
    print(f"result: K_certain={scan['k_certain']:,} "
          f"U_amb={scan['u_ambiguous']} of N={scan['n_total']:,} -> "
          f"{verdict['verdict']} "
          f"(C=[{verdict['c_lo_volume']:.6f}, "
          f"{verdict['c_hi_volume']:.6f}], band={verdict['band']:.6f})")
    print(f"artifact: {_ARTIFACT}")


# ------------------------------------------------ the manifest

MANIFEST_NOTE = (
    "raw sha256; all paths .gitattributes eol=lf pinned. The campaign "
    "inherits the O4b instrument surface (predicate, wrapper "
    "preflight, budget/meter) and the box measure from o4_sizing; the "
    "reservation follows the corrected o5-pilot copy; the O3' "
    "artifact, the O3 artifact o4_sizing imports, and the ambiguity "
    "pilot's FEASIBLE result are pinned because the frozen endpoints, "
    "the import surface and the U ceiling quote them. This manifest "
    "cannot certify itself: the run also demands the exact freeze "
    "branch head via --freeze-rev."
)

PROTOCOL_SURFACE = (
    "experiments/oracle/o5_count_campaign.py",
    "experiments/oracle/o5_reservation.py",
    "experiments/oracle/o4b_reservation.py",
    "experiments/oracle/o4b_budget.py",
    "experiments/oracle/o4b_meter.py",
    "experiments/oracle/o4b_g3a.py",
    "experiments/oracle/o4_g3_redesign.py",
    "experiments/oracle/o4_sizing.py",
    "experiments/oracle/o4_volume_audit.py",
    "experiments/positive_control/s1_schwarzschild_cost.py",
    "experiments/positive_control/probe_seed_ledger.py",
    "docs/prereg/p14_o5_count.md",
    "docs/prereg/p14_o3p_volume.json",
    "docs/prereg/p14_o3_volume.json",
    "docs/prereg/p14_o5_amb_pilot.json",
    "pyproject.toml",
)


def build_manifest() -> dict:
    return {
        "stage": ("O5 Poisson-count campaign freeze manifest "
                  "(content-addressed protocol surface)"),
        "note": MANIFEST_NOTE,
        "files": {rel: _sha256(_REPO / rel)
                  for rel in PROTOCOL_SURFACE},
        "environment": _environment(),
    }


def write_manifest() -> Path:
    _MANIFEST.write_text(
        json.dumps(build_manifest(), ensure_ascii=False, indent=1)
        + "\n", encoding="utf-8", newline="\n")
    return _MANIFEST


if __name__ == "__main__":
    main()
