"""S6 M=1.4 ambiguity pilot -- the frozen runner.

EXECUTION FOLLOWS THE QUOTA-DELEGATED S6 RULING: after this freeze
merges and the exact-checkout preflight passes, the pilot executes
ONCE from the freeze branch head under the frozen gates -- and NEVER
twice. No rerun, no seed reuse, no sample extension, no cap raise.

WHAT THE PILOT MEASURES. One fresh stream (`s6_m14_pilot`) draws
n = 6,024,777 events from the SAME rung box measure the M = 1.4
count campaign will sprinkle (r^2 dr on the rung box, cos(psi)
uniform to the rung cap, phi, t in [0, dt_rung]). For each event the
full tri-state `causal_relation` at m = 1.4 is asked about both legs
(p -> x, then x -> q, early exit); the event is AMBIGUOUS iff any
asked leg returns None. The ambiguous count k out of the fixed n is
the whole observation.

THE FROZEN DECISION RULE IS GENERAL IN k -- the O5 pilot's rule with
this rung's constants (n sized at 96-bit from the ACTUAL certified
endpoints so the k = 1 boundary keeps a positive margin):

    p_upper  = exact one-sided Clopper-Pearson upper bound at
               confidence 1 - alpha_pilot, from (k, n);
    lambda_U = A_provisional * SCALE_rung * p_upper;
    tail     = exact P(Poisson(lambda_U) > U_max), a finite sum;
    verdict  = FEASIBLE      if tail <= tail_budget
               INCONCLUSIVE  otherwise.

FIXED n; an INCONCLUSIVE pilot stops this RUNG for a PI ruling and
indicts nothing (the other rungs proceed independently -- verdict
separation).

BUDGET. max_calls = 2*n + 6,537: the hard per-event bound (two legs)
plus the deterministic full-wrapper G3a preflight, which runs BEFORE
the reservation and is charged to the same budget. The caps are
frozen: RAISING THEM DURING OR AFTER THE RUN IS FORBIDDEN.

ORDER. static preflight (digests, environment lock, clean tree at the
exact approved SHA, artifact absence, fresh-seed assertion, retained
refs, namespace probe) -> metered G3a wrapper preflight -> reservation
claim on `refs/o5pilot/reservation` (the claim RETURNS its object; the
attempt is nonce-unique; any uncertain push outcome is SEED POSSIBLY
SPENT, fail-closed) -> RNG construction -> fixed-n chunked scan with
an atomic non-verdict checkpoint per chunk -> reservation re-read ->
write-once publication. Every failure path files a write-once incident
that names its `failure_point` and preserves the partial tallies.

Run (after the freeze merges, from a clean exact checkout of the
freeze branch head):

    python experiments/oracle/s6_m14_amb_pilot.py --preflight \
        --freeze-rev <full 40-hex freeze branch head>
    python experiments/oracle/s6_m14_amb_pilot.py \
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
import o4b_budget  # noqa: E402
import o4b_g3a as g3a  # noqa: E402
import o4b_meter as meter  # noqa: E402
import s1_schwarzschild_cost as s1  # noqa: E402
import s6_m14_reservation as reservation  # noqa: E402
import s6_rungs as s6  # noqa: E402
from probe_seed_ledger import assert_fresh_scalar  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
_MANIFEST = (_REPO / "docs" / "prereg"
             / "p14_s6_m14_pilot_freeze_manifest.json")
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_s6_m14_pilot.json"
_INCIDENT = _REPO / "docs" / "prereg" / "p14_s6_m14_pilot_incident.json"
_CHECKPOINT = (_REPO / "docs" / "prereg"
               / "p14_s6_m14_pilot_checkpoint.json")

#: The deterministic full-wrapper preflight's exact metered call
#: count, frozen by the O4b sizing artifact and re-asserted at run
#: time after the preflight completes.
G3A_PREFLIGHT_CALLS = 6_517

#: The frozen configuration (delegation ruling). Caps are part of the
#: freeze: no auto-raise, ever.
FROZEN = {
    "n_events": 6_024_777,
    "u_max": 30,
    "alpha_pilot": 0.01,
    "tail_budget": 0.001,
    "a_provisional": 830.0,
    "seed_name": "s6_m14_pilot",
    "tol": s1.DEFAULT_TOL,
    "chunk": 65_536,
    "max_calls": 12_056_071,        # = 2*n_events + G3A_PREFLIGHT_CALLS
    "max_wall_s": 86_400.0,
}
assert FROZEN["max_calls"] == 2 * FROZEN["n_events"] + G3A_PREFLIGHT_CALLS

#: The rung geometry, frozen as literals and re-derived at import
#: through the one exact path -- a drift in either direction refuses
#: the import (the ladder discipline).
_M = 1.4
_RUNG = {
    "dt": 9.070565190742672,
    "r_lo": 11.367185829820095,
    "r_hi": 18.7053462469963,
    "psi_max": 0.623439673163968,
    "scale": 18141.08004658374,
}
_g = s6.rung_geometry(_M)
if (_g["dt"] != _RUNG["dt"] or _g["r_lo"] != _RUNG["r_lo"]
        or _g["r_hi"] != _RUNG["r_hi"]
        or _g["psi_max"] != _RUNG["psi_max"]
        or _g["scale"] != _RUNG["scale"]):
    raise SystemExit(
        f"s6_m14_pilot: frozen rung literals drifted from the one "
        f"exact path: {_RUNG} != derived")
_GEOMETRY = _g
del _g

#: The anchors, exactly the rung's.
_P = np.array([0.0, 12.0, 0.0, 0.0])
_Q = np.array([_RUNG["dt"], 18.0, 0.0, 0.0])


# ------------------------------------------------ the frozen rule
#
# THE DECIDING TAIL ARITHMETIC RUNS AT 96-BIT MPFR PRECISION (review
# PR #83 R1). The k = 1 boundary sits ~7.7e-10 of tail below the
# budget, while a double-precision binomial CDF at n ~ 1e7 carries an
# lgamma cancellation floor near 1e-8 relative -- enough, in
# principle, for the same frozen inputs to flip the verdict across
# hosts, because the environment lock pins gmpy2/MPFR but not the
# platform libm behind math.lgamma. gmpy2 at a fixed precision IS
# part of the locked instrument, so the deciding arithmetic lives
# there; the double-precision engine survives only in the contract
# tests, as one of the independent cross-checks. One step is NOT
# MPFR (post-merge adversarial review): `decide()` forms lambda_U by
# a double multiply of cp_upper's float return before re-promoting to
# 96 bits -- verified bit-identical in `tail` against a never-
# downcast 96-bit chain at the tightest frozen boundary (k = 1), and
# left as-is because the executed pilot's published decision must
# stay reproducible by this exact code path.

_PREC = 96


def _ctx():
    ctx = gmpy2.get_context()
    ctx.precision = _PREC
    return ctx


def binom_cdf(k: int, n: int, p) -> gmpy2.mpfr:
    """P(Bin(n, p) <= k) at 96-bit precision (k small)."""

    _ctx()
    p = gmpy2.mpfr(p)
    if p <= 0:
        return gmpy2.mpfr(1)
    if p >= 1:
        return gmpy2.mpfr(0)
    lp, lq = gmpy2.log(p), gmpy2.log1p(-p)
    lg_n1 = gmpy2.lgamma(gmpy2.mpfr(n + 1))[0]
    s = gmpy2.mpfr(0)
    for i in range(k + 1):
        s += gmpy2.exp(lg_n1
                       - gmpy2.lgamma(gmpy2.mpfr(i + 1))[0]
                       - gmpy2.lgamma(gmpy2.mpfr(n - i + 1))[0]
                       + i * lp + (n - i) * lq)
    return min(gmpy2.mpfr(1), s)


def cp_upper(k: int, n: int, alpha: float) -> float:
    """Exact one-sided Clopper-Pearson upper bound at 96-bit: the p
    with P(Bin(n, p) <= k) = alpha; closed form at k = 0."""

    if k < 0 or n <= 0 or not 0.0 < alpha < 1.0:
        raise ValueError(f"cp_upper({k}, {n}, {alpha})")
    if k >= n:
        return 1.0
    _ctx()
    a = gmpy2.mpfr(alpha)
    if k == 0:
        return float(1 - gmpy2.exp(gmpy2.log(a) / n))
    lo, hi = gmpy2.mpfr(k) / n, gmpy2.mpfr(1)
    for _ in range(300):
        mid = (lo + hi) / 2
        if binom_cdf(k, n, mid) > a:
            lo = mid
        else:
            hi = mid
        if hi - lo <= gmpy2.mpfr("1e-25") * hi:
            break
    return float((lo + hi) / 2)


def pois_tail_gt(u_max: int, lam: float) -> float:
    """Exact P(Poisson(lam) > u_max) at 96-bit: 1 minus a finite pmf
    sum."""

    if lam <= 0.0:
        return 0.0
    _ctx()
    lam_m = gmpy2.mpfr(lam)
    s = gmpy2.mpfr(0)
    for i in range(u_max + 1):
        s += gmpy2.exp(i * gmpy2.log(lam_m) - lam_m
                       - gmpy2.lgamma(gmpy2.mpfr(i + 1))[0])
    return float(max(gmpy2.mpfr(0), 1 - s))


def decide(k: int, cfg: dict = FROZEN) -> dict:
    """The frozen verdict, computed by the GENERAL rule for whatever
    k the scan produced."""

    p_up = cp_upper(k, cfg["n_events"], cfg["alpha_pilot"])
    lam_u = cfg["a_provisional"] * _RUNG["scale"] * p_up
    tail = pois_tail_gt(cfg["u_max"], lam_u)
    feasible = tail <= cfg["tail_budget"]
    return {
        "k_ambiguous": k,
        "n_events": cfg["n_events"],
        "p_upper": p_up,
        "p_upper_confidence": 1.0 - cfg["alpha_pilot"],
        "lambda_u": lam_u,
        "u_max": cfg["u_max"],
        "tail_p_u_gt_u_max": tail,
        "tail_budget": cfg["tail_budget"],
        "verdict": "FEASIBLE" if feasible else "INCONCLUSIVE",
        "rule": ("exact one-sided Clopper-Pearson upper at "
                 "1 - alpha_pilot from (k, n); lambda_U = "
                 "A_provisional * SCALE * p_upper; exact Poisson "
                 "tail P(U > U_max) <= tail_budget -- general in k, "
                 "fixed n, no extension"),
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
            f"pilot runs from a clean exact checkout only")
    if state["rev"] != want:
        raise SystemExit(
            f"preflight: HEAD is {state['rev']}, not the approved "
            f"freeze SHA {want} -- a clean tree with matching digests "
            f"is NOT enough; the pilot runs from the exact approved "
            f"commit only")
    for path in (_ARTIFACT, _INCIDENT):
        if path.exists():
            raise SystemExit(
                f"preflight: {path.name} already exists -- the pilot "
                f"is write-once and the first attempt stands")
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

def draw_events(rng, k: int) -> np.ndarray:
    """`k` events from the 4D box measure, exactly: r by the inverse
    CDF of r^2 dr on [R_LO, R_HI]; cos(psi) uniform on
    [cos(PSI_MAX), 1]; phi uniform on [0, 2pi) -- drawn to consume
    the stream by the literal measure even though the predicate is
    axisymmetric in phi; t uniform on [0, DT]. Returns (k, 4) rows
    (t, r, psi, phi)."""

    u = rng.random((k, 4))
    r3_lo, r3_hi = _RUNG["r_lo"] ** 3, _RUNG["r_hi"] ** 3
    r = np.cbrt(r3_lo + u[:, 0] * (r3_hi - r3_lo))
    cmin = math.cos(_RUNG["psi_max"])
    psi = np.arccos(cmin + u[:, 1] * (1.0 - cmin))
    phi = u[:, 2] * (2.0 * math.pi)
    t = u[:, 3] * _RUNG["dt"]
    return np.column_stack([t, r, psi, phi])


def event_ambiguous(t: float, r: float, psi: float,
                    tol: float) -> tuple[bool, int]:
    """The tri-state indicator with early exit: (ambiguous?, calls).
    None on ANY asked leg makes the event ambiguous; the second leg
    is only asked if the first decided."""

    x = np.array([t, r, psi, 0.0])
    leg1 = s1.causal_relation(_P, x, _M, tol)
    if leg1 is None:
        return True, 1
    leg2 = s1.causal_relation(x, _Q, _M, tol)
    return leg2 is None, 2


# ------------------------------------------------ checkpoints

_CK_REQUIRED = ("freeze_sha", "manifest_digest", "seed", "rng_stream",
                "rng_position", "events_done", "k_ambiguous", "calls",
                "chunk_provenance", "budget")


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
        "why": ("progress record, not a result: the fixed-n scan has "
                "not finished, so no count here is a verdict"),
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
    (review PR #83 R1, the o4b R26-R29 pattern): whether THIS call
    published is a fact the incident boundary must be able to read,
    because an interrupt one line later must not file an incident
    beside a published result. Cleanup cannot raise."""

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
    file, not from a flag (the o4b R29 rule): the artifact's
    reservation object is this attempt's nonce-unique claim commit,
    so a file carrying it can only have been written by this run --
    and unlike a flag, the answer survives any interrupt because it
    is recomputed from what is on disk."""

    if claim_object is None or not _ARTIFACT.exists():
        return False
    try:
        rec = json.loads(_ARTIFACT.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return rec.get("reservation", {}).get("object") == claim_object


# ------------------------------------------------ the run

def run_scan(rng, budget, cfg: dict, on_chunk=None) -> dict:
    """The fixed-n chunked scan. Every chunk commits an atomic
    checkpoint before the next begins; `on_chunk` receives the
    running tallies (tests use it)."""

    n, chunk = cfg["n_events"], cfg["chunk"]
    k_amb, done, calls = 0, 0, 0
    with meter.metered(budget):
        while done < n:
            size = min(chunk, n - done)
            state_before = rng.bit_generator.state
            events = draw_events(rng, size)
            consumed = 0
            for t, r, psi, _phi in events:
                amb, c = event_ambiguous(float(t), float(r),
                                         float(psi), cfg["tol"])
                calls += c
                if amb:
                    k_amb += 1
                done += 1
                consumed += 1
            yield_state = {
                "freeze_sha": cfg["_freeze_sha"],
                "manifest_digest": cfg["_digest"],
                "seed": cfg["_seed"],
                "rng_stream": cfg["seed_name"],
                "rng_position": rng.bit_generator.state,
                "events_done": done,
                "k_ambiguous": k_amb,
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
    return {"k_ambiguous": k_amb, "events_done": done, "calls": calls}


def file_incident(failure_point: str, reason: str, preserved: dict,
                  context: dict) -> None:
    record = {
        "kind": "incident",
        "run_kind": "s6_m14_pilot",
        "failure_point": failure_point,
        "termination_reason": reason,
        "verdict": None,
        "why_no_verdict": ("fail-closed: the fixed-n scan or its "
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
              f" seed {checks['seed']} fresh, retained refs verified, "
              "namespace writable, no artifacts. Nothing was "
              "observed.")
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
            # distinct from claim-time uncertainty (review PR #85 R1):
            # a confirmed claim whose pre-publication ownership
            # re-check could not complete stays a confirmed claim
            "publish_reverify_uncertainty": claimed.get(
                "reverify_uncertain"),
            "seed_spent": (claimed["object"] is not None
                           or bool(claimed.get("uncertain"))),
            "environment": _environment(),
        }

    tallies = {"k_ambiguous": 0, "events_done": 0, "calls": 0}
    try:
        # metered deterministic full-wrapper preflight, BEFORE the
        # reservation: a failure here has spent nothing
        budget.enter("g3a_preflight")
        with meter.metered(budget):
            wrapper = g3a.run_preflight(FROZEN["tol"], m=_M,
                                        geometry=_GEOMETRY)
        if not wrapper["passed"]:
            raise RuntimeError(
                f"full-wrapper preflight failed "
                f"{wrapper['failed_conditions']}")
        if budget.reserved != G3A_PREFLIGHT_CALLS:
            raise RuntimeError(
                f"wrapper preflight spent {budget.reserved} metered "
                f"calls, not the frozen {G3A_PREFLIGHT_CALLS}")

        # the point of no return -- the claim RETURNS its object.
        # claim()'s exception contract (review PR #83 R2): every
        # SystemExit it raises is PRE-PUSH by construction (retained
        # refs, authority reads, an already-held ref, commit-tree) --
        # this run attempted no write and spent nothing, so it may
        # not file the pilot's GLOBAL incident: in the two-checkout
        # race, the loser's incident would sit beside the winner's
        # legitimate result and poison the write-once provenance.
        # Everything AFTER the push is ClaimUncertain, which does file.
        budget.enter("reservation")
        try:
            claimed["object"] = reservation.claim({
                "campaign": "s6_m14_pilot",
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

        def on_chunk(state: dict) -> None:
            tallies.update(k_ambiguous=state["k_ambiguous"],
                           events_done=state["events_done"],
                           calls=state["calls"])

        scan = run_scan(rng, budget, cfg, on_chunk)
        tallies.update(scan)

        verdict = decide(scan["k_ambiguous"], FROZEN)

        budget.enter("publish")
        try:
            reservation.verify_still_held(claimed["object"])
        except reservation.ClaimUncertain as uncertain:
            # the publish-time re-verify can be uncertain too (post-
            # merge adversarial review): the scan is complete and the
            # seed long spent, so the incident must carry the
            # uncertainty record -- under its OWN key (review PR #85
            # R1): the claim itself already succeeded, and writing
            # this as claim-time uncertainty would make the incident
            # deny a claim that is in fact confirmed
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
            "run_kind": "s6_m14_pilot",
            "stage": ("O5 ambiguity pilot: fixed-n tri-state "
                      "indicator scan on the campaign box measure"),
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
        # THIS RUN's publication commit is decided two ways (review
        # PR #83 R1, the o4b R26-R29 pattern): the receipt stamped the
        # instant os.link returned, and the artifact asked whether it
        # names this run's nonce-unique claim object. An interrupt one
        # line after the link -- or a SystemExit out of a write-once
        # refusal where the existing file turns out to be OURS -- must
        # never file an incident beside a published result.
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
        # failure this program has already outlawed. A write-once
        # refusal against a FOREIGN artifact lands here too, so a
        # fully scanned run that could not publish still leaves its
        # incident.
        # the STAGE names the failure point (post-merge adversarial
        # review): a ClaimUncertain out of the publish-time re-verify
        # comes after the whole scan, and labelling it
        # `reservation_claim` would read as if no event was processed
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
    print(f"result: k={scan['k_ambiguous']} of "
          f"{scan['events_done']:,} -> {verdict['verdict']} "
          f"(p_upper={verdict['p_upper']:.4e}, "
          f"lambda_U={verdict['lambda_u']:.2f}, "
          f"tail={verdict['tail_p_u_gt_u_max']:.4e})")
    print(f"artifact: {_ARTIFACT}")


# ------------------------------------------------ the manifest

MANIFEST_NOTE = (
    "raw sha256; all paths .gitattributes eol=lf pinned. The rung "
    "pilot inherits the O5 pilot's instrument surface (predicate, "
    "wrapper preflight at the rung mass+geometry, budget/meter, "
    "reservation mechanics) with this rung's one-exact-path "
    "constants; the rung oracle artifact is pinned because n and A "
    "were sized from its actual endpoints, and p14_o3_volume.json "
    "because o4_sizing imports it (the O5 maintenance lesson). This "
    "manifest cannot certify itself: the run also demands the exact "
    "freeze branch head via --freeze-rev."
)

PROTOCOL_SURFACE = (
    "experiments/oracle/s6_m14_amb_pilot.py",
    "experiments/oracle/s6_m14_reservation.py",
    "experiments/oracle/s6_rungs.py",
    "experiments/oracle/o4b_reservation.py",
    "experiments/oracle/o4b_budget.py",
    "experiments/oracle/o4b_meter.py",
    "experiments/oracle/o4b_g3a.py",
    "experiments/oracle/o4_g3_redesign.py",
    "experiments/oracle/o4_sizing.py",
    "experiments/oracle/certified_flight_time.py",
    "experiments/oracle/certified_interval.py",
    "experiments/positive_control/s1_schwarzschild_cost.py",
    "experiments/positive_control/probe_seed_ledger.py",
    "docs/prereg/p14_s6_m14_pilot.md",
    "docs/prereg/p14_s6_m14_volume.json",
    "docs/prereg/p14_o3_volume.json",
    "pyproject.toml",
)


def build_manifest() -> dict:
    return {
        "stage": ("O5 ambiguity-pilot freeze manifest "
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
