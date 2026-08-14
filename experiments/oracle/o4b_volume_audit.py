"""O4b: the S1 volume-response audit, re-run under the redesigned G3.

EXECUTION IS NOT APPROVED. This module is part of the freeze. The
campaign needs a separate approval naming an exact freeze SHA, as
every campaign in this program does, and the freeze approval is never
the execution approval.

The estimand, the oracle and tau are O4's and are inherited, not
re-frozen: `o4_sizing` is digest-pinned and this module reads it. The
sampler is inherited the same way -- `o4_volume_audit._draw` and
`._ell` are IMPORTED rather than copied, because the sampler is part
of what the audit compares and a second copy of it could drift into
being a different instrument while still passing every digest check.

WHAT IS NEW IS THE ORDER AND WHAT SURVIVES.

    static preflight  ->  metered G3a  ->  G3a PASS
         ->  reservation claim  ->  the generator
         ->  G3b on the G1 stream's fixed prefix  ->  the rest of G1
         ->  G2

The claim sits AFTER G3a, not before it. The frozen rule is that a
G3a failure has spent no fresh seed, and the reservation ref is the
authority on whether a stream has been opened -- so a claim made and
then abandoned retires both seeds by policy, however carefully the
generator was left unconstructed. G3a is entirely deterministic (a
frozen case table and a solver), so it needs no stream, and its 6,537
calls are charged to the same budget object the campaign then
continues to spend.

O4 ran `G3a -> G1/G2 -> G3b`, spent twelve hours, found the
instrumentation broken at the end, and kept nothing. Three things
follow from that, and they are the substance of this module:

  * G3b is a precondition, so it runs before the bulk of G1.
  * G3b consumes the front of the G1 stream, so EVERY point it sees
    enters the G1 accumulator first and unconditionally -- see
    `o4b_stages.Prefix`, which is where that rule lives.
  * every stage boundary and every G1 chunk writes an atomic
    `partial` / `non_verdict` checkpoint, and every abort path writes
    a write-once incident carrying what was already paid for.

Run, after execution approval, from a clean checkout of the approved
commit:

    python experiments/oracle/o4b_volume_audit.py --preflight \
        --freeze-rev <40-hex>
    python experiments/oracle/o4b_volume_audit.py --freeze-rev <40-hex>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import empirical_bernstein as eb  # noqa: E402
import exact_binomial as xb  # noqa: E402
import gmpy2  # noqa: E402
import o4_g3_redesign as g3  # noqa: E402
import o4_sizing as sz  # noqa: E402
import o4_volume_audit as o4  # noqa: E402
import o4b_budget  # noqa: E402
import o4b_g3a as g3a  # noqa: E402
import o4b_g3b as g3b  # noqa: E402
import o4b_meter as meter  # noqa: E402
import o4b_reservation as reservation  # noqa: E402
import o4b_sizing as sizing  # noqa: E402
import o4b_stages as stages  # noqa: E402
import s1_schwarzschild_cost as s1  # noqa: E402
from probe_seed_ledger import assert_fresh_scalar  # noqa: E402

_MANIFEST = _REPO / "docs" / "prereg" / "p14_o4b_freeze_manifest.json"
_ARTIFACT = _REPO / "docs" / "prereg" / "p14_o4b_results.json"
_INCIDENT = _REPO / "docs" / "prereg" / "p14_o4b_incident.json"
_CHECKPOINT = _REPO / "docs" / "prereg" / "p14_o4b_checkpoint.json"

#: Write-once outputs the campaign itself creates. They are results,
#: not protocol drift, so the clean-tree check ignores exactly these.
#: The CHECKPOINT is not among them: it is overwritten by design and
#: is not a result, which is why it is a separate file from both.
_WRITE_ONCE = (_ARTIFACT, _INCIDENT)

#: The frozen configuration. Sizes and the oracle come from the pinned
#: `o4_sizing`; the G3 quantities come from the redesign module. They
#: are named here so a desynchronised edit fails loudly rather than
#: quietly changing what was frozen.
FROZEN = {
    "tau": sz.TAU,
    "delta_g1_per_side": sz.DELTA_G1,
    "alpha_g2_per_end": sz.ALPHA_G2,
    "leak_budget_frac": sz.LEAK_BUDGET,
    "n_g1": sz.N_G1,
    "n_g2": sz.g2_points(),
    "n_avail": g3.N_AVAIL,
    "k_g3b": g3.K_G3B,
    "eta": g3.ETA,
    "tol": s1.DEFAULT_TOL,
    "max_calls": sizing.MAX_CALLS,
    "max_wall_s": sizing.MAX_WALL_S,
    "chunk": 4096,
}


# ---------------------------------------------------------- the freeze

def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def environment() -> dict:
    return {
        "python": ".".join(platform.python_version_tuple()[:2]),
        "gmpy2": gmpy2.version(),
        "mpfr": gmpy2.mpfr_version(),
        "gmp": gmpy2.mp_version(),
        "numpy": np.__version__,
    }


#: The content-addressed protocol surface. Listed here rather than in
#: a script that was run once, so the manifest is REGENERABLE and a
#: contract test can compare the committed file against this list --
#: a surface file that quietly stopped being pinned would otherwise be
#: invisible.
#:
#: `o4_sizing` and `o4_volume_audit` are INHERITED: O4b reads the
#: estimand, the oracle, tau and the sampler from them instead of
#: re-freezing them, so an edit to either is a different campaign and
#: is pinned as such.
PROTOCOL_SURFACE = (
    "experiments/oracle/certified_interval.py",
    "experiments/oracle/certified_flight_time.py",
    "experiments/oracle/empirical_bernstein.py",
    "experiments/oracle/exact_binomial.py",
    "experiments/oracle/o4_sizing.py",
    "experiments/oracle/o4_volume_audit.py",
    "experiments/oracle/o4_g3_redesign.py",
    "experiments/oracle/o4b_sizing.py",
    "experiments/oracle/o4b_budget.py",
    "experiments/oracle/o4b_meter.py",
    "experiments/oracle/o4b_g3a.py",
    "experiments/oracle/o4b_g3b.py",
    "experiments/oracle/o4b_checkpoint.py",
    "experiments/oracle/o4b_stages.py",
    "experiments/oracle/o4b_reservation.py",
    "experiments/oracle/o4b_volume_audit.py",
    "experiments/positive_control/s1_schwarzschild_cost.py",
    "experiments/positive_control/probe_seed_ledger.py",
    "docs/prereg/p14_o4_volume_audit.md",
    "docs/prereg/p14_o4_g3_redesign.md",
    "docs/prereg/p14_o4_g3_prereg_reopen.md",
    "docs/prereg/p14_o4b_sizing.json",
    "docs/prereg/p14_o3_volume.json",
    "docs/prereg/p14_o4_replay_diagnostic.json",
    "pyproject.toml",
)

MANIFEST_NOTE = (
    "raw sha256; all paths are .gitattributes eol=lf pinned. "
    "s1_schwarzschild_cost.py is the INSTRUMENT UNDER AUDIT and is "
    "pinned so the audited artefact cannot drift. o4_sizing.py and "
    "o4_volume_audit.py are INHERITED from the O4 freeze -- O4b reads "
    "the estimand, the oracle, tau and the sampler from them rather "
    "than re-freezing them, so any edit to either is a different "
    "campaign. p14_o4_replay_diagnostic.json is the committed census "
    "the cost projection reads. This manifest cannot certify itself: "
    "a commit that edits a protocol file and re-pins the manifest in "
    "the same commit passes every digest check, which is why the "
    "runner also demands the exact approved freeze SHA."
)


def build_manifest() -> dict:
    """The manifest this tree would produce. Compared against the
    committed one by a contract test."""

    return {
        "stage": ("O4b freeze manifest (content-addressed protocol "
                  "surface)"),
        "note": MANIFEST_NOTE,
        "files": {rel: _sha256(_REPO / rel)
                  for rel in PROTOCOL_SURFACE},
        "environment": environment(),
    }


def write_manifest() -> Path:
    _MANIFEST.write_text(
        json.dumps(build_manifest(), ensure_ascii=False, indent=1)
        + "\n", encoding="utf-8", newline="\n")
    return _MANIFEST


def verify_freeze(stage: str) -> dict:
    manifest = json.loads(_MANIFEST.read_text(encoding="utf-8"))
    bad = [rel for rel, want in manifest["files"].items()
           if _sha256(_REPO / rel) != want]
    if bad:
        raise SystemExit(
            f"freeze verification ({stage}): digest mismatch on {bad} "
            f"-- refusing to run on a drifted protocol surface")
    env, locked = environment(), manifest["environment"]
    drift = {k: (locked[k], env.get(k)) for k in locked
             if env.get(k) != locked[k]}
    if drift:
        raise SystemExit(
            f"freeze verification ({stage}): environment drift "
            f"{drift} -- the solver and the rounding instrument are "
            f"part of the frozen apparatus")
    return manifest


def git_state() -> dict:
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=_REPO,
                         check=True, capture_output=True,
                         text=True).stdout.strip()
    porcelain = subprocess.run(["git", "status", "--porcelain"],
                               cwd=_REPO, check=True,
                               capture_output=True,
                               text=True).stdout.splitlines()
    mine = {p.relative_to(_REPO).as_posix()
            for p in (*_WRITE_ONCE, _CHECKPOINT)}
    dirt = [ln for ln in porcelain
            if ln.strip() and ln[3:].strip().strip('"') not in mine]
    return {"rev": rev, "dirty": bool(dirt), "dirt": dirt}


def verify_rev(stage: str, freeze_rev: str, state: dict) -> None:
    """The manifest cannot certify itself. A commit after the approved
    freeze that edits a protocol file and re-pins the manifest in the
    SAME commit passes every digest check, so the campaign also
    demands the exact 40-hex SHA named in the execution approval."""

    want = (freeze_rev or "").strip().lower()
    if len(want) != 40 or any(c not in "0123456789abcdef" for c in want):
        raise SystemExit(
            f"{stage}: --freeze-rev must be the full 40-hex commit "
            f"named in the execution approval, got {freeze_rev!r}")
    if state["rev"].lower() != want:
        raise SystemExit(
            f"{stage}: HEAD is {state['rev']} but the approved freeze "
            f"is {want} -- refusing to spend campaign seeds on a "
            f"commit the approval does not name")


def preflight(freeze_rev: str) -> dict:
    """Everything that must hold before a single seed is touched."""

    manifest = verify_freeze("preflight")
    state = git_state()
    verify_rev("preflight", freeze_rev, state)
    if state["dirty"]:
        raise SystemExit(
            f"preflight: working tree is dirty {state['dirt']} -- the "
            f"campaign runs from a clean exact checkout only")
    for path in _WRITE_ONCE:
        if path.exists():
            raise SystemExit(
                f"preflight: {path.name} already exists -- a campaign "
                f"was already started on these streams; the audit is "
                f"write-once and the first attempt stands")
    # ...and so must the checkpoint (review R22). It is not write-once
    # -- the run overwrites it deliberately, and `git_state` excuses it
    # from the clean-tree check for exactly that reason -- but a
    # checkpoint present BEFORE the first draw is a leftover from some
    # earlier attempt. Left alone, an early G3a failure would file a
    # fresh incident beside another run's partial statistics, and a
    # reader has nothing in either file to tell them apart.
    if _CHECKPOINT.exists():
        raise SystemExit(
            f"preflight: {_CHECKPOINT.name} already exists -- it is a "
            f"partial from an earlier attempt, and a new run must not "
            f"start next to one; move it aside deliberately or "
            f"explain it in the incident record")
    reservation.verify_o4_ref_retained()
    if (already := reservation.held()) is not None:
        raise SystemExit(
            f"preflight: {reservation.REF} is already held by "
            f"{already} -- these streams were opened by some checkout "
            f"and are spent regardless of whether a result was "
            f"published")
    reservation.probe_namespace()
    seeds = {name: assert_fresh_scalar(name)
             for name in ("o4b_g1_audit", "o4b_g2_leakage")}
    return {"manifest": manifest, "git": state, "seeds": seeds}


def publish_write_once(path: Path, payload: dict,
                       receipt: dict | None = None) -> bool:
    """Atomic no-clobber publication: fsynced temp in the same
    directory, then `os.link`, which fails atomically if the
    destination exists.

    `receipt` is stamped `committed` the instant `os.link` returns,
    and it is how the caller learns whether THIS CALL published
    (review R28). The file merely existing does not mean that: the
    destination can be created by something else during the eleven
    hours the campaign runs, in which case `os.link` fails, this
    function reports "the first observation stands", and a caller
    reading `path.exists()` would conclude its own result had been
    published -- crediting someone else's file to this run and filing
    no incident for a campaign whose seeds are already spent.

    Returns True on a commit, so a caller that can use the return
    value does not need the receipt at all. The receipt exists for the
    caller that cannot: an interrupt delivered between the return and
    the assignment would leave a plain flag unset."""

    text = json.dumps(payload, ensure_ascii=False, indent=1) + "\n"
    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())
        try:
            os.link(tmp, path)
        except FileExistsError:
            raise SystemExit(
                f"publish: {path.name} already exists -- the audit is "
                f"write-once and the first observation stands"
            ) from None
        if receipt is not None:
            receipt["committed"] = True     # the commit, recorded here
    finally:
        # THE COMMIT IS `os.link`, AND NOTHING AFTER IT MAY UNDO THE
        # CLAIM THAT IT HAPPENED (review R27). Removing the temporary
        # is housekeeping; a failure or an interrupt here would
        # otherwise propagate as a publication failure over a result
        # that IS published, and the caller would file an incident
        # beside it -- a run with both a verdict and an incident,
        # which no path is allowed to produce.
        try:
            tmp.unlink(missing_ok=True)
        except BaseException:                   # noqa: BLE001
            pass
    return True


# ------------------------------------------------------- the G3 state

def solver_state(r: float, theta: float, tol: float) -> dict:
    """The PREDICATE's view of one cluster, in the predicate's own
    angle.

    THIS IS NOT THE G1 ESTIMAND'S VIEW, and merging the two was the
    third of the abort's four non-adaptive margins (review R9). G1
    integrates `L_S1 = [dt - T1(theta) - T2(theta)]_+` at the drawn
    `theta`; `causal_relation` never sees `theta`, it recovers
    `dpsi = acos(cos theta)` from the coordinates and calls the solver
    with THAT. The two agree to within an ulp of angle, and the abort
    happened because an ulp of angle is not an ulp of time -- the
    recovery error is unbounded in ulps near zero and reaches 4.14e-13
    at theta = 1e-5, which is the same order as eta.

    So a probe placed from `T(theta)` and judged at `t_min(dpsi)` is
    placed in one coordinate and read in another. The window is built
    here, in the coordinate that reads it.

    `err` is kept, unlike in `_ell`, because the probes are PLACED
    with it. It still never enters the inference: it is a
    Gauss-Legendre stopping heuristic, not a certified enclosure, and
    an interval built from it would not be conservative."""

    dpsi = g3.wrapper_dpsi(theta)
    d1: dict = {}
    d2: dict = {}
    t1, err1 = s1.flight_time(sz.R_IN, r, dpsi, s1.M, tol, d1)
    t2, err2 = s1.flight_time(r, sz.R_OUT, dpsi, s1.M, tol, d2)
    return {"t1": t1, "err1": err1, "t2": t2, "err2": err2,
            "theta": theta, "dpsi": dpsi,
            "recovery_shift": dpsi - theta,
            "family1": d1.get("family"), "family2": d2.get("family")}


def judge_cluster(state: dict, eta: float) -> dict:
    """Eligibility and the three probes, for one candidate.

    `fully_testable` is the conjunction the redesign froze: a strictly
    positive robust width, a lower probe that stays at `dt >= 0`, and
    all three probes reaching the realized margin."""

    lo, hi = state["t1"], sz.DT - state["t2"]
    w = g3.w_robust(hi - lo, state["err1"], state["err2"])
    if not w > 0.0:
        return {"eligible": False, "fully_testable": False,
                "reason": "w_robust <= 0", "w_robust": w}

    # BOTH outside probes place their leg's `dt` BELOW its `t_min`:
    # `above`/`below` name where the intermediate event sits relative
    # to the window, and starving either leg means giving it less time
    # than it needs. The predicate must answer False at both.
    inside = g3b.place_inside(0.0, sz.DT, state["t1"], state["err1"],
                              state["t2"], state["err2"], eta)
    above = g3a.place_row(state["t2"], state["err2"], eta, above=False)
    below = g3a.place_row(state["t1"], state["err1"], eta, above=False)
    reached = {"inside": inside["reached"], "outside_above":
               above["reached"], "outside_below": below["reached"]}
    if not below["reached"] and "short-circuit" in (
            below.get("why") or ""):
        reason = "lower probe reaches dt < 0"
    elif not all(reached.values()):
        reason = "construction-unavailable"
    else:
        reason = None
    return {"eligible": True, "fully_testable": reason is None,
            "reason": reason, "w_robust": w, "reached": reached,
            "probes": {"inside": inside, "outside_above": above,
                       "outside_below": below}}


# ----------------------------------------------------------- the run

def _statistics(acc) -> dict:
    """The running accumulator, safe at zero samples.

    A cap can fire before the first point is accumulated -- during
    G3a, or on the first G3b chunk -- and the accumulator raises on an
    empty mean. Letting that propagate would take out the incident
    that the cap path exists to write, which is the failure mode in
    miniature: an error while recording an error, and nothing kept."""

    if acc.n == 0:
        return {"n": 0, "mean_z": None, "var_z": None,
                "why_no_moments": "no samples accumulated yet"}
    return {"n": acc.n, "mean_z": acc.mean, "var_z": acc.var}


class Campaign:
    """The four stages in their frozen order, with the checkpoint and
    the preservation rules attached to the transitions."""

    def __init__(self, cfg: dict, seeds: dict, freeze_sha: str,
                 digest: str, checkpoint_path: Path = _CHECKPOINT,
                 draw=None, state_of=None, ell=None) -> None:
        self.cfg = cfg
        self.seeds = seeds
        self.freeze_sha = freeze_sha
        self.digest = digest
        self.checkpoint_path = Path(checkpoint_path)
        self.budget = o4b_budget.Budget(cfg["max_calls"],
                                        cfg["max_wall_s"])
        self.draw = draw if draw is not None else o4._draw
        # `ell` is G1's estimand at the drawn theta; `state_of` is the
        # predicate's view at the recovered dpsi. Two coordinates, two
        # functions, and the runner never substitutes one for the other.
        self.ell = ell if ell is not None else o4._ell
        self.state_of = (state_of if state_of is not None
                         else solver_state)
        self.g1 = eb.Accumulator()
        self.prefix = stages.Prefix(self.g1, cfg["n_avail"])
        self.rng = None
        self.reservation_object: str | None = None
        self.g1_result: dict | None = None
        self.g2_result: dict | None = None
        self.mismatches: list[dict] = []
        self.scanned = 0
        # The draw is BATCHED: one `draw` call consumes the generator
        # for a whole chunk before any point is processed. So the
        # generator's position and the number of points actually
        # committed move at different times, and a stop in the middle
        # leaves them disagreeing (review R13). Both ends of the chunk
        # are therefore recorded: the state BEFORE the draw, which a
        # resume can restore and redraw from, and how many of the
        # drawn points were consumed.
        self.chunk_start_rng = None
        self.chunk_consumed = 0
        self.chunk_size = 0
        #: The part of the last chunk G3b drew but did not consume.
        #: It is HANDED TO G1, not discarded (review R16): the RNG has
        #: already advanced past it, so throwing it away would make G1
        #: skip that stretch of the frozen stream and integrate a
        #: different sample.
        self.pending: list[tuple[float, float]] = []

    # -- checkpoints ---------------------------------------------

    def checkpoint(self, point: str, extra: dict | None = None) -> None:
        stages.write_checkpoint(
            self.checkpoint_path, point,
            freeze_sha=self.freeze_sha, digest=self.digest,
            seed=self.seeds["o4b_g1_audit"], rng=self.rng,
            samples=self.g1.n,
            statistics={**_statistics(self.g1), **(extra or {}),
                        "chunk": self._chunk_position()},
            budget=self.budget)

    def _chunk_position(self) -> dict:
        """Where the batched draw stands.

        `rng_position` alone is the END of the last chunk drawn, which
        is not where the committed points end. A resume restores
        `state_before_draw`, redraws the same `size` points and skips
        the first `consumed` -- so the stopping condition and the
        generator move on the same boundary."""

        return {
            "state_before_draw": self.chunk_start_rng,
            "size": self.chunk_size,
            "consumed": self.chunk_consumed,
            "why": ("the draw is batched, so rng_position is the end "
                    "of the chunk while the accumulator holds only "
                    "the consumed prefix of it; resuming from "
                    "rng_position alone would skip the unconsumed "
                    "tail and change the frozen stream"),
        }

    def preserved(self) -> dict:
        """What a failure hands to the incident. Always the G1 prefix:
        it is paid for and is valid independently of anything G3b
        found."""

        finished = {}
        if self.g1_result is not None:
            finished["g1"] = self.g1_result
        if self.g2_result is not None:
            finished["g2"] = self.g2_result
        return {
            "g1_partial": _statistics(self.g1),
            "availability": self.prefix.report(),
            "scanned_points": self.scanned,
            "budget": self.budget.state(),
            # gates that DID complete before the failure. Kept because
            # the sample is paid for and valid, and marked because a
            # completed gate in an incident is still not a verdict --
            # the run that would have carried it did not finish.
            "completed_gates": finished,
            "completed_gates_are_not_a_verdict": (
                "these gates ran to their frozen sample size, but the "
                "run did not publish, so no sentence rests on them"),
            "is_not_a_verdict": (
                "a preserved partial sample: the stopping rule did "
                "not fire, so no gate has a status"),
        }

    def _chunk(self, k: int):
        """One batched draw, with both ends of it recorded."""

        self.chunk_start_rng = self.rng.bit_generator.state
        self.chunk_consumed = 0
        self.chunk_size = k
        return self.draw(self.rng, k, sz.R_LO, sz.R_HI, sz.PSI_MAX)

    def _g3b_done(self) -> bool:
        return (self.prefix.complete
                and self.prefix.total_fully_testable
                >= self.cfg["k_g3b"])

    # -- stages ---------------------------------------------------

    def run_g3a(self) -> dict:
        """First, and before the generator exists. A G3a failure has
        spent nothing, and that is a property of this ordering rather
        than a promise about it."""

        assert self.rng is None, "the stream was opened before G3a"
        with meter.metered(self.budget):
            result = g3a.run_preflight(self.cfg["tol"], self.cfg["eta"])
        if not result["passed"]:
            raise stages.StageFailure(
                "g3a", "INVALID",
                f"wrapper contract failed "
                f"{result['failed_conditions']}",
                preserved={"fresh_seed_touched": False},
                detail={"preflight": result})
        return result

    def run_g3b(self) -> dict:
        """The fixed prefix, then the scan, on the G1 stream.

        Every drawn point goes into the G1 accumulator through
        `Prefix.observe`, whether or not it is a candidate and whether
        or not it is eligible."""

        self.rng = np.random.default_rng(self.seeds["o4b_g1_audit"])
        cap = min(g3.scan_cap(), self.cfg["n_g1"])
        with meter.metered(self.budget):
            while self.scanned < cap and not self._g3b_done():
                k = min(self.cfg["chunk"], cap - self.scanned)
                points = list(zip(*self._chunk(k), strict=True))
                for i, (r, th) in enumerate(points):
                    # STOP AT THE STOPPING TIME, not at the chunk end.
                    # The draw is batched; running on to the end of the
                    # buffer would contract-check points that lie after
                    # the sample the rule defined, and a mismatch there
                    # could turn an already complete G3b INVALID.
                    if self._g3b_done():
                        # the RNG is already past these; G1 takes them
                        self.pending = [(float(a), float(b))
                                        for a, b in points[i:]]
                        break
                    self._one(float(r), float(th))
        report = self.prefix.report()
        if self.mismatches:
            raise stages.StageFailure(
                "g3b", "INVALID",
                f"{len(self.mismatches)} wrapper contract mismatches; "
                f"the redesign requires zero",
                preserved=self.preserved(),
                detail={"mismatches": self.mismatches[:16],
                        "availability": report})
        if self.prefix.total_fully_testable < self.cfg["k_g3b"]:
            raise stages.StageFailure(
                "g3b", "INCONCLUSIVE",
                f"scan reached the end of the frozen G1 sample with "
                f"{self.prefix.total_fully_testable} fully-testable "
                f"clusters, short of {self.cfg['k_g3b']}",
                preserved=self.preserved(),
                detail={"availability": report})
        self.checkpoint("g3b", {"availability": report})
        return report

    def _consume(self, points, target: int) -> None:
        """G1 points, in the order they were drawn."""

        for r, th in points:
            if self.g1.n >= target:
                break
            # the estimand, at the drawn theta -- the predicate's
            # coordinate has no business here
            ell, _, _ = self.ell(float(r), float(th), self.cfg["tol"])
            self.g1.add(ell / sz.DT)
            self._advance()

    def _advance(self) -> None:
        """The cursor, moved in the same breath as the accumulation.

        Not after the judgement: a cap inside `judge` would then leave
        an incident whose G1 sample includes a point the cursor says
        was never reached (review R17)."""

        self.scanned += 1
        self.chunk_consumed += 1

    def _one(self, r: float, theta: float) -> None:
        """One drawn point, in BOTH coordinates and in that order.

        `z` comes from the inherited `_ell`, at the drawn `theta`:
        that is the frozen estimand, and it goes into G1 whatever
        happens next. Only a candidate gets a second, separate state
        computed in the predicate's recovered `dpsi`, which is what
        G3b places its probes from. Computing one state and using it
        for both would put the estimand in the predicate's coordinate
        or the probes in the estimand's; the design separates them
        deliberately."""

        ell, _, _ = self.ell(r, theta, self.cfg["tol"])
        z = ell / sz.DT

        def judge() -> dict:
            state = self.state_of(r, theta, self.cfg["tol"])
            verdict = judge_cluster(state, self.cfg["eta"])
            if verdict.get("fully_testable"):
                self._check_contract(r, theta, state, verdict)
            return verdict

        self.prefix.observe(z, judge, on_accumulated=self._advance)

    #: What each probe must produce, LEG BY LEG. The design fixes
    #: both legs of all three probes, and the cost model froze six
    #: predicate calls per fully-testable cluster on exactly that
    #: basis (review R10). Checking only the starved leg would let a
    #: wrapper defect on the other one -- `None`, or a wrong `False`
    #: -- pass silently and a verdict be published over it, and it
    #: would spend four calls where the budget charged six.
    CONTRACT = {
        "inside": {"p_to_x": True, "x_to_q": True},
        "outside_above": {"p_to_x": True, "x_to_q": False},
        "outside_below": {"p_to_x": False, "x_to_q": True},
    }

    def _check_contract(self, r, theta, state, verdict) -> None:
        """Six answers per cluster, and every one of them decided.

        `None` never satisfies a row: an undecided answer at a probe
        placed a realized margin outside the error band is the
        contract failing, which is what the abort was."""

        for name, wanted in self.CONTRACT.items():
            t_x = self._probe_time(name, verdict["probes"][name])
            if t_x is None:
                continue                  # not constructible; tallied
            for leg, expected in wanted.items():
                # ONE call, then judge it, then the next call (review
                # R19). Making both calls and comparing afterwards
                # loses the first leg's disagreement when the second
                # call trips the cap -- and then the promotion rule
                # from R18 sees an empty list and settles for
                # INCONCLUSIVE, which is the very masking that rule
                # was added to prevent, one level further in.
                got = self._leg(r, theta, t_x, leg)
                if got is not expected:
                    self.mismatches.append({
                        "r": r, "theta": theta, "dpsi": state["dpsi"],
                        "probe": name, "leg": leg,
                        "want": expected, "got": repr(got),
                        "t_x": t_x,
                        "t1": state["t1"], "err1": state["err1"],
                        "t2": state["t2"], "err2": state["err2"]})

    @staticmethod
    def _probe_time(name: str, probe: dict) -> float | None:
        """The intermediate event's TIME, from whichever placement the
        probe produced. The outside probes place a leg's `dt`; the
        time is where that puts `x`."""

        if not probe.get("reached"):
            return None
        if name == "inside":
            return probe["t_x"]
        if name == "outside_above":
            return sz.DT - probe["dt"]     # starves x -> q
        return probe["dt"]                 # starves p -> x

    def _leg(self, r, theta, t_x, leg: str):
        """ONE predicate call, at one intermediate event.

        One call per invocation on purpose: the caller commits each
        answer before asking for the next, so an interruption cannot
        take an already-observed disagreement with it. Six of these
        per fully-testable cluster, which is what the budget charged."""

        x = np.array([t_x, r, theta, 0.0])
        if leg == "p_to_x":
            p = np.array([0.0, sz.R_IN, 0.0, 0.0])
            return s1.causal_relation(p, x, s1.M, self.cfg["tol"])
        q = np.array([sz.DT, sz.R_OUT, 0.0, 0.0])
        return s1.causal_relation(x, q, s1.M, self.cfg["tol"])

    def run_g1(self) -> dict:
        """The rest of the frozen sample, continuing the SAME stream.

        The prefix already consumed part of it; resuming a fresh
        generator here would re-draw points the estimator has
        already taken."""

        target = self.cfg["n_g1"]
        with meter.metered(self.budget):
            # whatever G3b drew and did not use comes first, in order
            self._consume(self.pending, target)
            self.pending = []
            while self.g1.n < target:
                k = min(self.cfg["chunk"], target - self.g1.n)
                self._consume(list(zip(*self._chunk(k), strict=True)),
                              target)
                self.checkpoint("g1_chunk")
        interval = eb.interval(self.g1, self.cfg["delta_g1_per_side"])
        lo, hi = interval.rescaled(sz.SCALE)
        gap = (lo - sz.V_HI, hi - sz.V_LO)
        band = self.cfg["tau"] * sz.V_REF
        if gap[0] >= -band and gap[1] <= band:
            status = "concordant"
        elif gap[1] < -band or gap[0] > band:
            status = "discordant"
        else:
            status = "inconclusive"
        result = {"status": status, "n": self.g1.n,
                  "v_s1_lo": lo, "v_s1_hi": hi,
                  "identified_discrepancy": list(gap),
                  "band_abs": band}
        self.checkpoint("g1_complete", {"g1": result})
        return result

    def run_g2(self) -> dict:
        """Leakage outside the certified box, on its own stream."""

        rng = np.random.default_rng(self.seeds["o4b_g2_leakage"])
        acc, leaks, target = eb.Accumulator(), 0, self.cfg["n_g2"]
        with meter.metered(self.budget):
            while acc.n < target:
                k = min(self.cfg["chunk"], target - acc.n)
                rs, ths = self.draw(rng, 2 * k, sz.PATCH_R[0],
                                    sz.PATCH_R[1], sz.PATCH_CAP)
                for r, th in zip(rs, ths, strict=True):
                    if acc.n >= target:
                        break
                    if (sz.R_LO <= r <= sz.R_HI) and th <= sz.PSI_MAX:
                        continue
                    ell, _, _ = self.ell(float(r), float(th),
                                         self.cfg["tol"])
                    acc.add(ell / sz.DT)
                    if ell > 0.0:
                        leaks += 1
        alpha = self.cfg["alpha_g2_per_end"]
        upper = sz.DT * sz.B_OUT * xb.cp_upper(leaks, acc.n, alpha)
        lower = sz.B_OUT * sz.DT * max(0.0, eb.lower_bound(acc, alpha))
        allowed = self.cfg["leak_budget_frac"] * sz.V_REF
        if upper <= allowed:
            status = "concordant"
        elif lower > allowed:
            status = "discordant"
        else:
            status = "inconclusive"
        result = {"status": status, "n": acc.n, "leaking_points": leaks,
                  "leak_upper_abs": upper, "leak_lower_abs": lower,
                  "budget_abs": allowed}
        self.checkpoint("g2_complete", {"g2": result})
        return result

    def _staged(self, stage: str, work):
        """Run one stage, turning a fired cap into that stage's
        `INCONCLUSIVE` failure.

        Without this the cap escapes as `CapReached` and the entry
        point, which catches `StageFailure`, lets it out (review R11).
        The cap can fire DURING the G3b scan, before the first `g3b`
        checkpoint exists -- and by then the seeds are reserved and
        spent, and the reservation forbids re-running them. Everything
        accumulated would vanish with the process, which is precisely
        the loss O4 suffered and this design exists to prevent.

        A cap is not a contract failure: the instrument said nothing
        wrong, the run simply ran out. Hence `INCONCLUSIVE`, and the
        preserved sample goes with it."""

        try:
            return work()
        except stages.StageFailure:
            raise                       # already the right shape
        except o4b_budget.CapReached as cap:
            self.checkpoint_on_cap(stage)
            if self.mismatches:
                raise self._contract_failure(stage, cap) from None
            raise stages.StageFailure(
                stage, "INCONCLUSIVE",
                f"completion budget exhausted ({cap.reason}) during "
                f"{stage}: {cap.reserved:,} calls reserved, "
                f"{cap.wall_s:.1f}s elapsed",
                preserved=self.preserved(),
                detail={"cap": {"reason": cap.reason,
                                "reserved": cap.reserved,
                                "completed": cap.completed,
                                "wall_s": cap.wall_s,
                                "refused_call_count": cap.wanted}},
            ) from None
        except BaseException as exc:
            # ANY other exit through a stage (review R14). The
            # fail-closed solver paths raise `SystemExit`, which is a
            # BaseException and would slip past `except Exception`
            # entirely; a numerical failure in `_ell` or
            # `causal_relation` raises an ordinary one. Either way the
            # seeds are already reserved and spent, the reservation
            # forbids re-running them, and without this the process
            # dies with no incident at all -- which is the observability
            # defect this freeze exists to close.
            # KeyboardInterrupt included (review R21). It is a
            # deliberate stop, but the seeds are spent all the same
            # once the reservation is claimed, and "the operator meant
            # it" is not a reason to leave no record of a run that
            # cannot be repeated. It is still fatal; it just leaves an
            # incident on the way out.
            self.checkpoint_on_cap(stage)
            if self.mismatches:
                raise self._contract_failure(stage, exc) from exc
            raise stages.StageFailure(
                stage, "ABORT",
                f"unhandled {type(exc).__name__} during {stage}: {exc}",
                preserved=self.preserved(),
                detail={"exception": {"type": type(exc).__name__,
                                      "message": str(exc)}},
            ) from exc

    def _contract_failure(self, stage: str, cause) -> stages.StageFailure:
        """An observed mismatch outranks whatever stopped the run
        afterwards (review R18).

        The frozen contract is that ONE probe disagreement is a
        contract failure. Mismatches are collected rather than raised
        on the spot so the record can carry several, but that meant a
        cap firing on the very next solver call inside the same
        `_check_contract` ended the run `INCONCLUSIVE` -- an
        instrumentation failure already observed, reported as "we ran
        out of budget", with the mismatches nowhere in the incident.

        `INCONCLUSIVE` indicts nothing; `INVALID` indicts the
        instrument. Having seen the instrument disagree, the second is
        the true sentence and the interruption is a detail of it."""

        return stages.StageFailure(
            stage, "INVALID",
            f"{len(self.mismatches)} wrapper contract mismatch(es) "
            f"observed before the run stopped; the redesign requires "
            f"zero",
            preserved=self.preserved(),
            detail={"mismatches": self.mismatches[:16],
                    "availability": self.prefix.report(),
                    "stopped_by": {"type": type(cause).__name__,
                                   "message": str(cause)},
                    "why_not_inconclusive": (
                        "the contract failure was observed first and "
                        "is what the run is evidence of; the "
                        "interruption only decided when it stopped")})

    def checkpoint_on_cap(self, stage: str) -> None:
        """Record progress on the way out, WITHOUT destroying a record
        of something already finished (review R15).

        Only G1 has a truthful point to write here: `g1_chunk`. A stop
        inside G3a or G3b has reached no checkpoint stage at all, and
        writing `g3b` would claim G3b finished. A stop inside G2 is
        the dangerous one: the last checkpoint is `g1_complete`, which
        carries G1's status and interval, and this exit publishes no
        results artifact -- overwriting it with a `g1_chunk` record
        would delete a G1 result the run had already established.

        In every other case the incident is the record: it carries
        `preserved()`, which holds the G1 partial, the availability
        report and the budget."""

        if stage == "g1":
            self.checkpoint("g1_chunk",
                            {"availability": self.prefix.report()})

    def run(self, on_g3a_passed=None, on_before_publish=None) -> dict:
        """The frozen order. Each stage's failure is raised, never
        swallowed, and the caller turns it into an incident.

        `on_g3a_passed` runs BETWEEN G3a and G3b, and in the campaign
        it is the reservation claim (review R20). The frozen rule is
        that a G3a failure has spent no fresh seed, and claiming the
        ref before G3a broke it in the only way that matters: the ref
        is the authority on whether a stream has been opened, so a
        claim made and then abandoned retires both seeds by policy
        even though no generator was ever constructed. Not
        constructing the generator is necessary and was never
        sufficient.

        The claim sits inside this method rather than in `main` so the
        order it belongs to is the order that states it, and so a test
        can watch the sequence."""

        g3a_result = self._staged("g3a", self.run_g3a)
        if on_g3a_passed is not None:
            # INSIDE the incident boundary (review R23). The claim is
            # the point of no return, so a failure there is the one
            # failure that most needs a record -- and the ref can be
            # created by a push whose reply is then lost, spending the
            # seeds while the call never returns. Attributed to `g3b`,
            # the stage it was entering.
            self.reservation_object = self._staged("g3b",
                                                   on_g3a_passed)
        availability = self._staged("g3b", self.run_g3b)
        g1 = self.g1_result = self._staged("g1", self.run_g1)
        g2 = self.g2_result = self._staged("g2", self.run_g2)
        # THE CLAIM IS RE-READ BEFORE ANYTHING IS PUBLISHED (review
        # R25). Eleven hours separate the claim from this point, and
        # nothing so far would notice the ref being deleted or
        # replaced in between -- so the unique attempt object the
        # nonce exists to create was never actually tied to the
        # result. A result published over someone else's claim is a
        # result whose provenance says nothing.
        if on_before_publish is not None:
            self._staged("g2", lambda: on_before_publish(
                self.reservation_object))
        return {
            "kind": "results", "run_kind": "campaign",
            "freeze_sha": self.freeze_sha,
            "manifest_digest": self.digest,
            "seeds": self.seeds,
            "reservation": {
                "ref": reservation.REF,
                "authority": reservation.CANONICAL_AUTHORITY,
                "object": self.reservation_object,
                "seeds_spent": True,
                "verified_at_exit": on_before_publish is not None,
                "why": ("the streams were opened by this attempt and "
                        "are retired; the object is the one the ref "
                        "held both when the claim was made and when "
                        "this result was published"),
            },
            "order": list(stages.STAGES),
            "g3a": {"passed": g3a_result["passed"],
                    "conditions": g3a_result["conditions"]},
            "availability": availability,
            "g1": g1, "g2": g2,
            "budget": self.budget.state(),
            "environment": environment(),
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--write-manifest", action="store_true",
                        help="regenerate the freeze manifest; a freeze "
                             "action, never part of a campaign run")
    parser.add_argument("--freeze-rev", default="")
    args = parser.parse_args()

    if args.write_manifest:
        print(f"wrote {write_manifest()}")
        return

    checks = preflight(args.freeze_rev)
    if args.preflight:
        print(f"preflight ok at {checks['git']['rev']}")
        print(f"  seeds {checks['seeds']}")
        print(f"  {reservation.RETAINED_REF} retained at "
              f"{reservation.RETAINED_OBJECT[:7]}")
        return

    digest = _sha256(_MANIFEST)
    campaign = Campaign(FROZEN, checks["seeds"], args.freeze_rev,
                        digest)
    claimed: dict = {"object": None}
    receipt: dict = {}

    def claim() -> None:
        payload = {
            "campaign": "o4b", "freeze_rev": args.freeze_rev,
            "manifest_sha256": digest, "seeds": checks["seeds"],
        }
        try:
            claimed["object"] = reservation.claim(payload)
        except reservation.ClaimUncertain as uncertain:
            # record before re-raising: `_staged` turns this into a
            # StageFailure and the incident must be able to say the
            # streams may already be open
            claimed["uncertain"] = uncertain.as_record()
            raise

    def verify_before_publish(obj: str | None) -> None:
        if obj is None:                          # pragma: no cover
            raise SystemExit("internal: no claim object to verify")
        reservation.verify_still_held(obj)

    def file_incident(failure: stages.StageFailure) -> None:
        publish_write_once(_INCIDENT, stages.incident(failure, {
            "freeze_sha": args.freeze_rev,
            "manifest_digest": digest,
            "seeds": checks["seeds"],
            # what the ledger has to be able to read: the seeds are
            # spent if and only if the ref was claimed
            "reservation_claimed": (
                "uncertain" if claimed.get("uncertain")
                else claimed["object"] is not None),
            "reservation_object": claimed["object"],
            "reservation_uncertainty": claimed.get("uncertain"),
            # fail-closed toward SPENT: a push whose reply was lost
            # still created the ref, and under-recording a claim is
            # the unsafe direction
            "seeds_spent": (claimed["object"] is not None
                            or bool(claimed.get("uncertain"))),
            "environment": environment(),
        }))

    # THE PUBLICATION IS INSIDE THE BOUNDARY TOO (review R26). By this
    # point the seeds are permanently spent and G1 and G2 have run to
    # their frozen sample sizes, so a failure at the very last step --
    # the artifact already existing, a full disk, an interrupt -- would
    # end the process with no record of a campaign that cannot be
    # repeated. The `g2_complete` checkpoint does not stand in for it:
    # it carries running statistics, not G1's final status and
    # interval.
    try:
        results = campaign.run(on_g3a_passed=claim,
                               on_before_publish=verify_before_publish)
        publish_write_once(_ARTIFACT, results, receipt=receipt)
    except stages.StageFailure as failure:
        file_incident(failure)
        raise SystemExit(
            f"{failure.stage}: {failure.outcome} -- "
            f"{failure.reason}; incident written to "
            f"{_INCIDENT.name}, no verdict published") from None
    except BaseException as exc:
        # THIS CALL's commit, not the file's existence (review R28).
        # `os.link` is the commit and the receipt is stamped the
        # instant it returns, so an interrupt after it still reports a
        # published result -- which is R27. But a file that appeared
        # from somewhere else during the run is NOT this run's result:
        # `os.link` would have failed, and reading `exists()` there
        # would credit another file to this campaign and file no
        # incident for seeds that are already spent.
        if receipt.get("committed"):
            raise SystemExit(
                f"{_ARTIFACT.name} was published; the failure came "
                f"after the commit ({type(exc).__name__}: {exc}) and "
                f"no incident is filed over a published result"
            ) from exc
        file_incident(stages.StageFailure(
            "g2", "ABORT",
            f"the campaign completed but the result could not be "
            f"published: {type(exc).__name__}: {exc}",
            preserved=campaign.preserved(),
            detail={"exception": {"type": type(exc).__name__,
                                  "message": str(exc)},
                    "at": "artifact publication",
                    "note": ("the gates ran to their frozen sample "
                             "sizes and are preserved here; the "
                             "result artifact was never written, so "
                             "no verdict is published")}))
        raise SystemExit(
            f"publication failed ({type(exc).__name__}); incident "
            f"written to {_INCIDENT.name}, no verdict published"
        ) from exc
    print(f"published {_ARTIFACT.name}: G1 {results['g1']['status']}, "
          f"G2 {results['g2']['status']}")


if __name__ == "__main__":
    main()
