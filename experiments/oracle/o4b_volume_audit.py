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

    G3a  ->  G3b on the G1 stream's fixed prefix  ->  the rest of G1
         ->  G2

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


def publish_write_once(path: Path, payload: dict) -> None:
    """Atomic no-clobber publication: fsynced temp in the same
    directory, then `os.link`, which fails atomically if the
    destination exists."""

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
    finally:
        tmp.unlink(missing_ok=True)


# ------------------------------------------------------- the G3 state

def solver_state(r: float, theta: float, tol: float) -> dict:
    """The predicate's own view of one cluster.

    `err` is kept, unlike in `_ell`, because the probes are PLACED
    with it. It still never enters the inference: it is a
    Gauss-Legendre stopping heuristic, not a certified enclosure, and
    an interval built from it would not be conservative."""

    d1: dict = {}
    d2: dict = {}
    t1, err1 = s1.flight_time(sz.R_IN, r, theta, s1.M, tol, d1)
    t2, err2 = s1.flight_time(r, sz.R_OUT, theta, s1.M, tol, d2)
    return {"t1": t1, "err1": err1, "t2": t2, "err2": err2,
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
    return {"n": acc.n, "mean_z": acc.mean, "var_z": acc.var}


class Campaign:
    """The four stages in their frozen order, with the checkpoint and
    the preservation rules attached to the transitions."""

    def __init__(self, cfg: dict, seeds: dict, freeze_sha: str,
                 digest: str, checkpoint_path: Path = _CHECKPOINT,
                 draw=None, state_of=None) -> None:
        self.cfg = cfg
        self.seeds = seeds
        self.freeze_sha = freeze_sha
        self.digest = digest
        self.checkpoint_path = Path(checkpoint_path)
        self.budget = o4b_budget.Budget(cfg["max_calls"],
                                        cfg["max_wall_s"])
        self.draw = draw if draw is not None else o4._draw
        self.state_of = (state_of if state_of is not None
                         else solver_state)
        self.g1 = eb.Accumulator()
        self.prefix = stages.Prefix(self.g1, cfg["n_avail"])
        self.rng = None
        self.mismatches: list[dict] = []
        self.scanned = 0

    # -- checkpoints ---------------------------------------------

    def checkpoint(self, point: str, extra: dict | None = None) -> None:
        stages.write_checkpoint(
            self.checkpoint_path, point,
            freeze_sha=self.freeze_sha, digest=self.digest,
            seed=self.seeds["o4b_g1_audit"], rng=self.rng,
            samples=self.g1.n,
            statistics={**_statistics(self.g1), **(extra or {})},
            budget=self.budget)

    def preserved(self) -> dict:
        """What a failure hands to the incident. Always the G1 prefix:
        it is paid for and is valid independently of anything G3b
        found."""

        return {
            "g1_partial": _statistics(self.g1),
            "availability": self.prefix.report(),
            "scanned_points": self.scanned,
            "budget": self.budget.state(),
            "is_not_a_verdict": (
                "a preserved partial sample: the stopping rule did "
                "not fire, so no gate has a status"),
        }

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
            while (self.scanned < cap
                   and (not self.prefix.complete
                        or self.prefix.total_fully_testable
                        < self.cfg["k_g3b"])):
                k = min(self.cfg["chunk"], cap - self.scanned)
                rs, ths = self.draw(self.rng, k, sz.R_LO, sz.R_HI,
                                    sz.PSI_MAX)
                for r, th in zip(rs, ths, strict=True):
                    self._one(float(r), float(th))
                    self.scanned += 1
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

    def _one(self, r: float, theta: float) -> None:
        state = self.state_of(r, theta, self.cfg["tol"])
        ell = sz.DT - state["t1"] - state["t2"]
        z = (ell if ell > 0.0 else 0.0) / sz.DT

        def judge() -> dict:
            verdict = judge_cluster(state, self.cfg["eta"])
            if verdict.get("fully_testable"):
                self._check_contract(r, theta, state, verdict)
            return verdict

        self.prefix.observe(z, judge)

    def _check_contract(self, r, theta, state, verdict) -> None:
        """The three probes must answer as the window says. One
        disagreement anywhere is a contract failure."""

        want = {"inside": True, "outside_above": False,
                "outside_below": False}
        for name, expected in want.items():
            probe = verdict["probes"][name]
            got = self._ask(r, theta, state, name, probe)
            if got is not expected:
                self.mismatches.append({
                    "r": r, "theta": theta, "probe": name,
                    "want": expected, "got": repr(got),
                    "t1": state["t1"], "err1": state["err1"],
                    "t2": state["t2"], "err2": state["err2"]})

    def _ask(self, r, theta, state, name, probe):
        """One predicate answer, at the time the probe placed.

        The inside probe asks about BOTH legs, because "the window is
        open at this time" is the conjunction; the outside probes ask
        about the single leg they starve."""

        if name == "inside":
            return self._pair(r, theta, probe["t_x"])
        return self._leg(r, theta, name, probe["dt"])

    def _pair(self, r, theta, t_x):
        """Both legs of the inside probe: connected means BOTH."""

        p = np.array([0.0, sz.R_IN, 0.0, 0.0])
        x = np.array([t_x, r, theta, 0.0])
        q = np.array([sz.DT, sz.R_OUT, 0.0, 0.0])
        first = s1.causal_relation(p, x, s1.M, self.cfg["tol"])
        second = s1.causal_relation(x, q, s1.M, self.cfg["tol"])
        if first is None or second is None:
            return None
        return bool(first and second)

    def _leg(self, r, theta, name, dt):
        if name == "outside_above":
            x = np.array([sz.DT - dt, r, theta, 0.0])
            q = np.array([sz.DT, sz.R_OUT, 0.0, 0.0])
            return s1.causal_relation(x, q, s1.M, self.cfg["tol"])
        p = np.array([0.0, sz.R_IN, 0.0, 0.0])
        x = np.array([dt, r, theta, 0.0])
        return s1.causal_relation(p, x, s1.M, self.cfg["tol"])

    def run_g1(self) -> dict:
        """The rest of the frozen sample, continuing the SAME stream.

        The prefix already consumed part of it; resuming a fresh
        generator here would re-draw points the estimator has
        already taken."""

        target = self.cfg["n_g1"]
        with meter.metered(self.budget):
            while self.g1.n < target:
                k = min(self.cfg["chunk"], target - self.g1.n)
                rs, ths = self.draw(self.rng, k, sz.R_LO, sz.R_HI,
                                    sz.PSI_MAX)
                for r, th in zip(rs, ths, strict=True):
                    state = self.state_of(r, th, self.cfg["tol"])
                    ell = sz.DT - state["t1"] - state["t2"]
                    self.g1.add((ell if ell > 0.0 else 0.0) / sz.DT)
                self.scanned += k
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
                    state = self.state_of(float(r), float(th),
                                          self.cfg["tol"])
                    ell = sz.DT - state["t1"] - state["t2"]
                    acc.add((ell if ell > 0.0 else 0.0) / sz.DT)
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

    def run(self) -> dict:
        """The frozen order. Each stage's failure is raised, never
        swallowed, and the caller turns it into an incident."""

        g3a_result = self.run_g3a()
        availability = self.run_g3b()
        g1 = self.run_g1()
        g2 = self.run_g2()
        return {
            "kind": "results", "run_kind": "campaign",
            "freeze_sha": self.freeze_sha,
            "manifest_digest": self.digest,
            "seeds": self.seeds,
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
    reservation.claim({
        "campaign": "o4b", "freeze_rev": args.freeze_rev,
        "manifest_sha256": digest, "seeds": checks["seeds"],
    })
    campaign = Campaign(FROZEN, checks["seeds"], args.freeze_rev,
                        digest)
    try:
        results = campaign.run()
    except stages.StageFailure as failure:
        publish_write_once(_INCIDENT, stages.incident(failure, {
            "freeze_sha": args.freeze_rev,
            "manifest_digest": digest,
            "seeds": checks["seeds"],
            "environment": environment(),
        }))
        raise SystemExit(
            f"{failure.stage}: {failure.outcome} -- "
            f"{failure.reason}; incident written to "
            f"{_INCIDENT.name}, no verdict published") from None
    publish_write_once(_ARTIFACT, results)
    print(f"published {_ARTIFACT.name}: G1 {results['g1']['status']}, "
          f"G2 {results['g2']['status']}")


if __name__ == "__main__":
    main()
