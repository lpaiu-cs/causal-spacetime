"""The order the O4b stages run in, and what survives each failure.

O4's order was `G3a -> G1/G2 -> G3b`. That is the order in which you
spend twelve hours and then find out the instrumentation was broken,
and it is what happened. G3b is a precondition for reading G1 at all:
until `causal_relation` is known to accept the window `L_S1` implies,
G1's numbers do not mean what they are taken to mean. So it runs on
the front of the G1 stream, and the rest of G1 follows.

    G3a  ->  G3b on the fixed prefix  ->  the rest of G1  ->  G2

THE PREFIX IS SHARED, WHICH IS WHERE THE BIAS WOULD COME FROM. G3b
does not draw its own sample; it consumes the beginning of G1's. Two
different filters could leak into the G1 accumulator through that, and
neither may:

  * G3b ELIGIBILITY. Accumulating only points with `W_robust > 0`
    would make `V_S1` an estimator of a conditional distribution.

  * `L_S1 > 0`. Dropping the points where S1 opened no window at all
    -- the `Z = 0` points -- biases `V_S1` upward. The fixed prefix's
    denominator is defined over `L_S1 > 0` candidates because it is
    the denominator of an AVAILABILITY report, and that is a different
    tally over the same stream, not the estimator's sample definition.

The rule that keeps both out is one line of control flow, and it is
worth stating as a rule rather than trusting to the reading: EVERY
DRAWN POINT ENTERS THE G1 ACCUMULATOR BEFORE ANY PREDICATE IS ASKED
ABOUT IT. Judging first and accumulating second is the same code with
the bias in it.

WHAT SURVIVES EACH FAILURE.

  G3a fails   incident; NO FRESH SEED HAS BEEN TOUCHED. The generator
              is not constructed until G3a passes, so this is a
              property of the code, not a promise about it.

  G3b fails   incident, and THE ACCUMULATED G1 PREFIX IS KEPT. That
              statistic is already paid for and is valid independently
              of the instrumentation check -- discarding it is what
              O4 did, and it cost twelve hours. Kept as `partial`, so
              it is a preserved sample and not a verdict.

  cap fires   incident, plus the last checkpoint. `INCONCLUSIVE`.

No path returns a verdict and an incident together. Fail-closed here
means no scientific sentence is published, NOT that nothing is
recorded -- the O4 abort's observability defect was the second
reading, and every stage below writes before it stops.
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

import o4b_checkpoint as ck  # noqa: E402

#: Ordered, and the order is part of the freeze.
STAGES = ("g3a", "g3b", "g1", "g2")

#: Terminal outcomes. `INVALID` means the instrument contract failed;
#: `INCONCLUSIVE` means the run ran out of budget or sample before the
#: contract sample was complete. Neither is a scientific verdict, and
#: they are not interchangeable: the first indicts the instrument, the
#: second indicts nothing.
OUTCOMES = ("VERDICT", "INVALID", "INCONCLUSIVE", "ABORT")


class StageFailure(Exception):
    """A stage that cannot continue. Carries what has to be preserved,
    because the alternative -- raising bare and reconstructing the
    context at the catch site -- is how O4 lost its statistics."""

    def __init__(self, stage: str, outcome: str, reason: str,
                 preserved: dict | None = None,
                 detail: dict | None = None,
                 failure_point: str | None = None) -> None:
        super().__init__(f"{stage}: {reason}")
        if stage not in STAGES:
            raise ValueError(f"{stage!r} is not a frozen stage")
        if outcome not in OUTCOMES or outcome == "VERDICT":
            raise ValueError(
                f"{outcome!r} is not a failure outcome; a stage that "
                f"failed does not produce a verdict")
        self.stage = stage
        self.outcome = outcome
        self.reason = reason
        self.preserved = preserved or {}
        self.detail = detail or {}
        # The STEP that failed, when it is not the stage itself. The
        # pre-publish re-verify and the reservation claim both run under
        # a stage label (`g2`, `g3b`) that would otherwise read as a
        # claim about that gate; `failure_point` records the actual step
        # so a reader does not mistake a post-gate stop for an
        # unfinished gate (incident defect 5b).
        self.failure_point = failure_point or stage


#: What an incident owns. A caller's run context may not supply these
#: (review R5): the record is write-once and is the provenance an audit
#: or a recovery reads, so a reused context carrying `stage: "g1"` or
#: `verdict: "concordant"` would file a G3b failure as a G1 result --
#: and break the invariant one line above it, that no path returns a
#: verdict and an incident together. Same boundary as the checkpoint's.
INCIDENT_KEYS = (
    "kind", "run_kind", "stage", "failure_point", "termination_reason",
    "outcome", "verdict", "why_no_verdict", "preserved", "unavailable",
    "detail",
)


def incident(failure: StageFailure, context: dict) -> dict:
    """The record every abort path writes. Write-once at the caller;
    this builds it so no path can compose a different shape."""

    clashing = [k for k in INCIDENT_KEYS if k in context]
    if clashing:
        raise ValueError(
            f"incident context supplies {clashing}, which the incident "
            f"itself owns -- the failure decides its stage, outcome "
            f"and that there is no verdict, and a run context cannot "
            f"overrule the failure it is recording")
    return {
        "kind": "incident",
        "run_kind": "campaign",
        "stage": failure.stage,
        "failure_point": failure.failure_point,
        "termination_reason": failure.reason,
        "outcome": failure.outcome,
        "verdict": None,
        "why_no_verdict": (
            "fail-closed: the stage that stopped is a precondition "
            "for the sentence, so no gate has a status"),
        "preserved": failure.preserved,
        "unavailable": failure.detail.get("unavailable", []),
        "detail": failure.detail,
        **context,
    }


def fresh_seed_touched(stage: str) -> bool:
    """Whether reaching this stage means a campaign stream was drawn
    from. G3a is entirely deterministic -- a frozen case table and a
    solver -- so a G3a failure spends nothing."""

    if stage not in STAGES:
        raise ValueError(f"{stage!r} is not a frozen stage")
    return STAGES.index(stage) >= STAGES.index("g3b")


class Prefix:
    """The shared front of the G1 stream.

    Holds THREE tallies that must not be confused.

      * The G1 estimator accumulates EVERY drawn point, for the whole
        run. It never stops and it never filters.

      * The availability report counts `L_S1 > 0` candidates and what
        happened to them -- and FREEZES at the `N_avail`-th candidate.

      * The scan carries on past that point, because G3b needs
        `K_G3B` fully-testable clusters and the prefix will not have
        supplied them all. Its counts are recorded separately.

    THE FREEZE IS THE WHOLE POINT OF THE FIXED PREFIX (review R6). The
    scan's stopping time depends on what it finds, so its denominator
    is chosen by the results; the fixed prefix exists precisely so one
    denominator is not. Letting the counters run on past `N_avail`
    would put the scan's points into the fixed report's numerator
    while the denominator stayed `N_avail` -- an `eligible_rate` above
    1 is the visible form of that, and a slightly inflated rate is the
    invisible one. Any prefix with even a single unavailable candidate
    reaches this, because the scan then has to keep going.
    """

    def __init__(self, accumulator, n_avail: int) -> None:
        self.acc = accumulator          # G1's estimator, unfiltered
        self.n_avail = n_avail
        # the fixed report, frozen at the N_avail-th candidate
        self.candidates = 0             # L_S1 > 0, the report's base
        self.eligible = 0
        self.fully_testable = 0
        self.zero_window = 0            # Z = 0: accumulated, not a
        self.reasons: dict[str, int] = {}   # candidate
        # the scan beyond it, counted apart
        self.scan_candidates = 0
        self.scan_fully_testable = 0
        self.scan_zero_window = 0

    @property
    def complete(self) -> bool:
        return self.candidates >= self.n_avail

    @property
    def total_fully_testable(self) -> int:
        """What the `K_G3B` completion rule counts: the contract sample
        is every fully-testable cluster, prefix and scan alike. Only
        the AVAILABILITY report is confined to the prefix."""

        return self.fully_testable + self.scan_fully_testable

    def observe(self, z: float, judge=None,
                on_accumulated=None) -> dict | None:
        """One drawn point.

        `z` enters the estimator FIRST and unconditionally -- before
        the prefix is complete, after it is complete, candidate or
        not. Only then is the point examined, and only if it is a
        candidate. `judge` returns a dict with `eligible` and
        `fully_testable`; it is not called for non-candidates, because
        they are not part of the availability question.

        TWO COMMIT BOUNDARIES, AND A FAILURE FALLS BETWEEN THEM
        CLEANLY (review R17). `judge` runs the solver, so it can raise
        -- a fired cap, a fail-closed path. Everything before that
        point commits together: the estimator takes `z` and
        `on_accumulated` advances the caller's cursor in the same
        breath, so an incident can never report a G1 sample that
        includes a point the cursor says was never reached. Everything
        the judgement decides commits only if the judgement RETURNS.
        An earlier version counted the candidate first, so a cap at
        the `N_avail`-th one could report `complete = True` and
        publish availability rates over a candidate that was never
        judged."""

        self.acc.add(z)                 # <- before any predicate
        if on_accumulated is not None:
            on_accumulated()            # <- the cursor, same boundary
        frozen = self.complete          # BEFORE this point is counted
        if z <= 0.0:
            if frozen:
                self.scan_zero_window += 1
            else:
                self.zero_window += 1
            return None
        if judge is None:
            if frozen:
                self.scan_candidates += 1
            else:
                self.candidates += 1
            return None
        verdict = judge()               # may raise: nothing counted
        if frozen:
            self.scan_candidates += 1
            # the cluster still counts toward K_G3B; it just does not
            # enter a report whose denominator was fixed before it
            if verdict.get("fully_testable"):
                self.scan_fully_testable += 1
            return verdict
        self.candidates += 1
        if verdict.get("eligible"):
            self.eligible += 1
        if verdict.get("fully_testable"):
            self.fully_testable += 1
        else:
            reason = verdict.get("reason", "unspecified")
            self.reasons[reason] = self.reasons.get(reason, 0) + 1
        return verdict

    def report(self) -> dict:
        """Availability, on the fixed denominator. Reported as rates
        only once the prefix is complete: a rate taken from a prefix
        that stopped early has a denominator chosen by the stopping,
        which is the very thing the fixed prefix exists to avoid."""

        base = {
            "n_avail_target": self.n_avail,
            "candidates": self.candidates,
            "complete": self.complete,
            "eligible": self.eligible,
            "fully_testable": self.fully_testable,
            "zero_window_points": self.zero_window,
            "accumulated_points": self.acc.n,
            "accumulates_every_drawn_point": True,
            "not_fully_testable_by_reason": dict(self.reasons),
            "beyond_the_prefix": {
                "candidates": self.scan_candidates,
                "fully_testable": self.scan_fully_testable,
                "zero_window_points": self.scan_zero_window,
                "excluded_from_the_rates_because": (
                    "the scan's stopping time depends on what it "
                    "finds, so these points would put a results-chosen "
                    "numerator over a fixed denominator"),
            },
            "total_fully_testable_for_k_g3b": self.total_fully_testable,
        }
        if not self.complete:
            base["rates_withheld"] = (
                "the fixed prefix did not complete, so its denominator "
                "was chosen by the stopping and no rate is reported")
            return base
        base["eligible_rate"] = self.eligible / self.n_avail
        base["fully_testable_rate"] = self.fully_testable / self.n_avail
        base["rate_denominator"] = "the fixed prefix, not the scan"
        return base


def checkpoint_payload(freeze_sha: str, digest: str,
                       seed: int, rng, rng_stream: str, samples: int,
                       statistics: dict, budget) -> dict:
    """Everything a checkpoint has to carry, assembled once.

    `rng_position` is the bit generator's own state, not a count of
    draws: resuming has to continue the same stream, and a count only
    reproduces it if every consumer drew exactly as before.

    `rng_stream` NAMES which stream `seed`/`rng_position` belong to. A
    checkpoint stamped `g2_complete` still serialises the G1 stream --
    G2 runs on its own generator, never stored here -- so without the
    name a recovery could read the G1 seed as G2's and replay the wrong
    sample (incident defect 5c)."""

    return {
        "freeze_sha": freeze_sha,
        "manifest_digest": digest,
        "seed": seed,
        "rng_stream": rng_stream,
        "rng_position": (rng.bit_generator.state
                         if rng is not None else None),
        "samples": samples,
        "statistics": statistics,
        "budget": budget.state() if hasattr(budget, "state") else budget,
    }


def write_checkpoint(path: Path, point: str, **payload) -> Path:
    """Thin, so the stages do not each learn the checkpoint format.

    `point` is one of the four checkpoint STAGES in `o4b_checkpoint`,
    which are finer than this module's four run stages -- G1 alone
    writes at every chunk and again at completion."""

    return ck.write(path, point, checkpoint_payload(**payload))
