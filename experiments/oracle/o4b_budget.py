"""The O4b completion budget: charged BEFORE spending, and honest
about what it can and cannot bound.

The O4 runner counted calls after making them (`spend(k)` incremented,
then compared), so a batch could carry the total past the cap and the
cap would only notice afterwards. A cap that is checked after the fact
is not a cap. Here every call, or every bounded batch, must be
RESERVED first: `reserved + k <= max_calls` is decided before the
solver is touched, and the reservation only stands if that succeeds.

RESERVED IS NOT COMPLETED. If a solver call inside a reserved batch
raises, the batch never finishes, and a counter that only tracked
reservations would tell the incident record a number the run did not
actually spend. Both are kept: the cap is charged against `reserved`,
which never under-charges, and the incident reports `completed`
alongside it. `completed` may never exceed `reserved` -- a completion
without a reservation is a call that was never charged, which is the
same defect the pre-charge rule exists to prevent.

WHAT THE TWO CAPS ACTUALLY BOUND. They are not symmetric, and saying
so precisely matters more than making them look alike:

  * CALLS. Bounded rigorously. No batch begins unless it fits, and a
    batch is at most `max_batch_calls`, so the overshoot past
    `max_calls` is zero -- the cap is never crossed at all.
  * WALL CLOCK. A START RULE only: no batch BEGINS at or after
    `max_wall_s`. A batch already in flight runs to its end, so the
    finish time can exceed the cap by one batch's duration -- and this
    module has no rigorous bound on that duration. The 768 us/pair in
    `p14_s1_cost.md` is a measured average ON ONE HOST, not a
    certified per-call maximum, so multiplying it by the batch size
    gives a NOMINAL PROJECTION and nothing stronger. A rigorous time
    bound would need a per-call timeout or a certified call-time
    ceiling; neither exists today.

So the run records what actually happened: the real finish time and
the real overrun go into the result and the incident, rather than a
bound being asserted that the evidence does not support.

The clock and the counters start before G3a, because G3a's calls are
inside the budget. A budget that began at G1 would be spending
unmetered.

Reaching a cap raises `CapReached`. It is not an error in the
instrument: the completion budget ran out, which the freeze records as
an incident and `INCONCLUSIVE`. Raising a cap during or after a run is
forbidden -- there is no method for it, and every limit is read-only.
"""

from __future__ import annotations

import time

#: The measured per-call cost from `p14_s1_cost.md` at tol 1e-8: 768
#: us/pair, ON THAT HOST. An average, not a maximum. Used ONLY to
#: project a nominal overrun, never to bound one.
NOMINAL_SECONDS_PER_CALL = 0.768e-3

#: The largest batch a single reservation may take. It bounds the
#: call overshoot exactly (to zero, since a batch that does not fit is
#: refused) and bounds the wall overrun only in units of CALLS -- one
#: batch -- not in seconds.
MAX_BATCH_CALLS = 8_192


class CapReached(Exception):
    """A completion budget ran out. Carries the counters so the
    incident record and the console cannot disagree."""

    def __init__(self, reason: str, reserved: int, completed: int,
                 wall_s: float, wanted: int) -> None:
        super().__init__(
            f"{reason}: refused a reservation of {wanted} call(s) at "
            f"{reserved:,} reserved / {completed:,} completed / "
            f"{wall_s:.1f}s")
        self.reason = reason
        self.reserved = reserved
        self.completed = completed
        self.wall_s = wall_s
        self.wanted = wanted


class Budget:
    """A pre-charged call and wall-clock budget.

    `stage` is advisory bookkeeping: it records where the budget was
    when it ran out, which the incident needs and the caps do not."""

    def __init__(self, max_calls: int, max_wall_s: float,
                 max_batch_calls: int = MAX_BATCH_CALLS,
                 clock=time.perf_counter) -> None:
        self._max_calls = int(max_calls)
        self._max_wall_s = float(max_wall_s)
        self._max_batch_calls = int(max_batch_calls)
        self._clock = clock
        self._t0 = clock()
        self._reserved = 0
        self._completed = 0
        self.stage = "g3a"
        self.per_stage: dict[str, int] = {}

    # every limit and counter is read-only; they move only through
    # `reserve` and `complete`, which enforce the invariants
    @property
    def max_calls(self) -> int:
        return self._max_calls

    @property
    def max_wall_s(self) -> float:
        return self._max_wall_s

    @property
    def max_batch_calls(self) -> int:
        return self._max_batch_calls

    @property
    def reserved(self) -> int:
        return self._reserved

    @property
    def completed(self) -> int:
        return self._completed

    @property
    def unfinished(self) -> int:
        return self._reserved - self._completed

    @property
    def nominal_overrun_s(self) -> float:
        """A PROJECTION of how far past `max_wall_s` a run may finish:
        one batch at the measured average call cost. Not a bound --
        the average is host-dependent and is not a per-call maximum."""

        return self._max_batch_calls * NOMINAL_SECONDS_PER_CALL

    @property
    def wall_s(self) -> float:
        return self._clock() - self._t0

    @property
    def wall_overrun_s(self) -> float:
        """How far past the wall cap the run actually is, right now."""

        return max(0.0, self.wall_s - self._max_wall_s)

    def enter(self, stage: str) -> None:
        self.stage = stage

    def reserve(self, k: int = 1) -> None:
        """Charge `k` calls, or refuse. Decided BEFORE the calls are
        made, so the call cap is never crossed rather than merely
        detected afterwards."""

        if k < 0:
            raise ValueError(f"cannot reserve {k} calls")
        if k > self._max_batch_calls:
            raise ValueError(
                f"batch of {k} exceeds max_batch_calls "
                f"{self._max_batch_calls}: the overrun statement "
                f"depends on batches being bounded")
        if self._reserved + k > self._max_calls:
            raise CapReached("max-calls", self._reserved,
                             self._completed, self.wall_s, k)
        wall = self.wall_s
        if wall >= self._max_wall_s:
            raise CapReached("max-wall", self._reserved,
                             self._completed, wall, k)
        self._reserved += k
        self.per_stage[self.stage] = (
            self.per_stage.get(self.stage, 0) + k)

    def complete(self, k: int = 1) -> None:
        """Record calls that actually finished.

        A batch that raises part way through completes fewer than it
        reserved, and the incident must be able to say so. Completing
        more than was reserved is refused: that would be a call the
        budget never charged for, which is exactly what the pre-charge
        rule exists to prevent."""

        if k < 0:
            raise ValueError(f"cannot complete {k} calls")
        if self._completed + k > self._reserved:
            raise ValueError(
                f"completing {k} would put completed at "
                f"{self._completed + k} against {self._reserved} "
                f"reserved: a completed call must have been charged")
        self._completed += k

    def state(self) -> dict:
        """What both the result and the incident report."""

        return {
            "reserved": self._reserved,
            "completed": self._completed,
            "unfinished": self.unfinished,
            "wall_s": self.wall_s,
            "wall_overrun_s": self.wall_overrun_s,
            "max_calls": self._max_calls,
            "max_wall_s": self._max_wall_s,
            "max_batch_calls": self._max_batch_calls,
            "reserved_by_stage": dict(self.per_stage),
            "stage": self.stage,
            "clock_started_at": "before G3a",
            "cap_is_charged_against": ("reserved, which never "
                                       "under-charges"),
            "call_cap_overshoot": (
                "none: a batch that does not fit is refused, so "
                "max_calls is never crossed"),
            "wall_cap_is": (
                "a start rule: no batch begins at or after max_wall_s, "
                "so the run can finish one batch late -- at most "
                f"{self._max_batch_calls:,} calls, but NO rigorous "
                "bound in seconds, because the 768 us/pair cost is a "
                "host-dependent average and not a per-call maximum"),
            "nominal_overrun_s": self.nominal_overrun_s,
            "nominal_overrun_is": ("a projection from the measured "
                                   "average; the recorded "
                                   "wall_overrun_s is the fact"),
        }
