"""One counting wrapper, on the one function that costs.

Every solver path in O4b is metered by patching `flight_time` itself
for the duration of a run. Nothing has to remember to route through a
helper: `_ell`, `causal_relation`, and the G3a table's bisections all
reach the solver through the same module attribute, so patching it
catches them all. A path that escaped the meter would be a call the
budget never charged for, and the whole point of pre-charging is that
no such call exists.

WHY `flight_time` AND NOT `causal_relation`. Patching both would
double-count: one `causal_relation` makes exactly one `flight_time`
call. Metering the solver is also the honest unit -- it is what the
sizing model counts, and it charges the `dt < 0` short circuit
correctly at zero, because that path never consults the solver.

WHY PER CALL AND NOT PER BATCH. `Budget` permits a bounded batch, and
the runner could reserve 8,192 at a time to match G1's chunk. It does
not. Reserving one call at a time makes the wall cap's overrun ONE
CALL rather than one batch, which is the strongest form of a rule that
has no rigorous bound in seconds.

The price is one clock read and one counter update per call, tens of
millions of times. That overhead is UNMEASURED, and the wall-time
projection does not include it.

The order is fixed: reserve, then call, then complete only if the call
returned. A solver exception leaves the reservation charged and the
completion unrecorded, which is what the incident needs to be able to
say.
"""

from __future__ import annotations

import sys
from contextlib import contextmanager
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                       / "positive_control"))

import s1_schwarzschild_cost as s1  # noqa: E402


@contextmanager
def metered(budget):
    """Charge `budget` for every `flight_time` call made inside.

    Restores the original on the way out, including when the body
    raises -- a half-installed meter on a shared module attribute
    would corrupt every later run in the same process."""

    original = s1.flight_time

    def counting(*args, **kwargs):
        budget.reserve(1)
        result = original(*args, **kwargs)
        budget.complete(1)
        return result

    s1.flight_time = counting
    try:
        yield budget
    finally:
        s1.flight_time = original


def unmetered_solver():
    """The real `flight_time`, for code that must not be charged --
    there is deliberately no such caller in the runner, and this
    exists so tests can prove the meter is installed rather than
    assumed."""

    return s1.flight_time
