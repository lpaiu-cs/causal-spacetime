"""Frozen constants and arithmetic for the G3 redesign.

The design is in `docs/prereg/p14_o4_g3_redesign.md`; this module is
its numeric half, kept apart from both the frozen campaign runner and
the replay diagnostic so that a constant has ONE home and the document
can be pinned against it.

There is no runner here. G3a and G3b are not implemented yet -- their
sizes and thresholds are frozen only after the census -- but every
value that must NOT be chosen after seeing the census lives here now,
which is the whole point of the split:

  * `ETA` is a float-separation margin (design section 4).
  * `MAX_NUDGES` bounds the placement search (design section 4.3).

Both are determined by rounding arguments alone. Neither may be
revisited in light of a frequency the census reports.
"""

from __future__ import annotations

#: The separation margin, in the same units as the flight times.
#: Chosen so that it dominates the rounding in forming a probe time
#: (all magnitudes <= dt = 8.5, `ulp(8.5) = 1.78e-15`, three or four
#: operations) while sitting six decades below the `err` scale, so it
#: cannot drive eligibility. It is NOT chosen from the frequency of
#: `err2 >= 1e-6`, and it is not a bare addend: the runner must CHECK
#: that the realized margin reaches it (design section 4.2).
ETA = 1e-12

#: The exploratory grid the census scans. Eligibility is strict:
#: `L - err1 - err2 - 2*eta > 0`.
ETA_GRID = (0.0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9,
            1e-8, 1e-7, 1e-6)

#: How many `nextafter` steps the probe placement may take before the
#: construction is declared unavailable, frozen HERE rather than after
#: the census -- the cap decides whether a given cluster yields a valid
#: probe or an unavailable one, so choosing it later would let the
#: census pick which clusters count.
#:
#: The bound comes from rounding alone. A candidate is placed at a
#: nominal margin of exactly `ETA`, and the arithmetic that forms it
#: (`t_x` in two or three operations, then `dt`, then the difference)
#: rounds by at most `0.5 * ulp(8.5) = 8.9e-16` each, so the realized
#: margin can fall short by at most about `3.6e-15`. One step moves
#: `dt` by one ulp of the probe time, worst case `ulp(1.0) = 2.2e-16`,
#: so about 17 steps suffice. 64 leaves a factor of ~3.7.
MAX_NUDGES = 64


def w_robust(length: float, err1: float, err2: float,
             eta: float = ETA) -> float:
    """`L - err1 - err2 - 2*eta`, in the predicate's coordinates."""

    return length - err1 - err2 - 2.0 * eta


def is_eligible(length: float, err1: float, err2: float,
                eta: float = ETA) -> bool:
    """Strictly positive robust window. Strict, not `>=`: a zero window
    admits no probe with a margin to spare."""

    return w_robust(length, err1, err2, eta) > 0.0


def lower_probe_time(lo: float, err1: float,
                     eta: float = ETA) -> float:
    """Probe (iii): strictly below the window's lower edge.

    Split out because eligibility needs it on its own -- whether this
    time is non-negative decides between a lower-boundary solver probe
    and the predicate's negative-`dt` short circuit, which are
    different checks and are counted apart (design section 3.4)."""

    return lo - err1 - eta


def probe_times(lo: float, hi: float, err1: float, err2: float,
                eta: float = ETA) -> dict:
    """The three nominal probe times of design section 3.3.

    Nominal: each is where the algebra puts the probe. Whether the
    REALIZED margin reaches `eta` is a separate check the runner makes
    against the predicate itself, and is not decided here."""

    return {
        "inside": 0.5 * ((lo + err1 + eta) + (hi - err2 - eta)),
        "outside_above": hi + err2 + eta,
        "outside_below": lower_probe_time(lo, err1, eta),
    }
