"""Cross-slice predicate reports for order-first diagnostics."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CrossSlicePredicateReport:
    """Report whether a cross-slice predicate is defined by available protocol."""

    predicate_name: str
    defined: bool
    reason: str
    value: float | bool | None


def same_position_without_transport() -> CrossSlicePredicateReport:
    """Report that same-position is undefined without cross-slice transport."""

    return CrossSlicePredicateReport(
        predicate_name="same_position",
        defined=False,
        reason=(
            "Cross-slice same-position is undefined without a supplied "
            "transport, anchor, persistence, or calibration rule."
        ),
        value=None,
    )


def same_direction_without_transport() -> CrossSlicePredicateReport:
    """Report that same-direction is undefined without cross-slice transport."""

    return CrossSlicePredicateReport(
        predicate_name="same_direction",
        defined=False,
        reason=(
            "Cross-slice same-direction is undefined without an orientation "
            "transport or equivalent reference protocol."
        ),
        value=None,
    )


def velocity_without_transport() -> CrossSlicePredicateReport:
    """Report that velocity is undefined without persistence and transport."""

    return CrossSlicePredicateReport(
        predicate_name="velocity",
        defined=False,
        reason=(
            "Velocity requires object persistence, time calibration, and a "
            "cross-slice spatial identification rule; without those structures "
            "it is undefined rather than false."
        ),
        value=None,
    )


def spatial_evolution_without_transport() -> CrossSlicePredicateReport:
    """Report that spatial evolution is undefined without cross-slice transport."""

    return CrossSlicePredicateReport(
        predicate_name="spatial_evolution",
        defined=False,
        reason=(
            "Comparing spatial structures across slices requires transport, "
            "anchors, persistence, calibration, or relational invariant "
            "structure. Without such structure, spatial evolution is undefined."
        ),
        value=None,
    )
