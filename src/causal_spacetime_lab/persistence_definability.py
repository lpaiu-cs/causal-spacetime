"""Definability reports for persistence and identity predicates."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PersistencePredicateReport:
    """Report whether a persistence-dependent predicate is defined."""

    predicate_name: str
    defined: bool
    reason: str
    value: float | bool | None


def object_identity_without_persistence() -> PersistencePredicateReport:
    """Report that cross-slice object identity is undefined without persistence."""

    return PersistencePredicateReport(
        predicate_name="object_identity",
        defined=False,
        reason=(
            "Cross-slice object identity is undefined without supplied "
            "persistence labels, a matching hypothesis, anchors, or dynamics."
        ),
        value=None,
    )


def pair_distance_history_without_persistence() -> PersistencePredicateReport:
    """Report that pair-distance histories require persistence or a hypothesis."""

    return PersistencePredicateReport(
        predicate_name="pair_distance_history",
        defined=False,
        reason=(
            "Pair-distance order history requires persistent object labels or "
            "an explicit persistence hypothesis linking objects across slices."
        ),
        value=None,
    )


def persistence_hypothesis_report(
    hypothesis_name: str,
    supplied: bool,
    inferred: bool,
    ambiguous: bool,
) -> PersistencePredicateReport:
    """Report that relational history is defined relative to a hypothesis."""

    status = "supplied" if supplied else "inferred" if inferred else "unspecified"
    suffix = " The hypothesis is ambiguous." if ambiguous else ""
    return PersistencePredicateReport(
        predicate_name="persistence_hypothesis",
        defined=bool(supplied or inferred),
        reason=(
            f"Relational history is defined relative to the {status} "
            f"persistence hypothesis '{hypothesis_name}', not absolutely."
            f"{suffix}"
        ),
        value=None if ambiguous else bool(supplied or inferred),
    )
