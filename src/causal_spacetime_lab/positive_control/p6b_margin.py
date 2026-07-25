"""Corrected P6b margins induced by each source's declared frozen gates."""

from __future__ import annotations

from collections.abc import Mapping

from causal_spacetime_lab.positive_control.geometry_score import minimum_gate_margin

BLOCKED_MARGIN = -5.0
SUPPORTED_SOURCES = frozenset({"P1", "P3", "P5", "P6"})


def _metric(row: Mapping[str, object], *names: str) -> float | None:
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return float(value)
    return None


def p6b_instrument_margin(
    source: str,
    row: Mapping[str, object],
    *,
    include_truth: bool = True,
) -> float:
    """Return the P6b margin with all declared order-only gates retained.

    P1's preregistered decision is the conjunction of held-out violation,
    truth-order error, and restart-order disagreement.  The historical P6b
    implementation omitted restart stability.  The descriptive
    ``include_truth=False`` variant removes only the truth-coordinate term;
    it deliberately keeps restart stability because stability is computed
    from repeated fits to the order data.
    """

    if source not in SUPPORTED_SOURCES:
        raise ValueError(f"unsupported P6b source: {source}")
    status = row.get("status")
    if not isinstance(status, str) or not status.strip():
        raise ValueError(f"{source} row has no valid status")
    if status != "ok":
        return BLOCKED_MARGIN

    heldout = _metric(row, "heldout", "heldout_violation")
    if heldout is None:
        raise ValueError(f"{source} row has no heldout value")

    truth = _metric(row, "truth", "truth_order_error")
    if source == "P1":
        stability = _metric(
            row,
            "restart_order_disagreement",
            "stability",
        )
        if stability is None:
            raise ValueError("P1 row has no restart-stability value")
        if include_truth and truth is None:
            raise ValueError("P1 row has no truth-order value")
        return minimum_gate_margin(
            heldout=heldout,
            heldout_max=0.05,
            truth_error=truth if include_truth else None,
            truth_error_max=0.15 if include_truth else None,
            stability_error=stability,
            stability_error_max=0.15,
        )

    if source in {"P5", "P6"} and include_truth and truth is None:
        raise ValueError(f"{source} row has no truth-order value")
    source_truth = truth if include_truth and source in {"P5", "P6"} else None
    null_gap = _metric(row, "null_gap")
    if null_gap is None:
        raise ValueError(f"{source} row has no null-gap value")
    return minimum_gate_margin(
        heldout=heldout,
        heldout_max=0.10,
        null_gap=null_gap,
        null_gap_min=0.10,
        truth_error=source_truth,
        truth_error_max=0.40 if source_truth is not None else None,
    )
