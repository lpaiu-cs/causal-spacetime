"""Continuous score induced by the frozen reconstruction gates."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np


def minimum_gate_margin(
    *,
    heldout: float,
    heldout_max: float,
    null_gap: float | None = None,
    null_gap_min: float | None = None,
    truth_error: float | None = None,
    truth_error_max: float | None = None,
) -> float:
    """Return the limiting normalized margin for the supplied frozen gates."""

    pairs = (
        (null_gap, null_gap_min, "null-gap"),
        (truth_error, truth_error_max, "truth-error"),
    )
    if heldout_max <= 0.0:
        raise ValueError("gate thresholds must be positive")
    margins = [(heldout_max - float(heldout)) / heldout_max]
    for value, threshold, name in pairs:
        if (value is None) != (threshold is None):
            raise ValueError(f"{name} value and threshold must be supplied together")
        if value is None:
            continue
        if threshold <= 0.0:
            raise ValueError("gate thresholds must be positive")
        if name == "null-gap":
            margins.append((float(value) - threshold) / threshold)
        else:
            margins.append((threshold - float(value)) / threshold)
    values = np.asarray(margins, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("gate quantities must be finite")
    return float(np.min(values))


def clipped_gate_margin_score(margins: Iterable[float]) -> float:
    """Map the limiting normalized gate margin to a score in ``[0, 1]``.

    A margin is zero on its frozen gate, positive on the passing side, and
    normalized by that gate's threshold. The minimum retains the logical AND
    of the frozen gates. The fixed offset places the decision boundary at
    0.5 and saturates once the limiting margin reaches +/- 0.5.
    """

    values = np.asarray(list(margins), dtype=float)
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("margins must be a non-empty finite one-dimensional list")
    return float(np.clip(0.5 + np.min(values), 0.0, 1.0))


def geometry_order_parameter(
    *,
    status: str,
    heldout: float | None = None,
    null_gap: float | None = None,
    truth_error: float | None = None,
    heldout_max: float = 0.10,
    null_gap_min: float = 0.10,
    truth_error_max: float = 0.40,
) -> float:
    """Return P7's total geometry score for one reconstruction attempt.

    Structural failures map to zero. For a numerically valid reconstruction,
    ``G >= 0.5`` is exactly equivalent to passing all three frozen P3/P5
    gates. Truth coordinates evaluate the order-intrinsic reconstruction but
    are never inputs to the fit itself.
    """

    if status != "ok":
        return 0.0
    if heldout is None or null_gap is None or truth_error is None:
        raise ValueError("valid reconstructions require all three gate quantities")
    margin = minimum_gate_margin(
        heldout=heldout,
        heldout_max=heldout_max,
        null_gap=null_gap,
        null_gap_min=null_gap_min,
        truth_error=truth_error,
        truth_error_max=truth_error_max,
    )
    return clipped_gate_margin_score((margin,))
