"""Continuous score induced by the frozen reconstruction gates."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np


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
    if heldout_max <= 0.0 or null_gap_min <= 0.0 or truth_error_max <= 0.0:
        raise ValueError("gate thresholds must be positive")
    return clipped_gate_margin_score(
        (
            (heldout_max - float(heldout)) / heldout_max,
            (float(null_gap) - null_gap_min) / null_gap_min,
            (truth_error_max - float(truth_error)) / truth_error_max,
        )
    )
