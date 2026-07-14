"""Focused guards for the frozen P6 hard-negative protocol."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p6_hard_negatives import (  # noqa: E402
    GATE_HELDOUT,
    GATE_NULL_GAP,
    GATE_TRUTH,
    LAYER_JITTER_WINDOW,
    N_ELEMENTS,
    _validate_frozen_constants,
)


def _frozen_constants() -> dict:
    return {
        "n_elements": N_ELEMENTS,
        "gate_heldout": GATE_HELDOUT,
        "gate_nullgap": GATE_NULL_GAP,
        "gate_truth": GATE_TRUTH,
        "layer_jitter_window": LAYER_JITTER_WINDOW,
    }


def test_module_constants_match_frozen_protocol():
    frozen = _frozen_constants()
    assert _validate_frozen_constants(frozen) is frozen


def test_frozen_protocol_drift_is_rejected():
    frozen = _frozen_constants()
    frozen["gate_truth"] = GATE_TRUTH + 0.01
    with pytest.raises(RuntimeError, match="gate_truth"):
        _validate_frozen_constants(frozen)
