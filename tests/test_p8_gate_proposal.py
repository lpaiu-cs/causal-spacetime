"""Regressions for the P8 Stage A gate proposal.

Both cases below are review findings on the first version of
`propose_gates`, and both would have been silent in production: the
Stage A that actually ran overlapped on both axes at once, which is the
one configuration that hides the first bug, and its held-out maximum sat
0.001 under the ceiling, which is what hid the second.
"""

from __future__ import annotations

import sys
from pathlib import Path

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p8_3plus1d import propose_gates  # noqa: E402


def _row(d2t, d3t, d3h, ctrl):
    return {
        "d2_truth": d2t, "d3_truth": d3t, "d3_heldout": d3h,
        "control_d3_heldout": ctrl, "control_status": "ok",
    }


def test_both_axes_overlapping_places_no_gate_on_either():
    """The configuration Stage A actually produced."""

    rows = [_row(0.133, 0.140, 0.099, 0.095), _row(0.190, 0.110, 0.030, 0.200)]
    out = propose_gates(rows)
    assert out["clusters_overlap_truth"]
    assert out["clusters_overlap_heldout"]
    assert out["proposed_gate_truth"] is None
    assert out["proposed_gate_heldout"] is None


def test_heldout_overlap_suppresses_only_the_heldout_gate():
    """Review finding 1. Truth separates, held-out overlaps: the first
    version keyed both gates on the truth overlap alone and would have
    emitted a held-out gate from the midpoint of two OVERLAPPING
    clusters. Now the axes are independent."""

    rows = [_row(0.170, 0.120, 0.099, 0.095), _row(0.200, 0.110, 0.050, 0.150)]
    out = propose_gates(rows)
    assert not out["clusters_overlap_truth"]
    assert out["clusters_overlap_heldout"]
    assert out["proposed_gate_truth"] is not None
    assert out["proposed_gate_heldout"] is None
    assert "held-out clusters overlap" in out["note"]


def test_ceiling_capping_below_the_pass_cluster_is_flagged():
    """Review finding 2. When the 0.10 ceiling pulls the gate below the
    pass cluster's own maximum, every seed fails H-SENS-3D. Permitted --
    the ceiling is a standard, not a knob -- but never silent."""

    rows = [_row(0.200, 0.120, 0.115, 0.180), _row(0.210, 0.110, 0.060, 0.190)]
    out = propose_gates(rows)
    assert not out["clusters_overlap_heldout"]
    assert out["proposed_gate_heldout"] == 0.10
    assert out["heldout_ceiling_excludes_pass_cluster"]
    assert "every seed would fail" in out["note"]


def test_clean_separation_emits_both_gates_without_flags():
    rows = [_row(0.180, 0.100, 0.060, 0.150), _row(0.200, 0.090, 0.040, 0.160)]
    out = propose_gates(rows)
    assert out["proposed_gate_truth"] == 0.14
    assert out["proposed_gate_heldout"] == 0.10
    assert "heldout_ceiling_excludes_pass_cluster" not in out
    assert "note" not in out


def test_the_recorded_stage_a_reproduces_no_gate_from_the_frozen_table():
    """The shipped Stage A artifact, pushed through the fixed logic,
    still places no gate -- the fix changes bookkeeping and guards, not
    this outcome."""

    import csv

    table = Path(__file__).resolve().parents[1] / "docs" / "prereg" / \
        "frozen" / "p8_stage_a_calibration.csv"
    rows = [
        {
            "d2_truth": float(r["d2_truth"]),
            "d3_truth": float(r["d3_truth"]),
            "d3_heldout": float(r["d3_heldout"]),
            "control_d3_heldout": float(r["control_d3_heldout"]),
            "control_status": r["control_status"],
        }
        for r in csv.DictReader(open(table, encoding="utf-8"))
        if r["status"] == "ok"
    ]
    assert len(rows) == 20
    out = propose_gates(rows)
    assert out["clusters_overlap_truth"]
    assert out["clusters_overlap_heldout"]
    assert out["proposed_gate_truth"] is None
    assert out["proposed_gate_heldout"] is None
