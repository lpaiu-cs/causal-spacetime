"""Regressions for the P9 paired decision rule.

All synthetic. The confirmatory seeds are used exactly once, so the
decision function has to be right BEFORE that run; these tests are what
that requirement rests on.
"""

from __future__ import annotations

import sys
from pathlib import Path

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p9_paired_rule import CEILING, PASS_MIN, SATURATION_TOL, decide  # noqa: E402


def _row(d2t=0.17, d3t=0.12, d4t=0.12, h3=0.07, ctrl=0.16, status="ok"):
    return {
        "d2_truth": d2t, "d3_truth": d3t, "d4_truth": d4t,
        "d3_heldout": h3, "control_d3_heldout": ctrl,
        "control_status": status,
    }


def test_clean_knee_supports_everything():
    rows = [_row() for _ in range(20)]
    out = decide(rows, 20)
    assert out["underfit_count"] == 20
    assert out["improvement_count"] == 0
    assert out["paired_spec_count"] == 20
    assert out["h_knee_supported"]
    assert out["h_spec_paired_supported"]
    assert out["p9_supported"]


def test_underfit_threshold_is_sharp_at_sixteen():
    """Sixteen underfits of twenty pass; fifteen of twenty do not. Both
    sides use exactly twenty rows so the count, not the denominator, is
    what moves."""

    inverted = [_row(d2t=0.10, d3t=0.12) for _ in range(4)]
    at_16 = [_row() for _ in range(16)] + inverted
    at_15 = [_row() for _ in range(15)] + inverted + [_row(d2t=0.10, d3t=0.12)]
    assert len(at_16) == len(at_15) == 20
    assert decide(at_16, 20)["p9_supported"]
    assert not decide(at_15, 20)["p9_supported"]


def test_five_improving_seeds_reject_the_knee():
    """t4 beating t3 by more than the tolerance on 5 of 20 means the
    error is still falling past d = 3 -- no knee."""

    improving = [_row(d4t=0.12 - SATURATION_TOL - 0.01) for _ in range(5)]
    steady = [_row() for _ in range(15)]
    out = decide(steady + improving, 20)
    assert out["improvement_count"] == 5
    assert not out["h_knee_no_improvement_not_rejected"]
    assert not out["p9_supported"]
    # exactly four is still a non-rejection
    assert decide(steady + improving[:4] + [_row()], 20)["p9_supported"]


def test_improvement_within_tolerance_does_not_count():
    """d4 slightly better than d3 is expected from extra parameters and
    must not count against the knee."""

    rows = [_row(d4t=0.12 - SATURATION_TOL + 0.001) for _ in range(20)]
    assert decide(rows, 20)["improvement_count"] == 0


def test_structural_block_counts_as_a_spec_pass():
    """A control that cannot even build is the strongest possible block."""

    rows = [_row(status="structural_block: pool failure") for _ in range(20)]
    assert decide(rows, 20)["paired_spec_count"] == 20


def test_nan_control_with_ok_status_is_not_a_spec_pass():
    rows = [_row(ctrl=float("nan")) for _ in range(20)]
    assert decide(rows, 20)["paired_spec_count"] == 0


def test_invalid_seeds_tighten_the_geq_legs():
    """Thresholds stay absolute while the denominator shrinks: 17 valid
    seeds must still produce 16 underfits, and 15 valid seeds cannot."""

    rows = [_row() for _ in range(17)]
    out = decide(rows, 20)
    assert out["valid_seed_count"] == 17
    assert out["p9_supported"]                # 17 >= 16 still reachable
    assert not decide(rows[:15], 20)["p9_supported"]   # 15 < 16 cannot pass


def test_invalid_seeds_loosen_the_leq_leg_and_that_asymmetry_is_pinned():
    """The absolute cap of 4 improving seeds does NOT get harder as the
    denominator shrinks -- fewer seeds mean fewer chances to observe an
    improving one, so the equivalence leg is statistically EASIER with
    losses. The rule's first docstring claimed 'strictly harder, never
    easier' across all legs; that was wrong for this one, and this test
    exists so the asymmetry stays visible rather than being re-assumed.

    Concretely: 4 improving seeds is a pass whether they sit among 20
    valid seeds or among 17 -- the cap does not scale down with the
    denominator.
    """

    improving = [_row(d4t=0.12 - SATURATION_TOL - 0.01) for _ in range(4)]
    among_20 = [_row() for _ in range(16)] + improving
    among_17 = [_row() for _ in range(13)] + improving
    assert decide(among_20, 20)["h_knee_no_improvement_not_rejected"]
    out = decide(among_17, 20)
    assert out["valid_seed_count"] == 17
    assert out["improvement_count"] == 4
    assert out["h_knee_no_improvement_not_rejected"]   # same cap, fewer seeds


def test_the_ceiling_is_reported_and_gates_nothing():
    """Every seed above the 0.10 ceiling: H-CEILING reports 0.0 and P9
    can still pass, because the ceiling is a diagnostic by prereg."""

    rows = [_row(h3=CEILING + 0.02) for _ in range(20)]
    out = decide(rows, 20)
    assert out["h_ceiling_fraction"] == 0.0
    assert out["p9_supported"]


def test_the_equivalence_leg_is_labelled_a_non_rejection():
    out = decide([_row() for _ in range(20)], 20)
    assert out["h_knee_no_improvement_is_a_non_rejection"]
    assert PASS_MIN == 16
