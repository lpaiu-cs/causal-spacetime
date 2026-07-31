"""Regressions for P13 (prereg as amended through design v3).

The load-bearing pins: the per-rung packing budget (P12's collapse
class), the m-matching that the whole contrast rests on, the flat
twin's calibration and normalization (both of which were wrong in the
first implementation), the inverted verdict table, and window privacy.

v3 adds two that exist because of 11.3: a significant but sub-margin
contrast must no longer earn the word DEGRADES, and the verdict and
the control must key to the SAME margin -- the invariant whose absence
was v2's defect.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p13_tau_ell import (  # noqa: E402
    DELTA_DETECT,
    DELTA_EQ,
    M_TARGET,
    M_TOLERANCE,
    N_CAP,
    N_EQ_COEFF,
    N_SUP_COEFF,
    PATCH,
    PILOT_BLOCKS,
    PILOT_SAMPLES,
    PILOT_TWIN_BLOCKS,
    RUNGS,
    STAGE_A_BLOCKS,
    STRIDE,
    TWIN_BAND_CENTRE,
    TWIN_BLOCKS,
    TWIN_RHO,
    VERIFY_BASE,
    control_result,
    eligible_pairs,
    patch_proper_volume,
    run_sample,
    sprinkle,
    tau_curved,
    verdict,
)


def test_every_rung_packs_six_disjoint_pairs():
    """The P12 collapse class: v1.0's constants completed 100/10/0/0
    percent. Each rung must now complete at design-check resolution."""

    for tau_c in RUNGS:
        done = sum(
            run_sample(tau_c, VERIFY_BASE[tau_c] + k,
                       single_stream=True)[1]
            for k in range(8)
        )
        assert done == 8, tau_c


def test_discreteness_is_held_fixed_across_rungs():
    """The design's whole claim rests on m being constant, so the
    realized counts must land on target at every rung -- otherwise the
    contrast mixes curvature with discreteness."""

    means = {}
    for tau_c in RUNGS:
        ms = [
            r["mean_m"]
            for k in range(8)
            if (r := run_sample(tau_c, VERIFY_BASE[tau_c] + k,
                                single_stream=True)[0]).get("mean_m")
        ]
        means[tau_c] = float(np.mean(ms))
    grand = float(np.mean(list(means.values())))
    assert abs(grand - M_TARGET) / M_TARGET < 0.12
    for tau_c, value in means.items():
        assert abs(value - grand) / grand < 0.10, (tau_c, value, grand)


def test_flat_twin_is_calibrated_and_normalized():
    """Two defects the first implementation had, both caught here.
    The twin's band centre and intensity are MEASURED (a guessed scale
    put its m at 172 and 872 against a target of 76), and its rho uses
    the FLAT proper volume -- the coordinate area -- where using the
    curved patch's volume mis-normalized tau_hat into a 39-70%
    relative error instead of ~18%."""

    for tau_c in RUNGS:
        assert TWIN_RHO[tau_c] == pytest.approx(
            2.0 * M_TARGET / TWIN_BAND_CENTRE[tau_c] ** 2
        )
    for tau_c, base in ((RUNGS[0], TWIN_BLOCKS[RUNGS[0]][0]),
                        (RUNGS[-1], TWIN_BLOCKS[RUNGS[-1]][0])):
        rows = [run_sample(tau_c, base + STRIDE * k, flat=True)
                for k in range(6)]
        ms = [r["mean_m"] for r, ok in rows if ok]
        errs = [r["median_relerr"] for r, ok in rows if ok]
        assert len(ms) == 6
        assert abs(np.mean(ms) - M_TARGET) / M_TARGET < 0.12
        # the twin is flat, so its error is the pure discreteness one
        assert 0.10 < float(np.mean(errs)) < 0.30


def test_flat_twin_proper_volume_is_the_coordinate_area():
    eta_lo, xhalf, _ = PATCH[1.00]
    curved = patch_proper_volume(eta_lo, xhalf)
    flat = (abs(eta_lo) - 1.0) * 2.0 * xhalf
    assert curved < flat            # Omega < 1 over the patch
    assert curved == pytest.approx((1.0 - 1.0 / abs(eta_lo)) * 2.0 * xhalf)


def test_eligibility_enforces_band_and_box_inside_patch():
    tau_c = 0.60
    eta_lo, xhalf, rho = PATCH[tau_c]
    rng = np.random.default_rng(7_400_000)
    eta, x = sprinkle(eta_lo, xhalf, rho, rng)
    band = (tau_c * 0.9, tau_c * 1.1)
    pool = eligible_pairs(eta, x, eta_lo, xhalf, band)
    assert pool.shape[0] > 0
    a, b = pool[:, 0], pool[:, 1]
    tau = tau_curved(eta[a], x[a], eta[b], x[b])
    assert np.all(tau >= band[0]) and np.all(tau <= band[1])
    u, v = eta + x, eta - x
    e1, x1 = (u[b] + v[a]) / 2.0, (u[b] - v[a]) / 2.0
    e2, x2 = (u[a] + v[b]) / 2.0, (u[a] - v[b]) / 2.0
    for e, xx in ((e1, x1), (e2, x2)):
        assert np.all(e >= eta_lo - 1e-12) and np.all(e <= -1.0 + 1e-12)
        assert np.all(np.abs(xx) <= xhalf + 1e-12)


def test_v3_verdict_table_keys_rows_1_and_2_to_the_margin():
    """Section 12.1's repair. P13's polarity is still the reverse of
    P11/P12 -- the null is no effect and the interesting outcome is
    POSITIVE -- but rows 1 and 2 now require clearing delta_eq instead
    of clearing zero."""

    assert verdict(0.03, 0.05, True) == "CURVATURE-DEGRADES"
    assert verdict(-0.05, -0.03, True) == "CURVATURE-HELPS"
    assert verdict(-0.01, 0.01, True) == "CURVATURE-ROBUST"
    assert verdict(-0.05, 0.01, True) == "INCONCLUSIVE"
    assert verdict(-0.01, 0.01, False) == "UNRESOLVED"
    assert DELTA_EQ == 0.02 and N_CAP == 3000


def test_a_significant_but_sub_margin_contrast_no_longer_says_degrades():
    """THE case v3 exists to catch, so it is tested directly rather
    than implied. Under v2 an interval excluding zero earned the word
    DEGRADES at any size, which is how a +0.024 contrast -- under half
    the margin, with its whole interval inside the margin -- came to
    carry a mechanism it could not support (11.2)."""

    # excludes zero, sits entirely inside the margin
    assert verdict(0.001, 0.004, True) == "CURVATURE-ROBUST"
    assert verdict(0.005, 0.015, True) == "CURVATURE-ROBUST"
    # excludes zero, straddles the margin edge: no verdict, not a word
    assert verdict(0.005, 0.030, True) == "INCONCLUSIVE"
    # and only a contrast wholly past the margin earns DEGRADES
    assert verdict(0.021, 0.030, True) == "CURVATURE-DEGRADES"

    # Campaign v2's frozen interval, for the record and NOT as a
    # re-scoring: it was judged under v2's rules and stands as
    # CURVATURE-DEGRADES (Section 12.1's closing paragraph). What this
    # line pins is that the repair targets the real defect.
    assert verdict(0.004444, 0.045607, True) == "INCONCLUSIVE"


def test_m_gate_covers_the_twin_arm():
    """Review C5. Section 5 gives the twin "the same eligibility
    conditions, K, rejection cap, fill rule and m-gate", and the first
    implementation gated only the curved rows -- so a twin rung could
    drift in discreteness and still have its contrast used as a
    control. Tested against the case it is meant to catch: a twin arm
    whose top rung drifts past tolerance must fail the gate even when
    the curved arm is perfectly matched."""

    def gate(curved, twin):
        """The Section 5 gate as run_stage_a applies it: each arm
        against its own grand mean."""
        out = []
        for means in (curved, twin):
            grand = sum(means) / len(means)
            out.append(all(abs(m - grand) / grand <= M_TOLERANCE
                           for m in means))
        return out[0] and out[1]

    flat = [76.0, 76.0, 76.0, 76.0]
    assert gate(flat, flat)
    # curved clean, twin drifting: must NOT pass
    assert not gate(flat, [76.0, 76.0, 76.0, 90.0])
    assert not gate(flat, [60.0, 76.0, 76.0, 76.0])
    # a pure LEVEL offset between the arms is not drift and must pass:
    # the gate exists to stop within-arm confounding, and the campaigns
    # ran with the twin about 1.7% below the curved arm
    assert gate(flat, [74.6, 74.6, 74.6, 74.6])
    # the campaigns' own realized figures, both arms
    assert gate([75.64, 75.79, 76.12, 76.00],
                [75.24, 74.66, 74.35, 74.05])       # v3
    assert gate([75.57, 75.58, 75.98, 75.89],
                [74.91, 74.72, 73.96, 74.16])       # v2


def test_equivalence_coefficient_is_two_sided():
    """Review C2, and the guard is written the way the rule demands --
    it derives the DELIVERED power back out of the constant instead of
    trusting the label beside it.

    ROBUST needs BOTH bounds inside the margin, so at Delta = 0 the
    probability is P(|Z| < z) = 2 Phi(z) - 1 with
    delta_eq / se = 1.960 + z. v3 first wrote z = Phi^-1(0.99) = 2.326
    into that slot, which delivers 98.0% while the design claimed 99%.
    The 90% slot's 1.645 was already Phi^-1(0.95), i.e. already
    two-sided, so the slip broke a convention that was right."""

    def phi(z):
        return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))

    z = math.sqrt(N_EQ_COEFF) * DELTA_EQ - 1.960
    assert z == pytest.approx(2.576, abs=1e-3)
    assert 2.0 * phi(z) - 1.0 == pytest.approx(0.99, abs=5e-4)
    # the one-sided quantile is what this test exists to reject
    assert 2.0 * phi(2.326) - 1.0 == pytest.approx(0.98, abs=5e-4)


def test_single_margin_invariant():
    """v2's resolution mismatch existed because the verdict keyed to 0
    while the control keyed to delta_eq. The invariant that would have
    caught it: on the positive side the verdict's row 1 and the
    control's CONFOUNDED branch must be the SAME predicate, and both
    sizings must key to the same margin. A guard is tested against the
    case it is meant to catch."""

    grid = np.linspace(-0.06, 0.06, 49)
    for lo in grid:
        for hi in grid:
            if hi < lo:
                continue
            degrades = verdict(lo, hi, True) == "CURVATURE-DEGRADES"
            confounded_high = (control_result(lo, hi) == "CONFOUNDED"
                               and lo > 0)
            assert degrades == confounded_high, (lo, hi)
    # the equivalence sizing reads the same margin, so lowering
    # delta_eq cannot silently leave the sample size behind
    assert N_EQ_COEFF == pytest.approx(
        (1.960 + 2.576) ** 2 / DELTA_EQ ** 2
    )
    # and the superiority requirement measures against the margin too,
    # because row 1 now clears it rather than clearing zero
    assert N_SUP_COEFF == pytest.approx(
        (1.960 + 1.282) ** 2 / (DELTA_DETECT - DELTA_EQ) ** 2
    )


def test_control_gate_is_an_equivalence_test_with_a_label_row():
    """v2 Section 10: the point threshold that failed a perfect
    control one campaign in five is replaced by a three-way
    equivalence test, evaluated by precedence."""

    # clean: the whole interval sits inside the margin
    assert control_result(-0.01, 0.01) == "CONTROL-CLEAN"
    # confounded: the interval lies wholly outside it
    assert control_result(0.03, 0.20) == "CONFOUNDED"
    assert control_result(-0.20, -0.03) == "CONFOUNDED"
    # the ROW 3 BOUNDARY CASE the review asked for: the interval
    # excludes zero but straddles the margin -- a point threshold
    # would call this a failure, the three-way test labels it and
    # lets the verdict stand
    assert control_result(0.01, 0.20) == "UNDERPOWERED-CONTROL"
    # merely wide but centred is also row 3, not clean
    assert control_result(-0.04, 0.04) == "UNDERPOWERED-CONTROL"
    # Campaign v2's own twin interval, [-0.0194, +0.0328] at n = 179,
    # was CONTROL-CLEAN against delta_eq = 0.05 and is only
    # UNDERPOWERED-CONTROL against 0.02: it no longer fits inside the
    # tighter margin, and it does not lie outside it either. That is
    # the whole argument for raising the power target alongside
    # lowering the margin (Section 12.3) -- tightening delta_eq without
    # buying precision converts clean controls into underpowered ones.
    # v2's record is judged at its own stamp and is not re-scored here.
    assert control_result(-0.019403, 0.032829) == "UNDERPOWERED-CONTROL"


def test_v3_cap_and_windows_moved_off_both_spent_campaigns():
    """Section 12.4. v1 spent 6000000-6475999, v2 spent
    8000000-8795999, and 11.3's diagnostic reached 14550000 in
    design-check space, so v3 starts above all three."""

    from p13_tau_ell import N_CAP as cap
    assert cap == 3000
    for base, _slots in (list(PILOT_BLOCKS.values())
                         + list(PILOT_TWIN_BLOCKS.values())
                         + list(STAGE_A_BLOCKS.values())
                         + list(TWIN_BLOCKS.values())):
        assert base >= 15_000_000
    for base in VERIFY_BASE.values():
        assert base >= 15_000_000
    # The reserve rule applies per block to what that block actually
    # draws: pilots draw PILOT_SAMPLES, Stage A and the twin draw up to
    # the cap. Raising the cap without widening the Stage A blocks is
    # the failure this catches.
    for _base, slots in (list(PILOT_BLOCKS.values())
                         + list(PILOT_TWIN_BLOCKS.values())):
        assert slots >= PILOT_SAMPLES + 20
    for _base, slots in (list(STAGE_A_BLOCKS.values())
                         + list(TWIN_BLOCKS.values())):
        assert slots >= cap + 20


def test_windows_are_private_and_clear_of_design_space():
    spans = []
    for base, slots in (list(PILOT_BLOCKS.values())
                        + list(PILOT_TWIN_BLOCKS.values())
                        + list(STAGE_A_BLOCKS.values())
                        + list(TWIN_BLOCKS.values())):
        span = set(range(base, base + STRIDE * slots))
        spans.append(span)
        assert min(span) >= 15_000_000
        # 22M+ is design-check space from v3 onward (Section 12.4)
        assert max(span) < 22_000_000
    for a in spans:
        for b in spans:
            assert a is b or not (a & b)
    for tau_c in RUNGS:
        vspan = set(range(VERIFY_BASE[tau_c], VERIFY_BASE[tau_c] + 2000))
        assert max(vspan) < min(min(s) for s in spans)
        for span in spans:
            assert not (vspan & span)


def test_v3_windows_are_clear_of_every_spent_range():
    """The freshness rule stated against the actual documented spans
    rather than against a floor, so a future edit that walks a block
    backwards into spent seeds fails here."""

    spent = [(6_000_000, 6_475_999),      # campaign v1 experimental
             (7_000_000, 7_999_999),      # v1 design-check
             (8_000_000, 8_795_999),      # campaign v2 experimental
             (9_000_000, 14_550_000)]     # 11.3's diagnostic
    blocks = [(b, STRIDE * s) for b, s in (
        list(PILOT_BLOCKS.values()) + list(PILOT_TWIN_BLOCKS.values())
        + list(STAGE_A_BLOCKS.values()) + list(TWIN_BLOCKS.values()))]
    blocks += [(b, 2000) for b in VERIFY_BASE.values()]
    for base, width in blocks:
        for lo, hi in spent:
            assert base > hi or base + width - 1 < lo, (base, lo, hi)
