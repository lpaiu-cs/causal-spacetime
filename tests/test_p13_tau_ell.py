"""Regressions for P13 (prereg v1.1, Section 8).

The load-bearing pins: the per-rung packing budget (P12's collapse
class), the m-matching that the whole contrast rests on, the flat
twin's calibration and normalization (both of which were wrong in the
first implementation), the inverted verdict table, and window privacy.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p13_tau_ell import (  # noqa: E402
    DELTA_EQ,
    M_TARGET,
    N_CAP,
    PATCH,
    PILOT_BLOCKS,
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


def test_inverted_verdict_table_is_single_valued():
    """P13's polarity is the reverse of P11/P12: the null is no effect
    and the interesting outcome is POSITIVE."""

    assert verdict(0.01, 0.04, True) == "CURVATURE-DEGRADES"
    assert verdict(-0.04, -0.01, True) == "CURVATURE-HELPS"
    assert verdict(-0.03, 0.03, True) == "CURVATURE-ROBUST"
    assert verdict(-0.09, 0.03, True) == "INCONCLUSIVE"
    assert verdict(-0.03, 0.03, False) == "UNRESOLVED"
    # rows 1 and 2 fire on significance, not size -- the pin the
    # prereg added under the table
    assert verdict(0.001, 0.004, True) == "CURVATURE-DEGRADES"
    assert DELTA_EQ == 0.05 and N_CAP == 300


def test_control_gate_is_an_equivalence_test_with_a_label_row():
    """v2 Section 10: the point threshold that failed a perfect
    control one campaign in five is replaced by a three-way
    equivalence test, evaluated by precedence."""

    # clean: the whole interval sits inside the margin
    assert control_result(-0.03, 0.03) == "CONTROL-CLEAN"
    # confounded: the interval lies wholly outside it
    assert control_result(0.06, 0.20) == "CONFOUNDED"
    assert control_result(-0.20, -0.06) == "CONFOUNDED"
    # the ROW 3 BOUNDARY CASE the review asked for: the interval
    # excludes zero but straddles the margin -- v1 would have called
    # this CONFOUNDED, v2 labels it and lets the verdict stand
    assert control_result(0.03, 0.20) == "UNDERPOWERED-CONTROL"
    # and the campaign v1 reading itself lands on row 3
    assert control_result(0.0318, 0.1355) == "UNDERPOWERED-CONTROL"
    # merely wide but centred is also row 3, not clean
    assert control_result(-0.09, 0.09) == "UNDERPOWERED-CONTROL"


def test_v2_cap_and_windows_moved_off_the_spent_campaign():
    from p13_tau_ell import N_CAP as cap
    assert cap == 300
    for base, _slots in (list(PILOT_BLOCKS.values())
                         + list(PILOT_TWIN_BLOCKS.values())
                         + list(STAGE_A_BLOCKS.values())
                         + list(TWIN_BLOCKS.values())):
        assert base >= 8_000_000      # campaign v1's seeds are spent


def test_windows_are_private_and_clear_of_design_space():
    spans = []
    for base, slots in (list(PILOT_BLOCKS.values())
                        + list(PILOT_TWIN_BLOCKS.values())
                        + list(STAGE_A_BLOCKS.values())
                        + list(TWIN_BLOCKS.values())):
        span = set(range(base, base + STRIDE * slots))
        spans.append(span)
        assert min(span) >= 8_000_000
        assert max(span) < 9_000_000     # 9M+ is design-check space
    for a in spans:
        for b in spans:
            assert a is b or not (a & b)
    for tau_c in RUNGS:
        vspan = set(range(VERIFY_BASE[tau_c], VERIFY_BASE[tau_c] + 2000))
        assert max(vspan) < min(min(s) for s in spans)
