"""Unit tests synchronizing the frozen S1 confirmatory runner and verifier."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(ROOT / "experiments" / "positive_control"))

import check_subspacing_confirmatory as checker  # noqa: E402
import subspacing_confirmatory as runner  # noqa: E402
from seed_windows import assert_windows_disjoint_and_fresh  # noqa: E402


def test_s1_seed_blocks_fresh_and_constant_parity():
    blocks = runner.check_seed_freshness()
    assert blocks == checker.FROZEN_BLOCKS
    assert runner.Z_TH == checker.FROZEN_Z_TH == 3.30
    assert runner.LOCK_HALF_WIDTH == checker.FROZEN_LOCK_WIDTH == 8.0
    np.testing.assert_allclose(runner.FREQUENCIES, checker.FROZEN_FREQUENCIES)
    jobs = runner.build_jobs()
    assert len(jobs) == 216
    # Injecting a spent exploratory or P14 seed must abort via SystemExit.
    bad_blocks = list(blocks) + [("collision", 750_000, 750_010)]
    with pytest.raises(SystemExit):
        assert_windows_disjoint_and_fresh(
            bad_blocks, runner.ALL_SPENT_RANGES, 0, 52_000_000, "test"
        )


def test_s1_ast_blindness_and_quad_oracle_and_stats_parity():
    ast_res = checker.audit_observer_ast()
    assert ast_res["passed"]
    oracle_res = checker.independent_rosen_order_check(seed=51_400_000, n_pts=40)
    assert oracle_res["passed"] and oracle_res["mismatches"] == 0
    # Synthetic raw dataset verifying identical Wilson/Welch/t gate evaluation
    rng = np.random.default_rng(1)
    runs = []
    for i in range(24):
        runs.append(
            dict(
                arm="N22",
                kind="subspacing",
                seed=i,
                whitened_peak_z=2.3,
                whitened_peak_frequency=24.0,
                q_hat=4.01,
                q_abs_err=0.01,
                calibration_rms=4e-5,
                c_ref_estimate=0.01,
                detected=False,
                freq_locked=False,
            )
        )
        runs.append(
            dict(
                arm="N23",
                kind="subspacing",
                seed=100 + i,
                whitened_peak_z=2.35,
                whitened_peak_frequency=28.0,
                q_hat=4.005,
                q_abs_err=0.005,
                calibration_rms=3e-5,
                c_ref_estimate=0.01,
                detected=False,
                freq_locked=False,
            )
        )
    for idx_arm, (arm, n, r_val, pk) in enumerate(
        (
            ("W1", 36, 64.0, 4.6),
            ("W1p", 36, 64.0, 4.5),
            ("W2", 32, 96.0, 4.4),
            ("W3", 20, 128.0, 3.5),
            ("X_sub_q1_r32", 16, 32.0, 5.5),
        )
    ):
        for i in range(n):
            z = pk + float(rng.normal(0, 0.2))
            runs.append(
                dict(
                    arm=arm,
                    kind="subspacing",
                    seed=10_000 * (idx_arm + 1) + i,
                    whitened_peak_z=z,
                    whitened_peak_frequency=r_val,
                    q_hat=4.01,
                    q_abs_err=0.01,
                    calibration_rms=4e-5,
                    c_ref_estimate=0.01,
                    detected=z >= 3.30,
                    freq_locked=True,
                )
            )
    for i in range(12):
        runs.append(
            dict(
                arm="X_null_q1",
                kind="subspacing",
                seed=9000 + i,
                whitened_peak_z=2.1,
                whitened_peak_frequency=40.0,
                q_hat=1.01,
                q_abs_err=0.01,
                calibration_rms=4e-5,
                c_ref_estimate=0.01,
                detected=False,
                freq_locked=False,
            )
        )
    for i in range(8):
        runs.append(
            dict(
                arm="X_weyl_r4",
                kind="weyl_jet",
                seed=8000 + i,
                electric_ratio=0.81 + 0.01 * (i % 2),
                c_norm=2.2,
            )
        )
        runs.append(
            dict(
                arm="X_weyl_r32",
                kind="weyl_jet",
                seed=8100 + i,
                electric_ratio=0.02 + 0.005 * (i % 2),
                c_norm=0.6,
            )
        )
    summary = runner.evaluate_gates(runs)
    assert summary["verdict"] == "CONFIRMED"
    rec = checker.recompute_and_verify_gates({"summary": summary, "runs": runs})
    assert rec["verdict"] == "CONFIRMED"
    assert rec["disagreements"] == []
    # Also ensure frozen summary JSON if present parses cleanly
    frozen_res = (
        ROOT / "docs" / "prereg" / "frozen" / "subspacing_confirmatory" / "results.json"
    )
    if frozen_res.exists():
        checker.recompute_and_verify_gates(json.loads(frozen_res.read_text("utf-8")))
