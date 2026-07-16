"""Regression tests pinning the P7 FSS re-scope reconnaissance verdict.

The window-disjointness conclusion (docs/p7_fss_rescope_recon.md) rests
on a few deterministic facts that must not silently drift: the frozen
protocol is envelope-blocked in the samplable range, it operates at
N = 500, the crystal negative stays structurally blocked, and the most
permissive re-scoped spec still fails calibration-grade operability at
the samplable boundary. Seeds are fixed; every assertion is a
deterministic replay of the probe's cells.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p7_fss_rescope_probe import (  # noqa: E402
    FROZEN_SPEC,
    SOFT_SPECS,
    evaluate,
    order_inputs,
)

from causal_spacetime_lab.positive_control.order_intrinsic import (  # noqa: E402
    select_disjoint_chains,
)
from causal_spacetime_lab.positive_control.two_orders import (  # noqa: E402
    bipartite_perm,
)


def test_frozen_protocol_is_envelope_blocked_in_the_samplable_range():
    """LIS ~ 2 sqrt(N) < 25 for N <= ~156: the frozen six-25-tick-chain
    protocol cannot even begin below the envelope, and stays blocked at
    N = 300 (measured, not just asymptotic)."""

    for n in (120, 300):
        for seed in (0, 1):
            rng = np.random.default_rng(9000 + seed)
            result = evaluate(rng.permutation(n), FROZEN_SPEC, seed)
            assert result["status"] == "block_chains", (n, seed, result)


def test_frozen_protocol_operates_at_n_500():
    """The positive control passes where the instrument was designed to
    live: at N = 500 at least one probe seed evaluates and passes."""

    outcomes = []
    for seed in (0, 1, 2):
        rng = np.random.default_rng(9000 + seed)
        outcomes.append(evaluate(rng.permutation(500), FROZEN_SPEC, seed))
    ok = [r for r in outcomes if r["status"] == "ok"]
    assert ok, outcomes
    assert any(r["G"] >= 0.5 for r in ok), ok


def test_crystal_stays_structurally_blocked_for_rescoped_specs():
    """Height-2 orders must never yield chains, at any probed size."""

    for n in (120, 500):
        causal, times, _ = order_inputs(bipartite_perm(n))
        assert len(select_disjoint_chains(causal, times, 3, 8)) == 0


def test_rescoped_spec_fails_calibration_grade_at_the_samplable_boundary():
    """The probe's headline negative, pinned: at N = 120 the most
    permissive spec evaluates but its beta = 0 positive control does not
    reach calibration-grade operability (G median >= 0.5 with >= 3/4
    passing) -- half of it fails, which disqualifies G as an order
    parameter at samplable N."""

    spec = SOFT_SPECS["3x8_t30"]
    results = []
    for seed in range(8):
        rng = np.random.default_rng(9000 + seed)
        results.append(evaluate(rng.permutation(120), spec, seed))
    ok = [r for r in results if r["status"] == "ok"]
    assert len(ok) >= 4, results
    g_values = [r["G"] for r in ok]
    passing = sum(1 for g in g_values if g >= 0.5)
    calibration_grade = (
        float(np.median(g_values)) >= 0.5 and passing * 4 >= len(ok) * 3
    )
    assert not calibration_grade, {
        "g_median": float(np.median(g_values)),
        "passing": passing,
        "ok": len(ok),
    }
