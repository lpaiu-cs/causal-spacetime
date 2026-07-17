"""Regression tests for the G2 density-coupled tick protocols (T1 v0.5).

Constructor invariants are exact -- a failure means the instrumentation
is broken, not noisy. The scaling assertions run a reduced grid and pin
only what Section 5 marks [PROVED] (the thinned protocol's exact error
prediction and its rho^{-1/2} law) plus generous sanity bands for the
measured harvested-chain characterization, whose sharp exponent the
document deliberately leaves open.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "theory"
sys.path.insert(0, str(EXPERIMENT_DIR))

import pytest  # noqa: E402
from t1_g2_density_scaling import (  # noqa: E402
    audit_harvested_chains,
    audit_order_only_chains,
    fit_exponent,
    run_arm,
)

from causal_spacetime_lab.density_coupled_clocks import (  # noqa: E402
    bracket_width_against_worldline,
    chain_is_causal,
    harvest_chain_from_sprinkling_1p1,
    harvest_order_only_chain_1p1,
    make_poisson_clock_chain_1p1,
    nearest_event_index,
)
from causal_spacetime_lab.observer import (  # noqa: E402
    make_stationary_observer_chain_1p1,
)
from causal_spacetime_lab.sprinkling import (  # noqa: E402
    sprinkle_1p1_causal_diamond,
)


def test_harvested_chain_is_a_causal_chain_inside_its_tube():
    """The audit invariants, exercised directly on one harvest: total
    causal order, tube and window containment, strictly increasing
    times, and bit-identical re-harvest (deterministic tie-breaking)."""

    bulk = sprinkle_1p1_causal_diamond(4000, T=2.0, seed=3)
    idx = harvest_chain_from_sprinkling_1p1(bulk, 0.2, 0.08, -0.6, 0.6)
    again = harvest_chain_from_sprinkling_1p1(bulk, 0.2, 0.08, -0.6, 0.6)
    chain = bulk[idx]

    assert idx.size >= 2
    assert np.array_equal(idx, again)
    assert chain_is_causal(chain)
    assert np.all(np.abs(chain[:, 1] - 0.2) <= 0.04)
    assert np.all((chain[:, 0] >= -0.6) & (chain[:, 0] <= 0.6))
    assert np.all(np.diff(chain[:, 0]) > 0)


def test_poisson_clock_is_reproducible_and_on_the_worldline():
    ticks_a = make_poisson_clock_chain_1p1(-0.5, 0.5, 80.0, 0.1, seed=5)
    ticks_b = make_poisson_clock_chain_1p1(-0.5, 0.5, 80.0, 0.1, seed=5)
    assert np.array_equal(ticks_a, ticks_b)
    assert np.all(ticks_a[:, 1] == 0.1)
    assert np.all(np.diff(ticks_a[:, 0]) > 0)
    assert ticks_a[0, 0] >= -0.5 and ticks_a[-1, 0] <= 0.5


def test_bracket_width_agrees_with_the_model_d_example():
    """On a uniform-grid worldline the general-worldline width must
    reproduce Lemma 2's null-aligned example exactly (W = 4 for target
    (0, 0.5) against ticks -0.75..0.75), and go NaN when unreachable."""

    chain_events, _ = make_stationary_observer_chain_1p1(1.5, 7, x=0.0)
    assert bracket_width_against_worldline(chain_events, 0.0, 0.5) == 4.0
    assert np.isnan(bracket_width_against_worldline(chain_events, 0.0, 5.0))


def test_constructor_audit_passes():
    assert audit_harvested_chains(seed=11, scenes=3)["passed"]


def test_thinned_protocol_obeys_the_proved_prediction():
    """The [PROVED] arm at reduced scale: rate exponent ~ 1, RMSE
    exponent ~ -1/2, and the RMSE tracks the exact sqrt(d / 2 lambda)
    prediction (loose band -- targets sharing a chain are correlated)."""

    rows = run_arm("thinned", rho_grid=(1000, 8000), seeds=6)
    ratios = [row["rmse"] / row["rmse_predicted"] for row in rows]
    assert 0.8 < fit_exponent(rows, "lam_mean") < 1.2
    assert -0.75 < fit_exponent(rows, "rmse") < -0.25
    assert 0.7 < float(np.mean(ratios)) < 1.3
    assert all(row["unreachable"] == 0 for row in rows)
    assert all(row["short_clocks"] == 0 for row in rows)


def test_harvested_rate_couples_at_the_discreteness_scale():
    """The [MEASURED] arm at reduced scale: the chain rate grows like
    sqrt(rho), NOT rho -- the fact that makes the rho^{-1/2} law
    protocol-dependent. Only sanity bands: the sharp error exponent is
    an open question the tests must not pretend to settle."""

    rows = run_arm("harvest_scaled", rho_grid=(1000, 8000), seeds=6)
    assert 0.3 < fit_exponent(rows, "lam_mean") < 0.7
    assert -0.45 < fit_exponent(rows, "rmse") < -0.05
    assert all(row["unreachable"] == 0 for row in rows)
    assert all(row["short_clocks"] == 0 for row in rows)


def test_order_only_chain_is_a_longest_anchor_to_anchor_chain():
    """Constructor invariants for the order-only harvest: total causal
    order, anchor endpoints, interval containment (in lightcone product
    order), strictly increasing times, and bit-identical re-harvest."""

    bulk = sprinkle_1p1_causal_diamond(4000, T=2.0, seed=7)
    bottom = nearest_event_index(bulk, -0.6, 0.3)
    top = nearest_event_index(bulk, 0.6, 0.3)
    idx = harvest_order_only_chain_1p1(bulk, bottom, top)
    again = harvest_order_only_chain_1p1(bulk, bottom, top)
    chain = bulk[idx]

    assert idx.size >= 4
    assert idx[0] == bottom and idx[-1] == top
    assert np.array_equal(idx, again)
    assert chain_is_causal(chain)
    assert np.all(np.diff(chain[:, 0]) > 0)
    u = bulk[:, 0] + bulk[:, 1]
    v = bulk[:, 0] - bulk[:, 1]
    interior = idx[1:-1]
    assert np.all((u[interior] > u[bottom]) & (u[interior] < u[top]))
    assert np.all((v[interior] > v[bottom]) & (v[interior] < v[top]))


def test_order_only_anchor_validation_rejects_non_causal_pairs():
    events = np.array([[0.0, 0.0], [0.1, 0.9], [1.0, 0.1]])
    with pytest.raises(ValueError):
        harvest_order_only_chain_1p1(events, 0, 0)  # identical anchors
    with pytest.raises(ValueError):
        harvest_order_only_chain_1p1(events, 0, 1)  # spacelike anchors
    idx = harvest_order_only_chain_1p1(events, 0, 2)  # timelike: fine
    assert list(idx) == [0, 2]


def test_order_only_harvest_is_invariant_under_order_preserving_maps():
    """The review's counterexample class, pinned: a spatial reflection
    (and a boost) preserves the labelled causal order and the anchors,
    so the harvested chain -- whose tie-break among maximum chains
    reads element labels, never coordinates -- must be identical
    index-by-index across these coordinate presentations."""

    for seed in (7, 21):
        bulk = sprinkle_1p1_causal_diamond(2000, T=2.0, seed=seed)
        bottom = nearest_event_index(bulk, -0.6, -0.3)
        top = nearest_event_index(bulk, 0.6, -0.3)
        idx = harvest_order_only_chain_1p1(bulk, bottom, top)

        reflected = bulk.copy()
        reflected[:, 1] = -reflected[:, 1]
        assert np.array_equal(
            idx, harvest_order_only_chain_1p1(reflected, bottom, top)
        )

        k = 1.7  # boost: (u, v) -> (k u, v / k) preserves the order
        u = bulk[:, 0] + bulk[:, 1]
        v = bulk[:, 0] - bulk[:, 1]
        boosted = np.column_stack(
            [(k * u + v / k) / 2.0, (k * u - v / k) / 2.0]
        )
        assert np.array_equal(
            idx, harvest_order_only_chain_1p1(boosted, bottom, top)
        )


def test_order_only_audit_passes_with_dp_cross_validation():
    """The audit includes an independent causal-matrix DP length check:
    patience sorting must return an order-theoretic LONGEST chain."""

    audit = audit_order_only_chains(seed=13, scenes=3)
    assert audit["dp_length_checks"] == 3
    assert audit["passed"], audit


def test_order_only_rate_couples_at_the_discreteness_scale():
    """Reduced-grid sanity for the order-only arm: discreteness-scale
    rate, no clock failures, and a recorded (positive) wandering RMS.
    The sharp error and wandering exponents belong to the full grid."""

    rows = run_arm("harvest_order_only", rho_grid=(1000, 8000), seeds=6)
    assert 0.3 < fit_exponent(rows, "lam_mean") < 0.7
    assert -0.50 < fit_exponent(rows, "rmse") < 0.05
    assert all(row["unreachable"] == 0 for row in rows)
    assert all(row["short_clocks"] == 0 for row in rows)
    assert all(row["transverse_rms"] > 0 for row in rows)
