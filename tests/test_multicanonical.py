"""Pin the multicanonical sampler to exact enumeration at small N.

At N <= 7 every permutation is enumerable, so the density of states g(S),
the canonical averages, and the free energy are all known analytically. The
Wang-Landau estimate is checked against those, not against another sampler --
an agreement between two Monte Carlo runs would only prove they share a bug.

The bar matters: this sampler exists to make a claim at N=600 that local
Metropolis cannot make. If it cannot reproduce N=6 exactly, nothing it says
at N=600 is worth reading.
"""

from __future__ import annotations

from itertools import permutations

import numpy as np
import pytest

from causal_spacetime_lab.positive_control.multicanonical import (
    _bin_index,
    action_range,
    effective_sample_size,
    multicanonical_2d_order,
    reweight_to_beta,
    wang_landau_2d_order,
)
from causal_spacetime_lab.positive_control.two_orders import (
    _IncrementalState,
    bipartite_perm,
    chain_observables,
    perm_to_causal_matrix,
)


def _exact_states(n: int, eps: float) -> list[tuple[float, dict[str, float]]]:
    """Every permutation of n elements with its action and observables."""

    states = []
    for perm in permutations(range(n)):
        pi = np.array(perm, dtype=int)
        action = _IncrementalState(pi, eps).action()
        observables = chain_observables(perm_to_causal_matrix(pi))
        observables["S"] = action
        states.append((action, observables))
    return states


def _exact_ln_g(
    n: int, eps: float, bin_edges: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Exact log density of states, binned onto the same grid as the estimate."""

    counts = np.zeros(bin_edges.size - 1, dtype=float)
    for action, _ in _exact_states(n, eps):
        index = _bin_index(action, bin_edges)
        assert index >= 0, "enumeration escaped the action window"
        counts[index] += 1.0
    occupied = counts > 0
    ln_g = np.zeros_like(counts)
    ln_g[occupied] = np.log(counts[occupied])
    ln_g[occupied] -= ln_g[occupied].min()
    return ln_g, occupied


def _exact_canonical_mean(n: int, eps: float, beta: float, key: str) -> float:
    """Exact <key> under the Gibbs measure exp(-beta S) over all permutations."""

    states = _exact_states(n, eps)
    actions = np.array([action for action, _ in states])
    values = np.array([observables[key] for _, observables in states])
    log_weights = -beta * actions
    log_weights -= log_weights.max()
    weights = np.exp(log_weights)
    return float((weights * values).sum() / weights.sum())


N = 6
EPS = 0.3


@pytest.fixture(scope="module")
def window() -> tuple[float, float]:
    """An action window guaranteed to contain every N=6 state."""

    actions = [action for action, _ in _exact_states(N, EPS)]
    span = max(actions) - min(actions)
    return min(actions) - 0.05 * span, max(actions) + 0.05 * span


@pytest.fixture(scope="module")
def wl(window: tuple[float, float]):
    s_min, s_max = window
    return wang_landau_2d_order(
        pi0=np.arange(N),
        eps=EPS,
        s_min=s_min,
        s_max=s_max,
        n_bins=24,
        seed=7,
        sweep_steps=4000,
        ln_f_final=1e-4,
        max_sweeps=4000,
    )


def test_wang_landau_converges_and_crosses_the_barrier(wl):
    assert wl.converged, f"ln_f stalled at {wl.final_ln_f}"
    # A flat-in-S walker that never traverses the window has not sampled the
    # barrier, which is the entire point of running this instead of Metropolis.
    assert wl.round_trips > 0
    assert 0.0 < wl.acceptance < 1.0


def test_wang_landau_escapes_the_flatness_stall_via_one_over_t():
    """Regression: a flatness-only schedule stalls on a rugged density of states.

    In the eps*N=12 regime the N=60 pilot managed a single ln_f halving in 24M
    moves -- the histogram over a many-decade g(S) simply never came flat, so
    ln_f froze at 0.5 and the run could never converge however long it ran.

    Forcing an unreachable flatness (1.0 demands a perfectly uniform histogram)
    reproduces that stall exactly. The 1/t schedule must take over and still
    drive ln_f to the target.
    """

    actions = [action for action, _ in _exact_states(N, EPS)]
    span = max(actions) - min(actions)

    result = wang_landau_2d_order(
        pi0=np.arange(N),
        eps=EPS,
        s_min=min(actions) - 0.05 * span,
        s_max=max(actions) + 0.05 * span,
        n_bins=24,
        seed=5,
        sweep_steps=2000,
        flatness=1.0,  # unsatisfiable: only the stage budget can drive ln_f down
        max_sweeps_per_stage=5,
        # The normalized Belardinelli-Pereyra tail decays as n_bins/updates, so
        # reaching ln_f_final needs >= n_bins/ln_f_final updates: 240k here,
        # comfortably inside the 2.4M-update budget.
        ln_f_final=1e-4,
        max_sweeps=1200,
    )

    # 1/t engages late by design -- it handles the fine-convergence tail, after
    # halving has already brought ln_f down to the 1/t curve.
    assert result.entered_one_over_t, "the 1/t schedule never engaged"
    assert result.converged, (
        f"1/t failed to drive ln_f to target: stalled at {result.final_ln_f:.2e}"
    )
    assert result.final_ln_f <= 1e-4


def test_wang_landau_ln_g_matches_exact_enumeration(wl, window):
    exact, occupied = _exact_ln_g(N, EPS, wl.bin_edges)

    # Only bins that actually hold states are meaningful; WL should have found
    # every one of them, and should not have invented mass in empty bins.
    assert np.array_equal(wl.visited, occupied), (
        "Wang-Landau visited a different bin set than exact enumeration"
    )

    # ln_g is defined up to an additive constant; both are shifted to min 0.
    error = np.abs(wl.ln_g[occupied] - exact[occupied])
    assert error.max() < 0.35, f"max |ln g error| = {error.max():.3f}"


@pytest.mark.parametrize("beta", [0.0, 0.4, 1.0, 2.0])
def test_reweighted_canonical_means_match_exact_gibbs(wl, beta):
    """The payoff: one flat-S run reproduces every beta, checked against exact."""

    production = multicanonical_2d_order(
        pi0=np.arange(N),
        eps=EPS,
        ln_g=wl.ln_g,
        bin_edges=wl.bin_edges,
        steps=400_000,
        seed=11,
        sample_every=20,
        burn_frac=0.2,
        visited=wl.visited,
    )
    assert production.round_trips > 0

    for key, tolerance in (("S", 0.02), ("n0", 0.35), ("height", 0.06)):
        estimate, _ = reweight_to_beta(production.samples, beta, key)
        exact = _exact_canonical_mean(N, EPS, beta, key)
        assert abs(estimate - exact) < tolerance, (
            f"beta={beta} {key}: reweighted {estimate:.4f} vs exact {exact:.4f}"
        )

    # Reweighting is only licensed where it retains statistical support.
    assert effective_sample_size(production.samples, beta) > 10.0


def test_ln_f_does_not_decay_before_the_walker_traverses_the_window():
    """Regression: ln_f reaching its target says nothing about ln_g being right.

    A stage budget that halves ln_f on a sweep count alone will drive ln_f onto
    the 1/t curve while the walker is still stuck in one basin -- and then the
    run reports converged=True on a ln_g that is simply wrong. That happened:
    an N=60 pilot "converged" with a multicanonical chain covering 4.7% of the
    window at 1.3% acceptance.

    Demanding an impossible number of round trips per stage must therefore
    freeze the schedule outright rather than let it decay on schedule.
    """

    actions = [action for action, _ in _exact_states(N, EPS)]
    span = max(actions) - min(actions)

    result = wang_landau_2d_order(
        pi0=np.arange(N),
        eps=EPS,
        s_min=min(actions) - 0.05 * span,
        s_max=max(actions) + 0.05 * span,
        n_bins=24,
        seed=5,
        sweep_steps=2000,
        flatness=1.0,
        max_sweeps_per_stage=1,
        min_round_trips_per_stage=10_000_000,  # unreachable
        ln_f_final=1e-6,
        max_sweeps=60,
    )

    assert not result.converged
    assert result.final_ln_f == pytest.approx(1.0), (
        "ln_f decayed despite the walker never meeting the round-trip bar"
    )
    assert not result.entered_one_over_t


def test_reweighting_ess_decays_away_from_flat_support(wl):
    """ESS must fall as beta grows -- silence here would hide extrapolation."""

    production = multicanonical_2d_order(
        pi0=np.arange(N),
        eps=EPS,
        ln_g=wl.ln_g,
        bin_edges=wl.bin_edges,
        steps=100_000,
        seed=13,
        sample_every=20,
        burn_frac=0.2,
        visited=wl.visited,
    )
    near = effective_sample_size(production.samples, 0.0)
    far = effective_sample_size(production.samples, 6.0)
    assert far < near


def test_production_rejects_bins_the_wang_landau_run_never_visited(wl):
    """Regression: an unvisited bin has NO ln_g estimate, not a low one.

    After the output shift its ln_g entry is 0.0 -- indistinguishable from the
    lowest-density visited bins -- so without the visited mask a production
    move into it is accepted freely and reweighted as if its density were
    known. The mask must make such bins behave like out-of-window states.
    """

    visited = wl.visited.copy()
    initial_bin = _bin_index(
        float(_IncrementalState(np.arange(N), EPS).action()), wl.bin_edges
    )
    # Poke a hole in some visited interior bin that is not the starting bin.
    candidates = [
        b for b in np.flatnonzero(visited) if b != initial_bin
    ]
    hole = candidates[len(candidates) // 2]
    visited[hole] = False

    production = multicanonical_2d_order(
        pi0=np.arange(N),
        eps=EPS,
        ln_g=wl.ln_g,
        bin_edges=wl.bin_edges,
        steps=150_000,
        seed=17,
        sample_every=20,
        burn_frac=0.1,
        visited=visited,
    )
    assert production.samples, "production produced no samples"
    for row in production.samples:
        assert _bin_index(row["S"], wl.bin_edges) != hole

    # And a start inside a hole must refuse outright rather than walk on a
    # weight that does not exist.
    bad_visited = wl.visited.copy()
    bad_visited[initial_bin] = False
    with pytest.raises(ValueError, match="never visited"):
        multicanonical_2d_order(
            pi0=np.arange(N),
            eps=EPS,
            ln_g=wl.ln_g,
            bin_edges=wl.bin_edges,
            steps=1000,
            seed=19,
            visited=bad_visited,
        )


def test_action_range_brackets_the_extreme_orders_not_just_random_ones():
    """Regression: a random-only bracket truncates the chain extreme.

    At N=20, eps=0.1 the identity permutation (the total chain) has a larger
    action than any of 40 random draws. A window fitted to random probes would
    exclude it, and Wang-Landau would reject every move into that region in
    silence. All three extreme orders must land inside the bracket.
    """

    n, eps = 20, 0.1
    low, high = action_range(n_elements=n, eps=eps, seed=3, probes=40)
    assert low < high
    for extreme in (
        np.arange(n),
        np.arange(n)[::-1],
        bipartite_perm(n),
    ):
        action = _IncrementalState(np.array(extreme, dtype=int), eps).action()
        assert low <= action <= high, f"extreme order at S={action:.3f} escaped"


def test_wang_landau_rejects_a_window_it_cannot_start_in():
    with pytest.raises(ValueError, match="outside"):
        wang_landau_2d_order(
            pi0=np.arange(N),
            eps=EPS,
            s_min=1e6,
            s_max=2e6,
            n_bins=10,
            seed=1,
        )
