"""Wang-Landau / multicanonical sampling over 2D orders.

Motivation (P7). At N=600 the local Metropolis chain of :mod:`.two_orders`
is not ergodic on any feasible budget for beta >= 16: random and bipartite
starts stay in different basins with start-separated mean-action gaps of
77--96, so no equilibrium beta_c can be extracted. A naive replica ladder
does not fix this -- matching a target exchange acceptance of 0.2 across an
action gap of ~77 needs a local beta spacing of order ``-log(0.2)/77 ~ 0.02``,
i.e. hundreds of replicas, and still guarantees nothing about barrier
crossing.

This module attacks the barrier directly. Instead of sampling
``exp(-beta S)`` at fixed beta, it estimates the density of states ``g(S)``
of the smeared action and then samples the *multicanonical* weight
``W(S) = 1/g(S)``, which is flat in S by construction. A flat-in-S walker
must traverse the action barrier to cover its range, so basin-to-basin
round trips become a measurable, reportable quantity rather than a hope.

Two consequences matter for P7:

1. One ln g(S) estimate reweights to *every* beta -- the beta grid stops
   being a set of independent, separately-unconverged experiments.
2. Convergence is falsifiable. Histogram flatness and the round-trip count
   are direct overlap evidence between the continuum and crystal basins;
   zero round trips means the run failed, and says so.

Correctness is pinned to exact enumeration at small N (see
``tests/test_multicanonical.py``): at N <= 7 every permutation is
enumerable, so the estimated ln g(S) and every reweighted canonical
observable are checked against the analytic answer, not against another
sampler.

Nothing here is preregistered. This is method development; the N=600
production protocol and any beta_c hypothesis must be frozen separately,
after this sampler is shown to converge at production scale.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .two_orders import _IncrementalState, chain_observables


@dataclass
class WangLandauResult:
    """Estimated log density of states and the diagnostics that qualify it."""

    ln_g: np.ndarray
    """Log density of states per bin, shifted so ``min(ln_g[visited]) == 0``."""

    bin_edges: np.ndarray
    visited: np.ndarray
    """True for bins the walker actually entered; ln_g is meaningless elsewhere."""

    final_ln_f: float
    """Modification factor reached. The run converged only if this hit ln_f_final."""

    sweeps: int
    round_trips: int
    """Completed low-S <-> high-S traversals. Zero means the barrier never fell."""

    acceptance: float
    converged: bool

    entered_one_over_t: bool = False
    """True once the run left the flatness regime for the 1/t schedule."""

    moves: int = 0

    def bin_centers(self) -> np.ndarray:
        return 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])


@dataclass
class MulticanonicalResult:
    """Flat-in-S production samples, carrying the weights needed to reweight."""

    samples: list[dict[str, float]] = field(default_factory=list)
    acceptance: float = 0.0
    round_trips: int = 0
    ln_g: np.ndarray = field(default_factory=lambda: np.zeros(0))
    bin_edges: np.ndarray = field(default_factory=lambda: np.zeros(0))


def _bin_index(action: float, bin_edges: np.ndarray) -> int:
    """Index of the bin holding ``action``, or -1 if outside the window."""

    if action < bin_edges[0] or action > bin_edges[-1]:
        return -1
    index = int(np.searchsorted(bin_edges, action, side="right") - 1)
    return min(max(index, 0), bin_edges.size - 2)


def action_range(
    n_elements: int, eps: float, seed: int = 0, probes: int = 200
) -> tuple[float, float]:
    """Heuristic bracket of the action window, from the known extreme orders.

    Probes four families: random permutations (the continuum bulk), the
    complete bipartite order (the crystal extreme), and both the identity and
    reversed permutations (the total chain and antichain -- the relation-count
    extremes). Random draws alone are *not* enough: at N=20, eps=0.1 the
    identity sits above the largest of 40 random actions, so a window fitted to
    random probes would silently truncate reachable states, and every move into
    the truncated region would be rejected without any error.

    This is a bracket over probed configurations, not a proven global bound.
    Padding is 10% of the observed span on each side. If a Wang-Landau run
    piles visits against either edge bin, the window is still too narrow and
    must be widened by hand -- flatness alone will not reveal it.
    """

    from .two_orders import bipartite_perm

    rng = np.random.default_rng(seed)
    actions = [
        _IncrementalState(rng.permutation(n_elements), eps).action()
        for _ in range(probes)
    ]
    for extreme in (
        bipartite_perm(n_elements),
        np.arange(n_elements),
        np.arange(n_elements)[::-1],
    ):
        actions.append(_IncrementalState(np.array(extreme, dtype=int), eps).action())

    low, high = min(actions), max(actions)
    pad = 0.1 * (high - low) if high > low else 1.0
    return low - pad, high + pad


def wang_landau_2d_order(
    pi0: np.ndarray,
    eps: float,
    s_min: float,
    s_max: float,
    n_bins: int,
    seed: int,
    sweep_steps: int = 10_000,
    flatness: float = 0.8,
    ln_f_init: float = 1.0,
    ln_f_final: float = 1e-6,
    max_sweeps: int = 100_000,
    max_sweeps_per_stage: int = 50,
    min_round_trips_per_stage: int = 1,
) -> WangLandauResult:
    """Estimate ln g(S) by Wang-Landau random walk over the action.

    The walker proposes the same random transpositions as the canonical
    sampler, but accepts with ``min(1, g(S_old)/g(S_new))``, which drives it
    toward *under*-visited actions instead of low ones. Moves leaving
    ``[s_min, s_max]`` are rejected, so the window must actually contain both
    basins -- use :func:`action_range` to set it.

    Each visit does ``ln_g[bin] += ln_f`` and ``H[bin] += 1``. ln_f is halved
    when the visit histogram over reached bins comes flat to ``flatness`` -- or,
    failing that, after ``max_sweeps_per_stage`` sweeps. The budget matters: on
    this ensemble g(S) spans many decades, flatness almost never passes, and a
    flatness-only schedule freezes ln_f (the N=60 pilot managed one halving in
    24M moves). Once ln_f has fallen to ``1/t``, the run switches permanently to
    the Belardinelli-Pereyra ``ln_f = 1/t`` schedule, which has no flatness test
    to stall on and still lets ln_g accumulate, since the harmonic sum diverges.

    The run converges when ln_f falls below ``ln_f_final``; if ``max_sweeps``
    is hit first, ``converged`` is False and the result must not be used for
    production reweighting.
    """

    if s_max <= s_min:
        raise ValueError("s_max must exceed s_min")
    if n_bins < 2:
        raise ValueError("n_bins must be at least 2")

    rng = np.random.default_rng(seed)
    bin_edges = np.linspace(s_min, s_max, n_bins + 1)
    ln_g = np.zeros(n_bins, dtype=float)
    histogram = np.zeros(n_bins, dtype=np.int64)
    visited = np.zeros(n_bins, dtype=bool)

    state = _IncrementalState(pi0, eps)
    action = state.action()
    current = _bin_index(action, bin_edges)
    if current < 0:
        raise ValueError(
            f"initial action {action:.4f} lies outside [{s_min:.4f}, {s_max:.4f}]"
        )

    # Round trips are counted between the outer fifths of the window: the
    # walker must reach the crystal end and come back for one to register.
    low_zone = n_bins // 5
    high_zone = n_bins - 1 - n_bins // 5
    last_end: int | None = None
    round_trips = 0

    ln_f = ln_f_init
    proposed_moves = 0
    accepted = 0
    sweeps = 0
    stage_sweeps = 0
    round_trips_at_stage_start = 0
    updates = 0
    one_over_t = False

    while ln_f > ln_f_final and sweeps < max_sweeps:
        for _ in range(sweep_steps):
            i, j = rng.integers(0, state.pi.size, 2)
            if i == j:
                continue
            proposed_moves += 1
            state.swap(i, j)
            proposed_action = state.action()
            candidate = _bin_index(proposed_action, bin_edges)
            if candidate < 0:
                state.swap(i, j)  # outside the window: reject, stay put
            elif np.log(rng.uniform()) < ln_g[current] - ln_g[candidate]:
                action = proposed_action
                current = candidate
                accepted += 1
            else:
                state.swap(i, j)

            updates += 1
            ln_g[current] += ln_f
            histogram[current] += 1
            visited[current] = True

            if one_over_t:
                # Belardinelli-Pereyra: once the flatness schedule has driven
                # ln_f below 1/t, keep it pinned there. The 1/t decay is slow
                # enough to preserve detailed balance asymptotically and fast
                # enough to terminate, and -- the reason it is here -- it does
                # not depend on a flatness test that a rugged density of states
                # may never satisfy.
                ln_f = 1.0 / updates

            if current <= low_zone:
                if last_end == 1:
                    round_trips += 1
                last_end = 0
            elif current >= high_zone:
                if last_end == 0:
                    round_trips += 1
                last_end = 1

        sweeps += 1
        if one_over_t:
            continue

        stage_sweeps += 1
        reached = histogram[visited]
        flat = bool(reached.size and reached.min() >= flatness * reached.mean())

        # Halve on flatness, or on an exhausted stage budget -- but never before
        # the walker has traversed the window at least once in the current stage.
        #
        # Both halves of that rule were forced by an observed failure:
        #
        # - Without the budget, ln_f freezes. g(S) here spans many decades, the
        #   flatness test almost never passes, and the N=60 pilot managed ONE
        #   halving in 24M moves. The textbook 1/t switch cannot rescue that
        #   either: a frozen ln_f compared against a shrinking 1/t never fires.
        #
        # - Without the round-trip requirement, the budget overshoots the other
        #   way. Halving every 50 sweeps regardless of exploration drove ln_f
        #   onto the 1/t curve while ln_g was still wrong, and the run then
        #   reported converged=True on a ln_g whose multicanonical chain covered
        #   4.7% of the window at 1.3% acceptance. ln_f reaching its target says
        #   nothing about ln_g being right.
        #
        # Tying the schedule to completed round trips makes the decay rate
        # answer to exploration rather than to the wall clock.
        explored = round_trips - round_trips_at_stage_start
        if (flat or stage_sweeps >= max_sweeps_per_stage) and (
            explored >= min_round_trips_per_stage
        ):
            ln_f *= 0.5
            histogram[:] = 0
            stage_sweeps = 0
            round_trips_at_stage_start = round_trips

        # Hand over to 1/t once it decays at least as fast as halving would.
        # Beyond this point there is no flatness test left to stall on, and the
        # harmonic sum still diverges, so ln_g keeps accumulating corrections.
        if ln_f <= 1.0 / updates:
            one_over_t = True

    if visited.any():
        ln_g[visited] -= ln_g[visited].min()
    ln_g[~visited] = 0.0

    return WangLandauResult(
        ln_g=ln_g,
        bin_edges=bin_edges,
        visited=visited,
        final_ln_f=ln_f,
        sweeps=sweeps,
        round_trips=round_trips,
        acceptance=accepted / proposed_moves if proposed_moves else 0.0,
        converged=ln_f <= ln_f_final,
        entered_one_over_t=one_over_t,
        moves=updates,
    )


def multicanonical_2d_order(
    pi0: np.ndarray,
    eps: float,
    ln_g: np.ndarray,
    bin_edges: np.ndarray,
    steps: int,
    seed: int,
    sample_every: int = 1000,
    burn_frac: float = 0.5,
) -> MulticanonicalResult:
    """Production run at fixed multicanonical weights ``W(S) = exp(-ln_g(S))``.

    Unlike the Wang-Landau stage, the weights no longer move, so the chain has
    a well-defined stationary distribution (flat in S over the window) and its
    samples are usable for :func:`reweight_to_beta`. Every sample carries the
    ln_g value at its own action, which is exactly the log-weight that has to
    be divided out when reweighting.
    """

    rng = np.random.default_rng(seed)
    state = _IncrementalState(pi0, eps)
    action = state.action()
    current = _bin_index(action, bin_edges)
    if current < 0:
        raise ValueError("initial action lies outside the multicanonical window")

    n_bins = bin_edges.size - 1
    low_zone = n_bins // 5
    high_zone = n_bins - 1 - n_bins // 5
    last_end: int | None = None
    round_trips = 0

    burn = int(steps * burn_frac)
    samples: list[dict[str, float]] = []
    proposed_moves = 0
    accepted = 0

    for step in range(steps):
        i, j = rng.integers(0, state.pi.size, 2)
        if i == j:
            continue
        proposed_moves += 1
        state.swap(i, j)
        proposed_action = state.action()
        candidate = _bin_index(proposed_action, bin_edges)
        if candidate < 0:
            state.swap(i, j)
        elif np.log(rng.uniform()) < ln_g[current] - ln_g[candidate]:
            action = proposed_action
            current = candidate
            accepted += 1
        else:
            state.swap(i, j)

        if current <= low_zone:
            if last_end == 1:
                round_trips += 1
            last_end = 0
        elif current >= high_zone:
            if last_end == 0:
                round_trips += 1
            last_end = 1

        if step >= burn and step % sample_every == 0:
            observables = chain_observables(state.C)
            observables["S"] = action
            observables["ln_weight"] = float(ln_g[current])
            samples.append(observables)

    return MulticanonicalResult(
        samples=samples,
        acceptance=accepted / proposed_moves if proposed_moves else 0.0,
        round_trips=round_trips,
        ln_g=ln_g,
        bin_edges=bin_edges,
    )


def reweight_to_beta(
    samples: list[dict[str, float]], beta: float, key: str
) -> tuple[float, float]:
    """Canonical mean and variance of ``key`` at ``beta`` from flat-S samples.

    A multicanonical sample at action S carries weight ``1/g(S)``, so the
    canonical average needs each sample re-weighted by
    ``g(S) * exp(-beta S)``, i.e. log-weight ``ln_g(S) - beta * S``. Weights
    are shifted by their maximum before exponentiating, which is a no-op on
    the normalized average but keeps exp() from overflowing at large beta.

    Returns (mean, variance). The variance is the reweighted second central
    moment, not a standard error -- it is what a susceptibility needs.
    """

    if not samples:
        raise ValueError("no samples to reweight")

    log_weights = np.array(
        [row["ln_weight"] - beta * row["S"] for row in samples], dtype=float
    )
    log_weights -= log_weights.max()
    weights = np.exp(log_weights)
    total = weights.sum()
    values = np.array([row[key] for row in samples], dtype=float)
    mean = float((weights * values).sum() / total)
    variance = float((weights * (values - mean) ** 2).sum() / total)
    return mean, variance


def effective_sample_size(samples: list[dict[str, float]], beta: float) -> float:
    """Kish ESS of the reweighting at ``beta``: ``(sum w)^2 / sum w^2``.

    This is the honest cost of reweighting far from where the flat-S run had
    support. If it collapses toward 1 at some beta, the ln g(S) estimate does
    not license a claim at that beta, however flat the histogram looked.
    """

    if not samples:
        raise ValueError("no samples to reweight")
    log_weights = np.array(
        [row["ln_weight"] - beta * row["S"] for row in samples], dtype=float
    )
    log_weights -= log_weights.max()
    weights = np.exp(log_weights)
    return float(weights.sum() ** 2 / (weights**2).sum())
