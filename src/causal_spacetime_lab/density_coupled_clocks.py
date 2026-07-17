"""Density-coupled tick protocols (T1 Section 5 / gap G2 instrumentation).

Three constructors whose tick statistics are coupled to the sprinkling
density, unlike ``make_stationary_observer_chain_1p1`` (whose grid never
sees ``rho``):

- ``make_poisson_clock_chain_1p1``: a synthetic worldline at exact ``x0``
  whose tick times form a Poisson process of a caller-chosen rate. With
  ``lam = rho * ell`` this realizes Model P with the ``lam ~ rho * ell``
  coupling that the T1 document's ``rho^{-1/2}`` law presumes.
- ``harvest_chain_from_sprinkling_1p1``: a longest causal chain among the
  *sprinkled events themselves* inside a spatial tube around ``x0`` — the
  coordinate-tube harvested clock. Its ticks are order elements, but the
  tube SELECTION uses embedded coordinates (``|x - x0|`` and a time
  window), so this is deliberately NOT called an order-intrinsic clock:
  an order-only causal set could not reproduce the harvest without extra
  geometric data. The chain's rate and fluctuation statistics are still
  properties of the sprinkling, not choices: measured, a maximal chain's
  rate grows like ``sqrt(rho)`` (the discreteness scale), NOT like
  ``rho``, so the protocols scale differently and any density-scaling
  claim must name its protocol.
- ``harvest_order_only_chain_1p1``: a longest causal chain between two
  DESIGNATED anchor events — the answer to the order-only design
  question the tube protocol left open (T1 Section 6, G2). Given the
  order and the two anchor labels, the selection reads order data
  alone; the chain is free to wander transversally, and that wandering
  is part of the protocol's physics, not an error to be corrected.

All are theory-track instrumentation: no frozen gate consumes them.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.causal import causal_matrix_1p1

_ATOL = 1e-12  # the causal_matrix_minkowski default interval tolerance


def make_poisson_clock_chain_1p1(
    t_min: float,
    t_max: float,
    lam: float,
    x: float,
    seed: int | np.random.Generator | None = None,
) -> NDArray[np.float64]:
    """Return worldline events ``(t_k, x)`` with Poisson tick times.

    Tick times are a homogeneous Poisson process of rate ``lam`` on
    ``[t_min, t_max]``, sorted. The worldline is exact (all ticks at
    spatial coordinate ``x``); only the clock statistics are stochastic.
    """

    if t_max <= t_min:
        raise ValueError("t_max must exceed t_min")
    if lam <= 0:
        raise ValueError("lam must be positive")
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    count = rng.poisson(lam * (t_max - t_min))
    times = np.sort(rng.uniform(t_min, t_max, size=count))
    return np.column_stack([times, np.full(count, float(x))])


def harvest_chain_from_sprinkling_1p1(
    events: NDArray[np.float64],
    x0: float,
    tube_width: float,
    t_min: float,
    t_max: float,
) -> NDArray[np.int_]:
    """Return bulk indices of one longest causal chain in a spatial tube.

    Candidate ticks are the sprinkled events with ``|x - x0| <=
    tube_width / 2`` and ``t in [t_min, t_max]``; among them a longest
    chain of the (null-inclusive) causal order is extracted by dynamic
    programming. Ties are broken deterministically (first maximum in
    time-sorted order), so the harvest is a pure function of its inputs.
    The returned indices are time-sorted; the chain's elements are
    *actual sprinkled events*, so tick positions wiggle inside the tube —
    that wiggle is part of the protocol, not an error to be corrected.
    """

    if tube_width <= 0:
        raise ValueError("tube_width must be positive")
    events = np.asarray(events, dtype=float)
    in_tube = (
        (np.abs(events[:, 1] - x0) <= tube_width / 2.0)
        & (events[:, 0] >= t_min)
        & (events[:, 0] <= t_max)
    )
    candidates = np.flatnonzero(in_tube)
    if candidates.size == 0:
        return np.empty(0, dtype=int)

    order = candidates[np.argsort(events[candidates, 0], kind="stable")]
    tube_events = events[order]
    causal = causal_matrix_1p1(tube_events)

    n = order.size
    best = np.ones(n, dtype=int)
    parent = np.full(n, -1, dtype=int)
    for j in range(n):
        preds = np.flatnonzero(causal[:, j])
        if preds.size:
            k = int(preds[np.argmax(best[preds])])
            if best[k] + 1 > best[j]:
                best[j] = best[k] + 1
                parent[j] = k
    j = int(np.argmax(best))
    chain: list[int] = []
    while j >= 0:
        chain.append(j)
        j = int(parent[j])
    return order[np.array(chain[::-1], dtype=int)]


def nearest_event_index(
    events: NDArray[np.float64],
    t: float,
    x: float,
) -> int:
    """Index of the sprinkled event nearest (Euclidean in ``(t, x)``) to a
    point.

    A SETUP helper for designating harvest anchors: the designation is
    coordinate-assisted by intent (like placing observers), and is not
    part of any order-only selection rule.
    """

    events = np.asarray(events, dtype=float)
    if events.shape[0] == 0:
        raise ValueError("no events to select an anchor from")
    d2 = (events[:, 0] - t) ** 2 + (events[:, 1] - x) ** 2
    return int(np.argmin(d2))


def harvest_order_only_chain_1p1(
    events: NDArray[np.float64],
    bottom_index: int,
    top_index: int,
) -> NDArray[np.int_]:
    """Return bulk indices of one longest causal chain between two
    designated anchor events -- the order-only harvest of T1's G2.

    Given the causal order and TWO DESIGNATED ELEMENTS (the anchors),
    the selection rule reads order data alone: candidates are the
    elements causally between the anchors, and the harvest is a longest
    chain from the bottom anchor to the top anchor. An order-only
    causal set with the two anchors labelled can reproduce this harvest
    exactly -- the upgrade over ``harvest_chain_from_sprinkling_1p1``,
    whose tube membership cannot be expressed without embedded
    coordinates. What remains non-order data is the anchor DESIGNATION
    itself (made once by the caller, e.g. via ``nearest_event_index``)
    and nothing else.

    Implementation notes. Every ingredient of the selection is a
    function of the labelled order alone -- where "the order" is the
    project's tolerant null-inclusive predicate, exactly as
    ``causal_matrix_minkowski`` evaluates it (``dt > 0`` and ``dt^2 >=
    dx^2 - atol`` with ``atol = 1e-12``): interval membership (causal
    comparability with both anchors, so events exactly on a null
    boundary are candidates), the per-element longest-chain lengths
    ``B(x)`` = length of a longest chain from ``x`` to the top anchor,
    and, among the typically many maximum chains, a greedy
    **minimal-label** choice: the chain is built bottom-up, at each
    step taking the candidate with ``B`` one less than the current
    element's that is causally above the current element and has the
    smallest bulk index. Ties therefore break by element LABELS,
    never by coordinates: relabelling can change the representative,
    but coordinate presentations of the same labelled order (e.g. a
    spatial reflection) cannot -- asserted as a reflection-invariance
    audit and regression.

    ``B`` is computed on the fast path in O(k log k) by patience
    sorting under the exact lightcone product order. That order can
    disagree with the tolerant predicate only on a *sliver pair*
    (``dt > 0`` with ``du * dv`` in ``[-atol, 0)``), and any such pair
    has ``|du| <= 1e-6`` or ``|dv| <= 1e-6`` (else ``|du * dv| >
    atol``), so slivers are detected by sorted-coordinate window
    scans; when none exist -- almost surely, for sprinkled data -- the
    exact and tolerant relations coincide on the candidate set and the
    fast path is provably the tolerant optimum. When a sliver IS
    present (near-null constructed inputs), ``B`` falls back to a
    direct dynamic program under the tolerant predicate itself.
    Reconstruction always uses the tolerant predicate. Events are
    assumed pairwise distinct (a.s. for sprinklings).
    """

    events = np.asarray(events, dtype=float)
    bottom = int(bottom_index)
    top = int(top_index)
    if bottom == top:
        raise ValueError("anchors must be distinct events")
    dt = events[top, 0] - events[bottom, 0]
    dx = events[top, 1] - events[bottom, 1]
    if not (dt > 0 and dt * dt - dx * dx >= -_ATOL):
        raise ValueError("anchors must satisfy bottom -> top causally")

    # null-inclusive interval membership, the causal-matrix predicate
    dt_b = events[:, 0] - events[bottom, 0]
    iv_b = dt_b * dt_b - (events[:, 1] - events[bottom, 1]) ** 2
    dt_t = events[top, 0] - events[:, 0]
    iv_t = dt_t * dt_t - (events[top, 1] - events[:, 1]) ** 2
    inside = (dt_b > 0) & (iv_b >= -_ATOL) & (dt_t > 0) & (iv_t >= -_ATOL)
    inside[bottom] = False
    inside[top] = False
    candidates = np.flatnonzero(inside)
    if candidates.size == 0:
        return np.array([bottom, top], dtype=int)

    ct = events[candidates, 0]
    cx = events[candidates, 1]
    su = ct + cx
    sv = ct - cx
    # B[i] = longest-chain length among candidates STARTING at i: an
    # order invariant (ends-at lengths of the reversed order)
    if _has_tolerance_sliver_pair(ct, cx):
        b_len = _tolerant_lengths_starting_at(ct, cx)
    else:
        b_len = _chain_lengths_ending_at(-su, -sv)
    m = int(b_len.max())

    # greedy minimal-label reconstruction of one maximum chain under
    # the tolerant predicate: at each step the feasible pool is
    # nonempty by the exchange argument (any successor of the current
    # element on a maximum chain has B exactly one less), and the pick
    # reads labels only
    chain_local: list[int] = []
    cur = -1
    for need in range(m, 0, -1):
        pool = np.flatnonzero(b_len == need)
        if cur >= 0:
            dt_p = ct[pool] - ct[cur]
            iv_p = dt_p * dt_p - (cx[pool] - cx[cur]) ** 2
            pool = pool[(dt_p > 0) & (iv_p >= -_ATOL)]
        pick = int(pool[np.argmin(candidates[pool])])
        chain_local.append(pick)
        cur = pick
    interior = candidates[np.array(chain_local, dtype=int)]
    return np.concatenate(([bottom], interior, [top])).astype(int)


def _has_tolerance_sliver_pair(
    t: NDArray[np.float64], x: NDArray[np.float64]
) -> bool:
    """Detect pairs where the exact lightcone product order disagrees
    with the tolerant causal predicate (``dt > 0``, ``du * dv`` in
    ``[-atol, 0)``). Any such pair has ``|du| <= 1e-6`` or ``|dv| <=
    1e-6``, so scanning sorted-coordinate windows of that width is
    exhaustive; window pairs are then tested exactly."""

    window = 1e-6
    u = t + x
    v = t - x
    for w, other in ((u, v), (v, u)):
        order = np.argsort(w, kind="stable")
        ws = w[order]
        os_ = other[order]
        start = 0
        for j in range(1, ws.size):
            while ws[j] - ws[start] > window:
                start += 1
            for i in range(start, j):
                du = ws[j] - ws[i]
                dv = os_[j] - os_[i]
                if du + dv < 0:  # orient so dt > 0
                    du, dv = -du, -dv
                if du + dv > 0 and du * dv >= -_ATOL and (du < 0 or dv < 0):
                    return True
    return False


def _tolerant_lengths_starting_at(
    t: NDArray[np.float64], x: NDArray[np.float64]
) -> NDArray[np.int_]:
    """Longest-path lengths to the future under the tolerant causal
    predicate, by direct dynamic programming in reverse time order.
    O(k^2); used only when a tolerance-sliver pair exists, which is
    measure-zero for sprinkled data."""

    order = np.argsort(-t, kind="stable")
    lengths = np.ones(t.size, dtype=int)
    for i in order:
        dt = t - t[i]
        iv = dt * dt - (x - x[i]) ** 2
        succ = (dt > 0) & (iv >= -_ATOL)
        if succ.any():
            lengths[i] = 1 + int(lengths[succ].max())
    return lengths


def _chain_lengths_ending_at(
    a: NDArray[np.float64], b: NDArray[np.float64]
) -> NDArray[np.int_]:
    """Length of the longest product-order chain ending at each element
    (``x <= y`` iff ``a_x <= a_y`` and ``b_x <= b_y``; distinct points
    with both coordinates non-decreasing are null-inclusively causal),
    by patience sorting in O(k log k). The returned lengths are order
    invariants of the input point set; the sort is internal machinery
    only."""

    order = np.lexsort((b, a))  # a ascending, then b ascending
    lengths = np.zeros(a.size, dtype=int)
    tails: list[float] = []
    for i in order:
        value = float(b[i])
        lo, hi = 0, len(tails)
        while lo < hi:  # bisect_right: first tail > value (chains may tie)
            mid = (lo + hi) // 2
            if tails[mid] <= value:
                lo = mid + 1
            else:
                hi = mid
        if lo == len(tails):
            tails.append(value)
        else:
            tails[lo] = value
        lengths[i] = lo + 1
    return lengths


def chain_is_causal(chain_events: NDArray[np.float64]) -> bool:
    """Audit: every consecutive pair is causally related (hence, by
    transitivity of the Minkowski causal order, the set is a chain)."""

    chain_events = np.asarray(chain_events, dtype=float)
    if chain_events.shape[0] < 2:
        return True
    dt = np.diff(chain_events[:, 0])
    dx = np.diff(chain_events[:, 1])
    return bool(np.all((dt > 0) & (dt * dt - dx * dx >= -1e-12)))


def bracket_width_against_worldline(
    tick_events: NDArray[np.float64],
    t: float,
    x: float,
) -> float:
    """Rank-unit bracket width of a target against arbitrary tick events.

    Unlike the uniform-grid path, harvested ticks are not at a single
    spatial coordinate, so the causal relation is evaluated per tick.
    Returns NaN when the target is unreachable (no predecessor or no
    successor among the ticks). Rank is the index in time-sorted order —
    the only clock, as everywhere in T1.
    """

    tick_events = np.asarray(tick_events, dtype=float)
    if tick_events.shape[0] == 0:
        return float("nan")
    dt = t - tick_events[:, 0]
    dx = x - tick_events[:, 1]
    interval = dt * dt - dx * dx
    precedes = (dt > 0) & (interval >= -1e-12)
    succeeds = (dt < 0) & (interval >= -1e-12)
    if not precedes.any() or not succeeds.any():
        return float("nan")
    return float(np.min(np.flatnonzero(succeeds)) - np.max(np.flatnonzero(precedes)))
