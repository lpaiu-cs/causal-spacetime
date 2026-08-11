"""Certified diamond-volume integrator (PR-O2) -- the L6 contract of
docs/theory/schwarzschild_volume_oracle_certification.md.

    V = 2 pi INT INT r^2 f(r) sin psi [dt - T1 - T2]_+ drho dpsi

over the L4 containment box, as a finite directed-rounding interval
sum of per-cell enclosures. Two certified cell modes:

- first-order (Lipschitz): center flight-time calls + the L1
  coordinate bounds |dS/drho| <= 2, |dS/dpsi| <= 2 w. Valid at ANY
  distance from the anchors (no curvature constant enters), so
  anchor neighborhoods need no ball excision -- strictly tighter
  than L6b's [0, U_ball] terms, and the fallback whenever the
  tangent model's preconditions fail.
- tangent model (L6a): center value + eikonal gradient enclosures
  (b hulls from the flight-time contract, radial signs from the arc
  geometry) + the L3 remainder with PATH-WISE curvature constants
  (kappa coth(kappa d) <= 1/d + kappa^2 d / 3, so no coth is ever
  evaluated). The affine-with-interval-coefficients hinge is
  integrated by an interval SUBCELL COVERING (L6a option ii) -- no
  point sampling anywhere.

The angle-cost constant of the pruning/distance functionals is
3 sqrt(3) M for M > 0 (the global optical-radius minimum, L4) and
2 R_MIN / pi for M = 0 -- the flat chord obeys
chord >= 2 sqrt(r1 r2) sin(psi/2) >= (2 R_MIN / pi) psi on the
shell, whereas a naive "w >= R_MIN" bound would be UNSOUND (flat
space has no photon-sphere floor).

Pruning (L6c) uses only closed-form lower bounds, so it can never
discard a cell that meets the diamond. Refinement intersects each
recomputed global interval with the previous one (L6d); an empty
intersection hard-stops the run. Failure to reach the frozen target
(V_hi - V_lo)/(V_hi + V_lo) <= 0.01 under the cost caps returns the
certified interval with status `target-not-met` -- never a silent
acceptance. Monte Carlo lives in `mc_diagnostic`, outside every
certified path, and is never a gate.
"""

from __future__ import annotations

import heapq
import time
from dataclasses import dataclass, field

from certified_flight_time import (
    M_DEFAULT,
    R_MAX,
    R_MIN,
    CertifiedFlightTime,
    angle_cost_iv,
    containment_certificate,
    flight_time_certified,
    tortoise_iv,
    w_iv,
)
from certified_interval import (
    CertificationError,
    Iv,
    iv_max,
    iv_min,
    iv_sum,
)

_TWO = Iv(2)
_RECOMPUTE_EVERY = 32   # splits between certified total recomputes
_MODES = ("pruned", "tangent", "first-order", "uncosted")


@dataclass(frozen=True)
class OracleConfig:
    """One anchor configuration plus the frozen target and the cost
    caps. `n_sub` is the flight-time quadrature; `k_micro` the
    subcell-covering resolution of the tangent-model hinge."""

    r_in: float
    r_out: float
    dt: float
    m: float = M_DEFAULT
    target_ratio: float = 0.01
    n_sub: int = 128
    k_micro: int = 4
    d_switch: float = 0.25      # optical distance below which cells
    #                             use first-order mode (no 1/d)
    max_calls: int = 200_000
    max_wall_s: float = 3600.0
    max_depth: int = 12
    init_rho: int = 12          # initial grid resolution
    init_psi: int = 12


@dataclass
class _Cell:
    rho0: float
    rho1: float
    psi0: float
    psi1: float
    depth: int
    contrib: Iv | None = None
    mode: str = "unevaluated"
    dead: bool = False


@dataclass
class _State:
    cfg: OracleConfig
    m: Iv
    dt: Iv
    rho_p: Iv
    rho_q: Iv
    wg: Iv
    rho_smin: Iv = None
    w_shell: Iv = None
    calls: int = 0
    t0: float = 0.0
    r_cache: dict[float, Iv] = field(default_factory=dict)


def _r_of_rho(st: _State, rho: float) -> Iv:
    """Certified r with tortoise(r) = rho, cached per grid point
    (strictly increasing map; bisection with certified signs)."""

    hit = st.r_cache.get(rho)
    if hit is not None:
        return hit
    target = Iv(rho)
    lo, hi = R_MIN, R_MAX
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        t_mid = tortoise_iv(Iv(mid), st.m)
        if t_mid.certainly_lt(target):
            lo = mid
        elif t_mid.certainly_gt(target):
            hi = mid
        else:
            break
    out = Iv(lo, hi)
    st.r_cache[rho] = out
    return out


def _cell_r(st: _State, c: _Cell) -> Iv:
    return _r_of_rho(st, c.rho0).hull(_r_of_rho(st, c.rho1))


def _weight(st: _State, rho_box: Iv, psi_box: Iv, r_box: Iv) -> Iv:
    """r^2 f sin psi = r (r - 2M) sin psi: monotone in r and in psi
    on the cap, so its exact interval sits at the corners (L6a
    step 5, with the (rho, psi) Jacobian dr = f drho)."""

    lo_r, hi_r = Iv(r_box.lo), Iv(r_box.hi)
    w_lo = lo_r * (lo_r - _TWO * st.m)
    w_hi = hi_r * (hi_r - _TWO * st.m)
    s = psi_box.sin()
    return Iv((w_lo * Iv(s.lo)).lo, (w_hi * Iv(s.hi)).hi)


def _hinge(x: Iv) -> Iv:
    zero = Iv(0)
    lo = x.lo if x.lo > 0 else zero.lo
    hi = x.hi if x.hi > 0 else zero.hi
    return Iv(lo, hi)


def _angle_lower(st: _State, c: _Cell, rho_a: Iv) -> Iv:
    """Certified angular-cost lower bound from anchor to cell: any
    path either stays in the shell -- length >= w(R_MIN) dpsi (w is
    monotone above the photon sphere) -- or reaches r < R_MIN and
    pays the tortoise exit cost from both endpoints. The minimum of
    the two branches is a valid bound and beats the global 3 sqrt(3)
    constant by ~2x on the shell; wg remains the floor."""

    in_shell = st.w_shell * Iv(c.psi0)
    exit_cost = ((rho_a - st.rho_smin)
                 + (Iv(c.rho0) - st.rho_smin))
    branch = iv_min(in_shell, iv_max(exit_cost, Iv(0)))
    return iv_max(branch, st.wg * Iv(c.psi0))


def _closed_form_lower(st: _State, c: _Cell) -> Iv:
    """Certified closed-form lower bound of T1 + T2 over the cell
    from the L4 functionals only (no solver calls)."""

    rho_box = Iv(c.rho0, c.rho1)
    t1_lo = iv_max(Iv((rho_box - st.rho_p).abs().lo),
                   _angle_lower(st, c, st.rho_p))
    t2_lo = iv_max(Iv((rho_box - st.rho_q).abs().lo),
                   _angle_lower(st, c, st.rho_q))
    return t1_lo + t2_lo


def _anchor_distance_lower(st: _State, c: _Cell, rho_a: Iv) -> Iv:
    rho_box = Iv(c.rho0, c.rho1)
    d_rho = (rho_box - rho_a).abs()
    return iv_max(Iv(d_rho.lo), _angle_lower(st, c, rho_a))


def _center_calls(st: _State, c: _Cell,
                  ) -> tuple[CertifiedFlightTime, CertifiedFlightTime,
                             Iv, float]:
    """Two certified flight-time calls at a cell-center point. The
    expansion center's rho coordinate is returned as the CERTIFIED
    interval tortoise(r_pt) -- the model must expand around the point
    actually evaluated, not around the float grid midpoint."""

    rho_c = 0.5 * (c.rho0 + c.rho1)
    psi_c = 0.5 * (c.psi0 + c.psi1)
    r_pt = _r_of_rho(st, rho_c).mid()
    ft1 = flight_time_certified(st.cfg.r_in, r_pt, psi_c, st.cfg.m,
                                n_sub=st.cfg.n_sub)
    ft2 = flight_time_certified(st.cfg.r_out, r_pt, psi_c, st.cfg.m,
                                n_sub=st.cfg.n_sub)
    st.calls += 2
    rho_ctr = tortoise_iv(Iv(r_pt), st.m)
    return ft1, ft2, rho_ctr, psi_c


def _first_order(st: _State, c: _Cell, r_cell: Iv) -> Iv:
    """L1 Lipschitz enclosure: valid everywhere, including against
    the anchors (no curvature constant enters)."""

    ft1, ft2, rho_ctr, psi_c = _center_calls(st, c)
    s_c = ft1.t + ft2.t
    hw_rho = Iv((Iv(c.rho0, c.rho1) - rho_ctr).abs().hi)
    hw_psi = Iv(max(psi_c - c.psi0, c.psi1 - psi_c))
    w_hi = Iv(w_iv(Iv(r_cell.hi), st.m).hi)
    spread = _TWO * hw_rho + _TWO * w_hi * hw_psi
    s_cell = s_c.widen(spread)
    hinge = _hinge(st.dt - s_cell)
    area = (Iv(c.rho1) - Iv(c.rho0)) * (Iv(c.psi1) - Iv(c.psi0))
    weight = _weight(st, Iv(c.rho0, c.rho1), Iv(c.psi0, c.psi1),
                     r_cell)
    return area * weight * hinge


def _w_prime(r: Iv, m: Iv) -> Iv:
    """w'(rho) = (1 - 3M/r)/sqrt(f) (L1), interval over the cell."""

    f = Iv(1) - _TWO * m / r
    return (Iv(1) - Iv(3) * m / r) / f.sqrt()


def _grad_rho(ft: CertifiedFlightTime, r_anchor: float,
              r_ctr: Iv, w_ctr: Iv) -> Iv:
    """Certified dT/drho enclosure AT the evaluated center point via
    the eikonal (L6a step 2): magnitude sqrt(1 - (b/w)^2), sign from
    the arc geometry, hulled when not certifiable."""

    if ft.b is None:
        return Iv(-1.0, 1.0)
    ratio = (ft.b / w_ctr).sq().intersect(Iv(0.0, 1.0))
    mag = (Iv(1) - ratio).sqrt0()
    if ft.family == "one-turn":
        return mag  # the non-anchor endpoint sits on the outgoing side
    if ft.family == "no-turn":
        if r_ctr.lo > r_anchor:
            return mag
        if r_ctr.hi < r_anchor:
            return -mag
    return Iv((-mag).lo, mag.hi)


def _r_linear(st: _State, c: _Cell, r_cell: Iv,
              rho0: float, rho1: float) -> Iv:
    """Certified r-range on a rho subinterval of the cell via the
    linear enclosure r(rho) in r(c.rho0) + [f_lo, f_hi] (rho -
    c.rho0) (dr/drho = f, monotone bounds from the cell's r-range) --
    avoids a tortoise inversion per microcell."""

    f_box = Iv(1) - _TWO * st.m / r_cell
    base = st.r_cache[c.rho0]
    off = Iv(rho0, rho1) - Iv(c.rho0)
    return (base + f_box * off).intersect(r_cell)


def _tangent_model(st: _State, c: _Cell, r_cell: Iv,
                   d1: Iv, d2: Iv) -> Iv:
    """L6a tangent-plane enclosure with path-wise curvature
    constants and an interval subcell covering of the hinge."""

    ft1, ft2, rho_ctr, psi_c = _center_calls(st, c)
    r_ctr = _r_of_rho(st, 0.5 * (c.rho0 + c.rho1))
    w_ctr = w_iv(r_ctr, st.m)
    g_rho = (_grad_rho(ft1, st.cfg.r_in, r_ctr, w_ctr)
             + _grad_rho(ft2, st.cfg.r_out, r_ctr, w_ctr))
    g_psi = ft1.b + ft2.b

    # L3c remainder with PATH-WISE kappa (L3a): the geodesic from an
    # anchor may undercut the cell's radial range down to the L2a
    # perihelion floor of the pair; kappa coth(kappa d) <= 1/d +
    # kappa^2 d / 3 evaluated at the certified d_min
    cos_half = (Iv(c.psi1) / _TWO).cos()
    lam = Iv(0)
    for r_anchor, d_min in ((st.cfg.r_in, d1), (st.cfg.r_out, d2)):
        r_geo = Iv(min(r_anchor, float(r_cell.lo))) * cos_half
        k2 = _hinge(_TWO * st.m / (r_geo * r_geo * r_geo)
                    - Iv(3) * st.m.sq() / (r_geo * r_geo).sq())
        d = Iv(d_min.lo)
        lam = lam + Iv(1) / d + k2 * d / Iv(3)
    w_hi = Iv(w_iv(Iv(r_cell.hi), st.m).hi)
    wp_hi = Iv(_w_prime(r_cell, st.m).abs().hi)
    m_rp = lam * w_hi + _TWO * wp_hi
    m_pp = lam * w_hi.sq() + _TWO * w_hi * wp_hi
    hw_r = Iv((Iv(c.rho0, c.rho1) - rho_ctr).abs().hi)
    hw_p = Iv(max(psi_c - c.psi0, c.psi1 - psi_c))
    e_c = (lam * hw_r.sq() + _TWO * m_rp * hw_r * hw_p
           + m_pp * hw_p.sq()) / _TWO

    s_c = ft1.t + ft2.t
    k = st.cfg.k_micro
    terms = []
    dr = (c.rho1 - c.rho0) / k
    dp = (c.psi1 - c.psi0) / k
    for i in range(k):
        r0, r1 = c.rho0 + i * dr, c.rho0 + (i + 1) * dr
        r_sub = _r_linear(st, c, r_cell, r0, r1)
        for j in range(k):
            p0, p1 = c.psi0 + j * dp, c.psi0 + (j + 1) * dp
            s_model = (s_c + g_rho * (Iv(r0, r1) - rho_ctr)
                       + g_psi * (Iv(p0, p1) - Iv(psi_c))
                       ).widen(e_c)
            hinge = _hinge(st.dt - s_model)
            area = (Iv(r1) - Iv(r0)) * (Iv(p1) - Iv(p0))
            weight = _weight(st, Iv(r0, r1), Iv(p0, p1), r_sub)
            terms.append(area * weight * hinge)
    return iv_sum(terms)


def _trivial_bound(st: _State, c: _Cell, r_cell: Iv) -> Iv:
    """The certified enclosure that costs NO solver calls: T1, T2 are
    nonnegative, so the hinge lies in [0, dt] and the cell integral
    in area * weight * [0, dt]. Sound, merely wide -- what a cell
    gets when the cost caps are already spent."""

    area = (Iv(c.rho1) - Iv(c.rho0)) * (Iv(c.psi1) - Iv(c.psi0))
    weight = _weight(st, Iv(c.rho0, c.rho1), Iv(c.psi0, c.psi1),
                     r_cell)
    return area * weight * Iv(0).hull(st.dt)


def _over_budget(st: _State) -> bool:
    if st.calls >= st.cfg.max_calls:
        return True
    # t0 == 0.0 means "not started under a clock" (direct use
    # outside assemble); only assemble stamps it
    return (st.t0 > 0.0
            and time.perf_counter() - st.t0 >= st.cfg.max_wall_s)


def _evaluate(st: _State, c: _Cell) -> None:
    lower = _closed_form_lower(st, c)
    if lower.certainly_gt(st.dt):
        c.contrib, c.mode = Iv(0), "pruned"
        return
    r_cell = _cell_r(st, c)
    if _over_budget(st):
        # the caps bind HERE, not only in the refinement loop: the
        # initial grid is solver work too, and a small cap must not
        # be blown through before the first target check (review R1)
        c.contrib, c.mode = _trivial_bound(st, c, r_cell), "uncosted"
        return
    d1 = _anchor_distance_lower(st, c, st.rho_p)
    d2 = _anchor_distance_lower(st, c, st.rho_q)
    if (st.cfg.m > 0.0 and d1.lo > st.cfg.d_switch
            and d2.lo > st.cfg.d_switch):
        try:
            c.contrib = _tangent_model(st, c, r_cell, d1, d2)
            c.mode = "tangent"
            return
        except CertificationError:
            pass  # first-order is always available
    c.contrib, c.mode = _first_order(st, c, r_cell), "first-order"


def _total(st: _State, cells: list[_Cell],
           prev: Iv | None) -> Iv:
    tot = iv_sum(c.contrib for c in cells
                 if not c.dead) * _TWO * Iv.pi()
    # L6d: refinement must nest -- an empty intersection means an
    # enclosure was wrong and intersect() hard-fails
    return tot if prev is None else tot.intersect(prev)


def assemble(cfg: OracleConfig, targets: list[float] | None = None,
             progress=None) -> dict:
    """Run the certified assembly. Returns the certified interval,
    the status against the target(s), and the cost counters; the
    interval is certified REGARDLESS of the status. An empty diamond
    (dt certainly below T(P, Q)) returns the exact [0, 0].

    `targets` runs the whole price ladder in ONE adaptive pass:
    the cost is recorded the moment the running interval first
    crosses each ratio, instead of restarting the refinement per
    rung (which throws away every previous rung's work). `progress`,
    if given, is called with each curve sample as it is produced --
    the caller streams it, so a long run is observable while it
    runs rather than only at the end."""

    m_i, dt_i = Iv(cfg.m), Iv(cfg.dt)
    rho_p = tortoise_iv(Iv(cfg.r_in), m_i)
    rho_q = tortoise_iv(Iv(cfg.r_out), m_i)
    if dt_i.certainly_lt(rho_q - rho_p):
        return {"v": Iv(0), "status": "empty-diamond", "ratio": None,
                "target_ratio": cfg.target_ratio, "calls": 0,
                "cells": 0, "modes": {}, "wall_s": 0.0,
                "certificate": None, "crossings": {}, "curve": []}
    cert = containment_certificate(cfg.r_in, cfg.r_out, cfg.dt,
                                   cfg.m)
    st = _State(cfg, m_i, dt_i, rho_p, rho_q, angle_cost_iv(m_i))
    st.rho_smin = tortoise_iv(Iv(R_MIN), m_i)
    st.w_shell = (w_iv(Iv(R_MIN), m_i) if cfg.m > 0.0
                  else angle_cost_iv(m_i))
    st.t0 = time.perf_counter()
    rho_lo = ((rho_p + rho_q - dt_i) / _TWO).lo_float()
    rho_hi = ((rho_p + rho_q + dt_i) / _TWO).hi_float()
    psi_hi = cert["psi_max"].hi_float()

    cells: list[_Cell] = []
    drho = (rho_hi - rho_lo) / cfg.init_rho
    dpsi = psi_hi / cfg.init_psi
    for i in range(cfg.init_rho):
        for j in range(cfg.init_psi):
            cells.append(_Cell(rho_lo + i * drho,
                               rho_lo + (i + 1) * drho,
                               j * dpsi, (j + 1) * dpsi, 0))
    heap: list[tuple[float, int, _Cell]] = []
    tie = 0
    for c in cells:
        _evaluate(st, c)
        heapq.heappush(heap, (-c.contrib.width(), tie, c))
        tie += 1

    pending = sorted(targets if targets else [cfg.target_ratio],
                     reverse=True)
    crossings: dict[float, dict] = {}
    curve: list[dict] = []
    state: dict = {"prev": None, "done": False}

    def checkpoint() -> Iv:
        """Recompute the certified total, record the curve sample,
        and credit every target it crosses. Runs on the periodic
        cadence AND once more after the loop exits: a crossing in
        the last few splits before a cap would otherwise be lost to
        `target-not-met` and misfiled as an extrapolation, which is
        exactly what the price ladder must never do (review R1)."""

        total = _total(st, cells, state["prev"])
        state["prev"] = total
        denom = float(total.lo + total.hi)
        if not (total.certainly_gt(Iv(0)) and denom > 0):
            return total
        sample = {"calls": st.calls,
                  "cells": sum(1 for c in cells if not c.dead),
                  "ratio": total.width() / denom,
                  "v_lo": float(total.lo), "v_hi": float(total.hi),
                  "wall_s": time.perf_counter() - st.t0}
        if not curve or curve[-1]["calls"] != sample["calls"]:
            curve.append(sample)
            if progress is not None:
                progress(sample)
        while pending and sample["ratio"] <= pending[0]:
            crossings[pending.pop(0)] = dict(sample)
        if not pending:
            state["done"] = True
        return total

    since_recompute = _RECOMPUTE_EVERY  # force an initial check
    while True:
        if since_recompute >= _RECOMPUTE_EVERY:
            since_recompute = 0
            checkpoint()
            if state["done"]:
                break
        if _over_budget(st) or not heap:
            break
        _, _, worst = heapq.heappop(heap)
        if (worst.dead or worst.depth >= cfg.max_depth
                or worst.mode == "pruned"):
            continue
        worst.dead = True
        w_r = worst.rho1 - worst.rho0
        w_p = (worst.psi1 - worst.psi0) * float(
            w_iv(Iv(_cell_r(st, worst).hi), st.m).hi)
        if w_r >= w_p:
            mid = 0.5 * (worst.rho0 + worst.rho1)
            kids = [_Cell(worst.rho0, mid, worst.psi0, worst.psi1,
                          worst.depth + 1),
                    _Cell(mid, worst.rho1, worst.psi0, worst.psi1,
                          worst.depth + 1)]
        else:
            mid = 0.5 * (worst.psi0 + worst.psi1)
            kids = [_Cell(worst.rho0, worst.rho1, worst.psi0, mid,
                          worst.depth + 1),
                    _Cell(worst.rho0, worst.rho1, mid, worst.psi1,
                          worst.depth + 1)]
        for kid in kids:
            _evaluate(st, kid)
            cells.append(kid)
            heapq.heappush(heap, (-kid.contrib.width(), tie, kid))
            tie += 1
        since_recompute += 1

    total = checkpoint()
    status = "target-met" if state["done"] else "target-not-met"
    denom = float(total.lo + total.hi)
    ratio = (total.width() / denom) if denom > 0 else None
    alive = [c for c in cells if not c.dead]
    return {
        "v": total,
        "status": status,
        "ratio": ratio,
        "target_ratio": cfg.target_ratio,
        "calls": st.calls,
        "cells": len(alive),
        "modes": {mode: sum(1 for c in alive if c.mode == mode)
                  for mode in _MODES},
        # which cell MODE carries the remaining width, not just how
        # many cells each has: the two answers differ sharply when
        # the anchor neighbourhoods (first-order, O(h) per cell) sit
        # among many well-behaved tangent cells (O(h^2)), and it is
        # the width split that says where to spend the next lever
        "width_by_mode": {
            mode: sum(c.contrib.width() for c in alive
                      if c.mode == mode)
            for mode in _MODES},
        "wall_s": time.perf_counter() - st.t0,
        "certificate": cert,
        "crossings": crossings,
        "curve": curve,
    }


def quadrature_floor_diagnostic(cfg: OracleConfig, v_estimate: float,
                                n_subs=(16, 32, 64), grid: int = 48,
                                probes: int = 12) -> list[dict]:
    """DIAGNOSTIC (not certified): the part of the assembled width
    that the flight-time enclosures themselves force, and how it
    scales with the quadrature.

    Every cell strictly inside the diamond carries its center
    flight-time's uncertainty no matter how small the cell gets, so
    cell refinement cannot push the assembled width below

        W_floor(n_sub) ~ width(T1 + T2) * 2 pi INT INT r^2 sin psi
                          dr dpsi   over the diamond's support.

    Only a finer quadrature shrinks it. The support integral is
    estimated on a fixed grid with the FAST non-certified S1 solver;
    width(T1 + T2) is measured with the CERTIFIED solver at interior
    probe points. Both are estimates: this names which lever binds,
    it never enters a certified interval.

    Computed standalone rather than accumulated inside the
    integrator: a per-cell tally that ignores the hinge counts cells
    OUTSIDE the diamond (whose true contribution width is zero) and
    overstates the floor -- the error this replaces."""

    import math
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                           / "positive_control"))
    from s1_schwarzschild_cost import flight_time as s1_ft

    cert = containment_certificate(cfg.r_in, cfg.r_out, cfg.dt,
                                   cfg.m)
    r_lo = cert["r_box"][0].lo_float()
    r_hi = cert["r_box"][1].hi_float()
    psi_hi = cert["psi_max"].hi_float()
    dr = (r_hi - r_lo) / grid
    dp = psi_hi / grid
    support = 0.0
    interior: list[tuple[float, float]] = []
    for i in range(grid):
        r = r_lo + (i + 0.5) * dr
        for j in range(grid):
            psi = (j + 0.5) * dp
            t1, _ = s1_ft(cfg.r_in, r, psi, cfg.m)
            t2, _ = s1_ft(cfg.r_out, r, psi, cfg.m)
            if cfg.dt - t1 - t2 > 0.0:
                support += r * r * math.sin(psi) * dr * dp
                interior.append((r, psi))
    omega = 2.0 * math.pi * support
    if not interior:
        return []
    step = max(1, len(interior) // probes)
    picked = interior[::step][:probes]

    rows = []
    for n_sub in n_subs:
        widths = []
        for r, psi in picked:
            f1 = flight_time_certified(cfg.r_in, r, psi, cfg.m,
                                       n_sub=n_sub)
            f2 = flight_time_certified(cfg.r_out, r, psi, cfg.m,
                                       n_sub=n_sub)
            widths.append(f1.t.width() + f2.t.width())
        mean_w = sum(widths) / len(widths)
        floor_w = mean_w * omega
        rows.append({
            "n_sub": n_sub,
            "mean_flight_time_width": mean_w,
            "max_flight_time_width": max(widths),
            "support_measure": omega,
            "floor_width": floor_w,
            "floor_ratio": floor_w / (2.0 * v_estimate),
            "probes": len(picked),
        })
    return rows


def mc_diagnostic(cfg: OracleConfig, n_samples: int,
                  seed: int) -> dict:
    """DIAGNOSTIC ONLY (L6d): a plain Monte Carlo estimate of V with
    a CLT interval, using the NON-certified S1 heuristic solver for
    speed. Never a gate: a disjoint comparison raises an
    investigation flag in the caller, nothing more. The seed is a
    spent diagnostic stream in the probe ledger."""

    import math
    import sys
    from pathlib import Path

    import numpy as np

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                           / "positive_control"))
    from s1_schwarzschild_cost import flight_time as s1_ft

    cert = containment_certificate(cfg.r_in, cfg.r_out, cfg.dt,
                                   cfg.m)
    r_lo = cert["r_box"][0].lo_float()
    r_hi = cert["r_box"][1].hi_float()
    psi_hi = cert["psi_max"].hi_float()
    rng = np.random.default_rng(seed)
    rs = rng.uniform(r_lo, r_hi, n_samples)
    psis = rng.uniform(0.0, psi_hi, n_samples)
    vals = np.empty(n_samples)
    for i, (r, psi) in enumerate(zip(rs, psis, strict=True)):
        t1, _ = s1_ft(cfg.r_in, float(r), float(psi), cfg.m)
        t2, _ = s1_ft(cfg.r_out, float(r), float(psi), cfg.m)
        vals[i] = (r * r * math.sin(psi)
                   * max(cfg.dt - t1 - t2, 0.0))
    box = (r_hi - r_lo) * psi_hi
    mean = 2.0 * math.pi * box * float(vals.mean())
    se = 2.0 * math.pi * box * float(vals.std(ddof=1)
                                     / math.sqrt(n_samples))
    return {"estimate": mean, "se": se,
            "ci95": (mean - 1.959964 * se, mean + 1.959964 * se),
            "n": n_samples, "seed": seed}
