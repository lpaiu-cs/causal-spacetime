"""P14 §8 P1: the measurements the guard still owes (design §4.6.2/.3).

EXPLORATORY. No gate, no threshold, no verdict, no confirmation seed
window appears in this file, and none may be added to it (design §8.2).
This module MEASURES; deciding what the numbers mean is a review's job.

Four debts, in the order §8 P1 lists them:

1. **Eligible fraction.** §4.6.3 settled where the eligible set is
   non-empty; this measures how much of the box it is -- per element
   and per pair -- against the analytic volume fraction the guard's
   own operating point implies.
2. **The feasibility product.** The next design decision needs
   (eligible volume) x (elements a usable box holds) x (effect size)
   in one table. The effect at the fattest eligible axis diamond is
   `R(a) - 1` with `a = w du`, by quadrature (`axis_volume_ratio`).
   Two sizings are reported and they answer DIFFERENT questions
   (R1.1): `n90_detect` sizes rejecting a no-shift null with the raw
   same-points difference `N_A - N_0` (`Var = rho V_dis`, exact,
   since shared points cancel where the predicates agree) -- can the
   effect be seen at all; and `sd_Z` is the per-sprinkling standard
   deviation of P2's preregistered residual `Z = N_A - r N_0`
   (`Var(Z) = rho (V_A + r^2 V_0 - 2 r V_int)`, §8 P2), under which
   `E[Z] = 0` -- so it sizes the PRECISION with which n sprinklings
   verify the prediction, not a detection. Neither is presented as
   the other.
3. **The order-invariant candidate (§4.6.2).** Eligibility by
   coordinates costs Class C its pure-order standing. The candidate
   proxy -- an element is interior iff it has at least `k` elements
   below AND above it -- is measured against the coordinate guard:
   agreement, and each side's exclusive admissions. Candidate-only
   elements are the DANGEROUS direction (their diamonds may leave the
   box); guard-only elements are lost sensitivity.
4. **Ambiguity and its price.** `ambiguous_fraction` per sprinkling
   as §5.1 defines it -- a pair is ambiguous when EITHER arm is
   undecided, so both arms are censused and the union taken (R1.2) --
   plus escalation counts per arm and a micro-benchmark of what one
   escalated decision costs against a generic one.

Run:  python experiments/positive_control/p14_probe_p1.py [seed]
"""

from __future__ import annotations

import math
import sys
import time
from dataclasses import dataclass

import numpy as np
from p14_plane_wave import (
    PlaneWaveGeometry,
    Slab,
    arms,
    causal_relation,
    class_c_eligible,
    conjugate_du,
    guard_insets,
    sprinkle,
    transverse_cost,
)

#: Two-sided 95% plus 90% power -- the house convention (§6.1 sizes
#: designs at 90%, never at break-even). Break-even keeps `z_beta = 0`.
_Z_ALPHA = 1.959964
_Z_BETA = 1.281552


# --------------------------------------------------------------------
# 1. The eligible region, read off the guard
# --------------------------------------------------------------------

def element_eligible(geometry: PlaneWaveGeometry, point: np.ndarray) -> bool:
    """Class C eligibility of ONE element.

    The pair rule factorizes -- `class_c_eligible` tests each endpoint
    against the same inset sub-box -- so a single element's label is
    the pair predicate applied to `(p, p)`. Calling the predicate
    itself, rather than reimplementing the box test, means this cannot
    drift from the pipeline (the R12.3 lesson); the factorization is
    pinned by a test, not assumed here.
    """

    return class_c_eligible(geometry, point, point)


def eligible_volume_fraction(geometry: PlaneWaveGeometry) -> float:
    """Analytic eligible coordinate volume over box volume.

    Zero when the guard blocks. This is the number the guard's own
    operating point maximizes (§4.6.1), so measuring the element
    fraction against it is a consistency check on the whole chain.
    """

    slab = geometry.slab
    insets = guard_insets(geometry)
    x_lim = slab.dx / 2.0 - insets.x
    y_lim = slab.dy / 2.0 - insets.y
    v_window = slab.dv - 2.0 * insets.v
    if min(x_lim, y_lim, v_window) <= 0.0:
        return 0.0
    return (2.0 * x_lim * 2.0 * y_lim * v_window) / (
        slab.dx * slab.dy * slab.dv)


# --------------------------------------------------------------------
# 2. Relations: the full matrix, its cost, and its ambiguity
# --------------------------------------------------------------------

@dataclass(frozen=True)
class RelationCensus:
    """One sprinkling's relation matrix and what it cost.

    `related[i, j]` is True iff element `i` precedes element `j`, over
    elements SORTED BY `u`, so only `i < j` is ever populated.
    `ambiguous` counts pairs undecided even after escalation -- §5.1's
    undecided set -- and `escalated` counts pairs that needed exact
    arithmetic at all.
    """

    related: np.ndarray
    ambiguous: int
    escalated: int
    seconds: float
    #: `(i, j)` indices (u-sorted) of the undecided pairs, so a paired
    #: census can take the UNION across arms -- §5.1 calls a pair
    #: ambiguous when EITHER arm is undecided (R1.2).
    ambiguous_pairs: tuple = ()

    @property
    def pairs(self) -> int:
        n = self.related.shape[0]
        return n * (n - 1) // 2

    @property
    def ambiguous_fraction(self) -> float:
        return self.ambiguous / self.pairs if self.pairs else 0.0


def relation_census(geometry: PlaneWaveGeometry,
                    points: np.ndarray) -> RelationCensus:
    """Every unordered pair, decided once, with the clock running."""

    order = np.argsort(points[:, 0], kind="stable")
    pts = points[order]
    n = len(pts)
    related = np.zeros((n, n), dtype=bool)
    undecided = []
    escalated = 0
    start = time.perf_counter()
    for i in range(n):
        for j in range(i + 1, n):
            rel = causal_relation(geometry, pts[i], pts[j])
            if rel.escalated:
                escalated += 1
            if rel.related is None:
                undecided.append((i, j))
            elif rel.related:
                related[i, j] = True
    return RelationCensus(related=related, ambiguous=len(undecided),
                          escalated=escalated,
                          seconds=time.perf_counter() - start,
                          ambiguous_pairs=tuple(undecided))


def escalation_cost_microbench(geometry: PlaneWaveGeometry,
                               reps: int = 400) -> tuple[float, float]:
    """(generic_us, escalated_us) per `causal_relation` call.

    The escalated pair is CONSTRUCTED, not sampled: `Dv` is set to the
    double-precision cost itself, so the margin interval contains zero
    and the Decimal path must run. Uniform sampling cannot find the
    cone (measure zero), which is why the cost of escalation has to be
    measured on a pair built to sit there.
    """

    slab, w = geometry.slab, geometry.w
    du = slab.du * 0.8
    p = np.array([0.0, 0.6 * slab.dv, 0.05, 0.1])
    cost = transverse_cost(du, 0.05, -0.1, 0.1, 0.15, w)
    q_on = np.array([du, float(p[1]) - cost, -0.1, 0.15])
    q_off = np.array([du, float(p[1]) - cost - 0.25 * slab.dv, -0.1, 0.15])

    start = time.perf_counter()
    for _ in range(reps):
        causal_relation(geometry, p, q_off)
    generic = (time.perf_counter() - start) / reps

    start = time.perf_counter()
    for _ in range(reps):
        rel = causal_relation(geometry, p, q_on)
    escalated = (time.perf_counter() - start) / reps
    assert rel.escalated, "the constructed pair no longer escalates"
    return generic * 1e6, escalated * 1e6


# --------------------------------------------------------------------
# 3. The order-invariant candidate (§4.6.2)
# --------------------------------------------------------------------

def order_interiority(related: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(below, above) element counts, from the relation matrix alone.

    Pure order: nothing here reads a coordinate. `below[i]` is the
    number of elements preceding `i`, `above[i]` the number it
    precedes.
    """

    return related.sum(axis=0), related.sum(axis=1)


def candidate_mask(below: np.ndarray, above: np.ndarray,
                   k: int) -> np.ndarray:
    """Interior iff at least `k` elements below AND above (§4.6.2).

    `k = 0` admits everything, and the mask shrinks monotonically in
    `k`. Deliberately NOT a function of interval cardinality --
    selecting on that is forbidden, since truncation biases the very
    quantity Class C counts.
    """

    return (below >= k) & (above >= k)


def _choose2(n: int) -> int:
    return n * (n - 1) // 2


@dataclass(frozen=True)
class Agreement:
    """Confusion between candidate and coordinate guard.

    Element-level counts, with PAIR-level admissions derived from them
    -- §8 P1 asks for the pairs each rule admits that the other
    rejects, and element rates do not equal pair rates: both rules
    factorize over endpoints, so admissions compound. A pair is
    candidate-admitted iff both endpoints are candidate-eligible
    (categories `both` or `candidate_only`), guard-admitted iff both
    are guard-eligible (`both` or `guard_only`), admitted by BOTH
    rules iff both endpoints are in `both` (R1.3).
    """

    both: int
    candidate_only: int
    guard_only: int
    neither: int

    @property
    def total(self) -> int:
        return self.both + self.candidate_only + self.guard_only + self.neither

    @property
    def rate(self) -> float:
        return (self.both + self.neither) / self.total if self.total else 1.0

    @property
    def pair_total(self) -> int:
        return _choose2(self.total)

    @property
    def pair_both(self) -> int:
        return _choose2(self.both)

    @property
    def pair_candidate_only(self) -> int:
        """Pairs the candidate admits that the guard rejects -- the
        dangerous direction: their diamonds carry no containment
        certificate."""

        return _choose2(self.both + self.candidate_only) - self.pair_both

    @property
    def pair_guard_only(self) -> int:
        """Pairs the guard admits that the candidate rejects -- lost
        sensitivity, not lost soundness."""

        return _choose2(self.both + self.guard_only) - self.pair_both

    @property
    def pair_rate(self) -> float:
        if not self.pair_total:
            return 1.0
        disagree = self.pair_candidate_only + self.pair_guard_only
        return 1.0 - disagree / self.pair_total


def agreement(candidate: np.ndarray, guard: np.ndarray) -> Agreement:
    return Agreement(
        both=int(np.sum(candidate & guard)),
        candidate_only=int(np.sum(candidate & ~guard)),
        guard_only=int(np.sum(~candidate & guard)),
        neither=int(np.sum(~candidate & ~guard)))


# --------------------------------------------------------------------
# 4. The feasibility product
# --------------------------------------------------------------------

def axis_volume_ratio(a: float, nodes: int = 400) -> float:
    """`V_A / V_0` for an on-axis pair, exactly, as a function of `a`.

    Slicing the diamond at `u' = t s`: the two legs' costs sum to a
    quadratic form `A_x x^2 + A_y y^2` with

        A_x = (w/2) sinh(a) / (sinh(at) sinh(a(1-t)))
        A_y = (w/2) sin(a)  / (sin(at)  sin(a(1-t)))

    (the same `sin(a)/(sin sin) > 0` identity that carries the guard,
    R11), the transverse integral of the thickness is
    `pi Dv^2 / (2 sqrt(A_x A_y))`, and `Dv` FACTORS OUT -- the ratio
    depends on `a = w du` alone. Gauss-Legendre in `t`, like the
    design check script, whose pinned values this reproduces exactly:
    `1.00400047` at `a = 1`, `1.07300802` at `a = 2`
    (`docs/prereg/p14_checks/p14_interval_volume_constant_a.py`, whose
    `T` is the u-separation).

    **The first version of this probe used `(w tau)^4 / 252` with
    `tau^2 = 2 du dv` instead.** That imports `dv` into an effect `dv`
    cancels out of, overstating delta by 2.4e4x at the large-`dv`
    operating point and understating it 14x at `a = 1` -- and it was
    caught by the probe's own consistency bound, `delta <= V_dis/V_A`,
    not by reading. The series itself is also not usable at the large
    `a` the guard admits: at `a = 2.4` the true ratio is `1.1815`
    against the series' `1.1317`, and it diverges toward the conjugate
    point, so the probe uses the quadrature and never the series.
    """

    t, weight = np.polynomial.legendre.leggauss(nodes)
    t = 0.5 * (t + 1.0)
    weight = 0.5 * weight
    z1, z2 = a * t, a * (1.0 - t)
    curved = np.sqrt(np.sinh(z1) * np.sinh(z2) * np.sin(z1) * np.sin(z2))
    integral = float(np.sum(weight * curved))
    return integral / math.sqrt(math.sinh(a) * math.sin(a)) / (a / 6.0)


def fattest_axis_diamond(
        geometry: PlaneWaveGeometry) -> tuple[np.ndarray, np.ndarray]:
    """The largest eligible on-axis pair: full `u` extent, full window.

    On the transverse axis the cost is exactly zero, so the pair is
    causal whenever the window is open. Its volume effect is
    `axis_volume_ratio(w du) - 1`, a function of the u-extent alone.
    """

    slab = geometry.slab
    insets = guard_insets(geometry)
    p = np.array([0.0, slab.dv - insets.v, 0.0, 0.0])
    q = np.array([slab.du, insets.v, 0.0, 0.0])
    if not class_c_eligible(geometry, p, q):
        raise ValueError("the axis pair is not eligible; the operating "
                         "point admits nothing and has no feasibility row")
    assert causal_relation(geometry, p, q).related is True
    return p, q


def _between(geometry: PlaneWaveGeometry, p: np.ndarray, q: np.ndarray,
             r_u: float, r_v: float, r_x: float, r_y: float) -> bool:
    """Is `r` in the diamond -- BOTH memberships, not their sum.

    The first version tested `c1 + c2 <= v_p - v_q`, which never reads
    `r_v`: that is the condition for the `(u, x, y)` COLUMN to meet the
    diamond somewhere, so it measured the diamond's shadow times the
    full `v` depth of the box -- 3x the true volume in the test case
    that caught it, against the flat diamond's closed answer
    `pi s^2 Dv^2 / 6`. Membership is two inequalities, one per
    endpoint, and the sampled `v` must satisfy both.
    """

    c1 = transverse_cost(r_u - float(p[0]), float(p[2]), r_x,
                         float(p[3]), r_y, geometry.w)
    if float(p[1]) - r_v < c1:
        return False
    c2 = transverse_cost(float(q[0]) - r_u, r_x, float(q[2]),
                         r_y, float(q[3]), geometry.w)
    return r_v - float(q[1]) >= c2


def diamond_volumes_mc(
        curved: PlaneWaveGeometry, flat: PlaneWaveGeometry,
        p: np.ndarray, q: np.ndarray, samples: int,
        rng: np.random.Generator) -> tuple[float, float, float, float]:
    """(V_curved, V_flat, V_disagree, V_intersect), by shared-sample MC.

    The SAME sample decides all three, so `V_disagree` -- the volume
    where the two arms' membership predicates differ -- is measured
    directly rather than assembled from two noisy estimates. That is
    the variance that governs the same-points paired count: within one
    realization `Var(N_A - N_0) = rho * V_disagree`, because shared
    points cancel everywhere the predicates agree.

    The sampling region is the slab restricted to `u` between the
    endpoints: for an ELIGIBLE pair the guard certifies the diamond
    stays inside the box, so the region contains both diamonds and the
    estimate is unbiased.
    """

    slab = curved.slab
    u = rng.uniform(float(p[0]), float(q[0]), samples)
    v = rng.uniform(0.0, slab.dv, samples)
    x = rng.uniform(-slab.dx / 2.0, slab.dx / 2.0, samples)
    y = rng.uniform(-slab.dy / 2.0, slab.dy / 2.0, samples)
    region = (float(q[0]) - float(p[0])) * slab.dv * slab.dx * slab.dy

    in_curved = in_flat = disagree = overlap = 0
    for i in range(samples):
        a = _between(curved, p, q, u[i], v[i], x[i], y[i])
        b = _between(flat, p, q, u[i], v[i], x[i], y[i])
        in_curved += a
        in_flat += b
        disagree += a != b
        overlap += a and b
    return (region * in_curved / samples,
            region * in_flat / samples,
            region * disagree / samples,
            region * overlap / samples)


def sprinklings_needed(delta_abs: float, sd_per_sprinkling: float,
                       z_beta: float = _Z_BETA) -> float:
    """Sprinklings to see an absolute shift at 95% two-sided + power.

    `sd_per_sprinkling` and `delta_abs` in the same units; the
    sprinkling is the replication unit (§5.2), never the pair.
    """

    if delta_abs <= 0.0:
        return math.inf
    return ((_Z_ALPHA + z_beta) * sd_per_sprinkling / delta_abs) ** 2


# --------------------------------------------------------------------
# The probe itself
# --------------------------------------------------------------------

#: (label, w, du, dv, dx, dy). All must admit pairs -- the probe fails
#: loudly on one that does not, rather than skipping it. Chosen to
#: span the admitting census: the §4.6.3 slice, its anisotropic
#: neighbour, the roomiest cell, and the small boxes at large `a`.
OPERATING_POINTS = (
    ("slice-a0.3", 1.0, 0.3, 0.2, 6.0, 6.0),
    ("slice-a0.6", 1.0, 0.6, 0.2, 6.0, 6.0),
    ("slice-a1.0", 1.0, 1.0, 0.2, 6.0, 6.0),
    ("aniso-a1.0", 1.0, 1.0, 1.0, 2.0, 6.0),
    ("roomy-a0.2", 1.0, 0.2, 16.0, 6.0, 6.0),
    ("high-a2.0", 1.0, 2.0, 0.5, 2.0, 2.0),
    ("edge-a2.4", 1.0, 2.4, 0.2, 1.2, 0.8),
)

#: Elements per sprinkling, fixed ACROSS operating points so the table
#: compares instruments at equal compute (the predicate is O(N^2)),
#: not at equal density -- density is the free knob, `rho = N / V_box`.
TARGET_N = 300

K_SWEEP = (1, 2, 4, 8, 16)


def probe_point(label: str, w: float, du: float, dv: float,
                dx: float, dy: float, *, seed: int,
                sprinklings: int = 24,
                mc_samples: int = 200_000) -> dict:
    slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
    curved, flat = arms(slab, w)
    assert slab.du < conjugate_du(w)  # §4.3, structural but stated
    v_box = slab.coordinate_volume
    rho = TARGET_N / v_box
    rng = np.random.default_rng(seed)

    frac_analytic = eligible_volume_fraction(curved)
    p, q = fattest_axis_diamond(curved)
    delta = axis_volume_ratio(w * slab.du) - 1.0

    v_curved, v_flat, v_dis, v_int = diamond_volumes_mc(
        curved, flat, p, q, mc_samples, rng)
    lam = rho * v_curved
    lam_flat = rho * v_flat
    # marginal: relative shift delta on a Poisson(lam) count
    n_marginal = sprinklings_needed(delta, 1.0 / math.sqrt(lam))
    n_marg_be = sprinklings_needed(delta, 1.0 / math.sqrt(lam), z_beta=0.0)
    # DETECTION sizing (R1.1: this is not P2's estimand and is no
    # longer labeled as it). Null: the volumes do not differ;
    # alternative: they differ by the predicted rho V_0 delta.
    # Statistic: the raw same-points difference N_A - N_0, whose
    # per-realization variance is EXACTLY rho V_dis -- shared points
    # cancel wherever the two predicates agree. The signal is the
    # PREDICTION, never the MC difference of the volumes (fourth-order
    # in a, drowned by MC noise that is first-order in V_dis; an early
    # run returned exactly 0).
    assert delta <= v_dis / v_curved + 4.0 / math.sqrt(mc_samples), (
        f"{label}: delta {delta} exceeds the disagreement fraction "
        f"{v_dis / v_curved} -- |V_A - V_0| <= V_dis is an identity, so "
        "one of the two is being computed wrongly")
    signal_abs = rho * v_flat * delta
    n_detect = sprinklings_needed(signal_abs, math.sqrt(rho * v_dis))
    n_detect_be = sprinklings_needed(signal_abs, math.sqrt(rho * v_dis),
                                     z_beta=0.0)
    # P2's preregistered residual Z = N_A - r N_0, r predicted:
    # E[Z] = 0 under the prediction, so there is nothing to detect --
    # sd_z sizes the PRECISION of the verification. After n
    # sprinklings the measured ratio carries a standard error of
    # sd_z / (rho v_flat sqrt(n)); P2 sets its own tolerance against
    # that. This probe only reports the number (R1.1).
    r_pred = 1.0 + delta
    var_z = rho * (v_curved + r_pred ** 2 * v_flat
                   - 2.0 * r_pred * v_int)
    sd_z = math.sqrt(max(var_z, 0.0))

    elem_fracs, pair_fracs, amb_fracs = [], [], []
    per_pair_us = []
    escalations = escalations_flat = 0
    ks = {k: [] for k in K_SWEEP}
    for _ in range(sprinklings):
        pts = sprinkle(curved, rho, rng)
        n = len(pts)
        if n < 2:
            continue
        elig = np.array([element_eligible(curved, pt) for pt in pts])
        m = int(elig.sum())
        elem_fracs.append(m / n)
        pair_fracs.append((m * (m - 1)) / (n * (n - 1)))
        census = relation_census(curved, pts)
        flat_census = relation_census(flat, pts)
        # R1.2: a pair is ambiguous when EITHER arm is undecided
        union = set(census.ambiguous_pairs) | set(
            flat_census.ambiguous_pairs)
        amb_fracs.append(len(union) / census.pairs)
        escalations += census.escalated
        escalations_flat += flat_census.escalated
        per_pair_us.append(census.seconds / census.pairs * 1e6)
        below, above = order_interiority(census.related)
        order = np.argsort(pts[:, 0], kind="stable")
        guard_sorted = elig[order]
        for k in K_SWEEP:
            ks[k].append(agreement(candidate_mask(below, above, k),
                                   guard_sorted))

    return {
        "label": label, "w": w, "a": w * du, "slab": slab, "rho": rho,
        "frac_analytic": frac_analytic,
        "elem_frac": float(np.mean(elem_fracs)),
        "pair_frac": float(np.mean(pair_fracs)),
        "delta": delta,
        "lam": lam, "lam_flat": lam_flat,
        "v_curved": v_curved, "v_flat": v_flat, "v_dis": v_dis,
        "v_int": v_int, "sd_z": sd_z,
        "n_marginal": n_marginal, "n_marginal_be": n_marg_be,
        "n_detect": n_detect, "n_detect_be": n_detect_be,
        "ambiguous_fraction": float(np.mean(amb_fracs)),
        "escalations": escalations,
        "escalations_flat": escalations_flat,
        "per_pair_us": float(np.mean(per_pair_us)),
        "sprinklings": len(elem_fracs),
        "agreement": {k: ks[k] for k in K_SWEEP},
    }


def main(seed: int = 20260808) -> None:
    print(f"P14 probe P1 -- exploratory, seed {seed}, "
          f"N = {TARGET_N} per sprinkling\n")

    rows = []
    for point in OPERATING_POINTS:
        row = probe_point(*point, seed=seed)
        rows.append(row)
        print(f"  {row['label']} done ({row['sprinklings']} sprinklings, "
              f"{row['per_pair_us']:.1f} us/pair)")

    print("\n== eligibility ==")
    print(f"{'point':>12} {'a':>4} {'f_vol':>7} {'f_elem':>7} "
          f"{'f_pair':>7} {'pairs/spr':>9}")
    for r in rows:
        pairs = r["pair_frac"] * TARGET_N * (TARGET_N - 1) / 2
        print(f"{r['label']:>12} {r['a']:4.1f} {r['frac_analytic']:7.4f} "
              f"{r['elem_frac']:7.4f} {r['pair_frac']:7.4f} {pairs:9.1f}")

    print("\n== feasibility, fattest eligible axis diamond ==")
    print(f"{'point':>12} {'a':>4} {'delta':>9} {'lam_A':>8} "
          f"{'V_dis/V_A':>9} {'n90_marg':>9} {'n90_detect':>10} "
          f"{'sd_Z':>8}")
    for r in rows:
        print(f"{r['label']:>12} {r['a']:4.1f} {r['delta']:9.2e} "
              f"{r['lam']:8.2f} {r['v_dis'] / r['v_curved']:9.4f} "
              f"{r['n_marginal']:9.3g} {r['n_detect']:10.3g} "
              f"{r['sd_z']:8.3f}")

    print("\n== order-invariant candidate vs coordinate guard "
          "(PAIR level, mean over sprinklings; R1.3) ==")
    print(f"{'point':>12} " + " ".join(f"{'k=' + str(k):>15}"
                                       for k in K_SWEEP))
    print(f"{'':>12} " + " ".join(f"{'agree/c-only%':>15}"
                                  for _ in K_SWEEP))
    for r in rows:
        cells = []
        for k in K_SWEEP:
            ags = r["agreement"][k]
            rate = float(np.mean([a.pair_rate for a in ags]))
            conly = float(np.mean(
                [a.pair_candidate_only / a.pair_total for a in ags]))
            cells.append(f"{rate:7.3f}/{conly:6.3f}")
        print(f"{r['label']:>12} " + " ".join(f"{c:>15}" for c in cells))
    print("  (element-level tables: rerun with the same seed and read "
          "`.rate` / `.candidate_only` off the returned Agreements)")

    print("\n== ambiguity (union over BOTH arms, §5.1) and cost ==")
    print(f"{'point':>12} {'ambig_frac':>10} {'esc_curved':>10} "
          f"{'esc_flat':>8} {'us/pair':>8}")
    for r in rows:
        print(f"{r['label']:>12} {r['ambiguous_fraction']:10.2e} "
              f"{r['escalations']:10d} {r['escalations_flat']:8d} "
              f"{r['per_pair_us']:8.1f}")

    generic_us, escalated_us = escalation_cost_microbench(
        arms(rows[0]["slab"], rows[0]["w"])[0])
    print(f"\nescalation micro-bench: generic {generic_us:.1f} us, "
          f"escalated {escalated_us:.1f} us "
          f"({escalated_us / generic_us:.0f}x)")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20260808)
