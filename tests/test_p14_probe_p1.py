"""Regressions for the P14 §8 P1 probe machinery.

These pin the MEASURING instruments, not the measured numbers -- the
run's numbers live in `docs/prereg/p14_probe_p1_results.md` with their
seed, and re-measuring them is a script run, not a test. What a test
must catch is the instrument drifting: the volume-ratio quadrature, the
MC volume estimator, the order-interiority counts, and the eligibility
factorization the probe leans on.
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

import p14_probe_p1 as probe  # noqa: E402
from p14_plane_wave import (  # noqa: E402
    Slab,
    arms,
    class_c_eligible,
    sprinkle,
)


def test_the_axis_volume_ratio_reproduces_the_design_pinned_values():
    """The case: the probe's effect size delta comes entirely from this
    quadrature, and the FIRST version of the probe used
    `(w tau)^4 / 252` with `tau^2 = 2 du dv` instead -- importing `dv`
    into an effect `dv` cancels out of, wrong by 2.4e4x at one
    operating point and 14x at another. It was caught by the
    consistency bound `delta <= V_dis / V_A`, an identity the probe now
    asserts on every run.

    The pin is against the design's own committed check script
    (`p14_interval_volume_constant_a.py`), whose values this must
    reproduce to the digit -- agreement between two independently
    written quadratures is the cross-validation.
    """

    assert probe.axis_volume_ratio(1.0) == pytest.approx(
        1.00400047, abs=5e-8)
    assert probe.axis_volume_ratio(2.0) == pytest.approx(
        1.07300802, abs=5e-8)

    # small a: the (wT)^4/252 series with the design's own next
    # coefficient bounds the truncation
    for a in (0.3, 0.5):
        series = 1.0 + a ** 4 / 252.0
        next_term = 3.1951e-05 * a ** 8
        assert abs(probe.axis_volume_ratio(a) - series) < 2.0 * next_term

    # and the series is NOT usable where the guard admits: at a = 2.4
    # the quadrature exceeds it by five points of ratio, which is the
    # reason the probe never evaluates the series
    assert probe.axis_volume_ratio(2.4) - (1.0 + 2.4 ** 4 / 252.0) > 0.04


def test_the_mc_volume_estimator_matches_the_flat_alexandrov_volume():
    """The case: lambda and V_dis come from `diamond_volumes_mc`, so a
    biased sampler poisons the whole feasibility table. The flat axis
    diamond is the one case with a closed answer -- slicing at `u'`,
    the transverse integral of the thickness is
    `pi Dv^2 u'(s-u')/s`, and integrating over `u'` gives
    `V = pi s^2 Dv^2 / 6`. The estimator must land on it within MC
    error, and the curved volume must exceed the flat one (defocusing
    wins over focusing at fourth order -- that is the sign of delta).
    """

    slab = Slab(du=1.0, dv=1.0, dx=2.0, dy=6.0)
    curved, flat = arms(slab, 1.0)
    p, q = probe.fattest_axis_diamond(curved)
    dv_win = float(p[1]) - float(q[1])
    exact = math.pi * slab.du ** 2 * dv_win ** 2 / 6.0

    rng = np.random.default_rng(7)
    vols = probe.diamond_volumes_mc(curved, flat, p, q, 120_000, rng)
    v_curved, v_flat = vols.curved, vols.flat
    v_dis, v_int = vols.disagree, vols.intersect

    # binomial MC error: sd ~ sqrt(V (R - V)) / sqrt(n), R the region
    region = slab.du * slab.dv * slab.dx * slab.dy
    sd = math.sqrt(exact * (region - exact) / 120_000)
    assert v_flat == pytest.approx(exact, abs=4.0 * sd)
    assert v_curved > v_flat
    # the identity the probe asserts per run, seen once here directly
    assert v_curved - v_flat <= v_dis + 1e-12
    # and the four volumes are one partition: A + 0 = dis + 2*int,
    # exact per-sample, so exact for the estimates too
    assert v_curved + v_flat == pytest.approx(v_dis + 2.0 * v_int,
                                              rel=1e-12)
    # resolution is part of the return: every estimate is a count of
    # quanta, and here the disagreement was actually sampled
    region = slab.du * slab.dv * slab.dx * slab.dy
    assert vols.quantum == pytest.approx(region / 120_000)
    assert vols.disagree_resolved
    # R3.1/R5.1: the UCB is the Clopper-Pearson upper limit for the
    # observed count, ABOVE the point estimate but within the reporting
    # factor, which is what "resolved" means.
    assert vols.disagree_ucb > vols.disagree
    assert vols.disagree_ucb <= probe._MC_CI_FACTOR * vols.disagree
    assert vols.disagree_hits > 0 and vols.samples == 120_000


def test_clopper_pearson_is_exact_and_two_sided():
    """The case, review R5.1: the first bracket paired two one-sided
    95% limits, giving ~90% two-sided coverage, and used the Poisson
    approximation where the sampler is a fixed-n Bernoulli experiment.
    Both fixed: exact binomial Clopper-Pearson with a 2.5% tail on each
    end. These reference values are the standard exact CP intervals.
    """

    lo, hi = probe.clopper_pearson(0, 100)
    assert lo == 0.0 and hi == pytest.approx(0.03621, abs=1e-4)
    lo, hi = probe.clopper_pearson(2, 100)
    assert lo == pytest.approx(0.00243, abs=1e-4)
    assert hi == pytest.approx(0.07038, abs=1e-4)
    # each tail is 2.5%, verified against the binomial CDF directly
    lo, hi = probe.clopper_pearson(5, 1000)
    assert probe._binom_cdf(5, 1000, hi) == pytest.approx(0.025, abs=1e-4)
    assert 1.0 - probe._binom_cdf(4, 1000, lo) == pytest.approx(
        0.025, abs=1e-4)
    # all-hits and no-hits collapse the matching endpoint
    assert probe.clopper_pearson(0, 50)[0] == 0.0
    assert probe.clopper_pearson(50, 50)[1] == 1.0

    # a single hit is NOT resolved -- its 95% CP upper is many times
    # the point estimate
    one = probe.DiamondVolumes(curved=1.0, flat=1.0, disagree=1.0,
                               intersect=0.0, quantum=1.0,
                               disagree_hits=1, samples=1000)
    assert not one.disagree_resolved
    assert one.disagree_ucb > 3.0 * one.disagree


def test_the_cp_bracket_is_a_true_two_sided_interval():
    """The case, review R4.1/R5.1: the bracket's lower endpoint must be
    a genuine lower bound, not the MC point estimate, and the two
    endpoints together must be a real 95% interval -- 2.5% in each
    tail, not 5%. The point estimate lies strictly inside."""

    for hits, samples in ((1, 1000), (5, 1000), (85, 200_000)):
        lo, hi = probe.clopper_pearson(hits, samples)
        point = hits / samples
        assert 0.0 < lo < point < hi
        # 90% (two one-sided 95%) would be narrower than this 95%
        assert probe._binom_cdf(hits, samples, hi) == pytest.approx(
            0.025, abs=2e-3)


def test_sd_z_is_the_exact_r_rho_vdis_identity_not_a_recombination():
    """The case, review R4.2: with the predicted r = V_A/V_0, the
    residual variance Var(Z) = rho(V_A + r^2 V_0 - 2 r V_int) reduces
    ALGEBRAICALLY to r*rho*V_dis, because V_A - r V_0 = 0 kills the
    cross terms. So sd_Z is driven by the same V_dis as the detection
    sizing and must inherit its resolution, rather than recombining
    three noisy MC volumes. This checks the identity on the MC volumes
    themselves (they satisfy the partition exactly per sample) and
    that the reported sd_Z uses it.
    """

    slab = Slab(du=1.0, dv=1.0, dx=2.0, dy=6.0)
    curved, flat = arms(slab, 1.0)
    p, q = probe.fattest_axis_diamond(curved)
    rng = np.random.default_rng(7)
    vols = probe.diamond_volumes_mc(curved, flat, p, q, 120_000, rng)
    r = vols.curved / vols.flat
    recombined = vols.curved + r ** 2 * vols.flat - 2.0 * r * vols.intersect
    assert recombined == pytest.approx(r * vols.disagree, rel=1e-9)


def test_a_zero_hit_disagreement_is_a_bound_and_never_a_zero_sd():
    """The case, review R2.1: at `slice-a0.3` the 200k MC sample saw
    zero disagreeing points, the point estimate 0 flowed into
    `sqrt(rho * v_dis)`, and the detection sizing returned 0 -- a
    zero-noise instrument manufactured by finite MC resolution, while
    the analytic `delta > 0` guarantees `V_dis > 0`. The estimator now
    carries its quantum: a zero count means `<= 3` quanta at 95%
    (rule of three), the sizing uses that bound, and the row is
    flagged unresolved with a floor from `|V_A - V_0| = delta V_0`.
    """

    monkey = pytest.MonkeyPatch()
    try:
        monkey.setattr(probe, "TARGET_N", 60)
        row = probe.probe_point("zero-hit", 1.0, 0.3, 0.2, 6.0, 6.0,
                                seed=4, sprinklings=2, mc_samples=4_000)
    finally:
        monkey.undo()
    assert not row["v_dis_resolved"], (
        "the fixture resolved V_dis; shrink mc_samples so the zero-hit "
        "path is exercised")
    assert row["n_detect"] > 0.0 and math.isfinite(row["n_detect"])
    assert row["n_detect_floor"] > 0.0
    assert row["n_detect"] >= row["n_detect_floor"]
    # the unresolved sizing rests on the rule-of-three quantum, not 0
    assert row["v_dis_quantum"] > 0.0


def test_element_eligibility_factorizes_exactly_as_the_probe_assumes():
    """The case: the probe derives per-element labels by calling the
    pair predicate on `(p, p)` and then treats a pair as eligible iff
    both elements are. That is only sound while `class_c_eligible`
    factorizes per endpoint. If eligibility ever became pair-coupled
    -- a relative-separation term, say -- every fraction in the probe
    would silently mean something else. This pins the factorization on
    real sprinkled points, including ineligible ones.
    """

    slab = Slab(du=1.0, dv=0.2, dx=6.0, dy=6.0)
    curved, _ = arms(slab, 1.0)
    rng = np.random.default_rng(3)
    pts = sprinkle(curved, 200.0 / slab.coordinate_volume, rng)
    labels = [probe.element_eligible(curved, pt) for pt in pts]
    assert any(labels) and not all(labels), "nothing to compare"
    for i in range(0, len(pts), 7):
        for j in range(0, len(pts), 11):
            if i == j:
                continue
            assert class_c_eligible(curved, pts[i], pts[j]) == (
                labels[i] and labels[j])


def test_order_interiority_counts_agree_with_an_independent_recount():
    """The case: below/above counts feed the §4.6.2 candidate, and the
    probe computes them as axis sums of the relation matrix. The
    oracle here recounts them element by element from the same matrix
    with plain loops -- independent of the numpy axis convention,
    which is exactly the thing to get backwards silently."""

    rng = np.random.default_rng(11)
    n = 40
    upper = np.triu(rng.random((n, n)) < 0.3, k=1)
    below, above = probe.order_interiority(upper)
    for i in range(n):
        assert below[i] == sum(upper[j, i] for j in range(n))
        assert above[i] == sum(upper[i, j] for j in range(n))


def test_the_candidate_mask_shrinks_monotonically_and_k0_admits_all():
    rng = np.random.default_rng(5)
    below = rng.integers(0, 20, 200)
    above = rng.integers(0, 20, 200)
    assert probe.candidate_mask(below, above, 0).all()
    previous = None
    for k in (1, 2, 4, 8, 16, 32):
        mask = probe.candidate_mask(below, above, k)
        if previous is not None:
            assert not np.any(mask & ~previous), "mask grew with k"
        previous = mask
    assert not probe.candidate_mask(below, above, 32).any()


def test_the_agreement_census_partitions_the_elements():
    rng = np.random.default_rng(9)
    cand = rng.random(500) < 0.3
    guard = rng.random(500) < 0.2
    a = probe.agreement(cand, guard)
    assert a.total == 500
    assert a.both + a.candidate_only == int(cand.sum())
    assert a.both + a.guard_only == int(guard.sum())
    assert 0.0 <= a.rate <= 1.0


def test_pair_level_admissions_match_materialized_pair_masks():
    """The case, review R1.3: §8 P1 asks for the PAIRS each rule
    admits that the other rejects, and element rates do not equal pair
    rates -- admissions compound over two endpoints. The probe derives
    the pair counts combinatorially from the four element categories;
    the oracle here materializes the O(N^2) pair masks and counts
    directly. If the derivation dropped a cross term (a candidate-only
    endpoint paired with a `both` endpoint, say), only the oracle
    would see it."""

    rng = np.random.default_rng(21)
    for trial in range(5):
        cand = rng.random(60) < rng.uniform(0.1, 0.6)
        guard = rng.random(60) < rng.uniform(0.1, 0.6)
        a = probe.agreement(cand, guard)

        upper = np.triu(np.ones((60, 60), dtype=bool), k=1)
        cand_pairs = np.outer(cand, cand) & upper
        guard_pairs = np.outer(guard, guard) & upper
        assert a.pair_total == int(upper.sum()), trial
        assert a.pair_both == int((cand_pairs & guard_pairs).sum()), trial
        assert a.pair_candidate_only == int(
            (cand_pairs & ~guard_pairs).sum()), trial
        assert a.pair_guard_only == int(
            (~cand_pairs & guard_pairs).sum()), trial
        # both exclusive directions are reported, not just candidate
        # (R3.2): an empty candidate has zero candidate-only and large
        # guard-only, which the two together distinguish
        empty = probe.agreement(np.zeros(60, bool), guard)
        assert empty.pair_candidate_only == 0
        assert empty.pair_guard_only == int(
            np.triu(np.outer(guard, guard), 1).sum())
        disagree = a.pair_candidate_only + a.pair_guard_only
        assert a.pair_rate == pytest.approx(1.0 - disagree / a.pair_total)


def test_the_relation_census_is_upper_triangular_and_unambiguous_here():
    """The case: the census sorts by `u` and fills only `i < j`; an
    entry below the diagonal would double-count every relation the
    interiority counts consume. Ambiguity must be zero on a generic
    sprinkling -- pairs land on the cone with probability zero, and a
    nonzero count here would mean the error bound has widened into the
    generic population."""

    slab = Slab(du=1.0, dv=0.2, dx=6.0, dy=6.0)
    curved, _ = arms(slab, 1.0)
    rng = np.random.default_rng(13)
    pts = sprinkle(curved, 80.0 / slab.coordinate_volume, rng)
    census = probe.relation_census(curved, pts)
    assert not np.any(np.tril(census.related))
    assert census.ambiguous == 0
    assert census.ambiguous == len(census.ambiguous_pairs)
    assert census.escalated == 0
    assert census.related.sum() > 0
    assert census.pairs == len(pts) * (len(pts) - 1) // 2


def test_the_escalation_microbench_forces_the_decimal_path():
    """The constructed on-cone pair must actually escalate -- the
    function asserts it internally -- and cost measurably more than a
    generic decision. Uniform sampling cannot reach the cone (measure
    zero), which is why the probe constructs the pair instead of
    hunting for one; 7.5M sampled pairs in the committed run produced
    zero escalations, which is the point."""

    slab = Slab(du=1.0, dv=0.2, dx=6.0, dy=6.0)
    curved, _ = arms(slab, 1.0)
    generic_us, escalated_us = probe.escalation_cost_microbench(
        curved, reps=40)
    assert escalated_us > 3.0 * generic_us


def test_probe_point_holds_its_own_consistency_checks(monkeypatch):
    """One miniature end-to-end run: the measured element fraction must
    sit near the analytic volume fraction (they are estimates of the
    same number), the pair fraction near its square, and the
    feasibility outputs must be finite and positive where the point
    admits pairs. `TARGET_N` is patched down so this stays cheap."""

    monkeypatch.setattr(probe, "TARGET_N", 80)
    # full MC so V_dis resolves (the R5.2 check needs a resolved row);
    # the MC runs once per probe_point, so this stays cheap
    row = probe.probe_point("mini", 1.0, 1.0, 0.2, 6.0, 6.0,
                            seed=2, sprinklings=6, mc_samples=200_000)
    f = row["frac_analytic"]
    assert row["elem_frac"] == pytest.approx(f, abs=4.0 * math.sqrt(
        f * (1 - f) / (80 * 6)))
    assert row["pair_frac"] == pytest.approx(row["elem_frac"] ** 2,
                                             abs=0.01)
    assert row["delta"] == pytest.approx(
        probe.axis_volume_ratio(1.0) - 1.0)
    assert 0.0 < row["lam"] < 80.0
    assert math.isfinite(row["n_detect"]) and row["n_detect"] > 0.0
    assert row["n_detect"] >= row["n_detect_floor"] > 0.0
    # R4.3: the marginal sizing is delta*lam_0 signal against
    # sqrt(lam_A) noise -- exceeds the old delta-vs-1/sqrt(lam_A) form
    # by (1+delta)^2. Reconstruct the corrected value and check.
    z = 1.959964 + 1.281552
    lam_0 = row["lam"] / (1.0 + row["delta"])
    expected = (z * math.sqrt(row["lam"]) / (row["delta"] * lam_0)) ** 2
    assert row["n_marginal"] == pytest.approx(expected, rel=1e-9)
    # R4.2: sd_Z from the exact identity, and its bracket brackets it
    assert row["sd_z"] == pytest.approx(
        math.sqrt((1.0 + row["delta"]) * row["rho"] * row["v_dis"]))
    assert row["sd_z_floor"] <= row["sd_z_ucb"]
    # R5.2: n_detect and sd_z share one resolution semantics. This
    # fixture resolves V_dis, so the headline n_detect is the POINT
    # sizing (from v_dis), strictly below the UCB sizing -- not the UCB
    # silently kept, which is what the two columns disagreeing looked
    # like before.
    assert row["v_dis_resolved"], "fixture no longer resolves V_dis"
    # n scales linearly with V_dis, so the headline (point) sizing over
    # the UCB sizing is exactly v_dis / v_dis_ucb -- confirming the
    # headline uses the point, not the UCB silently kept.
    assert row["n_detect"] / row["n_detect_ucb"] == pytest.approx(
        row["v_dis"] / row["v_dis_ucb"], rel=1e-9)
    assert row["n_detect"] < row["n_detect_ucb"]
    assert row["v_dis_floor"] <= row["v_dis"] <= row["v_dis_ucb"]
    # Var(Z) with the predicted r must sit near rho*V_dis -- they
    # differ at O(delta) -- and never below the parallelogram floor
    assert row["sd_z"] ** 2 == pytest.approx(
        (80.0 / (1.0 * 0.2 * 6.0 * 6.0)) * row["v_dis"], rel=0.05)
    assert row["ambiguous_fraction"] == 0.0
    assert row["escalations"] == 0 and row["escalations_flat"] == 0
    for k, agreements in row["agreement"].items():
        assert len(agreements) == row["sprinklings"], k
