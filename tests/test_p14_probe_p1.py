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
import warnings
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
    operating point and 14x at another. It was flagged by the
    consistency diagnostic `delta*V_0 <= V_dis` (a warning, since V_dis
    is measured; the hard per-run assertion is the exact sample-count
    partition inside `diamond_volumes_mc` -- review R10).

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


def test_resolution_checks_both_ends_of_the_ci_not_just_the_upper():
    """The case, review R8.1: `disagree_resolved` promised the whole
    95% CI within `_MC_CI_FACTOR` of the point, but checked only
    `ucb <= factor*point`. The CP interval is asymmetric, so 8-11 hits
    at `samples = 2e5` pass the upper test (`ucb ~ 1.97x point`) while
    the lower endpoint is `~0.43x point` -- outside factor 2 on the
    low side. Those rows would switch `n90 detect`/`sd_Z` from a
    bracket to a point while half the interval is still out of range.
    Both endpoints are required now.
    """

    # 8 hits: upper passes, lower fails -> NOT resolved
    eight = probe.DiamondVolumes(curved=1.0, flat=1.0, disagree=8.0,
                                 intersect=0.0, quantum=1.0,
                                 disagree_hits=8, samples=200_000)
    assert eight.disagree_ucb <= probe._MC_CI_FACTOR * eight.disagree
    assert eight.disagree_lcb < eight.disagree / probe._MC_CI_FACTOR
    assert not eight.disagree_resolved

    # a count large enough that BOTH ends are inside the factor
    many = probe.DiamondVolumes(curved=1.0, flat=1.0, disagree=250.0,
                                intersect=0.0, quantum=1.0,
                                disagree_hits=250, samples=200_000)
    assert many.disagree_lcb >= many.disagree / probe._MC_CI_FACTOR
    assert many.disagree_ucb <= probe._MC_CI_FACTOR * many.disagree
    assert many.disagree_resolved


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
    # the aniso-a1.0 geometry: its diamond has a large enough V_dis
    # that 200k samples give ~150 hits, so V_dis resolves under the
    # two-sided rule (R8.1) and the R5.2 point-vs-UCB check can run.
    # The MC runs once per probe_point, so this stays cheap.
    row = probe.probe_point("mini", 1.0, 1.0, 1.0, 2.0, 6.0,
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
    # R6.1: the unpaired sizing is a detectability test of lam_A vs
    # lam_0 with the null at the OTHER arm's mean -- NOT §8 P2's
    # each-arm-against-its-own-mean check (zero displacement under the
    # prediction). Distinct SDs under null (sqrt lam_0) and alternative
    # (sqrt lam_A). Reconstruct and check.
    z_a, z_b = 1.959964, 1.281552
    lam_0 = row["lam"] / (1.0 + row["delta"])
    gap = row["lam"] - lam_0
    expected = ((z_a * math.sqrt(lam_0) + z_b * math.sqrt(row["lam"]))
                / gap) ** 2
    assert row["n_unpaired"] == pytest.approx(expected, rel=1e-9)
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
    # (sd_Z's exact identity is pinned above; no geometry-hardcoded
    # duplicate here.)
    assert row["ambiguous_fraction"] == 0.0
    assert row["escalations"] == 0 and row["escalations_flat"] == 0
    for k, agreements in row["agreement"].items():
        assert len(agreements) == row["sprinklings"], k


def test_a_gross_delta_error_warns_but_does_not_crash(monkeypatch):
    """The case, review R10: the analytic-vs-MC consistency check
    compares the analytic floor `delta*V_0` against `disagree_ucb`, a
    RANDOM 95% CP upper limit. Because a correct run can see that limit
    dip below the true V_dis on its ~2.5% upper tail -- and the script
    accepts arbitrary seeds -- this must NOT be a hard assertion that
    turns a valid draw into a deterministic crash. It is a diagnostic
    warning: a gross delta error (here a forced 100x ratio) blows past
    the bound and warns, while the run still completes and returns its
    row. The hard, DETERMINISTIC invariant is instead the exact
    sample-count partition asserted inside `diamond_volumes_mc`, whose
    effect on the derived volumes is already pinned by
    `test_the_mc_volume_estimator_matches_the_flat_alexandrov_volume`.
    """

    monkeypatch.setattr(probe, "TARGET_N", 60)
    monkeypatch.setattr(probe, "axis_volume_ratio", lambda a, **k: 100.0)
    with pytest.warns(UserWarning, match=r"analytic floor delta\*V_0"):
        row = probe.probe_point("gross", 1.0, 1.0, 1.0, 2.0, 6.0,
                                 seed=3, sprinklings=3, mc_samples=50_000)
    # completed instead of crashing, and delta reflects the patched ratio
    assert row["delta"] == pytest.approx(99.0)
    assert math.isfinite(row["n_detect"]) and row["n_detect"] > 0.0
    # R9.2: continuing past the warning must not emit a REVERSED bracket.
    # This gross error drives analytic_floor above disagree_ucb, exactly
    # the case that would invert [floor, ucb]; the upper endpoint is
    # reconciled to the deterministic floor, so every reported interval
    # stays ordered (main() prints them as [lower, upper]).
    assert row["v_dis_floor"] <= row["v_dis_ucb"]
    assert row["n_detect_floor"] <= row["n_detect_ucb"]
    assert row["sd_z_floor"] <= row["sd_z_ucb"]


def test_a_floor_conflict_forces_a_resolved_looking_row_unresolved(
        monkeypatch):
    """The case, review R9.3: R9.2 ordered the bracket, but the bracket
    is only REPORTED when the row is unresolved, and `resolved` came from
    the raw CP width alone. In the warning's conflict case the MC point
    is `v_dis <= disagree_ucb < analytic_floor`, i.e. below the
    deterministic floor delta*V_0 -- an impossible V_dis -- yet a tight
    (many-hit) CP interval can still read as resolved. The row would then
    be reported by that impossible point estimate (headline `n_detect`
    from `v_dis`, point `sd_Z`). The floor conflict now forces the row
    unresolved, so the constrained bracket is what gets reported.
    """

    monkeypatch.setattr(probe, "TARGET_N", 60)
    # a large forced ratio drives analytic_floor above disagree_ucb ...
    monkeypatch.setattr(probe, "axis_volume_ratio", lambda a, **k: 100.0)
    with pytest.warns(UserWarning, match=r"analytic floor delta\*V_0"):
        row = probe.probe_point("conflict", 1.0, 1.0, 1.0, 2.0, 6.0,
                                 seed=5, sprinklings=3, mc_samples=200_000)
    # ... while the raw CP interval on those hits is tight enough to read
    # as resolved on its own -- a RESOLVED-LOOKING conflict, not a
    # few-hit unresolved row
    raw = probe.DiamondVolumes(
        curved=1.0, flat=1.0,
        disagree=row["v_dis_hits"] * row["v_dis_quantum"],
        intersect=0.0, quantum=row["v_dis_quantum"],
        disagree_hits=row["v_dis_hits"], samples=200_000)
    assert raw.disagree_resolved, "fixture is not a resolved-looking conflict"
    # the row is nonetheless reported unresolved, so the headline is the
    # reconciled bracket (ucb), never the point below the floor
    assert row["v_dis_resolved"] is False
    assert row["n_detect"] == row["n_detect_ucb"]


def test_a_point_below_the_floor_cannot_carry_a_resolved_headline(
        monkeypatch):
    """The case, review R9.4: R9.3 gated resolution on
    `floor <= disagree_ucb`, which only proves the CONFIDENCE INTERVAL
    is compatible with the deterministic floor -- not the point. A
    modest downward fluctuation gives `v_dis < floor <= ucb`: the raw
    CP interval is tight enough to resolve on its own, NO warning fires
    (the floor does not exceed the ucb), yet the point sits below its
    exact minimum `delta*V_0`, and R9.3's gate would still quote it as
    the headline. Resolution now requires the RAW point itself to
    satisfy the floor (which strictly subsumes the interval check,
    since the CP interval contains its point), and every point-derived
    output reads the floor-constrained `max(v_dis, floor)` -- the class
    fix: one constrained view of V_dis, applied once, consumed
    everywhere.
    """

    slab = Slab(du=1.0, dv=1.0, dx=2.0, dy=6.0)
    curved, flat = arms(slab, 1.0)
    p, q = probe.fattest_axis_diamond(curved)
    # replicate probe_point's MC volumes: the rng's first consumer is
    # diamond_volumes_mc, so the same seed reproduces them exactly
    vols = probe.diamond_volumes_mc(curved, flat, p, q, 200_000,
                                    np.random.default_rng(2))
    assert vols.disagree_resolved  # raw CP interval resolves on its own
    dv_win = float(p[1]) - float(q[1])
    v0 = math.pi * slab.du ** 2 * dv_win ** 2 / 6.0
    # force the analytic floor 5% ABOVE the measured point but still
    # below the CP upper limit: `v_dis < floor <= ucb`, R9.4's case,
    # DISTINCT from the R9.3 fixture where the floor exceeds the ucb
    floor = 1.05 * vols.disagree
    assert floor <= vols.disagree_ucb, "fixture must not exceed the ucb"
    monkeypatch.setattr(probe, "axis_volume_ratio",
                        lambda a, **k: 1.0 + floor / v0)
    monkeypatch.setattr(probe, "TARGET_N", 60)
    with warnings.catch_warnings():
        # distinct from the R9.1/R9.3 fixtures: this case is NOT a
        # gross error, so the diagnostic must stay silent
        warnings.simplefilter("error", UserWarning)
        row = probe.probe_point("dip", 1.0, 1.0, 1.0, 2.0, 6.0,
                                seed=2, sprinklings=3,
                                mc_samples=200_000)
    # the point below its exact minimum may not carry a point headline
    assert row["v_dis_resolved"] is False
    assert row["n_detect"] == row["n_detect_ucb"]
    # and the point-derived sd_z reads the constrained point (the
    # floor), never the impossible raw point below it
    assert row["sd_z"] == pytest.approx(math.sqrt(
        (1.0 + row["delta"]) * row["rho"] * floor))


SIZING_JSON = (Path(__file__).resolve().parents[1]
               / "docs" / "prereg" / "p14_probe_p1_sizing.json")


#: Cross-platform float tolerance model for the artifact reproduction
#: test. The artifact is committed from one platform and recomputed on
#: another; libm and LAPACK differ at the ulp level between them, and a
#: single ulp (~2.2e-16) in the O(1) quadrature `R(a)` AMPLIFIES in
#: `delta = R - 1` by `1/delta` -- at roomy-a0.2 (delta = 6.35e-6) one
#: ulp is a 3.5e-11 RELATIVE shift, and `n ~ 1/delta^2` doubles it,
#: which is exactly the 7e-11 the first uniform-1e-12 version of this
#: test died on in CI (a Windows-committed artifact, a Linux rerun).
#:
#: The budget is DERIVED from that mechanism, not chosen round (PR #43
#: review: a `1e-11/delta` numerator -- five orders above ulp scale --
#: admits ~1.6e-6 relative drift at roomy, ~22,000x the observed
#: platform envelope, letting `n_detect` be hand-edited by ~55 while
#: the test's contract says edits are pinned). `_R_ROUNDOFF` is a
#: 100-ulp allowance on `R(a)` (observed: 1 ulp); `delta` is compared
#: ABSOLUTELY against it, and `_AMPLIFIED_FIELDS` -- everything
#: scaling with `delta^-1` or `delta^-2` -- get
#: `rel = 1e-11 + 2 * _R_ROUNDOFF / delta` (the worst power's factor
#: 2). At roomy that is ~7e-9: two orders above the observed drift,
#: 200x tighter than before, and it pins `n_detect` to +-0.25 counts.
#:
#: The relative branches pass `abs=0.0` explicitly: pytest.approx
#: otherwise keeps its DEFAULT abs=1e-12 and accepts when EITHER
#: tolerance holds, which for small fields (slice-a0.3's
#: `v_dis_floor` = 2.7e-8) widens the intended window ~27,000x
#: (PR #44 review). Classification re-audited before removing the
#: floor: the analytic floor is delta-LINEAR only where it binds
#: (zero-hit rows' `v_dis_floor`/`sd_z_floor`, in `_AMPLIFIED`);
#: every `n_*` divides by `(delta*lam0)^2` (in `_AMPLIFIED`);
#: `lam`/`lam_flat` carry delta roundoff only RELATIVELY (~2e-14,
#: inside 1e-11); `v_dis_ucb`/`sd_z_ucb` are CP-bound with the floor
#: orders below; `sd_z` is floor-bound on zero-hit rows (`v_dis = 0`
#: makes `v_dis_point = max(0, floor)` -- there `sd_z = sd_z_floor ~
#: sqrt((1+delta) delta)`, a half-power of delta, ~3.5e-10 relative
#: at slice-a0.3), so it sits in the amplified class too (PR #44
#: review R2 -- the first audit wrongly called it point-bound
#: everywhere); and exact-zero fields (`v_dis` at zero hits) compare
#: as equal under `abs=0.0`.
_R_ROUNDOFF = 100.0 * math.ulp(1.0)
_AMPLIFIED_FIELDS = frozenset({
    "n_unpaired", "n_unpaired_be",
    "n_detect", "n_detect_be", "n_detect_floor", "n_detect_ucb",
    "v_dis_floor", "sd_z", "sd_z_floor",
})


def test_the_sizing_artifact_reproduces_field_for_field():
    """The case, PR #41 review R1: pinning only four MC fields of the
    committed artifact leaves `TARGET_N`, the quadrature, every derived
    sizing number, and any stale hand-edit of the JSON free to drift
    while the tests stay green -- the same class the erratum closed for
    the Markdown table, moved into JSON. So the whole payload is
    regenerated from the recorded seed through `sizing_artifact` -- the
    SAME function that writes the file -- and compared field by field:
    ints, bools, and strings exactly; floats under the documented
    cross-platform tolerance model above (the uniform 1e-12 of the
    first version was refuted by a single quadrature ulp amplified
    through the smallest delta). A change to any input (`TARGET_N`,
    `axis_volume_ratio`, the guard, the MC) or any edit to the
    committed file still fails this test.
    """

    import json

    committed = json.loads(SIZING_JSON.read_text(encoding="utf-8"))
    fresh = probe.sizing_artifact(committed["seed"],
                                  committed["mc_samples"])
    assert fresh["target_n"] == committed["target_n"]
    assert fresh["script"] == committed["script"]
    assert len(fresh["points"]) == len(committed["points"])
    for f_rec, c_rec in zip(fresh["points"], committed["points"],
                            strict=True):
        assert set(f_rec) == set(c_rec), c_rec.get("label")
        delta = c_rec["delta"]
        for k, cv in c_rec.items():
            fv = f_rec[k]
            where = f"{c_rec['label']}.{k}"
            if isinstance(cv, (str, bool)) or isinstance(cv, int):
                assert fv == cv, where
            elif k == "delta":
                assert fv == pytest.approx(cv, abs=_R_ROUNDOFF), where
            elif k in _AMPLIFIED_FIELDS:
                assert fv == pytest.approx(
                    cv, rel=1e-11 + 2.0 * _R_ROUNDOFF / delta,
                    abs=0.0), where
            else:
                assert fv == pytest.approx(cv, rel=1e-11,
                                           abs=0.0), where


def test_the_results_doc_embeds_the_rendered_feasibility_table():
    """The doc's feasibility table is RENDERED from the artifact by
    `feasibility_table` and embedded verbatim -- so every cell of every
    column is pinned to the computation, not only the hits column the
    erratum corrected (PR #41 review R1: validating one column still
    let the other six drift).
    """

    import json

    art = json.loads(SIZING_JSON.read_text(encoding="utf-8"))
    doc = (SIZING_JSON.parent / "p14_probe_p1_results.md").read_text(
        encoding="utf-8")
    assert probe.feasibility_table(art) in doc
