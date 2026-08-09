"""Regressions for the P14 §8 P3-C preflight machinery.

Pins the frozen verdict rules and CI constructions, the seed layout,
and the preflight/P4 artifacts' reproducibility from the committed
P3-E samples. The campaign runner is committed but NOT exercised
end-to-end here (execution is held until preflight approval); its
shared pipeline (`metric_cis`, `classify`) is what these tests pin.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

import p14_probe_p3c as p3c  # noqa: E402
import p14_probe_p3e as p3e  # noqa: E402
from p14_probe_p1 import clopper_pearson  # noqa: E402

PREFLIGHT = (Path(__file__).resolve().parents[1]
             / "docs" / "prereg" / "p14_probe_p3c_preflight.json")
P3E_ART = (Path(__file__).resolve().parents[1]
           / "docs" / "prereg" / "p14_probe_p3e_results.json")


def test_the_seed_layout_is_frozen_and_fresh():
    """The confirmation block runs on 20260851/52; both downgraded
    campaign blocks' streams -- 20260831/32 (invalid boundary CI at
    run time, PR #48 R1) and 20260841/42 (rule and results in one
    commit, freeze ordering unprovable, PR #48 R2) -- are burned
    with everything before them, so no downgraded stream can ever
    feed a confirmation again."""

    assert p3c.CAMPAIGN_SEEDS == {"curved": 20260851, "flat": 20260852}
    for s in (20260808, 777, 778, 779, 780, 781,
              20260811, 20260812, 20260813, 20260814,
              20260821, 20260822, 20260823, 20260824,
              20260831, 20260832, 20260841, 20260842):
        assert s in p3c.BURNED_SEEDS, s
    p3c.assert_seed_layout()


def test_the_frozen_constants_and_margins():
    """n = 4800/arm (the equivalence branch binds), margins are
    p3e's frozen set BY IDENTITY (one source), and the preflight
    replicate counts and seed paths are review-frozen."""

    assert p3c.N_ARM == 4_800
    assert p3c.P3C_MARGINS is p3e.P3C_MARGINS
    assert p3c.PREFLIGHT_NULL_REPS == 20_000
    assert p3c.PREFLIGHT_EFFECT_REPS == 4_000
    assert p3c.PREFLIGHT_NULL_SEED == (781, 50)
    assert p3c.PREFLIGHT_EFFECT_SEED == (781, 51)


def test_the_frozen_ci_constructions():
    """s carries the FULL variance (the null form sqrt(2/n) is
    anti-conservative at large observed s); BA carries its binomial
    variance maximum. AUC is DeLong (separate test below)."""

    n = 400
    lo, hi = p3c.ci_s(2.0, n)
    se = math.sqrt(2.0 / n + 4.0 / (4.0 * (n - 1)))
    assert hi - lo == pytest.approx(2 * 1.959964 * se, rel=1e-12)
    lo0, hi0 = p3c.ci_s(0.0, n)
    assert hi0 - lo0 == pytest.approx(
        2 * 1.959964 * math.sqrt(2.0 / n), rel=1e-12)
    assert (hi - lo) > (hi0 - lo0)  # full variance widens at s != 0
    lo, hi = p3c.ci_ba(0.9, n)
    assert hi - lo == pytest.approx(
        2 * 1.959964 / (2.0 * math.sqrt(n)), rel=1e-12)


def test_the_auc_ci_is_delong_and_beats_the_null_form_off_null():
    """The case, PR #47 review R1: the Mann-Whitney null variance
    `1/(6n)`-form is NOT an upper bound under alternatives, so a CI
    built on it can manufacture `confirmed`. The reviewer's
    counterexample, pinned: flat constant, curved two-valued above/
    below with probability p = 0.55 -- AUC = 0.55 and the DeLong
    variance reproduces the asymptotic `p(1-p)/m`, ~48% ABOVE the
    null form. Plus an exact small hand case with ties, and the
    identity mean(V10) = midrank AUC."""

    n = 1000
    rng = np.random.default_rng(3)
    curved = np.where(rng.random(n) < 0.55, 1.0, -1.0)
    flat = np.zeros(n)
    auc, lo, hi = p3c.auc_delong(curved, flat)
    p_hat = (curved > 0).mean()
    assert auc == pytest.approx(p_hat, rel=1e-12)
    se = (hi - auc) / 1.959964
    want_var = curved.size / (curved.size - 1) * p_hat * (1 - p_hat) / n
    assert se ** 2 == pytest.approx(want_var, rel=1e-9)
    assert se ** 2 > (2 * n + 1) / (12.0 * n * n)  # beats null form

    # exact hand case with ties: a = [2,2,3,3], b = [1,1,2,0]
    a = np.array([2.0, 2.0, 3.0, 3.0])
    b = np.array([1.0, 1.0, 2.0, 0.0])
    auc, lo, hi = p3c.auc_delong(a, b)
    assert auc == pytest.approx(15.0 / 16.0)
    v10 = np.array([0.875, 0.875, 1.0, 1.0])   # placements of a in b
    v01 = np.array([1.0, 1.0, 0.75, 1.0])      # placements of b in a
    var = v10.var(ddof=1) / 4 + v01.var(ddof=1) / 4
    assert ((hi - auc) / 1.959964) ** 2 == pytest.approx(var, rel=1e-9)


def test_the_boundary_rule_replaces_the_degenerate_wald_interval():
    """The case, PR #48 review R1: at complete separation every
    placement is 1, the influence variance vanishes, and the Wald/
    DeLong interval degenerates to [1, 1] -- which is NOT a valid
    finite-sample CI for the POPULATION AUC (the sample separating
    does not make the population overlap exactly zero). The frozen
    boundary rule is the exact FIXED DISJOINT-PAIR bound (review R2:
    the placements share both arms, so they are never independent
    trials): index-pair a_i with b_i, i < k = min(m, n), fixed a
    priori; across unpaired iid streams the indicators 1[a_i > b_i]
    are iid Bernoulli(P(a > b)) with P(a > b) <= AUC, and complete
    separation makes all k successes, so 0.025^(1/k) is an exact
    97.5% one-sided lower bound. Symmetric at zero. At k = 4800 this
    is 0.999232 -- the review's independently computed number."""

    k = 4
    auc, lo, hi = p3c.auc_delong(np.array([3.0, 4.0, 5.0, 6.0]),
                                 np.array([1.0, 2.0, 2.5, 2.9]))
    assert auc == 1.0 and hi == 1.0
    assert lo == pytest.approx(0.025 ** (1.0 / k), rel=1e-12)
    assert lo < 1.0
    # unequal sizes: only min(m, n) disjoint pairs exist, so the six
    # extra curved values add no Bernoulli trials
    auc, lo, hi = p3c.auc_delong(
        np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
        np.array([1.0, 2.0]))
    assert auc == 1.0 and lo == pytest.approx(0.025 ** 0.5, rel=1e-12)
    # symmetric at zero, again with k = min(m, n) = 1
    auc, lo, hi = p3c.auc_delong(np.array([1.0]),
                                 np.array([3.0, 4.0, 5.0]))
    assert auc == 0.0 and lo == 0.0
    assert hi == pytest.approx(1.0 - 0.025, rel=1e-12)
    # the review's n = 4800 value
    assert 0.025 ** (1.0 / 4800) == pytest.approx(0.999232, abs=5e-7)


def test_the_three_branch_verdict_rules():
    """confirmed = every CI entirely outside its band in the FROZEN
    direction; equivalent = every CI entirely inside; everything
    else -- straddles, mixtures, and wrong-direction exceedances --
    is inconclusive."""

    m = p3c.P3C_MARGINS
    es, ea, eb = m["eps_s"], m["eps_auc"], m["eps_ba"]
    inside = ((-es / 2, es / 2), (0.5 - ea / 2, 0.5 + ea / 2),
              (0.5 - eb / 2, 0.5 + eb / 2))
    above = ((es * 2, es * 3), (0.5 + ea * 2, 0.5 + ea * 3),
             (0.5 + eb * 2, 0.5 + eb * 3))
    below = ((-es * 3, -es * 2), (0.5 - ea * 3, 0.5 - ea * 2),
             (0.5 - eb * 3, 0.5 - eb * 2))
    strad = ((0.0, es * 2), (0.5, 0.5 + ea * 2), (0.5, 0.5 + eb * 2))
    assert p3c.classify(*above) == "confirmed"
    assert p3c.classify(*inside) == "equivalent"
    assert p3c.classify(*strad) == "inconclusive"
    # wrong direction: fully outside BELOW the band is never confirmed
    assert p3c.classify(*below) == "inconclusive"
    # a mixture (two above, one inside) is inconclusive
    assert p3c.classify(above[0], above[1], inside[2]) == "inconclusive"


def test_metric_cis_on_hand_data():
    """The shared pipeline on tiny hand data: threshold and direction
    from the first halves, applied to the second halves; midrank AUC
    under ties; equal training means classify everything flat."""

    a = np.array([2.0, 2.0, 3.0, 3.0])
    b = np.array([1.0, 1.0, 2.0, 0.0])
    m = p3c.metric_cis(a, b)
    # training halves: mean a = 2.0, mean b = 1.0 -> thr 1.5, dir +
    # test: a = (3, 3) both > 1.5 -> acc_a = 1; b = (2, 0): <=1.5
    # only for 0 -> acc_b = 0.5; ba = 0.75
    assert m["ba"] == pytest.approx(0.75)
    # AUC with ties: pairs a>b: 2>1 x4, 3>1 x4, 3>2 x2, 2>0 x2, 3>0
    # x2 -> 14; ties 2=2 x2 -> +1 ; total 15/16
    assert m["auc"] == pytest.approx(15.0 / 16.0)
    # the CI in the pipeline IS the DeLong one
    assert m["ci_auc"] == pytest.approx(
        p3c.auc_delong(a, b)[1:], rel=1e-12)
    # equal training means -> everything classifies flat, ba = 0.5
    m2 = p3c.metric_cis(np.array([1.0, 1.0, 5.0, 5.0]),
                        np.array([1.0, 1.0, 0.0, 9.0]))
    assert m2["ba"] == pytest.approx(0.5)


def test_the_preflight_artifact_is_certified_raw_integers():
    """The committed preflight record: raw branch counts summing to
    reps, exact CP endpoints recomputing from the integers, both
    lower bounds >= 0.90, margins matching the frozen set, and the
    independent-resampling statement present (index sharing would
    certify a paired design)."""

    art = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    p = art["preflight"]
    assert p["n_arm"] == p3c.N_ARM
    assert p["margins"] == pytest.approx(p3c.P3C_MARGINS, rel=1e-15)
    assert "INDEPENDENTLY" in p["resampling"]
    for block, reps, key in ((p["null"], p3c.PREFLIGHT_NULL_REPS,
                              "equivalent"),
                             (p["effect"], p3c.PREFLIGHT_EFFECT_REPS,
                              "confirmed")):
        assert block["reps"] == reps
        c = block["counts"]
        assert c["confirmed"] + c["equivalent"] + c["inconclusive"] \
            == reps
        lo, hi = clopper_pearson(c[key], reps)
        assert block[f"{key}_ci95_exact"] == pytest.approx(
            [lo, hi], rel=1e-12)
        assert lo >= 0.90, key
    assert p["certified"] is True


def test_the_p4_block_recomputes_from_the_p3e_samples():
    """P4's headline quantities recompute exactly from the committed
    P3-E aniso samples, and the pair-level joint bootstrap (one index
    vector for both arms) reproduces its CIs from the recorded
    seed."""

    art = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    p4 = art["p4"]
    p3e_art = json.loads(P3E_ART.read_text(encoding="utf-8"))
    fresh = p3c.p4_block(p3e_art)
    assert p4["paired_variance"] == pytest.approx(
        fresh["paired_variance"], rel=1e-12)
    assert p4["variance_gain"] == pytest.approx(
        fresh["variance_gain"], rel=1e-12)
    assert p4["arm_correlation"] == pytest.approx(
        fresh["arm_correlation"], rel=1e-12)
    assert p4["bootstrap"]["gain_ci95"] == pytest.approx(
        fresh["bootstrap"]["gain_ci95"], rel=1e-9)
    assert p4["bootstrap"]["corr_ci95"] == pytest.approx(
        fresh["bootstrap"]["corr_ci95"], rel=1e-9)
    assert "pair-level joint" in p4["bootstrap"]["resampling"]


CAMPAIGN = (Path(__file__).resolve().parents[1]
            / "docs" / "prereg" / "p14_probe_p3c_results.json")
EXPLORATORY = (Path(__file__).resolve().parents[1]
               / "docs" / "prereg"
               / "p14_probe_p3c_results_exploratory.json")
EXPLORATORY2 = (Path(__file__).resolve().parents[1]
                / "docs" / "prereg"
                / "p14_probe_p3c_results_exploratory2.json")


def test_the_campaign_verdict_recomputes_from_the_raw_samples():
    """The committed CONFIRMATION record (fresh streams, corrected
    boundary rule, run from a clean checkout of the freeze commit):
    metrics and verdict must recompute from the stored
    per-sprinkling samples through the same frozen pipeline; seeds,
    sizes, margins, and the frozen positive-sentence scope must be
    the preflight's; ambiguity totals must be recorded and zero
    (undecided is never a silent False, PR #48 review R1); and the
    results doc embeds the rendered verdict table verbatim."""

    if not CAMPAIGN.exists():
        pytest.skip("confirmation block not yet run at this freeze "
                    "point -- its results land in a commit AFTER the "
                    "one freezing the rule and seeds (PR #48 R2)")
    art = json.loads(CAMPAIGN.read_text(encoding="utf-8"))
    assert art["seeds"] == p3c.CAMPAIGN_SEEDS
    assert art["n_arm"] == p3c.N_ARM and art["e_n"] == p3c.E_N
    assert art["margins"] == pytest.approx(p3c.P3C_MARGINS, rel=1e-15)
    assert "고정된 유한 상자·밀도" in art["positive_sentence_scope"]
    for arm_name in ("curved", "flat"):
        assert art["ambiguity"][arm_name]["ambiguous"] == 0
    fa = np.asarray(art["raw"]["f_curved"])
    f0 = np.asarray(art["raw"]["f_flat"])
    assert len(fa) == len(f0) == p3c.N_ARM
    m = p3c.metric_cis(fa, f0)
    for k, v in art["metrics"].items():
        assert v == pytest.approx(m[k], rel=1e-12), k
    assert p3c.classify(m["ci_s"], m["ci_auc"], m["ci_ba"]) \
        == art["verdict"]
    doc = (CAMPAIGN.parent / "p14_probe_p3c_results.md").read_text(
        encoding="utf-8")
    assert p3c.campaign_table(art) in doc


@pytest.mark.parametrize("path,seeds", [
    (EXPLORATORY, {"curved": 20260831, "flat": 20260832}),
    (EXPLORATORY2, {"curved": 20260841, "flat": 20260842}),
], ids=["first-block-R1", "second-block-R2"])
def test_the_downgraded_campaign_records_are_kept_not_erased(path, seeds):
    """The PR #48 review downgraded both earlier campaigns: the
    first ran under the degenerate [1,1] AUC interval (R1), the
    second's results entered history in the same commit as the
    corrected rule, leaving the freeze-before-execution ordering
    unprovable (R2). Each record survives as an exploratory
    artifact: grade stated, verdict downgraded, seeds burned -- a
    protocol correction is transparent, never a deletion."""

    art = json.loads(path.read_text(encoding="utf-8"))
    assert "EXPLORATORY" in art["grade"]
    assert "downgraded" in art["verdict"]
    assert art["seeds"] == seeds
    assert set(art["seeds"].values()) <= set(p3c.BURNED_SEEDS)


@pytest.mark.slow
def test_the_campaign_streams_reproduce_on_a_prefix():
    """The full 2 x 4800-sprinkling rerun is ~35 minutes; the stream
    reproducibility that matters is pinned on a prefix: the first 20
    values of each arm (and their zero ambiguity), recomputed from
    the frozen seeds through the same `arm_samples`, must equal the
    committed raw record exactly -- for the confirmation record
    (once it exists) AND both downgraded exploratory ones."""

    from p14_plane_wave import Slab, arms
    label, w, du, dv, dx, dy = p3c.POINT
    slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
    rho = p3c.E_N / slab.coordinate_volume
    curved, flat = arms(slab, w)
    k = 20
    for path in (CAMPAIGN, EXPLORATORY, EXPLORATORY2):
        if not path.exists():  # results land after the freeze commit
            continue
        art = json.loads(path.read_text(encoding="utf-8"))
        fresh_a, amb_a, _ = p3c.arm_samples(
            curved, rho, art["seeds"]["curved"], k)
        fresh_0, amb_0, _ = p3c.arm_samples(
            flat, rho, art["seeds"]["flat"], k)
        assert amb_a == 0 and amb_0 == 0
        assert fresh_a == pytest.approx(
            np.asarray(art["raw"]["f_curved"][:k]), rel=1e-12), path
        assert fresh_0 == pytest.approx(
            np.asarray(art["raw"]["f_flat"][:k]), rel=1e-12), path


@pytest.mark.slow
def test_the_preflight_effect_branch_reproduces_exactly():
    """The effect-branch certification (4000 replicates of the full
    pipeline at the frozen seed path) must reproduce its committed
    raw counts exactly."""

    art = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    p3e_art = json.loads(P3E_ART.read_text(encoding="utf-8"))
    rec = next(r for r in p3e_art["ladder_g"]
               if r["label"] == p3c.POINT[0])
    fa = np.asarray(rec["raw"]["f_curved"])
    f0 = np.asarray(rec["raw"]["f_flat"])
    fresh = p3c._branch_counts(fa, f0, p3c.N_ARM,
                               p3c.PREFLIGHT_EFFECT_REPS,
                               p3c.PREFLIGHT_EFFECT_SEED)
    assert fresh == art["preflight"]["effect"]["counts"]


@pytest.mark.slow
def test_the_preflight_null_branch_reproduces_exactly():
    """The case, PR #47 review R2: the null/equivalence branch is
    what actually DECIDES n = 4800, yet the first version reran only
    the effect branch -- a later change to `metric_cis` or
    `classify` could shift the null pass rate while the committed
    artifact stayed stale and every test stayed green (exactly the
    scenario of this PR's own CI-construction change). The full 20k
    replicates rerun here at the frozen seed path and must equal the
    committed counts integer for integer."""

    art = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    p3e_art = json.loads(P3E_ART.read_text(encoding="utf-8"))
    rec = next(r for r in p3e_art["ladder_g"]
               if r["label"] == p3c.POINT[0])
    f0 = np.asarray(rec["raw"]["f_flat"])
    fresh = p3c._branch_counts(f0, f0, p3c.N_ARM,
                               p3c.PREFLIGHT_NULL_REPS,
                               p3c.PREFLIGHT_NULL_SEED)
    assert fresh == art["preflight"]["null"]["counts"]
