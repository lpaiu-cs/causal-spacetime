"""Regressions for the P14 §8 P3-E probe machinery.

Same philosophy as the P1/P2 suites: pin the INSTRUMENTS (chain
enumeration, AUC, the fit machinery, the frozen seed layout and
margins, derived-from-raw arithmetic), not the measured numbers --
those live in the committed artifact with their seeds. The P3-C
preliminary power record is pinned two ways: its CP endpoints
recompute from the RAW integers on every run, and the full simulation
reproduces the raw pass count under the `slow` marker.
"""

from __future__ import annotations

import itertools
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

import p14_probe_p3e as p3  # noqa: E402
from p14_plane_wave import Slab, arms, sprinkle  # noqa: E402
from p14_probe_p1 import (  # noqa: E402
    clopper_pearson,
    element_eligible,
    relation_census,
)

ARTIFACT = (Path(__file__).resolve().parents[1]
            / "docs" / "prereg" / "p14_probe_p3e_results.json")


def test_the_seed_layout_and_burn_list_are_frozen():
    """P3-E component seeds and the burn list are review-frozen
    literals; the layout must check itself and never touch a burned
    seed (P1 campaign, P2 design + campaign, P3 pilots, independent
    review's 779/780)."""

    assert p3.SEEDS == {"ladder_e": 20260821, "low_a": 20260822,
                        "ladder_g": 20260823, "class_c": 20260824}
    for s in (20260808, 777, 778, 779, 780, 781,
              20260811, 20260812, 20260813, 20260814):
        assert s in p3.BURNED_SEEDS, s
    p3.assert_seed_layout()


def test_the_frozen_protocol_constants():
    """Ladder shapes and sizes from the five-round design review:
    9 rungs at n=400 with the three lowest carrying TOTAL 800; G at
    7x250; class-C descriptive at {1200, 2400} x 72 fixed; margins
    frozen at full precision and NEVER the inflated set (rejected in
    review: margins do not move to meet a power result)."""

    assert p3.W_LADDER == (0.125, 0.177, 0.25, 0.354, 0.5, 0.707,
                           1.0, 1.414, 2.0)
    assert p3.LOW_A_RUNGS == (0.125, 0.177, 0.25)
    assert p3.N_LADDER == 400 and p3.N_LOW_EXTRA == 400
    assert p3.N_G == 250
    assert p3.CLASSC_DENSITIES == (1200, 2400) and p3.N_CLASSC == 72
    assert p3.FIT_DR_ASYMPTOTIC_MAX_A == 0.0625
    assert p3.P3C_MARGINS["eps_s"] == pytest.approx(
        3.605 * math.sqrt(2.0 / 4000.0), rel=1e-15)
    assert p3.P3C_MARGINS["eps_auc"] == pytest.approx(
        3.605 * math.sqrt(8001.0 / (12.0 * 16_000_000.0)), rel=1e-12)
    assert p3.P3C_MARGINS["eps_ba"] == pytest.approx(
        3.605 / (2.0 * math.sqrt(4000.0)), rel=1e-15)
    assert p3.P3C_N_CANDIDATE == 4_800


def test_chain_and_opportunity_counts_match_brute_force():
    """The interval chain counts come from matrix powers of the
    interior relation; brute-force enumeration over element subsets
    must agree on a real (small) sprinkling, and the opportunity
    counts are exact binomials of the interior cardinalities."""

    slab = Slab(du=1.0, dv=1.0, dx=2.0, dy=6.0)
    curved, _flat = arms(slab, 1.0)
    rho = 120 / slab.coordinate_volume
    rng = np.random.default_rng(9)
    pts = sprinkle(curved, rho, rng)
    order = np.argsort(pts[:, 0], kind="stable")
    elig = np.array([element_eligible(curved, p) for p in pts[order]])
    cen = relation_census(curved, pts)
    got = p3._interval_chain_stats(cen, elig)

    rel = cen.related
    idx = np.where(elig)[0]
    brute = {k: 0 for k in (2, 3, 4)}
    opp = {k: 0 for k in (2, 3, 4)}
    n_iv = 0
    for i, j in itertools.combinations(idx, 2):
        if not rel[i, j]:
            continue
        n_iv += 1
        interior = [z for z in range(rel.shape[0])
                    if rel[i, z] and rel[z, j]]
        for k in (2, 3, 4):
            opp[k] += math.comb(len(interior), k)
            for combo in itertools.combinations(interior, k):
                if all(rel[a, b] for a, b in
                       itertools.combinations(combo, 2)):
                    brute[k] += 1
    assert got["n_intervals"] == n_iv
    assert got["chains"] == brute
    assert got["opportunities"] == opp


def test_ambiguous_pairs_never_count_as_definite_flips():
    """The case, PR #45 review R1: `related` leaves an undecided pair
    False, so a naive comparison counts a pair ambiguous in one arm
    and related in the other as a DEFINITE flip -- overstating
    §5.1's D_lower and double-counting the pair in the upper bound.
    Hand-built censuses: pair (0,1) related in flat, AMBIGUOUS in
    curved; pair (1,2) a genuine definite gain. The ambiguous pair
    must appear only in the ambiguity count -- lower bound 1 flip,
    upper bound 1 + 1, nothing twice."""

    from p14_probe_p1 import RelationCensus

    n = 3
    rel_flat = np.zeros((n, n), dtype=bool)
    rel_flat[0, 1] = True
    rel_curved = np.zeros((n, n), dtype=bool)
    rel_curved[1, 2] = True          # definite gain
    flat_c = RelationCensus(related=rel_flat, ambiguous=0, escalated=0,
                            seconds=0.0, ambiguous_pairs=())
    curved_c = RelationCensus(related=rel_curved, ambiguous=1,
                              escalated=0, seconds=0.0,
                              ambiguous_pairs=((0, 1),))
    gained, lost, amb = p3.pair_flips(flat_c, curved_c)
    assert (gained, lost, amb) == (1, 0, 1)
    # symmetric: ambiguous in the FLAT arm, related in curved -- the
    # would-be "gain" (1,2) is masked, leaving only the ambiguity
    flat_c2 = RelationCensus(related=np.zeros((n, n), dtype=bool),
                             ambiguous=1, escalated=0, seconds=0.0,
                             ambiguous_pairs=((1, 2),))
    curved_c2 = RelationCensus(related=rel_curved, ambiguous=0,
                               escalated=0, seconds=0.0,
                               ambiguous_pairs=())
    gained, lost, amb = p3.pair_flips(flat_c2, curved_c2)
    assert (gained, lost, amb) == (0, 0, 1)


def test_the_bootstrap_preserves_the_sprinkling_cluster():
    """The case, PR #45 review R1: the ladder shares point sets, so a
    replicate must apply ONE index vector per stratum to all of that
    stratum's rungs. Construct full-stratum rungs where rung B's rows
    equal rung A's rows plus a constant: under cluster resampling the
    resampled means must differ by EXACTLY that constant in every
    replicate; per-rung independent draws (the refuted first version)
    break the identity almost surely."""

    rungs_b = str(p3.W_LADDER[4])
    n = 40
    base = np.linspace(0.01, 0.05, n)
    raw = {"full": {"n": n, "rungs": {}}, "low": {"n": 4, "rungs": {}}}
    for w in p3.W_LADDER:
        key = str(w)
        vals = base if key != rungs_b else base + 0.5
        raw["full"]["rungs"][key] = [(float(v), 0.0, 0.0) for v in vals]
    for w in p3.LOW_A_RUNGS:
        raw["low"]["rungs"][str(w)] = [(0.02, 0.0, 0.0)] * 4
    rng = np.random.default_rng(5)
    for _ in range(25):
        m = p3._rung_means(raw, rng)
        assert (m[p3.W_LADDER[4]]["d"] - m[p3.W_LADDER[3]]["d"]
                == pytest.approx(0.5, abs=1e-12))


def test_the_auc_implementation_on_hand_values():
    """AUC = P(a > b) + 0.5 P(a = b): {1,2,3} vs {0,1,2} has six
    wins and two ties over nine pairs -> 7/9; degenerate all-equal
    arms give exactly 0.5."""

    a = np.array([1.0, 2.0, 3.0])
    b = np.array([0.0, 1.0, 2.0])
    assert p3._mann_whitney_auc(a, b) == pytest.approx(7.0 / 9.0)
    same = np.ones(5)
    assert p3._mann_whitney_auc(same, same.copy()) == pytest.approx(0.5)


def test_fit_machinery_recovers_exact_power_laws_and_records_failures():
    """A synthetic raw record with D exactly ~ A and dr exactly ~ A^2
    must fit slopes 1 and 2 to machine precision; a rung whose mean
    is non-positive is a FIT FAILURE, recorded, never discarded."""

    raw = {"full": {"rungs": {}, "swap": [], "n": 4},
           "low": {"rungs": {}, "n": 4}}
    for w in p3.W_LADDER:
        a = w * w
        g = (a + a * a) / 2.0
        lo = (a - a * a) / 2.0  # negative beyond A = 1: d=A, dr=A^2
        raw["full"]["rungs"][str(w)] = [(g, lo, 0.0)] * 4
        if w in p3.LOW_A_RUNGS:
            raw["low"]["rungs"][str(w)] = [(g, lo, 0.0)] * 4
    means = p3._rung_means(raw, None)
    assert p3._fit_slope(means, "d", 4.0) == pytest.approx(1.0, abs=1e-9)
    assert p3._fit_slope(means, "dr", 0.0625) == pytest.approx(
        2.0, abs=1e-9)
    # a non-positive rung mean in range -> None (a recorded failure)
    means[0.125] = {**means[0.125], "dr": 0.0}
    assert p3._fit_slope(means, "dr", 0.0625) is None


def test_the_p3c_record_is_raw_integers_with_recomputable_endpoints():
    """PR condition: the joint-power artifact stores the RAW pass
    count, reps, seed path, and full-precision margins; the CP
    endpoints must recompute exactly from the integers, the
    preliminary status must be stated, and the certification rule is
    the CP LOWER bound (never the point estimate)."""

    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    p = art["p3c_preliminary"]
    assert isinstance(p["passes_raw"], int)
    assert p["reps"] == p3.P3C_POWER_REPS
    assert p["n_arm"] == p3.P3C_N_CANDIDATE
    assert p["seed_path"] == [781, 4820]
    assert "PRELIMINARY" in p["status"]
    lo, hi = clopper_pearson(p["passes_raw"], p["reps"])
    assert p["joint_power_ci95_exact"] == pytest.approx([lo, hi],
                                                        rel=1e-12)
    assert p["joint_power"] == pytest.approx(
        p["passes_raw"] / p["reps"], rel=1e-15)
    assert p["margins"] == pytest.approx(p3.P3C_MARGINS, rel=1e-15)
    assert lo >= 0.90, "preliminary candidate lost its certification"


def test_ladder_g_derived_fields_recompute_from_the_raw_samples():
    """Every ladder-G aggregate recomputes from the stored
    per-sprinkling 3a samples and flip sums -- the PR #41/#42 lesson:
    a partially pinned artifact drifts where it is not pinned. The
    stored f arrays are also what P3-C's power-model refit consumes,
    so their consistency with the aggregates is load-bearing."""

    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert len(art["ladder_g"]) == len(p3.LADDER_G_POINTS)
    for rec in art["ladder_g"]:
        n = rec["n"]
        fa = np.asarray(rec["raw"]["f_curved"])
        f0 = np.asarray(rec["raw"]["f_flat"])
        assert len(fa) == len(f0) == n
        df = fa - f0
        assert rec["mean_f_curved"] == pytest.approx(fa.mean(), rel=1e-12)
        assert rec["mean_f_flat"] == pytest.approx(f0.mean(), rel=1e-12)
        assert rec["delta_f"] == pytest.approx(df.mean(), rel=1e-9)
        assert rec["se_delta_f"] == pytest.approx(
            df.std(ddof=1) / math.sqrt(n), rel=1e-9)
        sd_arm = math.sqrt(0.5 * (fa.var(ddof=1) + f0.var(ddof=1)))
        assert rec["s_realized"] == pytest.approx(df.mean() / sd_arm,
                                                 rel=1e-9)
        # the realized-precision formula, not the near-zero one
        assert rec["se_s_realized"] == pytest.approx(
            math.sqrt(2.0 / n + rec["s_realized"] ** 2 / (4 * (n - 1))),
            rel=1e-9)
        d_mean = (rec["raw"]["sum_g"] + rec["raw"]["sum_l"]) / n
        assert rec["mean_d"] == pytest.approx(d_mean, rel=1e-12)
        assert rec["d_interval"][0] == pytest.approx(d_mean, rel=1e-12)
        assert rec["d_interval"][1] == pytest.approx(
            d_mean + rec["raw"]["sum_amb"] / n, rel=1e-9)
        assert rec["auc"] == pytest.approx(
            p3._mann_whitney_auc(fa, f0), rel=1e-12)


def test_ladder_e_point_fits_recompute_from_the_raw_record():
    """The artifact's fit block must equal a recomputation from the
    committed raw rung records: point slopes exactly, the swap
    contrast exactly, and the rung means themselves (the bootstrap
    CIs are seeded and are pinned by the slow rerun of the campaign
    stage, not here -- recomputing 4000 resamples in the fast suite
    is wasteful)."""

    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    raw = art["ladder_e_raw"]
    means = p3._rung_means(raw, None)
    fits = art["ladder_e_fits"]
    for w in p3.W_LADDER:
        m = fits["rungs"][str(w)]
        assert m["d"] == pytest.approx(means[w]["d"], rel=1e-12)
        assert m["dr"] == pytest.approx(means[w]["dr"], rel=1e-9)
        want_n = p3.N_LADDER + (p3.N_LOW_EXTRA
                                if w in p3.LOW_A_RUNGS else 0)
        assert m["n"] == want_n, f"total n per rung, w={w}"
    for name, qty, max_a in p3._FITS:
        assert fits[name]["slope"] == pytest.approx(
            p3._fit_slope(means, qty, max_a), rel=1e-9)
    sw = np.asarray(raw["full"]["swap"])
    plus = np.asarray(raw["full"]["rungs"][str(p3.W_SWAP)])[:len(sw)]
    d_p = plus[:, 0] + plus[:, 1]
    d_m = sw[:, 0] + sw[:, 1]
    assert fits["swap_contrast"]["paired_diff"] == pytest.approx(
        float((d_p - d_m).mean()), rel=1e-9)


def test_the_results_doc_embeds_the_rendered_tables():
    """The doc's three tables are RENDERED from the artifact and
    embedded verbatim (the established structure since PR #41)."""

    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    doc = (ARTIFACT.parent / "p14_probe_p3e_results.md").read_text(
        encoding="utf-8")
    assert p3.ladder_e_table(art) in doc
    assert p3.ladder_g_table(art) in doc
    assert p3.class_c_table(art) in doc


@pytest.mark.slow
def test_the_full_fit_block_reproduces_from_the_committed_raw():
    """The bootstrap CIs are seeded, so the ENTIRE fits block --
    point slopes, cluster-bootstrap CIs, failure counts, refit
    diagnostics, swap contrast -- must reproduce from the committed
    raw record (the partial-pinning lesson: unpinned CI fields are
    hand-editable while the fast tests stay green)."""

    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    fresh = p3.ladder_e_fits(art["ladder_e_raw"])
    committed = art["ladder_e_fits"]
    for name, *_ in p3._FITS:
        f, c = fresh[name], committed[name]
        assert f["slope"] == pytest.approx(c["slope"], rel=1e-12), name
        assert f["ci95"] == pytest.approx(c["ci95"], rel=1e-9), name
        assert f["fit_failures"] == c["fit_failures"], name
        assert f["refit_diagnostic_drop_highest"] == pytest.approx(
            c["refit_diagnostic_drop_highest"], rel=1e-12), name
    assert fresh["swap_contrast"] == pytest.approx(
        committed["swap_contrast"], rel=1e-9)


@pytest.mark.slow
def test_the_p3c_power_simulation_reproduces_its_raw_count():
    """Full pipeline (threshold and direction re-learned per
    replicate) at the frozen seed path must reproduce the committed
    raw pass count exactly -- the one field the fast recomputation
    cannot reach from the stored integers."""

    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    fresh = p3.p3c_joint_power()
    assert fresh["passes_raw"] == art["p3c_preliminary"]["passes_raw"]
    assert fresh["joint_power_ci95_exact"] == pytest.approx(
        art["p3c_preliminary"]["joint_power_ci95_exact"], rel=1e-12)
