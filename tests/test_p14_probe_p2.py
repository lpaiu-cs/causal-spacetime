"""Regressions for the P14 §8 P2 probe machinery.

Same philosophy as the P1 suite: pin the INSTRUMENTS (the exact
distributions, the counting predicate, the seed layout, the derived-
from-raw arithmetic), not the measured numbers -- those live in the
committed results artifact with their seeds, and the artifact itself
is pinned by recomputation: fast fields on every run, the full
edge-a2.4 campaign under the `slow` marker.
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

import p14_probe_p2 as p2  # noqa: E402
from p14_plane_wave import Slab, arms  # noqa: E402
from p14_probe_p1 import (  # noqa: E402
    axis_volume_ratio,
    diamond_volumes_mc,
    fattest_axis_diamond,
)

ARTIFACT = (Path(__file__).resolve().parents[1]
            / "docs" / "prereg" / "p14_probe_p2_results.json")


def test_the_seed_layout_is_the_frozen_one():
    """The design review froze the consumption unit (one rng stream
    per operating point) and the literals; a silent edit to any seed
    must fail here, and the layout must stay collision-free and
    disjoint from the burned list (P1's campaign seed and the design
    checks' 777)."""

    assert p2.CAMPAIGN_SEEDS == {"slice-a1.0": 20260811,
                                 "high-a2.0": 20260812,
                                 "edge-a2.4": 20260813}
    assert p2.MC_SEED == 20260814
    assert 20260808 in p2.BURNED_SEEDS and 777 in p2.BURNED_SEEDS
    p2.assert_seed_layout()  # and the layout checks itself


def test_the_frozen_protocol_constants():
    """n = 21000/400/100 (primary TOST at the V_dis UCB, marginal
    tau_m = delta binding at the first two, variance floor at the
    third), E[N_box] = 30000, mc = 2e6 -- the design review's numbers,
    not tunables."""

    ns = {lbl: n for lbl, *_, n in p2.OPERATING_POINTS}
    assert ns == {"slice-a1.0": 21_000, "high-a2.0": 400,
                  "edge-a2.4": 100}
    assert p2.N_BOX == 30_000
    assert p2.MC_SAMPLES == 2_000_000


def test_student_t_criticals_match_published_tables():
    """The frozen primary CI is Student-t, 95% two-sided; the repo has
    no scipy, so the quantile is built from the incomplete beta and
    must match the published two-sided 97.5% points."""

    for df, want in ((1, 12.7062), (10, 2.2281), (99, 1.9842),
                     (399, 1.9659)):
        assert p2.student_t_crit(df) == pytest.approx(want, abs=6e-4)
    # large-df limit approaches the normal critical from above
    t_big = p2.student_t_crit(20_999)
    assert 1.9599 < t_big < 1.9605
    assert p2.student_t_crit(99) > p2.student_t_crit(399) > t_big


def test_garwood_poisson_ci_matches_published_values():
    """The frozen marginal CI is the exact Garwood interval. Standard
    values: zero counts -> [0, 3.6889]; ten counts ->
    [4.7954, 18.3904]. And the defining property: the interval's
    endpoints put exactly 2.5% of the fitted Poisson mass beyond the
    observation."""

    lo, hi = p2.poisson_mean_ci(0)
    assert lo == 0.0 and hi == pytest.approx(3.6889, abs=1e-3)
    lo, hi = p2.poisson_mean_ci(10)
    assert lo == pytest.approx(4.7954, abs=1e-3)
    assert hi == pytest.approx(18.3904, abs=1e-3)
    # defining property via the same regularized gamma the CI uses:
    # P(T) at the endpoints is the tail mass by the gamma-Poisson link
    assert p2._gamma_p(10.0, lo) == pytest.approx(0.025, abs=1e-9)
    assert p2._gamma_p(11.0, hi) == pytest.approx(0.975, abs=1e-9)


def test_count_diamond_excludes_the_anchors_by_construction():
    """The anchors are external and not counted (§8 P2, the R5
    approval condition). The count's `u` window is strict, so a point
    AT an anchor's coordinates is never counted, while the diamond's
    own midpoint is."""

    slab = Slab(du=2.4, dv=0.2, dx=1.2, dy=0.8)
    curved, flat = arms(slab, 1.0)
    p, q = fattest_axis_diamond(curved)
    mid = 0.5 * (np.asarray(p) + np.asarray(q))
    pts = np.vstack([np.asarray(p), np.asarray(q), mid])
    assert p2.count_diamond(curved, p, q, pts) == 1
    assert p2.count_diamond(flat, p, q, pts) == 1


def test_the_p1_sizing_link_holds():
    """P2's sizing inputs are P1's committed artifact; the quadrature
    the two scripts share must agree exactly at every operating
    point."""

    p2.check_p1_link()


#: Every key a point record carries. Adding a field to `run_point`
#: without extending the recomputation below fails here first -- the
#: gap PR #42's review found (the Bonferroni fields were added after
#: the test and went unpinned) cannot recur silently.
_POINT_KEYS = frozenset({
    "label", "w", "du", "dv", "dx", "dy", "n", "seed",
    "a", "r_pred", "delta", "rho", "lam0", "lam_a", "tau", "tau_m",
    "raw", "zbar", "var_emp", "se", "t_crit",
    "theta_ci", "equivalent", "discriminates",
    "theta_ci_bonf", "equivalent_bonf", "discriminates_bonf",
    "marginal", "var_diag", "seconds",
})


@pytest.mark.slow
def test_the_results_artifact_derived_fields_recompute_from_raw():
    """The PR #41 lesson applied to P2 from the start: every derived
    field in the committed artifact must recompute from the raw
    sufficient statistics and the frozen protocol -- zbar and Var from
    the sums, the t-CI from those, the Bonferroni CI and all four
    labels, the marginal CIs from pooled totals that are THEMSELVES
    tied back to the raw sums, and the analytic side of the variance
    diagnostic from a fresh MC at the recorded sub-seed (PR #42
    review R1: the first version skipped the Bonferroni fields --
    exactly the ones that publish the joint sentence -- and never tied
    the marginal totals to the raw record). The one exception is the
    bootstrap CI, which needs the full per-sprinkling sample; it is
    pinned by the slow edge-a2.4 rerun. `_POINT_KEYS` makes the
    coverage a schema contract: a new field fails here until pinned.
    """

    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert len(art["points"]) == len(p2.OPERATING_POINTS)
    for k, rec in enumerate(art["points"]):
        assert set(rec) == _POINT_KEYS, rec.get("label")
        label, w, du, dv, dx, dy, n = p2.OPERATING_POINTS[k]
        assert rec["label"] == label and rec["n"] == n
        assert rec["seed"] == p2.CAMPAIGN_SEEDS[label]
        slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
        curved, flat = arms(slab, w)
        p, q = fattest_axis_diamond(curved)
        r_pred = axis_volume_ratio(w * du)
        assert rec["r_pred"] == pytest.approx(r_pred, rel=1e-15)
        delta = r_pred - 1.0
        assert rec["tau"] == pytest.approx(delta, rel=1e-15)
        assert rec["tau_m"] == pytest.approx(delta, rel=1e-15)
        rho = p2.N_BOX / slab.coordinate_volume
        assert rec["rho"] == pytest.approx(rho, rel=1e-15)
        dv_win = float(p[1]) - float(q[1])
        lam0 = rho * math.pi * du ** 2 * dv_win ** 2 / 6.0
        assert rec["lam0"] == pytest.approx(lam0, rel=1e-12)
        assert rec["lam_a"] == pytest.approx(r_pred * lam0, rel=1e-12)

        raw = rec["raw"]
        zbar = (raw["sum_na"] - r_pred * raw["sum_n0"]) / n
        assert rec["zbar"] == pytest.approx(zbar, rel=1e-12)
        var_emp = (raw["sum_z2"] - n * zbar * zbar) / (n - 1)
        assert rec["var_emp"] == pytest.approx(var_emp, rel=1e-12)
        se = math.sqrt(var_emp / n)
        assert rec["se"] == pytest.approx(se, rel=1e-12)
        t_crit = p2.student_t_crit(n - 1)
        assert rec["t_crit"] == pytest.approx(t_crit, rel=1e-12)
        lo = (zbar - t_crit * se) / rec["lam0"]
        hi = (zbar + t_crit * se) / rec["lam0"]
        assert rec["theta_ci"] == pytest.approx([lo, hi], rel=1e-9)
        assert rec["equivalent"] == (-delta < lo and hi < delta)
        assert rec["discriminates"] == (lo > -delta or hi < -delta)
        # the Bonferroni fields publish the JOINT sentence via
        # primary_table -- they recompute from the same raw record
        t_bonf = p2.student_t_crit(n - 1, p2._T_LEVEL_BONF3)
        b_lo = (zbar - t_bonf * se) / rec["lam0"]
        b_hi = (zbar + t_bonf * se) / rec["lam0"]
        assert rec["theta_ci_bonf"] == pytest.approx([b_lo, b_hi],
                                                     rel=1e-9)
        assert rec["equivalent_bonf"] == (-delta < b_lo
                                          and b_hi < delta)
        assert rec["discriminates_bonf"] == (b_lo > -delta
                                             or b_hi < -delta)

        # the marginal totals ARE the raw sums -- one record, not two
        assert rec["marginal"]["curved"]["total"] == raw["sum_na"]
        assert rec["marginal"]["flat"]["total"] == raw["sum_n0"]
        for arm_name, lam in (("curved", rec["lam_a"]),
                              ("flat", rec["lam0"])):
            m = rec["marginal"][arm_name]
            g_lo, g_hi = p2.poisson_mean_ci(m["total"])
            want = [g_lo / (n * lam), g_hi / (n * lam)]
            assert m["ratio_ci"] == pytest.approx(want, rel=1e-9)
            assert m["equivalent"] == (
                1.0 - delta < want[0] and want[1] < 1.0 + delta)

        d = rec["var_diag"]
        mc_rng = np.random.default_rng([p2.MC_SEED, k])
        vols = diamond_volumes_mc(curved, flat, p, q,
                                  d["mc_samples"], mc_rng)
        assert d["v_dis_hits"] == vols.disagree_hits
        floor = delta * lam0 / rho
        want_lo = r_pred * rho * max(vols.disagree_lcb, floor)
        want_hi = r_pred * rho * max(vols.disagree_ucb, floor)
        assert d["analytic_ci"] == pytest.approx(
            [want_lo, want_hi], rel=1e-9)


@pytest.mark.slow
def test_the_frozen_n_have_executable_sizing_provenance():
    """The cases, PR #42 reviews R1 and R2. R1: minima quoted only in
    prose are not a frozen design -- `sizing_block` is executable
    provenance. R2 sharpened it twice: (a) the DESIGN provenance is
    the pre-campaign 2e6 MC at the burned design seed (hits
    120/79/737, primary minima 1167/78/6), reproduced here
    deterministically -- the campaign's own MC is carried only as an
    outcome diagnostic, since it can merely show the realized data
    would also have justified `n`; (b) the discrete Garwood power is
    NOT monotone in `n` (it dips back below target at slice's
    n = 20941), so the certification is the DIRECT power evaluation
    at the frozen `n`, never a crossing point -- the first crossing
    stays informational.

    Floats are compared under a documented tolerance: the power's
    T-range endpoints are integers found by comparing libm-computed
    gamma quantiles against a boundary, so a cross-platform ulp can
    shift an endpoint by one and move the power by ~one Poisson pmf
    term (~4e-4 at slice) -- abs 2e-3 covers it, while any
    meaningful edit (1e-2+) still fails; the same flip moves the
    first-crossing locator by a few integers, hence abs <= 3.
    """

    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert p2.POWER_TARGET == 0.90
    assert p2.DESIGN_MC_SEED == 777
    fresh = p2.sizing_block(art)
    blk = art["sizing"]
    assert blk["power_target"] == fresh["power_target"]
    assert blk["design_mc"] == fresh["design_mc"]
    for c_rec, f_rec in zip(blk["points"], fresh["points"],
                            strict=True):
        assert set(c_rec) == set(f_rec)
        for key in ("label", "n_frozen", "design_mc_hits",
                    "n_primary_min_design", "n_primary_min_campaign",
                    "variance_floor"):
            assert c_rec[key] == f_rec[key], f"{c_rec['label']}.{key}"
        assert c_rec["marginal_power_at_frozen_n"] == pytest.approx(
            f_rec["marginal_power_at_frozen_n"], abs=2e-3)
        assert abs(c_rec["n_marginal_first_crossing"]
                   - f_rec["n_marginal_first_crossing"]) <= 3

    by_label = {r["label"]: r for r in fresh["points"]}
    # the design-review numbers, independently cross-computed there
    assert [by_label[k]["design_mc_hits"] for k in
            ("slice-a1.0", "high-a2.0", "edge-a2.4")] == [120, 79, 737]
    assert [by_label[k]["n_primary_min_design"] for k in
            ("slice-a1.0", "high-a2.0", "edge-a2.4")] == [1167, 78, 6]
    assert by_label["slice-a1.0"]["n_marginal_first_crossing"] == 20_939
    # the certification itself: power AT the frozen n, direct
    for r in fresh["points"]:
        assert r["marginal_power_at_frozen_n"] >= p2.POWER_TARGET
        assert r["n_frozen"] >= max(r["n_primary_min_design"],
                                    r["variance_floor"]), r["label"]


@pytest.mark.slow
def test_the_edge_point_reruns_to_the_committed_record_exactly():
    """One full end-to-end pin: rerunning edge-a2.4's whole campaign
    (n = 100, the cheap point) from its frozen seed must reproduce the
    committed record exactly -- raw sums, bootstrap CI, labels, every
    field. This is the one test that also pins the bootstrap, which
    the fast recomputation cannot reach from the raw sums alone."""

    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    rec = art["points"][2]
    assert rec["label"] == "edge-a2.4"
    fresh = p2.run_point("edge-a2.4", 1.0, 2.4, 0.2, 1.2, 0.8, 100, 2)

    def deep_equal(a, b, path):
        if isinstance(b, dict):
            assert set(a) == set(b), path
            for kk in b:
                deep_equal(a[kk], b[kk], f"{path}.{kk}")
        elif isinstance(b, list):
            assert len(a) == len(b), path
            for i, (x, y) in enumerate(zip(a, b, strict=True)):
                deep_equal(x, y, f"{path}[{i}]")
        elif isinstance(b, float) and not isinstance(b, bool):
            assert a == pytest.approx(b, rel=1e-12), path
        else:
            assert a == b, path

    for key, cv in rec.items():
        if key == "seconds":
            continue
        deep_equal(fresh[key], cv, key)


def test_the_results_doc_embeds_the_rendered_tables():
    """The doc's three result tables are RENDERED from the artifact
    and embedded verbatim (the PR #41 structure): primary, marginal,
    and diagnostic -- every cell pinned to the computation."""

    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    doc = (ARTIFACT.parent / "p14_probe_p2_results.md").read_text(
        encoding="utf-8")
    assert p2.primary_table(art) in doc
    assert p2.marginal_table(art) in doc
    assert p2.diagnostic_table(art) in doc
