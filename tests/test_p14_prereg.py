"""Regressions for the P14 preregistration stage (preparation commit).

Pins the frozen constants and sentences, both kinds of seed checks
(scalar ledger and SeedSequence structural uniqueness), the C1
pipeline on hand data, the mechanical execution-order gates, and the
code_version/ancestry contract. Artifact-recompute tests skip until
the preflight and campaign artifacts exist (doc §8: they land in
LATER commits than this code).
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

import p12_stage_b  # noqa: E402
import p14_prereg as pr  # noqa: E402
import p14_probe_p3c as p3c  # noqa: E402
from p14_probe_p2 import student_t_crit  # noqa: E402
from seed_windows import (  # noqa: E402
    assert_point_seeds_fresh,
    assert_windows_disjoint_and_fresh,
)

REPO = Path(__file__).resolve().parents[1]
PREFLIGHT = REPO / "docs" / "prereg" / "p14_prereg_preflight.json"
RESULTS = REPO / "docs" / "prereg" / "p14_prereg_results.json"
DOC = REPO / "docs" / "prereg" / "p14_prereg.md"


def test_the_frozen_constants_and_sentences():
    """Sizes, margin, replicate counts, operating-point identity with
    the probe chain, and the frozen reporting sentences (C2's may
    never weaken the standing P3-C record)."""

    assert pr.N_C1 == 3_000
    assert pr.N_C2 == p3c.N_ARM == 4_800
    assert pr.POINT is p3c.POINT and pr.E_N == p3c.E_N
    assert pr.EPS_DELTA == 0.0003578762
    es, sd = pr.EPS_DELTA_ANCHOR
    assert es * sd == pytest.approx(pr.EPS_DELTA, rel=5e-7)
    assert pr.B_JOINT == 4_000
    assert pr.B_C1_NULL == 20_000 and pr.B_C2_NULL == 20_000
    assert pr.BOOT_CHILDREN == {"joint_effect": 3, "c1_null": 1,
                                "c2_null": 2}
    assert pr.C2_SENTENCES["confirmed"] == "P3-C 분리를 독립적으로 재현했다."
    assert "병기한다" in pr.C2_SENTENCES["inconclusive"]
    assert "충돌한다" in pr.C2_SENTENCES["equivalent"]
    for s in pr.C2_SENTENCES.values():
        assert "미확립" not in s  # the v0.2 phrasing the review killed


def test_the_scalar_seed_ledger_catches_the_p13_overlap():
    """The R-5 finding, pinned: the v0.3-approved 20260861/71/72 sit
    inside P13 campaign v3's LEDGER range [15000000..21503999], so
    the mechanical check must reject them. What that proves is a
    LEDGER overlap -- freshness cannot be certified -- not an actual
    integer stream reuse (the review's reconstruction of P13 v3's
    58880 possible RNG seeds intersects the date-styled P14 seeds
    nowhere), so the probe record is not retroactively downgraded.
    The stage seeds moved to 40000000+, past P12's full allocation
    decade, and the check treats that decade [30000000..39999999] as
    reserved rather than only P12's used ceiling."""

    pr.assert_seed_layout()  # the frozen 40000061/71/72 layout passes
    assert pr.CAMPAIGN_SEEDS == {"c1_paired": 40_000_061,
                                 "c2_curved": 40_000_071,
                                 "c2_flat": 40_000_072}
    with pytest.raises(SystemExit, match="P13 campaign v3"):
        assert_point_seeds_fresh(
            {"c1_paired": 20260861}, pr._SPENT_SCALARS,
            pr._SPENT_RANGES, "P14 prereg")
    # the documented allocation boundary binds, not the used ceiling:
    # 34000061 clears every used range yet is still reserved ground
    with pytest.raises(SystemExit, match="P12 allocation decade"):
        assert_point_seeds_fresh(
            {"c1_paired": 34_000_061}, pr._SPENT_SCALARS,
            pr._SPENT_RANGES, "P14 prereg")
    # the probe chain's own consumed streams are spent scalars
    assert 20260851 in pr._SPENT_SCALARS
    assert 20260852 in pr._SPENT_SCALARS
    with pytest.raises(SystemExit, match="already spent"):
        assert_point_seeds_fresh(
            {"c2_curved": 20260851}, pr._SPENT_SCALARS,
            (), "P14 prereg")


def test_the_extracted_helper_keeps_p12_behavior():
    """The shared pure function inherits P12's exact semantics: the
    P12 wrapper still passes on its own constants, and the pure
    function aborts on window escape, spent overlap, and pairwise
    overlap."""

    blocks = p12_stage_b.assert_windows_disjoint_and_fresh()
    assert blocks and all(len(b) == 3 for b in blocks)
    with pytest.raises(SystemExit, match="leaves T's window"):
        assert_windows_disjoint_and_fresh(
            [("a", 5, 15)], (), 10, 20, "T")
    with pytest.raises(SystemExit, match="already-spent old"):
        assert_windows_disjoint_and_fresh(
            [("a", 12, 14)], (("old", 13, 30),), 10, 20, "T")
    with pytest.raises(SystemExit, match="blocks overlap"):
        assert_windows_disjoint_and_fresh(
            [("a", 10, 14), ("b", 14, 18)], (), 10, 20, "T")


def test_the_bootstrap_roots_are_entropy_vectors_not_spawn_keys():
    """The v0.4 semantic freeze: `[781, 60]` is an ENTROPY VECTOR in
    the probe chain's convention, and the spawn-key construction is
    a DIFFERENT stream -- the review verified the initial states
    differ, pinned here by first draws."""

    a = int(np.random.default_rng([781, 60]).integers(0, 2 ** 63))
    b = int(np.random.default_rng(
        np.random.SeedSequence(781, spawn_key=(60,))).integers(0, 2 ** 63))
    assert a == 6625364951367084896
    assert b == 1519964771348177218
    assert a != b


def test_the_spawn_layout_identities_are_exhaustively_unique():
    """Collision freedom is STRUCTURAL KEY UNIQUENESS, enumerated in
    full (prereg v0.4): joint 4000*3 = 12000, C1 null 20000, C2 null
    20000*2 = 40000 identities all distinct, and the three new roots
    distinct from every spent root."""

    seen: set[tuple] = set()
    for branch, b in (("joint_effect", pr.B_JOINT),
                      ("c1_null", pr.B_C1_NULL),
                      ("c2_null", pr.B_C2_NULL)):
        n_children = pr.BOOT_CHILDREN[branch]
        ids = []
        for kids in pr.boot_layout(branch, b):
            assert len(kids) == n_children
            for ss in kids:
                ids.append((tuple(ss.entropy), ss.spawn_key))
        assert len(ids) == b * n_children
        assert len(set(ids)) == len(ids), branch
        seen |= set(ids)
    assert len(seen) == 12_000 + 20_000 + 40_000  # cross-branch too
    roots = {tuple(v) for v in pr.BOOT_ROOTS.values()}
    assert len(roots) == 3
    assert not roots & set(pr.SPENT_BOOT_ROOTS)


def test_the_first_replicate_streams_reproduce():
    """State pin (reproducibility, NOT the collision proof): the first
    replicate's child streams generate the recorded first draws."""

    pins = {"joint_effect": [8383312881007714273,
                             1481462205926588177,
                             8545546324759902268],
            "c1_null": [8755339445751806349],
            "c2_null": [2458872342616347232, 1834739337975821059]}
    for branch, want in pins.items():
        kids = pr.boot_layout(branch, 1)[0]
        got = [int(np.random.default_rng(ss).integers(0, 2 ** 63))
               for ss in kids]
        assert got == want, branch


def test_c1_metrics_and_classify_on_hand_data():
    """The paired pipeline on tiny hand data: Student-t critical
    values against published tables (t_{0.975,3} = 3.182446,
    t_{0.975,4} = 2.776445), the CI arithmetic, and the three
    branches with the frozen direction (fully below the band is
    never confirmed)."""

    assert student_t_crit(3) == pytest.approx(3.182446, abs=5e-7)
    assert student_t_crit(4) == pytest.approx(2.776445, abs=5e-7)
    delta = np.array([1.0, 2.0, 3.0, 4.0])
    m = pr.c1_metrics(delta)
    se = math.sqrt(5.0 / 3.0) / 2.0
    assert m["mean"] == pytest.approx(2.5, rel=1e-12)
    assert m["ci"][1] - m["ci"][0] == pytest.approx(
        2 * 3.1824463052837046 * se, rel=1e-9)
    e = pr.EPS_DELTA
    assert pr.c1_classify((2 * e, 3 * e)) == "confirmed"
    assert pr.c1_classify((-e / 2, e / 2)) == "equivalent"
    assert pr.c1_classify((e / 2, 2 * e)) == "inconclusive"
    assert pr.c1_classify((-3 * e, -2 * e)) == "inconclusive"


def test_paired_samples_reads_the_same_points_twice():
    """A tiny paired smoke at low density: delta is exactly
    f_curved - f_flat per sprinkling, ambiguity totals are recorded
    per reading, and the stream is deterministic."""

    from p14_plane_wave import Slab, arms
    label, w, du, dv, dx, dy = pr.POINT
    slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
    rho = 30 / slab.coordinate_volume  # tiny E[N] for speed
    curved, flat = arms(slab, w)
    d1, fa, f0, amb = pr.paired_samples(curved, flat, rho, 12345, 3)
    assert d1 == pytest.approx(fa - f0, rel=1e-15)
    assert set(amb) == {"curved", "flat"}
    assert all(set(v) == {"ambiguous", "escalated"}
               for v in amb.values())
    d2, _, _, _ = pr.paired_samples(curved, flat, rho, 12345, 3)
    assert d1 == pytest.approx(d2, rel=1e-15)


def test_the_campaign_gates_are_mechanical(tmp_path, monkeypatch):
    """Doc §8's order is enforced by three gates checked BEFORE any
    computation, exercised here against temporary paths so the tests
    survive the real manifest landing at the final freeze (PR #49
    review P2): (1) no manifest -> refuse; (2) manifest/artifact
    digest mismatch -> refuse; (3) execution sources changed since
    certification -> refuse (PR #49 review P1: without this, code
    edited between P and F runs an UNCERTIFIED pipeline while every
    other check passes)."""

    monkeypatch.setattr(pr, "_FREEZE_MANIFEST", tmp_path / "none.json")
    with pytest.raises(SystemExit, match="no freeze manifest"):
        pr.run_campaign()

    pre_path = tmp_path / "preflight.json"
    man_path = tmp_path / "freeze.json"
    monkeypatch.setattr(pr, "_PREFLIGHT_ARTIFACT", pre_path)
    monkeypatch.setattr(pr, "_FREEZE_MANIFEST", man_path)
    pre_path.write_text(json.dumps(
        {"certified": True, "source_digests": pr._source_digests()}),
        encoding="utf-8")
    man_path.write_text(json.dumps({"preflight_digest": "not-it"}),
                        encoding="utf-8")
    with pytest.raises(SystemExit, match="does not match the freeze"):
        pr.run_campaign()

    stale = dict(pr._source_digests())
    stale["p14_probe_p3c.py"] = "0" * 64
    pre_path.write_text(json.dumps(
        {"certified": True, "source_digests": stale}), encoding="utf-8")
    man_path.write_text(json.dumps(
        {"preflight_digest": pr._sha256(pre_path)}), encoding="utf-8")
    with pytest.raises(SystemExit,
                       match=r"sources changed.*p14_probe_p3c\.py"):
        pr.run_campaign()

    # and the digest table covers exactly the frozen source list
    assert set(pr._source_digests()) == set(pr._FROZEN_SOURCES)


def test_the_doc_carries_the_frozen_numbers():
    """The design document and the runner cannot diverge on the
    numbers a reader would copy."""

    doc = DOC.read_text(encoding="utf-8")
    for token in ("0.0003578762", "40000061", "40000071", "40000072",
                  "[781, 60]", "n_C1 = **3000**", "4,000", "20,000",
                  "1025c50"):
        assert token in doc, token


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True,
                          text=True, check=True).stdout.strip()


def test_the_code_version_contract():
    """`code_version` returns HEAD only on a clean checkout; a dirty
    worktree aborts (the recorded SHA must not lie about what ran)."""

    probe = REPO / "._cv_dirty_probe"
    if _git("status", "--porcelain"):
        pytest.skip("worktree already dirty; contract not testable here")
    assert pr.code_version() == _git("rev-parse", "HEAD")
    probe.write_text("x", encoding="utf-8")
    try:
        with pytest.raises(SystemExit, match="CLEAN checkout"):
            pr.code_version()
    finally:
        probe.unlink()


def test_the_ancestry_contract_once_artifacts_exist():
    """P (preflight's code_version) must be an ancestor of HEAD, and
    F (results') must descend from P -- constant identity alone
    cannot prove which commit executed (prereg v0.4). Skips until
    the artifacts land in their later commits."""

    if not PREFLIGHT.exists():
        pytest.skip("preflight artifact not yet produced -- it lands "
                    "in a commit after the preparation commit (doc §8)")
    p_sha = json.loads(PREFLIGHT.read_text(encoding="utf-8"))[
        "code_version"]
    assert subprocess.run(
        ["git", "merge-base", "--is-ancestor", p_sha, "HEAD"],
        cwd=REPO).returncode == 0
    if not RESULTS.exists():
        pytest.skip("campaign artifact not yet produced -- it lands "
                    "after the final freeze (doc §8)")
    f_sha = json.loads(RESULTS.read_text(encoding="utf-8"))[
        "code_version"]
    assert subprocess.run(
        ["git", "merge-base", "--is-ancestor", p_sha, f_sha],
        cwd=REPO).returncode == 0
    assert subprocess.run(
        ["git", "merge-base", "--is-ancestor", f_sha, "HEAD"],
        cwd=REPO).returncode == 0
