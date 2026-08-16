"""O3' freeze contract tests.

Same discipline as the O3 freeze (rules and caps pinned before any
execution; the run separately approved from a clean exact checkout),
plus what is new in O3': the immutable-O3 intersection base, the
certification-inconsistency rule, the V_ref' recommendation, and the
two-curve cost-projection range the caps are headroom over."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))

import o3p_frozen_volume as o3p  # noqa: E402
import o3p_projection as proj  # noqa: E402


def test_manifest_digests_match_the_working_tree():
    """Every pinned file matches its frozen digest, the runner pins
    ITSELF, and the pinned set adds the O3' surface AND the immutable
    O3 result the intersection reads."""

    manifest = o3p.verify_digests("test")
    assert set(manifest["files"]) == {
        "experiments/oracle/certified_interval.py",
        "experiments/oracle/certified_flight_time.py",
        "experiments/oracle/volume_oracle.py",
        "experiments/oracle/o3p_frozen_volume.py",
        "experiments/oracle/o3p_projection.py",
        "docs/theory/schwarzschild_volume_oracle_certification.md",
        "docs/prereg/p14_oracle_price.json",
        "docs/prereg/p14_oracle_mode_width.json",
        "docs/prereg/p14_o3_volume.json",
        "docs/prereg/p14_o3p_frozen_volume.md",
        "pyproject.toml",
    }


def test_manifest_locks_the_rounding_instrument():
    manifest = json.loads(o3p._MANIFEST.read_text(encoding="utf-8"))
    env = manifest["environment"]
    assert set(env) == {"python", "gmpy2", "mpfr", "gmp"}
    assert env["mpfr"].startswith("MPFR ")
    assert env["gmp"].startswith("GMP ")


def test_environment_verification_refuses_drift(monkeypatch):
    manifest = json.loads(o3p._MANIFEST.read_text(encoding="utf-8"))
    monkeypatch.setattr(o3p, "_environment",
                        lambda: dict(manifest["environment"]))
    o3p.verify_environment("test", manifest)
    drifted = dict(manifest["environment"], mpfr="MPFR 0.0.0")
    monkeypatch.setattr(o3p, "_environment", lambda: drifted)
    with pytest.raises(SystemExit, match="environment drift"):
        o3p.verify_environment("test", manifest)


def test_frozen_configuration_is_the_ruled_one():
    """The PI ruling's exact freeze: the adopted target 0.005, the
    inherited O3 knobs, and the 900k/24h caps -- never to be raised."""

    assert o3p.FROZEN == {
        "r_in": 12.0, "r_out": 18.0, "dt": 8.5, "m": 1.0,
        "target_ratio": 0.005,
        "n_sub": 16, "k_micro": 4, "d_switch": 0.25,
        "init_rho": 12, "init_psi": 12,
        "max_depth": 18,
        "max_calls": 900_000,
        "max_wall_s": 86_400.0,
    }


def test_the_freeze_commit_carried_no_result():
    """The freeze-commit-split discipline, as a permanent historical
    fact (the O3 pattern): the freeze merge commit must NOT contain
    the result artifact -- the separately approved execution wrote it
    later. The working tree DOES contain it once the results commit
    lands, so the assertion is against the freeze merge's tree."""

    probe = subprocess.run(
        ["git", "cat-file", "-e",
         "3b869c3a9cb70f6d1d1354c190260c3f8358892f:"
         "docs/prereg/p14_o3p_volume.json"],
        cwd=_REPO, capture_output=True)
    assert probe.returncode != 0


_SHA = "a" * 40


def test_preflight_refuses_when_a_result_exists(monkeypatch,
                                                tmp_path):
    fake = tmp_path / "p14_o3p_volume.json"
    fake.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(o3p, "_ARTIFACT", fake)
    monkeypatch.setattr(o3p, "verify_freeze",
                        lambda stage: {"files": {}})
    monkeypatch.setattr(o3p, "_git_state",
                        lambda: {"rev": _SHA, "dirty": False})
    with pytest.raises(SystemExit, match="write-once"):
        o3p.preflight(_SHA)


def test_preflight_refuses_a_dirty_tree(monkeypatch):
    monkeypatch.setattr(o3p, "verify_freeze",
                        lambda stage: {"files": {}})
    monkeypatch.setattr(o3p, "_git_state",
                        lambda: {"rev": _SHA, "dirty": True})
    with pytest.raises(SystemExit, match="clean exact checkout"):
        o3p.preflight(_SHA)


def test_preflight_refuses_a_missing_o3_base(monkeypatch, tmp_path):
    """The O3 interval is the intersection base; a checkout without
    it cannot run the campaign."""

    monkeypatch.setattr(o3p, "verify_freeze",
                        lambda stage: {"files": {}})
    monkeypatch.setattr(o3p, "_git_state",
                        lambda: {"rev": _SHA, "dirty": False})
    monkeypatch.setattr(o3p, "_ARTIFACT", tmp_path / "absent.json")
    monkeypatch.setattr(o3p, "_O3_ARTIFACT", tmp_path / "no_o3.json")
    with pytest.raises(SystemExit, match="intersection base"):
        o3p.preflight(_SHA)


# ------------------------------ the approved-SHA gate (PR #81 R2)

def test_a_missing_freeze_rev_is_refused_by_the_cli():
    """--freeze-rev is REQUIRED: argparse exits 2 before any freeze
    machinery runs."""

    script = str(_REPO / "experiments" / "oracle"
                 / "o3p_frozen_volume.py")
    out = subprocess.run([sys.executable, script, "--preflight"],
                         capture_output=True, text=True)
    assert out.returncode == 2
    assert "--freeze-rev" in out.stderr


@pytest.mark.parametrize("bad", ["21cff16", "a" * 39, "a" * 41,
                                 "g" * 40, "", "  "])
def test_a_short_or_malformed_rev_is_refused(bad):
    with pytest.raises(SystemExit, match="full 40-hex"):
        o3p.require_full_sha(bad)


def test_a_clean_matching_tree_at_another_sha_is_refused(monkeypatch):
    """The manifest cannot certify itself: clean tree + every digest
    matching is NOT enough, because a later commit can edit a protocol
    file and re-pin the manifest in the same commit. Only the exact
    approved full SHA runs."""

    monkeypatch.setattr(o3p, "verify_freeze",
                        lambda stage: {"files": {}})   # digests "match"
    monkeypatch.setattr(o3p, "_git_state",
                        lambda: {"rev": "b" * 40, "dirty": False})
    with pytest.raises(SystemExit, match="exact\\s+approved commit"):
        o3p.preflight(_SHA)


def _stub_res():
    from certified_interval import Iv
    return {"v": Iv(1) / Iv(3), "ratio": 0.004,
            "status": "target-met", "termination_reason": "target-met",
            "calls": 10, "cells": 3, "wall_s": 0.1,
            "max_depth_reached": 2, "cells_at_max_depth": 0,
            "uncosted_cells": 0, "modes": {}, "raw_width_by_mode": {},
            "raw_total_before_intersection": 0.0,
            "certified_total_after_intersection": 0.0,
            "intersection_active": False, "curve": []}


def test_a_mid_run_rebind_to_another_clean_commit_refuses(monkeypatch):
    """Review PR #81 R3: the baseline is the preflight's VERIFIED
    state and the exit rev is compared to the approved SHA directly.
    A worktree switched to a different clean commit between the
    approval check and the run -- even one whose re-pinned digests
    all match -- must refuse before publishing anything."""

    approved, other = "a" * 40, "b" * 40
    monkeypatch.setattr(o3p, "preflight", lambda rev: {
        "manifest": {}, "git": {"rev": approved, "dirty": False}})
    monkeypatch.setattr(o3p, "_git_state",
                        lambda: {"rev": other, "dirty": False})
    monkeypatch.setattr(o3p, "verify_freeze", lambda stage: {})
    monkeypatch.setattr(o3p, "assemble",
                        lambda cfg, progress=None: _stub_res())
    published = []
    monkeypatch.setattr(o3p, "_publish_write_once",
                        lambda p, s: published.append(p))
    monkeypatch.setattr(sys, "argv",
                        ["prog", "--freeze-rev", approved])
    with pytest.raises(SystemExit, match="approved freeze SHA"):
        o3p.main()
    assert not published


def test_main_publishes_on_an_unmoved_approved_checkout(monkeypatch,
                                                        tmp_path):
    """The real main() wiring end to end (the O4b lesson): approved
    SHA held from preflight through exit, artifact published with the
    gated V_ref' block and the intersection."""

    approved = "a" * 40
    monkeypatch.setattr(o3p, "preflight", lambda rev: {
        "manifest": {}, "git": {"rev": approved, "dirty": False}})
    monkeypatch.setattr(o3p, "_git_state",
                        lambda: {"rev": approved, "dirty": False})
    monkeypatch.setattr(o3p, "verify_freeze", lambda stage: {})
    monkeypatch.setattr(o3p, "assemble",
                        lambda cfg, progress=None: _stub_res())
    out = tmp_path / "p14_o3p_volume.json"
    monkeypatch.setattr(o3p, "_ARTIFACT", out)
    monkeypatch.setattr(sys, "argv",
                        ["prog", "--freeze-rev", approved])
    o3p.main()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert art["code"] == {"start": {"rev": approved, "dirty": False},
                           "end": {"rev": approved, "dirty": False}}
    # the stub interval [1/3-ish] is disjoint from the O3 base, so
    # the record must carry the inconsistency and a BLOCKED V_ref'
    assert art["intersection_with_o3"][
        "certification_inconsistency"] is True
    assert art["v_ref_prime_recommendation"]["value"] is None


def test_the_exact_approved_sha_passes(monkeypatch, tmp_path):
    o3_base = tmp_path / "o3.json"
    o3_base.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(o3p, "verify_freeze",
                        lambda stage: {"files": {}})
    monkeypatch.setattr(o3p, "_git_state",
                        lambda: {"rev": _SHA, "dirty": False})
    monkeypatch.setattr(o3p, "_ARTIFACT", tmp_path / "absent.json")
    monkeypatch.setattr(o3p, "_O3_ARTIFACT", o3_base)
    checks = o3p.preflight(_SHA)
    assert checks["git"]["rev"] == _SHA
    # case/whitespace normalization: the approval names one commit
    assert o3p.preflight("  " + _SHA.upper() + "  ")["git"]["rev"] \
        == _SHA


# ------------------------------------------- the intersection rule

def test_the_o3_base_is_the_committed_artifact():
    base = o3p.o3_interval()
    o3 = json.loads(
        (o3p._O3_ARTIFACT).read_text(encoding="utf-8"))["result"]
    assert base["v_lo"] == o3["v_lo"] == 56.21273686780051
    assert base["v_hi"] == o3["v_hi"] == 57.348018526434714
    assert base["artifact"] == "p14_o3_volume.json"


def test_a_contained_interval_is_consistent():
    x = o3p.intersect_with_o3(56.5, 57.0)
    assert x["consistent"] is True
    assert x["certification_inconsistency"] is False
    assert (x["v_lo"], x["v_hi"]) == (56.5, 57.0)
    assert "STANDALONE" in x["note"]


def test_a_partial_overlap_is_consistent_and_clipped():
    x = o3p.intersect_with_o3(57.0, 58.0)
    assert x["consistent"] is True
    assert x["v_lo"] == 57.0
    assert x["v_hi"] == 57.348018526434714


def test_an_empty_intersection_is_a_certification_inconsistency():
    """Disjointness does not pick a winner: it indicts at least one
    certification, machine-readably, and blocks downstream use."""

    x = o3p.intersect_with_o3(58.0, 59.0)
    assert x["consistent"] is False
    assert x["certification_inconsistency"] is True
    assert "v_lo" not in x
    assert "No downstream stage may consume" in x["note"]
    assert "does not adjudicate" in x["note"]


def test_v_ref_prime_is_blocked_on_inconsistency():
    """Review PR #81 R1: a failed record must carry no usable
    midpoint -- a direct-reading consumer (the o4_sizing pattern)
    would otherwise size against a certification inconsistency."""

    bad = o3p.intersect_with_o3(58.0, 59.0)
    rec = o3p.v_ref_prime(bad, 58.0, 59.0)
    assert rec["value"] is None
    assert rec["status"].startswith("BLOCKED")
    assert "no downstream stage may size" in rec["status"]
    assert "definition" not in rec


def test_v_ref_prime_is_the_standalone_midpoint_on_consistency():
    good = o3p.intersect_with_o3(56.5, 57.0)
    rec = o3p.v_ref_prime(good, 56.5, 57.0)
    assert rec["value"] == 0.5 * (56.5 + 57.0)
    assert rec["definition"] == (
        "midpoint of the STANDALONE O3' interval")
    assert "recommendation only" in rec["status"]


# ------------------------------------------ serialization/publish

def test_result_serialization_is_outward():
    from certified_interval import Iv
    third = Iv(1) / Iv(3)
    res = {"v": third, "ratio": 0.0, "status": "target-met",
           "termination_reason": "target-met", "calls": 1,
           "cells": 1, "wall_s": 0.0, "max_depth_reached": 0,
           "cells_at_max_depth": 0, "uncosted_cells": 0,
           "modes": {}, "raw_width_by_mode": {},
           "raw_total_before_intersection": 0.0,
           "certified_total_after_intersection": 0.0,
           "intersection_active": False}
    out = o3p._serialize_result(res)
    import gmpy2
    assert gmpy2.mpfr(out["v_lo"], 200) <= third.lo
    assert gmpy2.mpfr(out["v_hi"], 200) >= third.hi
    assert out["v_lo"] < out["v_hi"]


def test_target_not_met_serializes_with_full_provenance():
    """A fired cap publishes the interval reached at that moment plus
    the exact termination reason -- the same shape as success, never
    a suppressed record."""

    from certified_interval import Iv
    res = {"v": Iv(1) / Iv(3), "ratio": 0.25,
           "status": "target-not-met",
           "termination_reason": "max_calls", "calls": 900_000,
           "cells": 7, "wall_s": 1.0, "max_depth_reached": 9,
           "cells_at_max_depth": 2, "uncosted_cells": 0,
           "modes": {}, "raw_width_by_mode": {},
           "raw_total_before_intersection": 0.0,
           "certified_total_after_intersection": 0.0,
           "intersection_active": False}
    out = o3p._serialize_result(res)
    assert out["status"] == "target-not-met"
    assert out["termination_reason"] == "max_calls"
    assert out["v_lo"] < out["v_hi"]        # still a certified interval


def test_publish_is_atomic_and_write_once(tmp_path):
    dst = tmp_path / "result.json"
    o3p._publish_write_once(dst, '{"first": true}\n')
    with pytest.raises(SystemExit, match="write-once"):
        o3p._publish_write_once(dst, '{"second": true}\n')
    assert json.loads(dst.read_text(encoding="utf-8")) == {
        "first": True}
    assert list(tmp_path.iterdir()) == [dst]


def test_manifest_paths_are_lf_pinned():
    attrs = (_REPO / ".gitattributes").read_text(encoding="utf-8")
    manifest = json.loads(o3p._MANIFEST.read_text(encoding="utf-8"))
    for rel in manifest["files"]:
        name = rel.rsplit("/", 1)[-1]
        parent = rel.rsplit("/", 1)[0] if "/" in rel else ""
        covered = any(
            line.split()[0] in (rel, name,
                                f"{parent}/*{Path(name).suffix}")
            and "eol=lf" in line
            for line in attrs.splitlines()
            if line.strip() and not line.startswith("#"))
        assert covered, rel


def test_cli_is_fail_closed():
    script = str(_REPO / "experiments" / "oracle"
                 / "o3p_frozen_volume.py")
    bad = subprocess.run([sys.executable, script, "--preflght"],
                         capture_output=True, text=True)
    assert bad.returncode == 2
    helped = subprocess.run([sys.executable, script, "--help"],
                            capture_output=True, text=True)
    assert helped.returncode == 0
    assert "--preflight" in helped.stdout


def test_docstring_keeps_the_epistemic_grade():
    flat = " ".join(o3p.__doc__.split())
    assert "EXECUTION IS NOT YET APPROVED" in flat
    assert "FORBIDDEN" in flat
    assert "NO seed" in flat
    assert "IMMUTABLE" in flat
    assert "CERTIFICATION INCONSISTENCY" in flat
    assert "not a claim of optimality" in flat.lower()
    assert "Projections, not certifications" in flat


# ------------------------------------------ the projection contract

def test_the_projection_table_reproduces_from_the_artifacts():
    """The freeze document's cost table is executable: both curves
    fitted separately, the range and the central projection pinned,
    and the caps sit at least 20% above the worst fitted window."""

    t = proj.table(0.005)
    assert len(t["rows"]) == 8                     # 4 windows x 2 curves
    sources = {r["source"] for r in t["rows"]}
    assert sources == {"neighbor ladder", "O3 own curve"}
    lo, hi = t["range_calls"]
    assert round(lo) == 467_365                    # O3 full-curve fit
    assert round(hi) == 749_784                    # O3 ratio<=0.02 fit
    assert round(t["central_projection_calls"]) == 551_501
    assert o3p.FROZEN["max_calls"] >= 1.20 * hi
    assert t["grade"] == "projection, not certification"


def test_projection_exponents_are_fitted_not_assumed():
    """Every fitted exponent is a real number from the committed
    artifacts, spanning the p=2 central value rather than equalling
    it."""

    t = proj.table(0.005)
    ps = [r["p"] for r in t["rows"]]
    assert min(ps) < 2.0 < max(ps)
    assert all(1.5 < p < 2.6 for p in ps)
