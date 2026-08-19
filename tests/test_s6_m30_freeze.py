"""S6 rung M = 3.0 oracle-freeze contract tests.

The frozen configuration pinned to the ladder's one exact path; the
import gate's refusals; the SHA gates; the preflight refusals; the
write-once publication; the manifest-vs-tree identity; and the doc
discipline. NO RESULTS: this freeze publishes no artifact -- the
approved run happens from the merged freeze head, not from here."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import s6_m30_frozen_volume as run  # noqa: E402
import s6_rungs as s6  # noqa: E402

_SHA = "a" * 40


def test_the_frozen_configuration_is_the_ladder_rung():
    assert run.FROZEN == {
        "r_in": 12.0, "r_out": 18.0,
        "dt": 12.442423039673733, "m": 3.0,
        "mu": 0.4,
        "target_ratio": 0.005,
        "n_sub": 16, "k_micro": 4, "d_switch": 0.25,
        "init_rho": 12, "init_psi": 12,
        "max_depth": 18,
        "max_calls": 900_000,
        "max_wall_s": 86_400.0,
    }
    # the one exact path is the only source
    assert run.FROZEN["dt"] == s6.RUNG_CONSTANTS[3.0]["dt"]
    assert run.FROZEN["mu"] == s6.RUNG_CONSTANTS[3.0]["mu"]
    assert run.FROZEN["dt"] == s6.dt(3.0)


def test_the_import_gate_verified_the_whole_lemma_table():
    assert run.LEMMA_TABLE["all_pass"] is True
    assert set(run.LEMMA_TABLE["rows"]) >= {
        "exterior", "k_negative", "w_monotone", "q_positive",
        "perihelion_patch", "l5_winding", "l4_nonempty",
        "l4_inner_shell", "l4_outer_shell", "l4_polar_cap"}


def test_the_import_gate_refuses_a_drifted_constant(monkeypatch):
    monkeypatch.setitem(run.FROZEN, "dt", 9.73)
    with pytest.raises(SystemExit, match="one exact path"):
        run._import_gate()


def test_the_import_gate_refuses_a_failed_lemma_row(monkeypatch):
    def failed(m):
        return {"all_pass": False,
                "rows": {"l5_winding": {"certified": False,
                                        "margin": -1.0}}}

    monkeypatch.setattr(run.s6, "lemma_table", failed)
    with pytest.raises(SystemExit, match="lemma table FAILED"):
        run._import_gate()


@pytest.mark.parametrize("bad", ["0fd0328", "b" * 39, "g" * 40, ""])
def test_short_or_malformed_freeze_rev_is_refused(bad):
    with pytest.raises(SystemExit, match="full 40-hex"):
        run.require_full_sha(bad)


def _pass_static(monkeypatch, tmp_path, rev=_SHA):
    monkeypatch.setattr(run, "verify_freeze", lambda stage: {})
    monkeypatch.setattr(run, "_git_state",
                        lambda: {"rev": rev, "dirty": False})
    monkeypatch.setattr(run, "_ARTIFACT", tmp_path / "a.json")


def test_preflight_passes_only_on_the_exact_sha(monkeypatch, tmp_path):
    _pass_static(monkeypatch, tmp_path)
    assert run.preflight(_SHA)["git"]["rev"] == _SHA
    with pytest.raises(SystemExit, match="exact\\s+approved commit"):
        run.preflight("b" * 40)


def test_preflight_refuses_a_dirty_tree(monkeypatch, tmp_path):
    _pass_static(monkeypatch, tmp_path)
    monkeypatch.setattr(run, "_git_state",
                        lambda: {"rev": _SHA, "dirty": True})
    with pytest.raises(SystemExit, match="dirty"):
        run.preflight(_SHA)


def test_preflight_refuses_a_leftover_artifact(monkeypatch, tmp_path):
    _pass_static(monkeypatch, tmp_path)
    (tmp_path / "a.json").write_text("{}", encoding="utf-8")
    with pytest.raises(SystemExit, match="write-once"):
        run.preflight(_SHA)


def test_preflight_reruns_the_lemma_table(monkeypatch, tmp_path):
    _pass_static(monkeypatch, tmp_path)
    monkeypatch.setattr(run.s6, "lemma_table",
                        lambda m: {"all_pass": False, "rows": {}})
    with pytest.raises(SystemExit, match="no longer passes"):
        run.preflight(_SHA)


def test_publication_is_write_once(tmp_path):
    dest = tmp_path / "v.json"
    run._publish_write_once(dest, "{\"first\": true}\n")
    with pytest.raises(SystemExit, match="write-once"):
        run._publish_write_once(dest, "{\"second\": true}\n")
    assert json.loads(dest.read_text(encoding="utf-8")) == {
        "first": True}


def test_v_ref_rung_is_a_recommendation_only():
    v = run.v_ref_rung(60.0, 62.0)
    assert v["value"] == 61.0
    assert "recommendation only" in v["status"]
    assert "96-bit" in v["status"]


def test_the_committed_manifest_is_the_tree_it_describes():
    committed = json.loads(run._MANIFEST.read_text(encoding="utf-8"))
    assert committed == run.build_manifest()
    assert set(committed["files"]) == set(run.PROTOCOL_SURFACE)
    assert ("experiments/oracle/s6_rungs.py"
            in committed["files"])


def test_the_doc_freezes_the_discipline():
    doc = (_REPO / "docs" / "prereg"
           / "p14_s6_m30_volume.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    # this rung's execution IS approved -- but the run is still bound
    # to one clean checkout of the merged freeze head, and the seeds
    # belong to the later stages, not to this deterministic one
    assert "Execution of this rung **is approved**" in flat
    assert "never the merge commit" in flat
    assert "no seed, no sprinkling, no reservation" in flat
    assert "`40,000,301` remains unallocated and unspent" in flat
    assert "no auto-raise, ever" in flat
    assert "published as is" in flat
    assert "12.442423039673733" in flat
    assert "| **0.4** |" in flat
    assert "recommendation only" in flat
    # the design extrapolation is a cost projection, never a gate
    assert "are **NOT frozen anywhere**" in flat
    assert "`A ≈ 469`" in flat
    assert "no joint primary" in flat
