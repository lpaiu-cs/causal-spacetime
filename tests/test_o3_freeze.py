"""O3 freeze contract tests.

Freeze-commit discipline: the frozen configuration, the
content-addressed manifest, the write-once/no-results boundary, and
the fail-closed CLI are pinned BEFORE any execution -- the campaign
run itself is separately approved and happens from a clean exact
checkout. The digest checks are platform-independent (LF-pinned
paths); the environment lock binds the execution host only, so its
MECHANISM is tested here, never the CI host's own versions."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))

import o3_frozen_volume as o3  # noqa: E402


def test_manifest_digests_match_the_working_tree():
    """Every pinned file matches its frozen digest, the runner pins
    ITSELF, and the pinned set covers the whole certified surface
    plus the measurement basis of the configuration choice plus the
    dependency-lock surface."""

    manifest = o3.verify_digests("test")
    files = set(manifest["files"])
    assert files == {
        "experiments/oracle/certified_interval.py",
        "experiments/oracle/certified_flight_time.py",
        "experiments/oracle/volume_oracle.py",
        "experiments/oracle/o3_frozen_volume.py",
        "docs/theory/schwarzschild_volume_oracle_certification.md",
        "docs/prereg/p14_oracle_price.json",
        "docs/prereg/p14_oracle_mode_width.json",
        "pyproject.toml",
    }


def test_manifest_locks_the_rounding_instrument():
    manifest = json.loads(o3._MANIFEST.read_text(encoding="utf-8"))
    env = manifest["environment"]
    assert set(env) == {"python", "gmpy2", "mpfr", "gmp"}
    assert env["mpfr"].startswith("MPFR ")
    assert env["gmp"].startswith("GMP ")


def test_environment_verification_refuses_drift(monkeypatch):
    """The mechanism, not the CI host: a drifted MPFR must refuse,
    an exact match must pass."""

    manifest = json.loads(o3._MANIFEST.read_text(encoding="utf-8"))
    monkeypatch.setattr(o3, "_environment",
                        lambda: dict(manifest["environment"]))
    o3.verify_environment("test", manifest)  # exact match passes
    drifted = dict(manifest["environment"], mpfr="MPFR 0.0.0")
    monkeypatch.setattr(o3, "_environment", lambda: drifted)
    with pytest.raises(SystemExit, match="environment drift"):
        o3.verify_environment("test", manifest)


def test_frozen_configuration_is_the_ruled_one():
    """The ruling's exact freeze: anchors, target, algorithm knobs,
    and caps -- with the caps never to be raised."""

    assert o3.FROZEN == {
        "r_in": 12.0, "r_out": 18.0, "dt": 8.5, "m": 1.0,
        "target_ratio": 0.01,
        "n_sub": 16, "k_micro": 4, "d_switch": 0.25,
        "init_rho": 12, "init_psi": 12,
        "max_depth": 18,
        "max_calls": 600_000,
        "max_wall_s": 86_400.0,
    }


def test_freeze_commit_carried_no_result():
    """The freeze-commit-split discipline, as a permanent historical
    fact: the freeze merge commit must NOT contain the result
    artifact -- the separately approved execution wrote it later.
    (The working tree DOES contain it once the results commit lands,
    so the assertion is against the freeze commit's tree, not the
    checkout.)"""

    probe = subprocess.run(
        ["git", "cat-file", "-e",
         "785148ecf8be8b7b1baaa1f3866d5d827b8dfdf7:"
         "docs/prereg/p14_o3_volume.json"],
        cwd=_REPO, capture_output=True)
    assert probe.returncode != 0


def test_preflight_refuses_when_a_result_exists(monkeypatch,
                                                tmp_path):
    """Write-once: a present result artifact must abort preflight
    before anything else can run."""

    fake = tmp_path / "p14_o3_volume.json"
    fake.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(o3, "_ARTIFACT", fake)
    monkeypatch.setattr(o3, "verify_freeze",
                        lambda stage: {"files": {}})
    monkeypatch.setattr(o3, "_git_state",
                        lambda: {"rev": "x", "dirty": False})
    with pytest.raises(SystemExit, match="write-once"):
        o3.preflight()


def test_preflight_refuses_a_dirty_tree(monkeypatch):
    monkeypatch.setattr(o3, "verify_freeze",
                        lambda stage: {"files": {}})
    monkeypatch.setattr(o3, "_git_state",
                        lambda: {"rev": "x", "dirty": True})
    with pytest.raises(SystemExit, match="clean exact checkout"):
        o3.preflight()


def test_result_serialization_is_outward():
    """R1: the artifact's binary64 endpoints must still ENCLOSE the
    MPFR interval -- a nearest float() can shrink a certified
    interval inward, after which the stored numbers no longer
    contain the true value. 1/3 is not binary64-representable, so
    both directions must move strictly outward."""

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
    out = o3._serialize_result(res)
    import gmpy2
    assert gmpy2.mpfr(out["v_lo"], 200) <= third.lo
    assert gmpy2.mpfr(out["v_hi"], 200) >= third.hi
    assert out["v_lo"] < out["v_hi"]


def test_publish_is_atomic_and_write_once(tmp_path):
    """R1: publication must be no-clobber ATOMICALLY -- the exists()
    check in preflight is not linked to the write, and a plain 'w'
    open would let a second long run overwrite the first
    observation, while a crash mid-write would leave a truncated
    destination that preflight then treats as the result. The
    hard-link publish fails on an existing destination and never
    exposes a partial file."""

    dst = tmp_path / "result.json"
    o3._publish_write_once(dst, '{"first": true}\n')
    assert json.loads(dst.read_text(encoding="utf-8")) == {
        "first": True}
    with pytest.raises(SystemExit, match="write-once"):
        o3._publish_write_once(dst, '{"second": true}\n')
    # the first observation stands and no temp litter remains
    assert json.loads(dst.read_text(encoding="utf-8")) == {
        "first": True}
    assert list(tmp_path.iterdir()) == [dst]


def test_pyproject_is_lf_pinned_for_the_manifest():
    """R1: every digest-pinned path must be eol=lf so the raw sha256
    is checkout-independent -- pyproject.toml is in the manifest, so
    an autocrlf checkout without the attribute would materialize
    CRLF and fail preflight on an exact freeze checkout."""

    attrs = (_REPO / ".gitattributes").read_text(encoding="utf-8")
    manifest = json.loads(o3._MANIFEST.read_text(encoding="utf-8"))
    for rel in manifest["files"]:
        name = rel.rsplit("/", 1)[-1]
        parent = rel.rsplit("/", 1)[0] if "/" in rel else ""
        covered = any(
            line.split()[0] in (rel, name, f"{parent}/*{Path(name).suffix}")
            and "eol=lf" in line
            for line in attrs.splitlines()
            if line.strip() and not line.startswith("#"))
        assert covered, rel


def test_cli_is_fail_closed():
    """Unknown arguments exit 2 and --help exits 0, both before any
    freeze machinery runs -- same contract as the S4/S5 runners."""

    script = str(_REPO / "experiments" / "oracle"
                 / "o3_frozen_volume.py")
    bad = subprocess.run([sys.executable, script, "--preflght"],
                         capture_output=True, text=True)
    assert bad.returncode == 2
    helped = subprocess.run([sys.executable, script, "--help"],
                            capture_output=True, text=True)
    assert helped.returncode == 0
    assert "--preflight" in helped.stdout


def test_docstring_keeps_the_epistemic_grade():
    """The ruling's wording boundary: n_sub=16 is a freeze CHOICE
    based on current measurements, not an optimality claim; raising
    the caps is forbidden; execution needs separate approval."""

    doc = o3.__doc__
    assert "NOT a claim of\n  optimality" in doc or (
        "NOT a claim of optimality" in doc.replace("\n  ", " "))
    assert "EXECUTION IS NOT YET APPROVED" in doc
    assert "FORBIDDEN" in doc