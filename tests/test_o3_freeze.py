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


def test_freeze_commit_carries_no_result():
    """The freeze-commit-split discipline: the result artifact must
    NOT exist until the separately approved execution writes it."""

    assert not o3._ARTIFACT.exists()


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