from __future__ import annotations

# ruff: noqa: E402,I001

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.abspath("experiments"))


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    paths = [os.path.abspath("src"), os.path.abspath("experiments")]
    if env.get("PYTHONPATH"):
        paths.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


def _python_executable() -> str:
    venv_python = Path(".venv/bin/python")
    return str(venv_python) if venv_python.exists() else sys.executable


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture()
def m40_bundle(tmp_path: Path) -> Path:
    (tmp_path / "manifests_v2").mkdir(parents=True)
    (tmp_path / "manifests_v2" / "old_manifest.json").write_text(
        '{"manifest_id":"old","profile_label":"old_profile"}',
        encoding="utf-8",
    )
    _write_csv(
        tmp_path / "data" / "v3_manifest_family_design.csv",
        [
            {
                "family_name": "rank_gap_multi_manifest_v3",
                "family_kind": "structured",
                "planned_manifest_count": "5",
            }
        ],
    )
    return tmp_path


def _run(script: str, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [_python_executable(), f"experiments/{script}", *args],
        cwd=Path.cwd(),
        env=_env_with_src(),
        check=False,
        capture_output=True,
        text=True,
    )


def test_exp215_exact_sanity() -> None:
    from experiments.exp215_protocol_metadata_exact_sanity import run_experiment

    rows = run_experiment()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp216_exact_sanity() -> None:
    from experiments.exp216_response_profile_invariance_exact_sanity import (
        run_experiment,
    )

    rows = run_experiment()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp217_audit_smoke(m40_bundle: Path) -> None:
    result = _run(
        "exp217_current_profile_protocol_invariance_audit.py",
        "--output-dir",
        str(m40_bundle),
    )

    assert result.returncode == 0, result.stderr
    assert (
        m40_bundle / "data" / "response_profile_protocol_invariance_audit.csv"
    ).exists()


def test_exp218_patch_design(m40_bundle: Path) -> None:
    result = _run(
        "exp218_v3_protocol_invariant_patch_design.py",
        "--output-dir",
        str(m40_bundle),
    )

    assert result.returncode == 0, result.stderr
    assert (m40_bundle / "data" / "v3_protocol_invariant_family_patch.csv").exists()


def test_exp219_patched_preregistration_export(m40_bundle: Path) -> None:
    test_exp218_patch_design(m40_bundle)
    result = _run(
        "exp219_v3_protocol_patched_preregistration_export.py",
        "--output-dir",
        str(m40_bundle),
    )

    assert result.returncode == 0, result.stderr
    assert (
        m40_bundle / "remediation" / "v3_protocol_patched_preregistration_m40.json"
    ).exists()


def test_exp220_patch_audit(m40_bundle: Path) -> None:
    test_exp218_patch_design(m40_bundle)
    result = _run(
        "exp220_v3_protocol_patch_audit.py",
        "--output-dir",
        str(m40_bundle),
    )

    assert result.returncode == 0, result.stderr
    assert (m40_bundle / "data" / "v3_protocol_patch_audit.csv").exists()


def test_exp221_language_audit(m40_bundle: Path) -> None:
    result = _run(
        "exp221_protocol_invariance_language_audit.py",
        "--output-dir",
        str(m40_bundle),
    )

    assert result.returncode == 0, result.stderr
    assert (m40_bundle / "data" / "protocol_invariance_language_audit.csv").exists()


def test_exp222_no_execution_audit(m40_bundle: Path) -> None:
    test_exp219_patched_preregistration_export(m40_bundle)
    result = _run(
        "exp222_protocol_patch_no_execution_audit.py",
        "--output-dir",
        str(m40_bundle),
    )

    assert result.returncode == 0, result.stderr
    assert (m40_bundle / "data" / "protocol_patch_no_execution_audit.csv").exists()


def test_exp223_final_sanity() -> None:
    result = _run("exp223_protocol_invariance_final_sanity.py")

    assert result.returncode == 0, result.stderr
    assert Path("outputs/data/protocol_invariance_final_sanity.csv").exists()
