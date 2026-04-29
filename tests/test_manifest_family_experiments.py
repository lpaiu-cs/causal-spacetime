from __future__ import annotations

# ruff: noqa: E402,I001

import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.abspath("experiments"))

from experiments.exp148_manifest_family_exact_sanity import run_experiment as run_exp148
from experiments.exp153_manifest_family_no_retuning_audit import (
    run_experiment as run_exp153,
)
from experiments.exp154_manifest_family_report_card import run_experiment as run_exp154
from experiments.exp155_manifest_family_comparison_exact_sanity import (
    run_experiment as run_exp155,
)
from manifest_representation_experiment_helpers import build_exact_manifest


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    paths = [os.path.abspath("src"), os.path.abspath("experiments")]
    if env.get("PYTHONPATH"):
        paths.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


def _skip_without_matplotlib() -> None:
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")


def _prepare_manifest(tmp_path) -> None:  # type: ignore[no-untyped-def]
    build_exact_manifest(tmp_path, "manifest.json")


def test_exp148_exact_sanity(tmp_path) -> None:  # type: ignore[no-untyped-def]
    rows = run_exp148(tmp_path)

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp149_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    _prepare_manifest(tmp_path)
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp149_manifest_family_fit_comparison.py",
            "--manifest-dir",
            str(tmp_path / "manifests"),
            "--dims",
            "1",
            "2",
            "--steps",
            "20",
            "--restarts",
            "1",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert "Wrote manifest family fit comparison" in result.stdout


def test_exp150_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    _prepare_manifest(tmp_path)
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp150_manifest_family_null_taxonomy.py",
            "--manifest-dir",
            str(tmp_path / "manifests"),
            "--embedding-dim",
            "1",
            "--null-repetitions",
            "1",
            "--steps",
            "20",
            "--restarts",
            "1",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert "Wrote manifest family null taxonomy" in result.stdout


def test_exp151_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    _prepare_manifest(tmp_path)
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp151_manifest_family_stricter_criteria.py",
            "--manifest-dir",
            str(tmp_path / "manifests"),
            "--dims",
            "1",
            "--steps",
            "20",
            "--restarts",
            "1",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert "Wrote manifest family stricter criteria" in result.stdout


def test_exp152_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    _prepare_manifest(tmp_path)
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp152_manifest_family_failed_manifest_accounting.py",
            "--manifest-dir",
            str(tmp_path / "manifests"),
            "--generate-failed-controls",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert "Wrote manifest family failed-manifest accounting" in result.stdout


def test_exp153_audit(tmp_path) -> None:  # type: ignore[no-untyped-def]
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "manifest_family_failed_manifest_accounting.csv").write_text(
        "row_type\nfamily_count\n",
        encoding="utf-8",
    )
    rows = run_exp153(tmp_path)

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp154_summary() -> None:
    rows = run_exp154()

    assert isinstance(rows, list)


def test_exp155_exact_sanity(tmp_path) -> None:  # type: ignore[no-untyped-def]
    rows = run_exp155(tmp_path)

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
