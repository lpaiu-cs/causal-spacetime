from __future__ import annotations

import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.abspath("experiments"))

from manifest_representation_experiment_helpers import build_exact_manifest

from experiments.exp140_manifest_representability_exact_sanity import (
    run_experiment as run_exp140,
)
from experiments.exp146_frozen_manifest_representation_summary import (
    run_experiment as run_exp146,
)
from experiments.exp147_manifest_representation_no_metric_exact_sanity import (
    run_experiment as run_exp147,
)


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


def test_exp140_exact_sanity(tmp_path) -> None:  # type: ignore[no-untyped-def]
    rows = run_exp140(tmp_path)

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp141_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    _prepare_manifest(tmp_path)
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp141_frozen_manifest_ordinal_representation.py",
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
    assert "Wrote frozen manifest ordinal representation" in result.stdout


def test_exp142_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    _prepare_manifest(tmp_path)
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp142_frozen_manifest_representation_nulls.py",
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
    assert "Wrote frozen manifest representation nulls" in result.stdout


def test_exp143_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    _prepare_manifest(tmp_path)
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp143_frozen_manifest_fit_stability.py",
            "--manifest-dir",
            str(tmp_path / "manifests"),
            "--embedding-dim",
            "1",
            "--restart-count",
            "2",
            "--steps",
            "20",
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
    assert "Wrote frozen manifest fit stability" in result.stdout


def test_exp144_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    _prepare_manifest(tmp_path)
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp144_frozen_manifest_dimension_complexity_curve.py",
            "--manifest-dir",
            str(tmp_path / "manifests"),
            "--candidate-dims",
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
    assert "Wrote frozen manifest dimension-complexity curve" in result.stdout


def test_exp145_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _prepare_manifest(tmp_path)
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp145_failed_manifest_no_fit_controls.py",
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
    assert "Wrote failed manifest no-fit controls" in result.stdout


def test_exp146_summary_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    rows = run_exp146(tmp_path)

    assert isinstance(rows, list)


def test_exp147_exact_sanity(tmp_path) -> None:  # type: ignore[no-untyped-def]
    rows = run_exp147(tmp_path)

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
