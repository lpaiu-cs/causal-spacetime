from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp40_representation_stability_exact_sanity import run_experiment


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    src_path = os.path.abspath("src")
    env["PYTHONPATH"] = (
        src_path
        if not env.get("PYTHONPATH")
        else f"{src_path}{os.pathsep}{env['PYTHONPATH']}"
    )
    return env


def _skip_without_matplotlib() -> None:
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")


def test_exp36_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp36_heldout_ordinal_embedding_validation.py",
            "--n-points-values",
            "10",
            "--constraint-counts",
            "100",
            "--repetitions",
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
        timeout=90,
    )

    assert "Wrote held-out ordinal embedding validation data" in result.stdout
    assert (tmp_path / "data" / "heldout_ordinal_embedding_validation.csv").exists()


def test_exp37_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp37_embedding_stability_under_subsampling.py",
            "--n-points",
            "10",
            "--total-constraints",
            "120",
            "--subset-sizes",
            "30",
            "--num-subsets",
            "2",
            "--repetitions",
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
        timeout=90,
    )

    assert "Wrote embedding stability data" in result.stdout
    assert (tmp_path / "data" / "embedding_stability_under_subsampling.csv").exists()


def test_exp38_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp38_observer_order_null_baseline.py",
            "--N",
            "80",
            "--tick-values",
            "16",
            "--constraint-counts",
            "100",
            "--repetitions",
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
        timeout=90,
    )

    assert "Wrote observer-order null baseline data" in result.stdout
    assert (tmp_path / "data" / "observer_order_null_baseline.csv").exists()


def test_exp39_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp39_effective_metric_complexity_curve.py",
            "--candidate-dims",
            "1",
            "2",
            "--n-points",
            "10",
            "--num-constraints",
            "100",
            "--repetitions",
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
        timeout=90,
    )

    assert "Wrote effective metric complexity curve data" in result.stdout
    assert (tmp_path / "data" / "effective_metric_complexity_curve.csv").exists()


def test_exp40_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
