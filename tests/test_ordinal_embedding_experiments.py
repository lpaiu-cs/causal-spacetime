from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp35_ordinal_embedding_exact_sanity import run_experiment


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


def test_exp31_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp31_ordinal_embedding_recovery.py",
            "--true-dims",
            "1",
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

    assert "Wrote ordinal embedding recovery data" in result.stdout
    assert (tmp_path / "data" / "ordinal_embedding_recovery.csv").exists()


def test_exp32_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp32_ordinal_dimension_selection.py",
            "--true-dims",
            "1",
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

    assert "Wrote ordinal dimension selection data" in result.stdout
    assert (tmp_path / "data" / "ordinal_dimension_selection.csv").exists()


def test_exp33_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp33_noisy_incomplete_order_embedding.py",
            "--n-points",
            "10",
            "--constraint-counts",
            "100",
            "--flip-probabilities",
            "0.0",
            "0.05",
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

    assert "Wrote noisy incomplete order embedding data" in result.stdout
    assert (tmp_path / "data" / "noisy_incomplete_order_embedding.csv").exists()


def test_exp34_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp34_observer_distance_order_embedding.py",
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

    assert "Wrote observer distance-order embedding data" in result.stdout
    assert (tmp_path / "data" / "observer_distance_order_embedding.csv").exists()


def test_exp35_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
