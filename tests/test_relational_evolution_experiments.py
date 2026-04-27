from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp59_relational_evolution_exact_sanity import run_experiment


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


def test_exp54_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp54_predicate_definability_table.py",
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

    assert "Wrote predicate definability table" in result.stdout
    assert (tmp_path / "data" / "predicate_definability_table.csv").exists()


def test_exp55_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp55_relational_shape_history_without_transport.py",
            "--slice-count",
            "6",
            "--object-count",
            "4",
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

    assert "Wrote relational shape history data" in result.stdout
    output = tmp_path / "data" / "relational_shape_history_without_transport.csv"
    assert output.exists()


def test_exp56_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp56_relational_history_gauge_invariance.py",
            "--slice-count",
            "6",
            "--object-count",
            "4",
            "--repetitions",
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

    assert "Wrote relational history gauge data" in result.stdout
    assert (tmp_path / "data" / "relational_history_gauge_invariance.csv").exists()


def test_exp57_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp57_observer_slice_relational_evolution.py",
            "--slice-count",
            "6",
            "--object-count",
            "4",
            "--tick-count",
            "32",
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

    assert "Wrote observer-slice relational evolution data" in result.stdout
    assert (tmp_path / "data" / "observer_slice_relational_evolution.csv").exists()


def test_exp58_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp58_relational_invariants_vs_velocity.py",
            "--slice-count",
            "6",
            "--object-count",
            "4",
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

    assert "Wrote relational invariants versus velocity data" in result.stdout
    assert (tmp_path / "data" / "relational_invariants_vs_velocity.csv").exists()


def test_exp59_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
