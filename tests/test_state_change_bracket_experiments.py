from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp77_state_change_reference_bracket_exact_sanity import (
    run_experiment,
)


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


def test_exp77_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp74_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp74_state_change_reference_bracket_diagnostics.py",
            "--num-systems-values",
            "3",
            "--max-events-values",
            "30",
            "--trigger-probability-values",
            "0.2",
            "--repetitions",
            "1",
            "--random-candidate-count",
            "2",
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

    assert "Wrote state-change reference bracket diagnostics" in result.stdout
    output = tmp_path / "data" / "state_change_reference_bracket_diagnostics.csv"
    assert output.exists()


def test_exp75_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp75_state_change_bracket_rank_reference_dependence.py",
            "--num-systems",
            "3",
            "--max-events",
            "40",
            "--trigger-probability-values",
            "0.2",
            "--repetitions",
            "1",
            "--top-k",
            "3",
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

    assert "Wrote state-change bracket-rank reference dependence" in result.stdout
    output = tmp_path / "data" / "state_change_bracket_rank_reference_dependence.csv"
    assert output.exists()


def test_exp76_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp76_state_change_reference_bracket_coverage_vs_trigger_density.py",
            "--num-systems",
            "3",
            "--max-events",
            "40",
            "--trigger-probability-values",
            "0.1",
            "0.3",
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

    assert "Wrote state-change bracket coverage" in result.stdout
    output = (
        tmp_path
        / "data"
        / "state_change_reference_bracket_coverage_vs_trigger_density.csv"
    )
    assert output.exists()


def test_exp77_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp77_state_change_reference_bracket_exact_sanity.py",
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

    assert "Wrote state-change reference bracket exact sanity" in result.stdout
    output = tmp_path / "data" / "state_change_reference_bracket_exact_sanity.csv"
    assert output.exists()
