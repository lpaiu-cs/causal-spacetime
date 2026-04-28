from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp67_state_change_exact_sanity import run_experiment


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


def test_exp67_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp67_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp67_state_change_exact_sanity.py",
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

    assert "Wrote state-change exact sanity data" in result.stdout
    assert (tmp_path / "data" / "state_change_exact_sanity.csv").exists()


def test_exp68_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp68_state_change_toy_model.py",
            "--num-systems-values",
            "3",
            "--max-events-values",
            "20",
            "--trigger-probability-values",
            "0.2",
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

    assert "Wrote state-change toy-model summary" in result.stdout
    assert (tmp_path / "data" / "state_change_toy_model_summary.csv").exists()


def test_exp69_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp69_state_change_observer_chain_diagnostic.py",
            "--num-systems",
            "3",
            "--max-events",
            "20",
            "--trigger-probability",
            "0.2",
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

    assert "Wrote state-change observer chain diagnostic" in result.stdout
    output = tmp_path / "data" / "state_change_observer_chain_diagnostic.csv"
    assert output.exists()
