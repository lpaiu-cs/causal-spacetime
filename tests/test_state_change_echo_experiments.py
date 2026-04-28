from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp78_state_change_echo_exact_sanity import run_experiment


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


def test_exp78_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp78_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp78_state_change_echo_exact_sanity.py",
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

    assert "Wrote state-change echo exact sanity" in result.stdout
    assert (tmp_path / "data" / "state_change_echo_exact_sanity.csv").exists()


def test_exp79_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp79_state_change_echo_order_diagnostics.py",
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
            "--emission-count",
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

    assert "Wrote state-change echo-order diagnostics" in result.stdout
    assert (tmp_path / "data" / "state_change_echo_order_diagnostics.csv").exists()


def test_exp80_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp80_state_change_echo_reference_dependence.py",
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

    assert "Wrote state-change echo reference dependence" in result.stdout
    assert (tmp_path / "data" / "state_change_echo_reference_dependence.csv").exists()


def test_exp81_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp81_state_change_echo_emission_sensitivity.py",
            "--num-systems",
            "3",
            "--max-events",
            "40",
            "--trigger-probability",
            "0.2",
            "--repetitions",
            "1",
            "--emission-count",
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

    assert "Wrote state-change echo emission sensitivity" in result.stdout
    assert (tmp_path / "data" / "state_change_echo_emission_sensitivity.csv").exists()


def test_exp82_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp82_state_change_echo_coverage_vs_trigger_density.py",
            "--num-systems",
            "3",
            "--max-events",
            "40",
            "--trigger-probability-values",
            "0.1",
            "0.3",
            "--repetitions",
            "1",
            "--emission-count",
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

    assert "Wrote state-change echo coverage" in result.stdout
    output = tmp_path / "data" / "state_change_echo_coverage_vs_trigger_density.csv"
    assert output.exists()
