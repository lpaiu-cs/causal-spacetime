from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp70_observer_chain_exact_sanity import run_experiment


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


def test_exp70_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp70_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp70_observer_chain_exact_sanity.py",
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

    assert "Wrote observer-chain exact sanity data" in result.stdout
    assert (tmp_path / "data" / "observer_chain_exact_sanity.csv").exists()


def test_exp71_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp71_observer_chain_candidate_ranking.py",
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

    assert "Wrote observer-chain candidate ranking" in result.stdout
    assert (tmp_path / "data" / "observer_chain_candidate_ranking.csv").exists()


def test_exp72_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp72_observer_chain_coverage_vs_trigger_probability.py",
            "--num-systems",
            "3",
            "--max-events",
            "30",
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

    assert "Wrote observer-chain coverage data" in result.stdout
    output = tmp_path / "data" / "observer_chain_coverage_vs_trigger_probability.csv"
    assert output.exists()


def test_exp73_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp73_observer_chain_interval_profile.py",
            "--num-systems",
            "3",
            "--max-events",
            "30",
            "--trigger-probability",
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

    assert "Wrote observer-chain interval profile data" in result.stdout
    assert (tmp_path / "data" / "observer_chain_interval_profile.csv").exists()
