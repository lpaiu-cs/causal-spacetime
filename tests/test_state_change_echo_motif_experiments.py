from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp83_echo_motif_exact_sanity import run_experiment


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


def test_exp83_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp83_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp83_echo_motif_exact_sanity.py",
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

    assert "Wrote echo motif exact sanity" in result.stdout
    assert (tmp_path / "data" / "echo_motif_exact_sanity.csv").exists()


def test_exp84_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp84_planted_echo_motif_recovery.py",
            "--reference-lengths",
            "16",
            "--motif-counts",
            "10",
            "--delay-ranks",
            "2",
            "3",
            "5",
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

    assert "Wrote planted echo motif recovery" in result.stdout
    assert (tmp_path / "data" / "planted_echo_motif_recovery.csv").exists()


def test_exp85_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp85_echo_motif_background_interference.py",
            "--num-systems",
            "5",
            "--max-events",
            "80",
            "--trigger-probability-values",
            "0.05",
            "0.20",
            "--motif-count",
            "6",
            "--delay-ranks",
            "2",
            "3",
            "5",
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

    assert "Wrote echo motif background interference" in result.stdout
    assert (tmp_path / "data" / "echo_motif_background_interference.csv").exists()


def test_exp86_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp86_echo_motif_density_resolution.py",
            "--reference-length",
            "32",
            "--motif-counts",
            "10",
            "30",
            "--delay-rank-sets",
            "small",
            "medium",
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

    assert "Wrote echo motif density resolution" in result.stdout
    assert (tmp_path / "data" / "echo_motif_density_resolution.csv").exists()


def test_exp87_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp87_echo_motif_reference_choice_visibility.py",
            "--num-systems",
            "5",
            "--max-events",
            "80",
            "--trigger-probability",
            "0.20",
            "--motif-count",
            "6",
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

    assert "Wrote echo motif reference choice visibility" in result.stdout
    output = tmp_path / "data" / "echo_motif_reference_choice_visibility.csv"
    assert output.exists()
