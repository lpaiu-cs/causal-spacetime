from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp88_echo_shortcut_exact_sanity import run_experiment as run_exp88
from experiments.exp93_echo_interference_exact_sanity import run_experiment as run_exp93


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


def test_exp88_exact_sanity_rows_pass() -> None:
    rows = run_exp88()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp89_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp89_echo_shortcut_injection_sweep.py",
            "--reference-length",
            "16",
            "--motif-count",
            "6",
            "--delay-ranks",
            "3",
            "5",
            "--shortcut-probabilities",
            "0.0",
            "0.3",
            "--shortcut-modes",
            "target_to_early_reference",
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

    assert "Wrote echo shortcut injection sweep" in result.stdout
    assert (tmp_path / "data" / "echo_shortcut_injection_sweep.csv").exists()


def test_exp90_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp90_echo_background_edge_perturbation.py",
            "--reference-length",
            "16",
            "--motif-count",
            "6",
            "--delay-ranks",
            "3",
            "5",
            "--edge-probabilities",
            "0.0",
            "0.01",
            "--repetitions",
            "1",
            "--max-edges",
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

    assert "Wrote echo background edge perturbation" in result.stdout
    assert (tmp_path / "data" / "echo_background_edge_perturbation.csv").exists()


def test_exp91_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp91_echo_motif_path_length_robustness.py",
            "--reference-length",
            "16",
            "--motif-count",
            "6",
            "--delay-ranks",
            "5",
            "--outward-steps-values",
            "0",
            "1",
            "--return-steps-values",
            "0",
            "1",
            "--shortcut-probability",
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

    assert "Wrote echo motif path-length robustness" in result.stdout
    assert (tmp_path / "data" / "echo_motif_path_length_robustness.csv").exists()


def test_exp92_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp92_echo_shortcut_reference_dependence.py",
            "--num-systems",
            "5",
            "--max-events",
            "80",
            "--trigger-probability",
            "0.20",
            "--motif-count",
            "6",
            "--shortcut-probability",
            "0.3",
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

    assert "Wrote echo shortcut reference dependence" in result.stdout
    assert (tmp_path / "data" / "echo_shortcut_reference_dependence.csv").exists()


def test_exp93_exact_sanity_rows_pass() -> None:
    rows = run_exp93()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp93_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp93_echo_interference_exact_sanity.py",
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

    assert "Wrote echo interference exact sanity" in result.stdout
    assert (tmp_path / "data" / "echo_interference_exact_sanity.csv").exists()
