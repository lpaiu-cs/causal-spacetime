from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp94_echo_coarse_graining_exact_sanity import (
    run_experiment as run_exp94,
)
from experiments.exp99_echo_return_spectrum_stability_exact_sanity import (
    run_experiment as run_exp99,
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


def test_exp94_exact_sanity() -> None:
    rows = run_exp94()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp95_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp95_echo_event_thinning_stability.py",
            "--reference-length",
            "16",
            "--motif-count",
            "6",
            "--delay-ranks",
            "3",
            "5",
            "--keep-probabilities",
            "1.0",
            "0.4",
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

    assert "Wrote echo event thinning stability" in result.stdout
    assert (tmp_path / "data" / "echo_event_thinning_stability.csv").exists()


def test_exp96_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp96_echo_reference_subsampling_resolution.py",
            "--reference-length",
            "32",
            "--motif-count",
            "10",
            "--delay-ranks",
            "2",
            "3",
            "5",
            "--strides",
            "1",
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

    assert "Wrote echo reference subsampling resolution" in result.stdout
    assert (tmp_path / "data" / "echo_reference_subsampling_resolution.csv").exists()


def test_exp97_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp97_echo_edge_thinning_fragility.py",
            "--reference-length",
            "16",
            "--motif-count",
            "6",
            "--delay-ranks",
            "3",
            "5",
            "--removal-probabilities",
            "0.0",
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

    assert "Wrote echo edge thinning fragility" in result.stdout
    assert (tmp_path / "data" / "echo_edge_thinning_fragility.csv").exists()


def test_exp98_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp98_echo_shortcut_classification_under_coarse_graining.py",
            "--reference-length",
            "16",
            "--motif-count",
            "6",
            "--delay-ranks",
            "3",
            "5",
            "--shortcut-probability",
            "0.3",
            "--event-keep-probability",
            "0.4",
            "--edge-removal-probability",
            "0.15",
            "--reference-strides",
            "1",
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

    assert "Wrote echo shortcut classification under coarse-graining" in result.stdout
    output = tmp_path / "data"
    assert (output / "echo_shortcut_classification_under_coarse_graining.csv").exists()


def test_exp99_exact_sanity() -> None:
    rows = run_exp99()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
