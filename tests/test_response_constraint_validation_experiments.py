from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp126_response_constraint_pool_exact_sanity import (
    run_experiment as run_exp126,
)
from experiments.exp132_response_constraint_validation_exact_sanity import (
    run_experiment as run_exp132,
)


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    paths = [os.path.abspath("src"), os.path.abspath("experiments")]
    if env.get("PYTHONPATH"):
        paths.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


def _skip_without_matplotlib() -> None:
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")


def test_exp126_exact_sanity() -> None:
    rows = run_exp126()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp127_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp127_response_constraint_heldout_protocol_validation.py",
            "--reference-length",
            "48",
            "--emission-positions",
            "6",
            "12",
            "18",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--train-fractions",
            "0.5",
            "--max-constraints",
            "100",
            "--min-margins",
            "0.0",
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
    assert "Wrote response constraint held-out validation" in result.stdout


def test_exp128_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp128_response_constraint_bootstrap_stability.py",
            "--reference-length",
            "48",
            "--emission-positions",
            "6",
            "12",
            "18",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--max-constraints",
            "100",
            "--bootstrap-count",
            "3",
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
    assert "Wrote response constraint bootstrap stability" in result.stdout


def test_exp129_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp129_response_constraint_null_separation.py",
            "--reference-length",
            "48",
            "--emission-positions",
            "6",
            "12",
            "18",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--max-constraints",
            "100",
            "--null-repetitions",
            "2",
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
    assert "Wrote response constraint null separation" in result.stdout


def test_exp130_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp130_response_constraint_pool_coverage.py",
            "--target-counts",
            "12",
            "--protocol-counts",
            "3",
            "--unique-rank-count",
            "4",
            "--max-constraints-values",
            "50",
            "--min-margins",
            "0.0",
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
    assert "Wrote response constraint pool coverage" in result.stdout


def test_exp131_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp131_response_constraint_validation_gate_summary.py",
            "--reference-length",
            "48",
            "--emission-positions",
            "6",
            "12",
            "18",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--max-constraints",
            "100",
            "--min-margin-values",
            "0.0",
            "--repetitions",
            "1",
            "--bootstrap-count",
            "3",
            "--null-repetitions",
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
    assert "Wrote response constraint validation gate summary" in result.stdout


def test_exp132_exact_sanity() -> None:
    rows = run_exp132()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
