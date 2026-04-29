from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp113_response_order_underdetermination_exact_sanity import (
    run_experiment as run_exp113,
)
from experiments.exp116_response_representability_requirement_table import (
    run_experiment as run_exp116,
)
from experiments.exp118_response_profile_requirement_exact_sanity import (
    run_experiment as run_exp118,
)


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    src_path = os.path.abspath("src")
    exp_path = os.path.abspath("experiments")
    paths = [src_path, exp_path]
    if env.get("PYTHONPATH"):
        paths.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


def _skip_without_matplotlib() -> None:
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")


def test_exp113_exact_sanity() -> None:
    rows = run_exp113()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp114_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp114_single_reference_response_order_underdetermination.py",
            "--target-counts",
            "8",
            "--unique-rank-counts",
            "3",
            "--layout-count",
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

    assert "Wrote single-reference underdetermination" in result.stdout


def test_exp115_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp115_multi_reference_response_profile_diagnostics.py",
            "--reference-length",
            "48",
            "--emission-positions",
            "6",
            "12",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
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

    assert "Wrote multi-reference response profile diagnostics" in result.stdout


def test_exp116_table_experiment() -> None:
    rows = run_exp116()

    assert rows
    assert rows[0]["level"] == "scalar_response_rank"


def test_exp117_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp117_response_profile_stability_under_protocol_variants.py",
            "--reference-length",
            "48",
            "--emission-positions",
            "6",
            "12",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--shortcut-probability-values",
            "0.0",
            "0.3",
            "--reference-strides",
            "1",
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

    assert "Wrote response profile stability" in result.stdout


def test_exp118_exact_sanity() -> None:
    rows = run_exp118()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
