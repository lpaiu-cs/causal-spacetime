from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp119_pairwise_response_profile_exact_sanity import (
    run_experiment as run_exp119,
)
from experiments.exp125_pairwise_response_null_admissibility_exact_sanity import (
    run_experiment as run_exp125,
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


def test_exp119_exact_sanity() -> None:
    rows = run_exp119()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp120_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp120_pairwise_response_protocol_comparison.py",
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
    assert "Wrote pairwise response protocol comparison" in result.stdout


def test_exp121_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp121_pairwise_response_null_baselines.py",
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
            "--null-repetitions",
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
    assert "Wrote pairwise response null baselines" in result.stdout


def test_exp122_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp122_pairwise_response_protocol_variant_stability.py",
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
    assert "Wrote pairwise response protocol variant stability" in result.stdout


def test_exp123_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp123_pairwise_response_missing_data_sensitivity.py",
            "--target-count",
            "20",
            "--protocol-count",
            "3",
            "--reachable-probabilities",
            "0.5",
            "1.0",
            "--unique-rank-count",
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
    assert "Wrote pairwise response missing-data sensitivity" in result.stdout


def test_exp124_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp124_pairwise_response_protocol_choice_dependence.py",
            "--target-count",
            "20",
            "--protocol-count",
            "3",
            "--reachable-probability",
            "0.8",
            "--unique-rank-count",
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
    assert "Wrote pairwise response protocol-choice dependence" in result.stdout


def test_exp125_exact_sanity() -> None:
    rows = run_exp125()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
