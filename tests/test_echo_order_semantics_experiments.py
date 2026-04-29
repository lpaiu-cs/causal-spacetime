from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp107_echo_order_semantics_exact_sanity import (
    run_experiment as run_exp107,
)
from experiments.exp111_echo_terminology_audit import (
    run_experiment as run_exp111,
)
from experiments.exp112_response_representability_exact_sanity import (
    run_experiment as run_exp112,
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


def test_exp107_exact_sanity() -> None:
    rows = run_exp107()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp108_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp108_gated_echo_protocol_comparison.py",
            "--reference-length",
            "24",
            "--motif-count",
            "8",
            "--delay-ranks",
            "5",
            "8",
            "--shortcut-probability-values",
            "0.0",
            "0.5",
            "--gate-delay-ranks",
            "1",
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

    assert "Wrote gated echo protocol comparison" in result.stdout


def test_exp109_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp109_echo_delay_spacing_tie_resolution.py",
            "--reference-length-values",
            "32",
            "--emission-position",
            "4",
            "--layer-sets",
            "compact",
            "medium",
            "--targets-per-layer",
            "3",
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

    assert "Wrote echo delay spacing tie resolution" in result.stdout


def test_exp110_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp110_response_order_scalar_representability.py",
            "--reference-length",
            "32",
            "--emission-position",
            "4",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--shortcut-probabilities",
            "0.0",
            "0.5",
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

    assert "Wrote response-order scalar representability" in result.stdout


def test_exp111_audit() -> None:
    rows = run_exp111()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp112_exact_sanity() -> None:
    rows = run_exp112()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
