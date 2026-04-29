from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp100_echo_response_signature_exact_sanity import (
    run_experiment as run_exp100,
)
from experiments.exp106_echo_response_signature_stability_exact_sanity import (
    run_experiment as run_exp106,
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


def test_exp100_exact_sanity() -> None:
    rows = run_exp100()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp101_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp101_layered_echo_response_order_recovery.py",
            "--reference-length",
            "24",
            "--emission-position",
            "4",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer-values",
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

    assert "Wrote layered echo response order recovery" in result.stdout
    assert (tmp_path / "data" / "layered_echo_response_order_recovery.csv").exists()


def test_exp102_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp102_echo_response_signature_coarse_protocol_stability.py",
            "--reference-length",
            "24",
            "--emission-position",
            "4",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--event-keep-probabilities",
            "1.0",
            "0.5",
            "--reference-strides",
            "1",
            "2",
            "--edge-removal-probabilities",
            "0.0",
            "0.1",
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

    assert "Wrote echo response signature coarse protocol stability" in result.stdout


def test_exp103_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp103_echo_response_shortcut_robust_core.py",
            "--reference-length",
            "24",
            "--emission-position",
            "4",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--shortcut-probabilities",
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

    assert "Wrote echo response shortcut robust core" in result.stdout


def test_exp104_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp104_echo_response_reference_protocol_dependence.py",
            "--num-systems",
            "5",
            "--max-events",
            "120",
            "--trigger-probability",
            "0.20",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
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

    assert "Wrote echo response reference-protocol dependence" in result.stdout


def test_exp105_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp105_echo_response_order_precondition_diagnostics.py",
            "--reference-length",
            "24",
            "--emission-position",
            "4",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--shortcut-probabilities",
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

    assert "Wrote echo response order precondition diagnostics" in result.stdout


def test_exp106_exact_sanity() -> None:
    rows = run_exp106()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)

