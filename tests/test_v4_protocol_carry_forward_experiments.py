from __future__ import annotations

# ruff: noqa: E402,I001

import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.abspath("experiments"))

from experiments.exp279_v4_protocol_carry_forward_exact_sanity import (
    run_experiment as run_exp279,
)


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    paths = [os.path.abspath("src"), os.path.abspath("experiments")]
    if env.get("PYTHONPATH"):
        paths.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


def _python_executable() -> str:
    venv_python = Path(".venv/bin/python")
    return str(venv_python) if venv_python.exists() else sys.executable


def _run(script: str, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [_python_executable(), f"experiments/{script}", *args],
        cwd=Path.cwd(),
        env=_env_with_src(),
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.fixture(scope="module")
def m45_output_dir() -> Path:
    output_dir = Path("outputs")
    commands = [
        ("exp280_v4_protocol_bundle_input_audit.py", "--output-dir", "outputs"),
        (
            "exp281_v4_protocol_precondition_audit.py",
            "--output-dir",
            "outputs",
            "--manifest-dir",
            "outputs/manifests_v4",
        ),
        (
            "exp282_v4_protocol_carry_forward_decision.py",
            "--output-dir",
            "outputs",
            "--manifest-dir",
            "outputs/manifests_v4",
        ),
        (
            "exp283_v4_protocol_threshold_sensitivity.py",
            "--output-dir",
            "outputs",
            "--heldout-thresholds",
            "0.20",
            "0.25",
            "--null-gap-thresholds",
            "0.05",
            "0.10",
            "--stricter-pass-thresholds",
            "0.25",
            "0.50",
        ),
        (
            "exp284_v4_protocol_registry_export.py",
            "--output-dir",
            "outputs",
            "--manifest-dir",
            "outputs/manifests_v4",
        ),
        ("exp285_v4_protocol_stress_test_handoff_plan.py", "--output-dir", "outputs"),
        (
            "exp286_v4_protocol_failed_provisional_accounting.py",
            "--output-dir",
            "outputs",
        ),
        (
            "exp287_v4_protocol_failure_decomposition.py",
            "--output-dir",
            "outputs",
            "--manifest-dir",
            "outputs/manifests_v4",
        ),
        (
            "exp288_v4_protocol_no_retuning_no_refit_audit.py",
            "--output-dir",
            "outputs",
        ),
        ("exp289_v4_protocol_carry_forward_report_card.py", "--output-dir", "outputs"),
        ("exp290_v4_protocol_carry_forward_final_sanity.py", "--output-dir", "outputs"),
    ]
    for command in commands:
        result = _run(*command)
        assert result.returncode == 0, result.stderr
    return output_dir


def test_exp279_exact_sanity() -> None:
    rows = run_exp279()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp280_input_audit(m45_output_dir: Path) -> None:
    assert (m45_output_dir / "data" / "v4_protocol_bundle_input_audit.csv").exists()


def test_exp281_precondition_audit(m45_output_dir: Path) -> None:
    assert (m45_output_dir / "data" / "v4_protocol_precondition_audit.csv").exists()


def test_exp282_cli_smoke(m45_output_dir: Path) -> None:
    assert (
        m45_output_dir
        / "data"
        / "v4_protocol_cross_family_robustness_decisions.csv"
    ).exists()


def test_exp283_sensitivity(m45_output_dir: Path) -> None:
    assert (m45_output_dir / "data" / "v4_protocol_threshold_sensitivity.csv").exists()


def test_exp284_registry_export(m45_output_dir: Path) -> None:
    assert (
        m45_output_dir
        / "carry_forward_v4_protocol"
        / "carry_forward_registry_v4_protocol.json"
    ).exists()


def test_exp285_handoff_plan(m45_output_dir: Path) -> None:
    assert (
        m45_output_dir / "data" / "v4_protocol_stress_test_handoff_plan.csv"
    ).exists()


def test_exp286_accounting(m45_output_dir: Path) -> None:
    assert (
        m45_output_dir / "data" / "v4_protocol_failed_provisional_accounting.csv"
    ).exists()


def test_exp287_failure_decomposition(m45_output_dir: Path) -> None:
    assert (
        m45_output_dir
        / "data"
        / "v4_protocol_carry_forward_failure_decomposition.csv"
    ).exists()


def test_exp288_no_retuning_audit(m45_output_dir: Path) -> None:
    assert (
        m45_output_dir
        / "data"
        / "v4_protocol_carry_forward_no_retuning_audit.csv"
    ).exists()


def test_exp289_report_card(m45_output_dir: Path) -> None:
    assert (
        m45_output_dir / "data" / "v4_protocol_carry_forward_report_card.csv"
    ).exists()


def test_exp290_final_sanity(m45_output_dir: Path) -> None:
    assert (
        m45_output_dir / "data" / "v4_protocol_carry_forward_final_sanity.csv"
    ).exists()
