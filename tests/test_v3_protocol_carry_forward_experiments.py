from __future__ import annotations

# ruff: noqa: E402,I001

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("experiments"))

from experiments.exp240_v3_protocol_carry_forward_exact_sanity import (
    run_experiment as run_exp240,
)

from tests.v3_protocol_test_helpers import (
    family_from_manifest_dir,
    write_v3_protocol_bundle,
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


def _run_script(script: str, output_dir: Path, *args: str) -> None:
    result = subprocess.run(
        [_python_executable(), f"experiments/{script}", *args],
        cwd=Path.cwd(),
        env=_env_with_src(),
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert output_dir.exists()


def _prepare_bundle(tmp_path: Path, m41_manifest_dir: Path) -> Path:
    write_v3_protocol_bundle(tmp_path, family_from_manifest_dir(m41_manifest_dir))
    return tmp_path


def test_exp240_exact_sanity() -> None:
    rows = run_exp240()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp241_input_audit(tmp_path: Path, m41_manifest_dir: Path) -> None:
    bundle = _prepare_bundle(tmp_path, m41_manifest_dir)
    _run_script(
        "exp241_v3_protocol_bundle_input_audit.py",
        bundle,
        "--output-dir",
        str(bundle),
    )

    assert (bundle / "data" / "v3_protocol_bundle_input_audit.csv").exists()


def test_exp242_precondition_audit(tmp_path: Path, m41_manifest_dir: Path) -> None:
    bundle = _prepare_bundle(tmp_path, m41_manifest_dir)
    _run_script(
        "exp242_v3_protocol_precondition_audit.py",
        bundle,
        "--output-dir",
        str(bundle),
        "--manifest-dir",
        str(m41_manifest_dir),
    )

    assert (bundle / "data" / "v3_protocol_precondition_audit.csv").exists()


def test_exp243_cli_smoke(tmp_path: Path, m41_manifest_dir: Path) -> None:
    bundle = _prepare_bundle(tmp_path, m41_manifest_dir)
    _run_script(
        "exp243_v3_protocol_carry_forward_decision.py",
        bundle,
        "--output-dir",
        str(bundle),
        "--manifest-dir",
        str(m41_manifest_dir),
    )

    assert (
        bundle / "data" / "v3_protocol_cross_family_robustness_decisions.csv"
    ).exists()


def test_exp244_sensitivity(tmp_path: Path, m41_manifest_dir: Path) -> None:
    bundle = _prepare_bundle(tmp_path, m41_manifest_dir)
    # The script uses output_dir/manifests_v3 for preconditions.
    link_target = bundle / "manifests_v3"
    link_target.mkdir(exist_ok=True)
    if link_target != m41_manifest_dir:
        for path in m41_manifest_dir.glob("*.json"):
            (link_target / path.name).write_text(path.read_text(), encoding="utf-8")
    _run_script(
        "exp244_v3_protocol_threshold_sensitivity.py",
        bundle,
        "--output-dir",
        str(bundle),
        "--heldout-thresholds",
        "0.2",
        "--null-gap-thresholds",
        "0.1",
        "--stricter-pass-thresholds",
        "0.5",
    )

    assert (bundle / "data" / "v3_protocol_threshold_sensitivity.csv").exists()


def test_exp245_registry_export(tmp_path: Path, m41_manifest_dir: Path) -> None:
    test_exp243_cli_smoke(tmp_path, m41_manifest_dir)
    _run_script(
        "exp245_v3_protocol_registry_export.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
        "--manifest-dir",
        str(m41_manifest_dir),
    )

    assert (
        tmp_path
        / "carry_forward_v3_protocol"
        / "carry_forward_registry_v3_protocol.json"
    ).exists()


def test_exp246_handoff_plan(tmp_path: Path, m41_manifest_dir: Path) -> None:
    test_exp245_registry_export(tmp_path, m41_manifest_dir)
    _run_script(
        "exp246_v3_protocol_stress_test_handoff_plan.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (tmp_path / "data" / "v3_protocol_stress_test_handoff_plan.csv").exists()


def test_exp247_accounting(tmp_path: Path, m41_manifest_dir: Path) -> None:
    test_exp243_cli_smoke(tmp_path, m41_manifest_dir)
    _run_script(
        "exp247_v3_protocol_failed_provisional_accounting.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (
        tmp_path / "data" / "v3_protocol_failed_provisional_accounting.csv"
    ).exists()


def test_exp248_failure_decomposition(tmp_path: Path, m41_manifest_dir: Path) -> None:
    bundle = _prepare_bundle(tmp_path, m41_manifest_dir)
    _run_script(
        "exp248_v3_protocol_failure_decomposition.py",
        bundle,
        "--output-dir",
        str(bundle),
        "--manifest-dir",
        str(m41_manifest_dir),
    )

    assert (
        bundle / "data" / "v3_protocol_carry_forward_failure_decomposition.csv"
    ).exists()


def test_exp249_no_retuning_audit(tmp_path: Path, m41_manifest_dir: Path) -> None:
    test_exp244_sensitivity(tmp_path, m41_manifest_dir)
    _run_script(
        "exp249_v3_protocol_no_retuning_no_refit_audit.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (
        tmp_path / "data" / "v3_protocol_carry_forward_no_retuning_audit.csv"
    ).exists()


def test_exp250_report_card(tmp_path: Path, m41_manifest_dir: Path) -> None:
    test_exp246_handoff_plan(tmp_path, m41_manifest_dir)
    _run_script(
        "exp250_v3_protocol_carry_forward_report_card.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (tmp_path / "data" / "v3_protocol_carry_forward_report_card.csv").exists()


def test_exp251_final_sanity(tmp_path: Path, m41_manifest_dir: Path) -> None:
    test_exp250_report_card(tmp_path, m41_manifest_dir)
    _run_script(
        "exp251_v3_protocol_carry_forward_final_sanity.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (
        tmp_path / "data" / "v3_protocol_carry_forward_final_sanity.csv"
    ).exists()
