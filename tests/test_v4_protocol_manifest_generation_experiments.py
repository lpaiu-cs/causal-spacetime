from __future__ import annotations

# ruff: noqa: E402,I001

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.abspath("experiments"))

from causal_spacetime_lab.state_change_manifest_v4_design import (
    default_v4_protocol_family_designs,
)
from causal_spacetime_lab.state_change_manifest_v4_preregistration import (
    build_v4_protocol_preregistration_spec,
    write_v4_protocol_preregistration_spec,
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


def _write_m43_prereg(output_dir: Path) -> None:
    spec = build_v4_protocol_preregistration_spec(
        default_v4_protocol_family_designs()
    )
    write_v4_protocol_preregistration_spec(
        spec,
        output_dir / "remediation" / "v4_protocol_preregistration_spec_m43.json",
    )


def _write_m43_no_execution_audit(output_dir: Path) -> None:
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "v4_no_execution_audit.csv").write_text(
        "check,passed\nm43_no_execution_audit_not_retroactive,1.0\n",
        encoding="utf-8",
    )


def _copy_root_csv_to_output(output_dir: Path, name: str) -> None:
    source = Path("outputs/data") / name
    if source.exists():
        target = output_dir / "data" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)


@pytest.fixture(scope="module")
def m44_output_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    output_dir = tmp_path_factory.mktemp("m44_outputs")
    _write_m43_prereg(output_dir)
    _write_m43_no_execution_audit(output_dir)
    commands = [
        (
            "exp267_v4_protocol_manifest_generation.py",
            "--output-dir",
            str(output_dir),
            "--max-constraints",
            "80",
            "--bootstrap-count",
            "1",
            "--null-repetitions",
            "1",
        ),
        (
            "exp268_v4_protocol_manifest_metadata_audit.py",
            "--output-dir",
            str(output_dir),
        ),
        (
            "exp269_v4_protocol_family_fit_diagnostics.py",
            "--manifest-dir",
            str(output_dir / "manifests_v4"),
            "--steps",
            "40",
            "--restarts",
            "1",
        ),
        (
            "exp270_v4_protocol_null_taxonomy_diagnostics.py",
            "--manifest-dir",
            str(output_dir / "manifests_v4"),
            "--steps",
            "40",
            "--restarts",
            "1",
            "--null-repetitions",
            "1",
        ),
        (
            "exp271_v4_protocol_stricter_criteria_diagnostics.py",
            "--manifest-dir",
            str(output_dir / "manifests_v4"),
            "--steps",
            "40",
            "--restarts",
            "1",
        ),
        (
            "exp272_v4_protocol_failed_accounting.py",
            "--manifest-dir",
            str(output_dir / "manifests_v4"),
        ),
        (
            "exp273_v4_protocol_coverage_metrics.py",
            "--manifest-dir",
            str(output_dir / "manifests_v4"),
        ),
        (
            "exp274_v4_protocol_restart_latent_order_stability.py",
            "--manifest-dir",
            str(output_dir / "manifests_v4"),
            "--steps",
            "40",
            "--restart-count",
            "2",
        ),
    ]
    for command in commands:
        result = _run(*command)
        assert result.returncode == 0, result.stderr
    for name in [
        "v4_protocol_manifest_family_fit_summary.csv",
        "v4_protocol_manifest_family_null_taxonomy.csv",
        "v4_protocol_manifest_family_stricter_criteria.csv",
        "v4_protocol_manifest_family_failed_accounting.csv",
        "v4_protocol_manifest_family_coverage_metrics.csv",
        "v4_protocol_manifest_family_restart_stability.csv",
        "v4_protocol_manifest_family_latent_order_stability.csv",
    ]:
        _copy_root_csv_to_output(output_dir, name)
    for command in [
        ("exp275_v4_protocol_no_retuning_audit.py", "--output-dir", str(output_dir)),
        (
            "exp276_v4_protocol_required_metric_aggregation.py",
            "--output-dir",
            str(output_dir),
        ),
        (
            "exp277_v4_protocol_diagnostic_complete_bundle_report.py",
            "--output-dir",
            str(output_dir),
        ),
        (
            "exp278_v4_protocol_manifest_generation_final_sanity.py",
            "--output-dir",
            str(output_dir),
        ),
    ]:
        result = _run(*command)
        assert result.returncode == 0, result.stderr
    return output_dir


def test_exp265_exact_sanity() -> None:
    from experiments.exp265_v4_protocol_execution_spec_exact_sanity import (
        run_experiment,
    )

    rows = run_experiment()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp266_mapping_sanity() -> None:
    from experiments.exp266_v4_protocol_mapping_exact_sanity import run_experiment

    rows = run_experiment()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp267_cli_smoke(m44_output_dir: Path) -> None:
    assert (m44_output_dir / "data" / "v4_protocol_manifest_generation.csv").exists()
    assert list((m44_output_dir / "manifests_v4").glob("*.json"))


def test_exp268_metadata_audit(m44_output_dir: Path) -> None:
    assert (
        m44_output_dir / "data" / "v4_protocol_manifest_metadata_audit.csv"
    ).exists()


def test_exp269_cli_smoke(m44_output_dir: Path) -> None:
    assert (
        m44_output_dir / "data" / "v4_protocol_manifest_family_fit_summary.csv"
    ).exists()


def test_exp270_cli_smoke(m44_output_dir: Path) -> None:
    assert (
        m44_output_dir / "data" / "v4_protocol_manifest_family_null_taxonomy.csv"
    ).exists()


def test_exp271_cli_smoke(m44_output_dir: Path) -> None:
    assert (
        m44_output_dir / "data" / "v4_protocol_manifest_family_stricter_criteria.csv"
    ).exists()


def test_exp272_failed_accounting(m44_output_dir: Path) -> None:
    assert (
        m44_output_dir
        / "data"
        / "v4_protocol_manifest_family_failed_accounting.csv"
    ).exists()


def test_exp273_coverage_metrics(m44_output_dir: Path) -> None:
    assert (
        m44_output_dir / "data" / "v4_protocol_manifest_family_coverage_metrics.csv"
    ).exists()


def test_exp274_stability_metrics(m44_output_dir: Path) -> None:
    assert (
        m44_output_dir / "data" / "v4_protocol_manifest_family_restart_stability.csv"
    ).exists()
    assert (
        m44_output_dir
        / "data"
        / "v4_protocol_manifest_family_latent_order_stability.csv"
    ).exists()


def test_exp275_no_retuning_audit(m44_output_dir: Path) -> None:
    assert (m44_output_dir / "data" / "v4_protocol_no_retuning_audit.csv").exists()


def test_exp276_aggregation(m44_output_dir: Path) -> None:
    assert (
        m44_output_dir / "data" / "v4_protocol_cross_family_robustness_metrics.csv"
    ).exists()


def test_exp277_bundle_report(m44_output_dir: Path) -> None:
    assert (
        m44_output_dir / "data" / "v4_protocol_diagnostic_complete_bundle_report.csv"
    ).exists()


def test_exp278_final_sanity(m44_output_dir: Path) -> None:
    assert (
        m44_output_dir / "data" / "v4_protocol_manifest_generation_final_sanity.csv"
    ).exists()
