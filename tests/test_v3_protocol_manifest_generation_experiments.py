from __future__ import annotations

# ruff: noqa: E402,I001

import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.abspath("experiments"))

from causal_spacetime_lab.state_change_manifest_v3_protocol_patch import (
    default_v3_protocol_invariant_family_patches,
)
from causal_spacetime_lab.state_change_manifest_v3_protocol_preregistration import (
    build_v3_protocol_patched_preregistration,
    write_v3_protocol_patched_preregistration,
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


def _write_m40_prereg(output_dir: Path) -> None:
    spec = build_v3_protocol_patched_preregistration(
        output_dir / "remediation" / "v3_preregistration_spec_m39.json",
        default_v3_protocol_invariant_family_patches(),
    )
    write_v3_protocol_patched_preregistration(
        spec,
        output_dir / "remediation" / "v3_protocol_patched_preregistration_m40.json",
    )


@pytest.fixture(scope="module")
def m41_output_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    output_dir = tmp_path_factory.mktemp("m41_outputs")
    _write_m40_prereg(output_dir)
    commands = [
        (
            "exp228_v3_protocol_manifest_generation.py",
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
            "exp229_v3_protocol_manifest_metadata_audit.py",
            "--output-dir",
            str(output_dir),
        ),
        (
            "exp230_v3_protocol_family_fit_diagnostics.py",
            "--manifest-dir",
            str(output_dir / "manifests_v3"),
            "--steps",
            "50",
            "--restarts",
            "1",
        ),
        (
            "exp231_v3_protocol_null_taxonomy_diagnostics.py",
            "--manifest-dir",
            str(output_dir / "manifests_v3"),
            "--steps",
            "50",
            "--restarts",
            "1",
            "--null-repetitions",
            "1",
        ),
        (
            "exp232_v3_protocol_stricter_criteria_diagnostics.py",
            "--manifest-dir",
            str(output_dir / "manifests_v3"),
            "--steps",
            "50",
            "--restarts",
            "1",
        ),
        (
            "exp233_v3_protocol_failed_accounting.py",
            "--manifest-dir",
            str(output_dir / "manifests_v3"),
        ),
        (
            "exp234_v3_protocol_coverage_metrics.py",
            "--manifest-dir",
            str(output_dir / "manifests_v3"),
        ),
        (
            "exp235_v3_protocol_restart_latent_order_stability.py",
            "--manifest-dir",
            str(output_dir / "manifests_v3"),
            "--steps",
            "50",
            "--restart-count",
            "2",
        ),
        (
            "exp236_v3_protocol_no_retuning_audit.py",
            "--output-dir",
            str(output_dir),
        ),
        (
            "exp237_v3_protocol_required_metric_aggregation.py",
            "--output-dir",
            str(output_dir),
        ),
        (
            "exp238_v3_protocol_diagnostic_complete_bundle_report.py",
            "--output-dir",
            str(output_dir),
        ),
        (
            "exp239_v3_protocol_manifest_generation_final_sanity.py",
            "--output-dir",
            str(output_dir),
        ),
    ]
    for command in commands:
        result = _run(*command)
        assert result.returncode == 0, result.stderr
    return output_dir


def test_exp224_exact_sanity(m40_prereg_path) -> None:
    from experiments.exp224_v3_protocol_execution_spec_exact_sanity import (
        run_experiment,
    )

    rows = run_experiment(m40_prereg_path.parents[1])

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp225_parameter_sanity() -> None:
    from experiments.exp225_v3_protocol_parameter_completeness_sanity import (
        run_experiment,
    )

    rows = run_experiment()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp226_provenance_sanity() -> None:
    from experiments.exp226_handoff_provenance_exact_sanity import run_experiment

    rows = run_experiment()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp227_template_sanity(m40_prereg_path) -> None:
    from experiments.exp227_top_down_handoff_template_exact_sanity import (
        run_experiment,
    )

    rows = run_experiment()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp228_cli_smoke(m41_output_dir: Path) -> None:
    assert (m41_output_dir / "data" / "v3_protocol_manifest_generation.csv").exists()
    assert list((m41_output_dir / "manifests_v3").glob("*.json"))


def test_exp229_metadata_audit(m41_output_dir: Path) -> None:
    assert (
        m41_output_dir / "data" / "v3_protocol_manifest_metadata_audit.csv"
    ).exists()


def test_exp230_cli_smoke(m41_output_dir: Path) -> None:
    assert Path("outputs/data/v3_protocol_manifest_family_fit_summary.csv").exists()


def test_exp231_cli_smoke(m41_output_dir: Path) -> None:
    assert (
        Path("outputs/data/v3_protocol_manifest_family_null_taxonomy.csv")
    ).exists()


def test_exp232_cli_smoke(m41_output_dir: Path) -> None:
    assert (
        Path("outputs/data/v3_protocol_manifest_family_stricter_criteria.csv")
    ).exists()


def test_exp233_cli_smoke(m41_output_dir: Path) -> None:
    assert (
        Path("outputs/data/v3_protocol_manifest_family_failed_accounting.csv")
    ).exists()


def test_exp234_cli_smoke(m41_output_dir: Path) -> None:
    assert (
        Path("outputs/data/v3_protocol_manifest_family_coverage_metrics.csv")
    ).exists()


def test_exp235_cli_smoke(m41_output_dir: Path) -> None:
    assert (
        Path("outputs/data/v3_protocol_manifest_family_restart_stability.csv")
    ).exists()
    assert (
        Path("outputs/data/v3_protocol_manifest_family_latent_order_stability.csv")
    ).exists()


def test_exp236_no_retuning_audit(m41_output_dir: Path) -> None:
    assert (m41_output_dir / "data" / "v3_protocol_no_retuning_audit.csv").exists()


def test_exp237_aggregation(m41_output_dir: Path) -> None:
    assert (
        m41_output_dir / "data" / "v3_protocol_cross_family_robustness_metrics.csv"
    ).exists()


def test_exp238_bundle_report(m41_output_dir: Path) -> None:
    assert (
        m41_output_dir / "data" / "v3_protocol_diagnostic_complete_bundle_report.csv"
    ).exists()


def test_exp239_final_sanity(m41_output_dir: Path) -> None:
    assert (
        m41_output_dir / "data" / "v3_protocol_manifest_generation_final_sanity.csv"
    ).exists()
