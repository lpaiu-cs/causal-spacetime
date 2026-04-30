from __future__ import annotations

# ruff: noqa: E402,I001

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("experiments"))

import pytest

from experiments.exp183_v2_manifest_spec_exact_sanity import (
    run_experiment as run_exp183,
)
from experiments.v2_manifest_experiment_helpers import (
    ensure_m36_plan_and_spec_for_tests,
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


@pytest.fixture(scope="module")
def v2_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    output_dir = tmp_path_factory.mktemp("v2_bundle")
    ensure_m36_plan_and_spec_for_tests(output_dir)
    _run_script(
        "exp184_v2_manifest_generation.py",
        output_dir,
        "--output-dir",
        str(output_dir),
        "--max-constraints",
        "80",
        "--bootstrap-count",
        "1",
        "--null-repetitions",
        "1",
    )
    return output_dir


def test_exp183_exact_sanity() -> None:
    rows = run_exp183()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp184_cli_smoke(v2_bundle: Path) -> None:
    assert (v2_bundle / "data" / "v2_manifest_generation.csv").exists()
    assert list((v2_bundle / "manifests_v2").glob("*.json"))


def test_exp185_cli_smoke(v2_bundle: Path) -> None:
    _run_script(
        "exp185_v2_family_fit_diagnostics.py",
        v2_bundle,
        "--manifest-dir",
        str(v2_bundle / "manifests_v2"),
        "--output-dir",
        str(v2_bundle),
        "--dims",
        "1",
        "--steps",
        "20",
        "--restarts",
        "1",
    )

    assert (v2_bundle / "data" / "v2_manifest_family_fit_summary.csv").exists()


def test_exp186_cli_smoke(v2_bundle: Path) -> None:
    _run_script(
        "exp186_v2_null_taxonomy_diagnostics.py",
        v2_bundle,
        "--manifest-dir",
        str(v2_bundle / "manifests_v2"),
        "--output-dir",
        str(v2_bundle),
        "--embedding-dim",
        "1",
        "--null-repetitions",
        "1",
        "--steps",
        "20",
        "--restarts",
        "1",
    )

    assert (v2_bundle / "data" / "v2_manifest_family_null_taxonomy.csv").exists()


def test_exp187_cli_smoke(v2_bundle: Path) -> None:
    _run_script(
        "exp187_v2_stricter_criteria_diagnostics.py",
        v2_bundle,
        "--manifest-dir",
        str(v2_bundle / "manifests_v2"),
        "--output-dir",
        str(v2_bundle),
        "--dims",
        "1",
        "--steps",
        "20",
        "--restarts",
        "1",
    )

    assert (v2_bundle / "data" / "v2_manifest_family_stricter_criteria.csv").exists()


def test_exp188_cli_smoke(v2_bundle: Path) -> None:
    _run_script(
        "exp188_v2_failed_accounting.py",
        v2_bundle,
        "--manifest-dir",
        str(v2_bundle / "manifests_v2"),
        "--output-dir",
        str(v2_bundle),
    )

    assert (v2_bundle / "data" / "v2_manifest_family_failed_accounting.csv").exists()


def test_exp189_cli_smoke(v2_bundle: Path) -> None:
    _run_script(
        "exp189_v2_coverage_metrics.py",
        v2_bundle,
        "--manifest-dir",
        str(v2_bundle / "manifests_v2"),
        "--output-dir",
        str(v2_bundle),
    )

    assert (v2_bundle / "data" / "v2_manifest_family_coverage_metrics.csv").exists()


def test_exp190_cli_smoke(v2_bundle: Path) -> None:
    _run_script(
        "exp190_v2_restart_latent_order_stability.py",
        v2_bundle,
        "--manifest-dir",
        str(v2_bundle / "manifests_v2"),
        "--output-dir",
        str(v2_bundle),
        "--embedding-dim",
        "1",
        "--restart-count",
        "2",
        "--steps",
        "20",
    )

    assert (v2_bundle / "data" / "v2_manifest_family_restart_stability.csv").exists()
    assert (
        v2_bundle / "data" / "v2_manifest_family_latent_order_stability.csv"
    ).exists()


def test_exp191_no_retuning_audit(v2_bundle: Path) -> None:
    _run_script(
        "exp191_v2_no_retuning_audit.py",
        v2_bundle,
        "--output-dir",
        str(v2_bundle),
    )

    assert (v2_bundle / "data" / "v2_no_retuning_audit.csv").exists()


def test_exp192_aggregation(v2_bundle: Path) -> None:
    test_exp185_cli_smoke(v2_bundle)
    test_exp186_cli_smoke(v2_bundle)
    test_exp187_cli_smoke(v2_bundle)
    test_exp188_cli_smoke(v2_bundle)
    test_exp189_cli_smoke(v2_bundle)
    test_exp190_cli_smoke(v2_bundle)
    test_exp191_no_retuning_audit(v2_bundle)
    _run_script(
        "exp192_v2_required_metric_aggregation.py",
        v2_bundle,
        "--output-dir",
        str(v2_bundle),
    )

    assert (v2_bundle / "data" / "v2_cross_family_robustness_metrics.csv").exists()


def test_exp193_bundle_report(v2_bundle: Path) -> None:
    test_exp192_aggregation(v2_bundle)
    _run_script(
        "exp193_v2_diagnostic_complete_bundle_report.py",
        v2_bundle,
        "--output-dir",
        str(v2_bundle),
    )

    assert (
        v2_bundle / "data" / "v2_diagnostic_complete_bundle_report.csv"
    ).exists()


def test_exp194_final_sanity(v2_bundle: Path) -> None:
    test_exp193_bundle_report(v2_bundle)
    _run_script(
        "exp194_v2_manifest_generation_final_sanity.py",
        v2_bundle,
        "--output-dir",
        str(v2_bundle),
    )

    assert (
        v2_bundle / "data" / "v2_manifest_generation_final_sanity.csv"
    ).exists()
