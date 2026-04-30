from __future__ import annotations

# ruff: noqa: E402,I001

import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("experiments"))

import pytest

from experiments.exp195_v2_carry_forward_exact_sanity import (
    run_experiment as run_exp195,
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


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _metric_row(family_name: str, family_kind: str) -> dict[str, str]:
    return {
        "family_name": family_name,
        "family_kind": family_kind,
        "manifest_count": "3" if family_kind == "structured" else "1",
        "fitted_fraction": "1" if family_kind == "structured" else "0",
        "no_fit_fraction": "0" if family_kind == "structured" else "1",
        "mean_heldout_violation": "0.05" if family_kind == "structured" else "1",
        "mean_generalization_gap": "0.01",
        "stricter_threshold_pass_fraction": "1",
        "destructive_null_gap": "0.2",
        "symmetry_control_gap": "0.01",
        "target_coverage_fraction": "1",
        "pair_node_coverage_fraction": "0.8",
        "restart_std": "0.01",
        "latent_order_disagreement": "0.1",
        "no_retuning_audit_pass": "1",
        "failed_accounting_present": "1",
    }


def _write_v2_bundle(output_dir: Path) -> None:
    data = output_dir / "data"
    metrics = [
        _metric_row("strong_v2", "structured"),
        _metric_row("failed_v2", "failed_control"),
    ]
    _write_csv(data / "v2_cross_family_robustness_metrics.csv", metrics)
    _write_csv(
        data / "v2_diagnostic_completeness_check.csv",
        [
            {
                "family_name": row["family_name"],
                "family_kind": row["family_kind"],
                "required_metric_count": "14",
                "available_metric_count": "14",
                "missing_metric_count": "0",
                "completeness_fraction": "1",
                "missing_metrics": "",
                "not_carry_forward_evaluated": "1",
            }
            for row in metrics
        ],
    )
    for filename in [
        "v2_diagnostic_complete_bundle_report.csv",
        "v2_manifest_generation.csv",
        "v2_manifest_family_fit_summary.csv",
        "v2_manifest_family_null_taxonomy.csv",
        "v2_manifest_family_stricter_criteria.csv",
        "v2_manifest_family_failed_accounting.csv",
        "v2_manifest_family_coverage_metrics.csv",
        "v2_manifest_family_restart_stability.csv",
        "v2_manifest_family_latent_order_stability.csv",
        "v2_no_retuning_audit.csv",
    ]:
        _write_csv(data / filename, [{"family_name": "strong_v2", "ok": "1"}])
    _write_csv(
        data / "cross_family_robustness_criteria_table.csv",
        [
            {"criterion": "name", "value": "default_cross_family_robustness"},
            {"criterion": "min_manifest_count", "value": "3"},
            {"criterion": "min_fitted_fraction", "value": "0.5"},
            {"criterion": "max_no_fit_fraction", "value": "0.5"},
            {"criterion": "max_mean_heldout_violation", "value": "0.2"},
            {"criterion": "max_mean_generalization_gap", "value": "0.1"},
            {"criterion": "min_stricter_threshold_pass_fraction", "value": "0.5"},
            {"criterion": "min_destructive_null_gap", "value": "0.1"},
            {"criterion": "max_symmetry_control_gap", "value": "0.1"},
            {"criterion": "min_target_coverage_fraction", "value": "0.8"},
            {"criterion": "min_pair_node_coverage_fraction", "value": "0.5"},
            {"criterion": "max_restart_std", "value": "0.1"},
            {"criterion": "max_latent_order_disagreement", "value": "0.3"},
            {"criterion": "require_no_retuning_audit", "value": "True"},
            {"criterion": "require_failed_accounting", "value": "True"},
        ],
    )
    (output_dir / "manifests_v2").mkdir(parents=True, exist_ok=True)


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


@pytest.fixture()
def v2_eval_bundle(tmp_path: Path) -> Path:
    _write_v2_bundle(tmp_path)
    return tmp_path


def test_exp195_exact_sanity() -> None:
    rows = run_exp195()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp196_input_audit(v2_eval_bundle: Path) -> None:
    _run_script(
        "exp196_v2_bundle_input_audit.py",
        v2_eval_bundle,
        "--output-dir",
        str(v2_eval_bundle),
    )

    assert (v2_eval_bundle / "data" / "v2_bundle_input_audit.csv").exists()


def test_exp197_cli_smoke(v2_eval_bundle: Path) -> None:
    _run_script(
        "exp197_v2_carry_forward_decision.py",
        v2_eval_bundle,
        "--output-dir",
        str(v2_eval_bundle),
    )

    assert (
        v2_eval_bundle / "data" / "v2_cross_family_robustness_decisions.csv"
    ).exists()


def test_exp198_sensitivity(v2_eval_bundle: Path) -> None:
    _run_script(
        "exp198_v2_carry_forward_threshold_sensitivity.py",
        v2_eval_bundle,
        "--output-dir",
        str(v2_eval_bundle),
        "--heldout-thresholds",
        "0.2",
        "--null-gap-thresholds",
        "0.1",
        "--stricter-pass-thresholds",
        "0.5",
    )

    assert (
        v2_eval_bundle / "data" / "v2_carry_forward_threshold_sensitivity.csv"
    ).exists()


def test_exp199_registry_export(v2_eval_bundle: Path) -> None:
    test_exp197_cli_smoke(v2_eval_bundle)
    _run_script(
        "exp199_v2_carry_forward_registry_export.py",
        v2_eval_bundle,
        "--output-dir",
        str(v2_eval_bundle),
        "--manifest-dir",
        str(v2_eval_bundle / "manifests_v2"),
    )

    assert (
        v2_eval_bundle
        / "carry_forward_v2"
        / "carry_forward_registry_v2.json"
    ).exists()


def test_exp200_handoff_plan(v2_eval_bundle: Path) -> None:
    test_exp199_registry_export(v2_eval_bundle)
    _run_script(
        "exp200_v2_stress_test_handoff_plan.py",
        v2_eval_bundle,
        "--output-dir",
        str(v2_eval_bundle),
    )

    assert (v2_eval_bundle / "data" / "v2_stress_test_handoff_plan.csv").exists()


def test_exp201_accounting(v2_eval_bundle: Path) -> None:
    test_exp197_cli_smoke(v2_eval_bundle)
    _run_script(
        "exp201_v2_failed_provisional_accounting.py",
        v2_eval_bundle,
        "--output-dir",
        str(v2_eval_bundle),
    )

    assert (v2_eval_bundle / "data" / "v2_failed_provisional_accounting.csv").exists()


def test_exp202_failure_decomposition(v2_eval_bundle: Path) -> None:
    _run_script(
        "exp202_v2_carry_forward_failure_decomposition.py",
        v2_eval_bundle,
        "--output-dir",
        str(v2_eval_bundle),
    )

    assert (
        v2_eval_bundle / "data" / "v2_carry_forward_failure_decomposition.csv"
    ).exists()


def test_exp203_no_retuning_audit(v2_eval_bundle: Path) -> None:
    test_exp197_cli_smoke(v2_eval_bundle)
    test_exp198_sensitivity(v2_eval_bundle)
    (v2_eval_bundle / "manifests_v2" / "placeholder.json").write_text(
        "{}",
        encoding="utf-8",
    )
    _run_script(
        "exp203_v2_carry_forward_no_retuning_audit.py",
        v2_eval_bundle,
        "--output-dir",
        str(v2_eval_bundle),
    )

    assert (
        v2_eval_bundle / "data" / "v2_carry_forward_no_retuning_audit.csv"
    ).exists()


def test_exp204_report_card(v2_eval_bundle: Path) -> None:
    test_exp200_handoff_plan(v2_eval_bundle)
    test_exp202_failure_decomposition(v2_eval_bundle)
    _run_script(
        "exp204_v2_carry_forward_report_card.py",
        v2_eval_bundle,
        "--output-dir",
        str(v2_eval_bundle),
    )

    assert (v2_eval_bundle / "data" / "v2_carry_forward_report_card.csv").exists()


def test_exp205_final_sanity(v2_eval_bundle: Path) -> None:
    test_exp204_report_card(v2_eval_bundle)
    _run_script(
        "exp205_v2_carry_forward_final_sanity.py",
        v2_eval_bundle,
        "--output-dir",
        str(v2_eval_bundle),
    )

    assert (v2_eval_bundle / "data" / "v2_carry_forward_final_sanity.csv").exists()
