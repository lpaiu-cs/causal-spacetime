from __future__ import annotations

# ruff: noqa: E402,I001

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.abspath("experiments"))


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


def _metric_row(name: str, heldout: str = "0.5") -> dict[str, str]:
    return {
        "family_name": name,
        "family_kind": "structured",
        "manifest_count": "1",
        "fitted_fraction": "1",
        "no_fit_fraction": "0",
        "mean_heldout_violation": heldout,
        "mean_generalization_gap": "0.01",
        "stricter_threshold_pass_fraction": "1",
        "destructive_null_gap": "0.2",
        "symmetry_control_gap": "0.01",
        "target_coverage_fraction": "1",
        "pair_node_coverage_fraction": "1",
        "restart_std": "0.01",
        "latent_order_disagreement": "0.01",
        "no_retuning_audit_pass": "1",
        "failed_accounting_present": "1",
    }


@pytest.fixture()
def m39_bundle(tmp_path: Path) -> Path:
    data = tmp_path / "data"
    rows = [_metric_row("rank_gap_more_protocol_columns_v2")]
    _write_csv(data / "v2_cross_family_robustness_decision_metrics.csv", rows)
    _write_csv(
        data / "v2_cross_family_robustness_decisions.csv",
        [
            {
                **rows[0],
                "decision": "blocked",
                "passed": "0",
                "failed_reasons": "low_manifest_count;high_heldout_violation",
                "warning_reasons": "",
                "missing_inputs": "",
                "diagnostic_complete": "1",
            }
        ],
    )
    return tmp_path


def _run(script: str, output_dir: Path, *args: str) -> None:
    result = subprocess.run(
        [_python_executable(), f"experiments/{script}", *args],
        cwd=Path.cwd(),
        env=_env_with_src(),
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_exp206_exact_sanity() -> None:
    from experiments.exp206_v2_blocking_analysis_exact_sanity import run_experiment

    rows = run_experiment()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp207_cli_smoke(m39_bundle: Path) -> None:
    _run(
        "exp207_v2_blocked_root_cause_audit.py",
        m39_bundle,
        "--output-dir",
        str(m39_bundle),
    )

    assert (m39_bundle / "data" / "v2_blocked_root_cause_audit.csv").exists()


def test_exp208_cli_smoke(m39_bundle: Path) -> None:
    _run(
        "exp208_v2_criterion_margin_report.py",
        m39_bundle,
        "--output-dir",
        str(m39_bundle),
    )

    assert (m39_bundle / "data" / "v2_criterion_margin_report.csv").exists()


def test_exp209_cli_smoke(m39_bundle: Path) -> None:
    _run(
        "exp209_v2_structural_vs_measured_blocking.py",
        m39_bundle,
        "--output-dir",
        str(m39_bundle),
        "--hypothetical-manifest-count",
        "3",
    )

    assert (m39_bundle / "data" / "v2_structural_count_counterfactual.csv").exists()


def test_exp210_design_table() -> None:
    _run("exp210_v3_manifest_family_design.py", Path("outputs"))

    assert Path("outputs/data/v3_manifest_family_design.csv").exists()


def test_exp211_preregistration_export(m39_bundle: Path) -> None:
    _run("exp210_v3_manifest_family_design.py", m39_bundle)
    _run(
        "exp211_v3_preregistration_export.py",
        m39_bundle,
        "--output-dir",
        str(m39_bundle),
    )

    assert (m39_bundle / "remediation" / "v3_preregistration_spec_m39.json").exists()


def test_exp212_no_execution_audit(m39_bundle: Path) -> None:
    test_exp211_preregistration_export(m39_bundle)
    _run("exp212_v3_no_execution_audit.py", m39_bundle, "--output-dir", str(m39_bundle))

    assert (m39_bundle / "data" / "v3_no_execution_audit.csv").exists()


def test_exp213_report_card(m39_bundle: Path) -> None:
    test_exp207_cli_smoke(m39_bundle)
    test_exp208_cli_smoke(m39_bundle)
    test_exp209_cli_smoke(m39_bundle)
    test_exp211_preregistration_export(m39_bundle)
    test_exp212_no_execution_audit(m39_bundle)
    _run(
        "exp213_v2_blocked_decision_report_card.py",
        m39_bundle,
        "--output-dir",
        str(m39_bundle),
    )

    assert (m39_bundle / "data" / "v2_blocked_decision_report_card.csv").exists()


def test_exp214_final_sanity(m39_bundle: Path) -> None:
    test_exp213_report_card(m39_bundle)
    _run(
        "exp214_v2_blocked_v3_preregistration_final_sanity.py",
        m39_bundle,
    )

    assert Path("outputs/data/v2_blocked_v3_preregistration_final_sanity.csv").exists()
