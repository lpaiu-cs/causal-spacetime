from __future__ import annotations

# ruff: noqa: E402,I001

import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("experiments"))

from experiments.exp252_v3_protocol_blocking_exact_sanity import (
    run_experiment as run_exp252,
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


def _write(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _metric_rows() -> list[dict[str, str]]:
    base = {
        "manifest_count": "3",
        "fitted_fraction": "1",
        "no_fit_fraction": "0",
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
        "dominant_handoff_provenance_type": "hybrid_template_instantiated_from_profile",
    }
    structured = {
        **base,
        "family_name": "structured_family",
        "family_kind": "structured",
        "mean_heldout_violation": "0.35",
        "latent_order_disagreement": "0.6",
    }
    failed = {
        **base,
        "family_name": "failed_controls_v3",
        "family_kind": "failed_control",
        "mean_heldout_violation": "nan",
    }
    return [structured, failed]


def _decision_rows() -> list[dict[str, str]]:
    return [
        {
            "family_name": "structured_family",
            "family_kind": "structured",
            "decision": "blocked",
            "diagnostic_complete": "1",
            "preconditions_passed": "1",
            "failed_reasons": "high_heldout_violation;high_latent_order_disagreement",
        },
        {
            "family_name": "failed_controls_v3",
            "family_kind": "failed_control",
            "decision": "failed_control",
            "diagnostic_complete": "1",
            "preconditions_passed": "1",
            "failed_reasons": "failed_control_family",
        },
    ]


def _prepare_inputs(output_dir: Path) -> None:
    data = output_dir / "data"
    metrics = _metric_rows()
    _write(data / "v3_protocol_cross_family_robustness_decision_metrics.csv", metrics)
    _write(data / "v3_protocol_cross_family_robustness_metrics.csv", metrics)
    _write(data / "v3_protocol_cross_family_robustness_decisions.csv", _decision_rows())
    _write(
        data / "v3_protocol_precondition_audit.csv",
        [
            {
                "family_name": "structured_family",
                "family_kind": "structured",
                "preconditions_passed": "1",
                "failed_preconditions": "",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_family_fit_comparison.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "structured_family",
                "family_kind": "structured",
                "train_violation_rate": "0.1",
                "heldout_violation_rate": "0.35",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_family_coverage_metrics.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "structured_family",
                "target_coverage_fraction": "1",
                "pair_node_coverage_fraction": "1",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_family_restart_stability.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "structured_family",
                "restart_std": "0.01",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_family_latent_order_stability.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "structured_family",
                "latent_order_disagreement": "0.6",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_metadata_audit.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "structured_family",
                "measurement_protocol_id": "p",
                "profile_invariance_status": "protocol_invariant",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_family_failed_accounting.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "structured_family",
                "eligible": "1",
                "failed_reasons": "",
                "handoff_provenance_type": "hybrid_template_instantiated_from_profile",
            }
        ],
    )
    _write(
        data / "v3_protocol_manifest_family_null_taxonomy.csv",
        [
            {
                "family_name": "structured_family",
                "taxonomy_class": "destructive_null",
                "mean_heldout_violation_rate": "0.4",
                "structured_heldout_violation_rate": "0.35",
            },
            {
                "family_name": "structured_family",
                "taxonomy_class": "symmetry_control",
                "mean_heldout_violation_rate": "0.36",
                "structured_heldout_violation_rate": "0.35",
            },
        ],
    )


def _run_prereq_m43(output_dir: Path) -> None:
    _prepare_inputs(output_dir)
    for script in [
        "exp253_v3_protocol_blocked_root_cause_audit.py",
        "exp254_v3_protocol_criterion_margin_report.py",
        "exp258_v3_protocol_report_only_counterfactuals.py",
        "exp260_v4_protocol_preregistration_export.py",
    ]:
        _run_script(script, output_dir, "--output-dir", str(output_dir))


def test_exp252_exact_sanity() -> None:
    rows = run_exp252()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp253_cli_smoke(tmp_path: Path) -> None:
    _prepare_inputs(tmp_path)
    _run_script(
        "exp253_v3_protocol_blocked_root_cause_audit.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (tmp_path / "data" / "v3_protocol_blocked_root_cause_audit.csv").exists()


def test_exp254_cli_smoke(tmp_path: Path) -> None:
    _run_prereq_m43(tmp_path)

    assert (tmp_path / "data" / "v3_protocol_criterion_margin_report.csv").exists()


def test_exp255_cli_smoke(tmp_path: Path) -> None:
    _prepare_inputs(tmp_path)
    _run_script(
        "exp255_v3_protocol_manifest_level_drilldown.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (tmp_path / "data" / "v3_protocol_manifest_level_drilldown.csv").exists()


def test_exp256_cli_smoke(tmp_path: Path) -> None:
    _prepare_inputs(tmp_path)
    _run_script(
        "exp256_v3_protocol_null_failure_drilldown.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (tmp_path / "data" / "v3_protocol_null_failure_drilldown.csv").exists()


def test_exp257_cli_smoke(tmp_path: Path) -> None:
    _prepare_inputs(tmp_path)
    _run_script(
        "exp257_v3_protocol_stability_failure_drilldown.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (tmp_path / "data" / "v3_protocol_stability_failure_drilldown.csv").exists()


def test_exp258_counterfactual_report(tmp_path: Path) -> None:
    _prepare_inputs(tmp_path)
    _run_script(
        "exp258_v3_protocol_report_only_counterfactuals.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (tmp_path / "data" / "v3_protocol_single_fix_counterfactual.csv").exists()


def test_exp259_v4_design_table(tmp_path: Path) -> None:
    _run_script("exp259_v4_protocol_family_design.py", tmp_path)

    assert Path("outputs/data/v4_protocol_family_design.csv").exists()


def test_exp260_v4_preregistration_export(tmp_path: Path) -> None:
    _run_script(
        "exp260_v4_protocol_preregistration_export.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (
        tmp_path / "remediation" / "v4_protocol_preregistration_spec_m43.json"
    ).exists()


def test_exp261_v4_no_execution_audit(tmp_path: Path) -> None:
    _run_script(
        "exp260_v4_protocol_preregistration_export.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )
    _run_script(
        "exp261_v4_no_execution_audit.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (tmp_path / "data" / "v4_no_execution_audit.csv").exists()


def test_exp262_report_card(tmp_path: Path) -> None:
    _run_prereq_m43(tmp_path)
    _run_script(
        "exp262_v3_protocol_blocked_decision_report_card.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (
        tmp_path / "data" / "v3_protocol_blocked_decision_report_card.csv"
    ).exists()


def test_exp263_no_retuning_audit(tmp_path: Path) -> None:
    _run_prereq_m43(tmp_path)
    _run_script(
        "exp263_v3_protocol_no_retuning_audit.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (tmp_path / "data" / "v3_protocol_blocked_no_retuning_audit.csv").exists()


def test_exp264_final_sanity(tmp_path: Path) -> None:
    _run_prereq_m43(tmp_path)
    _run_script(
        "exp262_v3_protocol_blocked_decision_report_card.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )
    _run_script(
        "exp264_v3_protocol_blocked_v4_final_sanity.py",
        tmp_path,
        "--output-dir",
        str(tmp_path),
    )

    assert (
        tmp_path / "data" / "v3_protocol_blocked_v4_final_sanity.csv"
    ).exists()
