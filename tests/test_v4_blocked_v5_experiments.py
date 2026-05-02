from __future__ import annotations

# ruff: noqa: E402,I001

import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("experiments"))

from experiments.exp291_v4_blocking_exact_sanity import run_experiment as run_exp291

from causal_spacetime_lab.state_change_manifest_v5_design import (
    default_v5_protocol_family_designs,
    v5_protocol_family_design_table,
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


def _metric_rows(prefix: str = "v4") -> list[dict[str, str]]:
    structured_name = (
        "rank_gap_earliest_full_stability_v4"
        if prefix == "v4"
        else "rank_gap_earliest_full_reference_v3"
    )
    return [
        {
            "family_name": structured_name,
            "family_kind": "structured",
            "manifest_count": "5",
            "fitted_fraction": "1",
            "no_fit_fraction": "0",
            "mean_heldout_violation": "0.42",
            "mean_generalization_gap": "0.03",
            "stricter_threshold_pass_fraction": "0",
            "destructive_null_gap": "0.05",
            "symmetry_control_gap": "0.01",
            "target_coverage_fraction": "1",
            "pair_node_coverage_fraction": "0.4",
            "restart_std": "0.01",
            "latent_order_disagreement": "0.6",
            "no_retuning_audit_pass": "1",
            "failed_accounting_present": "1",
            "dominant_handoff_provenance_type": (
                "hybrid_template_instantiated_from_profile"
            ),
        },
        {
            "family_name": f"failed_controls_{prefix}",
            "family_kind": "failed_control",
            "manifest_count": "3",
            "fitted_fraction": "0",
            "no_fit_fraction": "1",
            "mean_heldout_violation": "1",
            "mean_generalization_gap": "0",
            "stricter_threshold_pass_fraction": "0",
            "destructive_null_gap": "-1",
            "symmetry_control_gap": "1",
            "target_coverage_fraction": "0.5",
            "pair_node_coverage_fraction": "0.2",
            "restart_std": "0",
            "latent_order_disagreement": "0",
            "no_retuning_audit_pass": "1",
            "failed_accounting_present": "1",
            "dominant_handoff_provenance_type": "report_only_control",
        },
    ]


def _prepare_inputs(output_dir: Path) -> None:
    data = output_dir / "data"
    v4_metrics = _metric_rows("v4")
    _write(
        data / "v4_protocol_cross_family_robustness_decision_metrics.csv",
        v4_metrics,
    )
    _write(data / "v4_protocol_cross_family_robustness_metrics.csv", v4_metrics)
    _write(data / "v3_protocol_cross_family_robustness_metrics.csv", _metric_rows("v3"))
    _write(
        data / "v4_protocol_cross_family_robustness_decisions.csv",
        [
            {
                "family_name": "rank_gap_earliest_full_stability_v4",
                "family_kind": "structured",
                "decision": "blocked",
                "diagnostic_complete": "1",
                "preconditions_passed": "1",
                "failed_reasons": (
                    "high_heldout_violation;high_latent_order_disagreement"
                ),
            },
            {
                "family_name": "failed_controls_v4",
                "family_kind": "failed_control",
                "decision": "failed_control",
                "diagnostic_complete": "1",
                "preconditions_passed": "1",
                "failed_reasons": "failed_control_family",
            },
        ],
    )
    _write(
        data / "v4_protocol_precondition_audit.csv",
        [
            {
                "family_name": "rank_gap_earliest_full_stability_v4",
                "family_kind": "structured",
                "preconditions_passed": "1",
                "failed_preconditions": "",
            }
        ],
    )
    _write(
        data / "v4_protocol_manifest_family_fit_comparison.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "rank_gap_earliest_full_stability_v4",
                "family_kind": "structured",
                "train_violation_rate": "0.1",
                "heldout_violation_rate": "0.42",
            }
        ],
    )
    _write(
        data / "v4_protocol_manifest_family_coverage_metrics.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "rank_gap_earliest_full_stability_v4",
                "target_coverage_fraction": "1",
                "pair_node_coverage_fraction": "0.4",
            }
        ],
    )
    _write(
        data / "v4_protocol_manifest_family_restart_stability.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "rank_gap_earliest_full_stability_v4",
                "restart_std": "0.01",
            }
        ],
    )
    _write(
        data / "v4_protocol_manifest_family_latent_order_stability.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "rank_gap_earliest_full_stability_v4",
                "latent_order_disagreement": "0.6",
            }
        ],
    )
    _write(
        data / "v4_protocol_manifest_metadata_audit.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "rank_gap_earliest_full_stability_v4",
                "measurement_protocol_id": "p",
                "profile_invariance_status": "protocol_invariant",
            }
        ],
    )
    _write(
        data / "v4_protocol_manifest_family_failed_accounting.csv",
        [
            {
                "manifest_id": "m1",
                "family_name": "rank_gap_earliest_full_stability_v4",
                "eligible": "1",
                "failed_reasons": "",
                "handoff_provenance_type": "hybrid_template_instantiated_from_profile",
            }
        ],
    )
    _write(
        data / "v4_protocol_manifest_family_null_taxonomy.csv",
        [
            {
                "family_name": "rank_gap_earliest_full_stability_v4",
                "taxonomy_class": "destructive_null",
                "mean_heldout_violation_rate": "0.47",
                "structured_heldout_violation_rate": "0.42",
            },
            {
                "family_name": "rank_gap_earliest_full_stability_v4",
                "taxonomy_class": "symmetry_control",
                "mean_heldout_violation_rate": "0.43",
                "structured_heldout_violation_rate": "0.42",
            },
        ],
    )
    _write(
        data / "v5_protocol_family_design.csv",
        [
            {key: str(value) for key, value in row.items()}
            for row in v5_protocol_family_design_table(
                default_v5_protocol_family_designs()
            )
        ],
    )


def _run_prereq_m46(output_dir: Path) -> None:
    _prepare_inputs(output_dir)
    for script in [
        "exp292_v4_blocked_root_cause_audit.py",
        "exp293_v4_criterion_margin_report.py",
        "exp298_v3_to_v4_delta_audit.py",
        "exp299_v4_report_only_counterfactuals.py",
        "exp301_v5_remediation_iteration_risk_audit.py",
        "exp302_v5_protocol_preregistration_export.py",
    ]:
        _run_script(script, output_dir, "--output-dir", str(output_dir))


def test_exp291_exact_sanity() -> None:
    rows = run_exp291()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp292_to_exp299_cli_smoke(tmp_path: Path) -> None:
    _prepare_inputs(tmp_path)
    for script in [
        "exp292_v4_blocked_root_cause_audit.py",
        "exp293_v4_criterion_margin_report.py",
        "exp294_v4_manifest_level_drilldown.py",
        "exp295_v4_null_failure_drilldown.py",
        "exp296_v4_stability_failure_drilldown.py",
        "exp297_v4_coverage_failure_drilldown.py",
        "exp298_v3_to_v4_delta_audit.py",
        "exp299_v4_report_only_counterfactuals.py",
    ]:
        _run_script(script, tmp_path, "--output-dir", str(tmp_path))

    assert (tmp_path / "data" / "v4_single_fix_counterfactual.csv").exists()


def test_exp300_v5_design_table() -> None:
    _run_script("exp300_v5_protocol_family_design.py", Path("outputs"))

    assert Path("outputs/data/v5_protocol_family_design.csv").exists()


def test_exp301_to_exp306_cli_smoke(tmp_path: Path) -> None:
    _run_prereq_m46(tmp_path)
    for script in [
        "exp303_v5_no_execution_audit.py",
        "exp304_v4_blocked_decision_report_card.py",
        "exp305_v4_blocked_no_retuning_audit.py",
        "exp306_v4_blocked_v5_final_sanity.py",
    ]:
        _run_script(script, tmp_path, "--output-dir", str(tmp_path))

    assert (tmp_path / "data" / "v4_blocked_v5_final_sanity.csv").exists()
