from __future__ import annotations

# ruff: noqa: E402,I001

import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("experiments"))

from experiments.exp165_carry_forward_failure_decomposition_exact_sanity import (
    run_experiment as run_exp165,
)
from experiments.exp166_carry_forward_failure_decomposition import (
    ExperimentConfig as Exp166Config,
)
from experiments.exp166_carry_forward_failure_decomposition import (
    run_experiment as run_exp166,
    write_outputs as write_exp166_outputs,
)
from experiments.exp167_cross_family_diagnostic_completeness_audit import (
    ExperimentConfig as Exp167Config,
)
from experiments.exp167_cross_family_diagnostic_completeness_audit import (
    run_experiment as run_exp167,
    write_outputs as write_exp167_outputs,
)
from experiments.exp168_stress_test_stop_condition_audit import (
    ExperimentConfig as Exp168Config,
)
from experiments.exp168_stress_test_stop_condition_audit import (
    run_experiment as run_exp168,
)
from experiments.exp169_upstream_remediation_design_table import (
    ExperimentConfig as Exp169Config,
)
from experiments.exp169_upstream_remediation_design_table import (
    run_experiment as run_exp169,
)
from experiments.exp170_missing_metric_impact_report import (
    ExperimentConfig as Exp170Config,
)
from experiments.exp170_missing_metric_impact_report import (
    run_experiment as run_exp170,
    write_outputs as write_exp170_outputs,
)
from experiments.exp171_failure_decomposition_no_retuning_audit import (
    run_experiment as run_exp171,
)
from experiments.exp172_carry_forward_failure_report_card import (
    ExperimentConfig as Exp172Config,
)
from experiments.exp172_carry_forward_failure_report_card import (
    run_experiment as run_exp172,
)
from experiments.exp173_carry_forward_failure_final_sanity import (
    run_experiment as run_exp173,
)
from experiments.exp157_cross_family_robustness_criteria_table import (
    run_experiment as run_exp157,
    write_outputs as write_exp157_outputs,
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


def _write_m34_outputs(output_dir: Path) -> None:
    data = output_dir / "data"
    _write_csv(
        data / "cross_family_robustness_metrics.csv",
        [
            {
                "family_name": "eligible_rank_gap",
                "family_kind": "structured",
                "manifest_count": "4",
                "fitted_fraction": "1.0",
                "no_fit_fraction": "0.0",
                "mean_heldout_violation": "0.3",
                "mean_generalization_gap": "0.12",
                "stricter_threshold_pass_fraction": "0.0",
                "destructive_null_gap": "0.2",
                "symmetry_control_gap": "0.03",
                "target_coverage_fraction": "nan",
                "pair_node_coverage_fraction": "nan",
                "restart_std": "nan",
                "latent_order_disagreement": "nan",
                "no_retuning_audit_pass": "1.0",
                "failed_accounting_present": "1.0",
            }
        ],
    )
    _write_csv(
        data / "cross_family_robustness_decisions.csv",
        [
            {
                "family_name": "eligible_rank_gap",
                "family_kind": "structured",
                "decision": "blocked",
                "passed": "0",
                "failed_reasons": "high_heldout_violation",
                "warning_reasons": "missing_low_target_coverage",
            }
        ],
    )
    write_exp157_outputs(run_exp157(), output_dir)


def test_exp165_exact_sanity() -> None:
    rows = run_exp165()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp166_cli_smoke(tmp_path: Path) -> None:
    _write_m34_outputs(tmp_path)
    result = subprocess.run(
        [
            _python_executable(),
            "experiments/exp166_carry_forward_failure_decomposition.py",
            "--output-dir",
            str(tmp_path),
        ],
        cwd=Path.cwd(),
        env=_env_with_src(),
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "data" / "carry_forward_failure_decomposition.csv").exists()


def test_exp167_cli_smoke(tmp_path: Path) -> None:
    _write_m34_outputs(tmp_path)
    result = subprocess.run(
        [
            _python_executable(),
            "experiments/exp167_cross_family_diagnostic_completeness_audit.py",
            "--output-dir",
            str(tmp_path),
        ],
        cwd=Path.cwd(),
        env=_env_with_src(),
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr


def test_exp168_stop_audit(tmp_path: Path) -> None:
    _write_m34_outputs(tmp_path)
    rows = run_exp168(Exp168Config(output_dir=tmp_path))

    assert float(rows[0]["stress_tests_allowed"]) == 0.0


def test_exp169_remediation_table(tmp_path: Path) -> None:
    _write_m34_outputs(tmp_path)
    records, summary = run_exp166(Exp166Config(output_dir=tmp_path))
    write_exp166_outputs(records, summary, tmp_path)

    rows = run_exp169(Exp169Config(output_dir=tmp_path))

    assert any(row["proposal_name"] == "add_protocol_columns" for row in rows)


def test_exp170_missing_metric_report(tmp_path: Path) -> None:
    _write_m34_outputs(tmp_path)
    records, summary = run_exp166(Exp166Config(output_dir=tmp_path))
    write_exp166_outputs(records, summary, tmp_path)
    completeness = run_exp167(Exp167Config(output_dir=tmp_path))
    write_exp167_outputs(completeness, tmp_path)

    rows = run_exp170(Exp170Config(output_dir=tmp_path))
    write_exp170_outputs(rows, tmp_path)

    assert rows[0]["decision_changed"] == 0.0


def test_exp171_no_retuning_audit(tmp_path: Path) -> None:
    _write_m34_outputs(tmp_path)

    rows = run_exp171(tmp_path)

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp172_report_card(tmp_path: Path) -> None:
    _write_m34_outputs(tmp_path)
    records, summary = run_exp166(Exp166Config(output_dir=tmp_path))
    write_exp166_outputs(records, summary, tmp_path)
    completeness = run_exp167(Exp167Config(output_dir=tmp_path))
    write_exp167_outputs(completeness, tmp_path)
    missing = run_exp170(Exp170Config(output_dir=tmp_path))
    write_exp170_outputs(missing, tmp_path)

    rows = run_exp172(Exp172Config(output_dir=tmp_path))

    assert rows[0]["decision"] == "blocked"


def test_exp173_final_sanity(tmp_path: Path) -> None:
    _write_csv(
        tmp_path / "data" / "carry_forward_failure_report_card.csv",
        [
            {
                "family_name": "eligible_rank_gap",
                "decision": "blocked",
                "hard_measured_failures": "mean_heldout_violation",
            }
        ],
    )

    rows = run_exp173(tmp_path)

    assert all(float(row["passed"]) == 1.0 for row in rows)
