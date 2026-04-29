from __future__ import annotations

# ruff: noqa: E402,I001

import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("experiments"))

from experiments.exp174_remediation_plan_exact_sanity import (
    run_experiment as run_exp174,
)
from experiments.exp175_diagnostic_complete_schema_export import (
    run_experiment as run_exp175,
)
from experiments.exp177_new_manifest_family_design_v2 import (
    run_experiment as run_exp177,
)
from experiments.exp178_preregistered_remediation_plan_export import (
    ExperimentConfig as Exp178Config,
)
from experiments.exp178_preregistered_remediation_plan_export import (
    run_experiment as run_exp178,
    write_outputs as write_exp178_outputs,
)
from experiments.exp179_future_manifest_run_spec import (
    ExperimentConfig as Exp179Config,
)
from experiments.exp179_future_manifest_run_spec import (
    run_experiment as run_exp179,
    write_outputs as write_exp179_outputs,
)
from experiments.exp180_remediation_no_execution_audit import (
    ExperimentConfig as Exp180Config,
)
from experiments.exp180_remediation_no_execution_audit import (
    run_experiment as run_exp180,
)
from experiments.exp181_remediation_plan_report_card import (
    ExperimentConfig as Exp181Config,
)
from experiments.exp181_remediation_plan_report_card import (
    run_experiment as run_exp181,
)
from experiments.exp182_remediation_plan_final_sanity import (
    run_experiment as run_exp182,
)
from experiments.exp176_failure_to_remediation_mapping import (
    ExperimentConfig as Exp176Config,
)
from experiments.exp176_failure_to_remediation_mapping import (
    run_experiment as run_exp176,
    write_outputs as write_exp176_outputs,
)
from experiments.exp177_new_manifest_family_design_v2 import (
    write_outputs as write_exp177_outputs,
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


def _write_m35_outputs(output_dir: Path) -> None:
    data = output_dir / "data"
    _write_csv(
        data / "carry_forward_failure_summary.csv",
        [
            {
                "family_name": "eligible_rank_gap",
                "family_kind": "structured",
                "root_cause_category": "heldout_failure",
                "status": "measured_failure",
                "count": "1",
            },
            {
                "family_name": "eligible_rank_gap",
                "family_kind": "structured",
                "root_cause_category": "missing_metric",
                "status": "missing_metric",
                "count": "1",
            },
        ],
    )
    _write_csv(
        data / "cross_family_diagnostic_completeness_audit.csv",
        [
            {
                "family_name": "eligible_rank_gap",
                "required_metric_count": "14",
                "available_metric_count": "10",
                "missing_metric_count": "4",
                "completeness_fraction": "0.714",
                "missing_metrics": "target_coverage_fraction;restart_std",
            }
        ],
    )


def test_exp174_exact_sanity() -> None:
    rows = run_exp174()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp175_schema_export() -> None:
    rows, summary = run_exp175()

    assert len(rows) == 14
    assert summary


def test_exp176_cli_smoke(tmp_path: Path) -> None:
    _write_m35_outputs(tmp_path)
    result = subprocess.run(
        [
            _python_executable(),
            "experiments/exp176_failure_to_remediation_mapping.py",
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
    assert (tmp_path / "data" / "failure_to_remediation_mapping.csv").exists()


def test_exp177_design_table() -> None:
    rows = run_exp177()

    assert all(row["execution_status"] == "planned_only" for row in rows)


def test_exp178_plan_export(tmp_path: Path) -> None:
    _write_m35_outputs(tmp_path)
    plan_path, rows = run_exp178(Exp178Config(output_dir=tmp_path))

    assert plan_path.exists()
    assert rows[0]["execution_allowed_in_current_milestone"] == 0.0


def test_exp179_future_spec(tmp_path: Path) -> None:
    _write_m35_outputs(tmp_path)
    plan_path, rows = run_exp178(Exp178Config(output_dir=tmp_path))
    write_exp178_outputs(rows, tmp_path)
    assert plan_path.exists()

    spec_path, spec_rows = run_exp179(Exp179Config(output_dir=tmp_path))

    assert spec_path.exists()
    assert spec_rows[0]["allowed_to_execute_now"] == 0.0


def test_exp180_no_execution_audit(tmp_path: Path) -> None:
    _write_m35_outputs(tmp_path)
    _, plan_rows = run_exp178(Exp178Config(output_dir=tmp_path))
    write_exp178_outputs(plan_rows, tmp_path)
    _, spec_rows = run_exp179(Exp179Config(output_dir=tmp_path))
    write_exp179_outputs(spec_rows, tmp_path)

    rows = run_exp180(Exp180Config(output_dir=tmp_path))

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp181_report_card(tmp_path: Path) -> None:
    _write_m35_outputs(tmp_path)
    mapping_rows = run_exp176(Exp176Config(output_dir=tmp_path))
    write_exp176_outputs(mapping_rows, tmp_path)
    family_rows = run_exp177()
    write_exp177_outputs(family_rows, tmp_path)
    _, plan_rows = run_exp178(Exp178Config(output_dir=tmp_path))
    write_exp178_outputs(plan_rows, tmp_path)
    _, spec_rows = run_exp179(Exp179Config(output_dir=tmp_path))
    write_exp179_outputs(spec_rows, tmp_path)

    rows = run_exp181(Exp181Config(output_dir=tmp_path))

    assert any(row["row_type"] == "global_stop_condition" for row in rows)


def test_exp182_final_sanity() -> None:
    rows = run_exp182()

    assert all(float(row["passed"]) == 1.0 for row in rows)
