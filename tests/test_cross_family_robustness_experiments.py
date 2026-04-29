from __future__ import annotations

# ruff: noqa: E402,I001

import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("experiments"))

from experiments.exp156_cross_family_robustness_exact_sanity import (
    run_experiment as run_exp156,
)
from experiments.exp157_cross_family_robustness_criteria_table import (
    run_experiment as run_exp157,
)
from experiments.exp160_carry_forward_registry_export import (
    ExperimentConfig as Exp160Config,
)
from experiments.exp160_carry_forward_registry_export import (
    run_experiment as run_exp160,
)
from experiments.exp161_cross_family_failed_provisional_accounting import (
    ExperimentConfig as Exp161Config,
)
from experiments.exp161_cross_family_failed_provisional_accounting import (
    run_experiment as run_exp161,
)
from experiments.exp162_stress_test_handoff_plan import (
    ExperimentConfig as Exp162Config,
)
from experiments.exp162_stress_test_handoff_plan import (
    run_experiment as run_exp162,
)
from experiments.exp163_cross_family_robustness_report_card import (
    ExperimentConfig as Exp163Config,
)
from experiments.exp163_cross_family_robustness_report_card import (
    run_experiment as run_exp163,
)
from experiments.exp164_cross_family_robustness_final_sanity import (
    run_experiment as run_exp164,
)
from manifest_representation_experiment_helpers import build_exact_manifest


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


def _write_family_outputs(output_dir: Path) -> None:
    data = output_dir / "data"
    build_exact_manifest(output_dir, "cross_family_test_manifest.json")
    _write_csv(
        data / "manifest_family_fit_summary.csv",
        [
            {
                "family_name": "eligible_rank_gap",
                "family_kind": "structured",
                "embedding_dim": "1",
                "manifest_count": "4",
                "fitted_count": "4",
                "no_fit_count": "0",
                "mean_train_violation": "0.05",
                "mean_heldout_violation": "0.10",
                "mean_generalization_gap": "0.05",
                "median_heldout_violation": "0.10",
                "best_heldout_violation": "0.05",
                "worst_heldout_violation": "0.15",
            }
        ],
    )
    _write_csv(
        data / "manifest_family_null_taxonomy.csv",
        [
            {
                "manifest_id": "m",
                "null_type": "shuffled_sides",
                "mean_heldout_violation_rate": "0.35",
                "taxonomy_class": "destructive_null",
            },
            {
                "manifest_id": "m",
                "null_type": "permuted_targets",
                "mean_heldout_violation_rate": "0.11",
                "taxonomy_class": "symmetry_control",
            },
        ],
    )
    _write_csv(
        data / "manifest_family_stricter_criteria.csv",
        [
            {
                "family_name": "eligible_rank_gap",
                "family_kind": "structured",
                "threshold_pass": "1",
            }
        ],
    )
    _write_csv(
        data / "manifest_family_failed_manifest_accounting.csv",
        [
            {
                "row_type": "family_count",
                "family_name": "eligible_rank_gap",
                "count": "4",
            }
        ],
    )
    _write_csv(
        data / "manifest_family_no_retuning_audit.csv",
        [{"check": "no_retuning_statement", "passed": "1", "value": "ok"}],
    )
    _write_csv(
        data / "manifest_family_report_card.csv",
        [
            {
                "family_name": "eligible_rank_gap",
                "family_kind": "structured",
                "fitted_count": "4",
                "no_fit_count": "0",
            }
        ],
    )


def test_exp156_exact_sanity() -> None:
    rows = run_exp156()

    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp157_criteria_table() -> None:
    rows = run_exp157()

    assert any(row["criterion"] == "max_mean_heldout_violation" for row in rows)


def test_exp158_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _write_family_outputs(tmp_path)
    result = subprocess.run(
        [
            _python_executable(),
            "experiments/exp158_cross_family_robustness_decision.py",
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
    assert (tmp_path / "data" / "cross_family_robustness_decisions.csv").exists()


def test_exp159_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _write_family_outputs(tmp_path)
    result = subprocess.run(
        [
            _python_executable(),
            "experiments/exp159_cross_family_robustness_threshold_sensitivity.py",
            "--output-dir",
            str(tmp_path),
            "--heldout-thresholds",
            "0.20",
            "--null-gap-thresholds",
            "0.05",
            "--stricter-pass-thresholds",
            "0.25",
        ],
        cwd=Path.cwd(),
        env=_env_with_src(),
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr


def test_exp160_registry_export(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _write_family_outputs(tmp_path)
    run_exp160(Exp160Config(output_dir=tmp_path))

    assert (tmp_path / "carry_forward" / "carry_forward_registry.json").exists()


def test_exp161_accounting(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _write_family_outputs(tmp_path)
    subprocess.run(
        [
            _python_executable(),
            "experiments/exp158_cross_family_robustness_decision.py",
            "--output-dir",
            str(tmp_path),
        ],
        cwd=Path.cwd(),
        env=_env_with_src(),
        check=True,
    )

    assert run_exp161(Exp161Config(output_dir=tmp_path))


def test_exp162_handoff_plan(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _write_family_outputs(tmp_path)
    run_exp160(Exp160Config(output_dir=tmp_path))

    assert run_exp162(Exp162Config(output_dir=tmp_path))


def test_exp163_report_card(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _write_family_outputs(tmp_path)
    subprocess.run(
        [
            _python_executable(),
            "experiments/exp158_cross_family_robustness_decision.py",
            "--output-dir",
            str(tmp_path),
        ],
        cwd=Path.cwd(),
        env=_env_with_src(),
        check=True,
    )
    run_exp160(Exp160Config(output_dir=tmp_path))
    run_exp162(Exp162Config(output_dir=tmp_path))

    assert run_exp163(Exp163Config(output_dir=tmp_path))


def test_exp164_final_sanity(tmp_path) -> None:  # type: ignore[no-untyped-def]
    data = tmp_path / "data"
    data.mkdir(parents=True, exist_ok=True)
    (data / "cross_family_robustness_report_card.csv").write_text(
        "family_name,decision\nblocked_family,blocked\n",
        encoding="utf-8",
    )
    rows = run_exp164(tmp_path)

    assert all(float(row["passed"]) == 1.0 for row in rows)
