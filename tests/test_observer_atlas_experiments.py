from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp15_exact_poincare_map_sanity import run_experiment


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    src_path = os.path.abspath("src")
    env["PYTHONPATH"] = (
        src_path
        if not env.get("PYTHONPATH")
        else f"{src_path}{os.pathsep}{env['PYTHONPATH']}"
    )
    return env


def test_exp14_observer_atlas_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp14_observer_atlas_consistency.py",
            "--n-values",
            "25",
            "--tick-values",
            "8",
            "--repetitions",
            "1",
            "--num-invariant-pairs",
            "10",
            "--seed",
            "23",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert "Wrote chart-event results" in result.stdout
    assert (tmp_path / "data" / "observer_atlas_chart_events.csv").exists()
    assert (tmp_path / "data" / "observer_atlas_transition_summary.csv").exists()
    assert (tmp_path / "data" / "observer_atlas_loop_summary.csv").exists()
    assert (
        tmp_path
        / "figures"
        / "observer_atlas_transition_beta_error_vs_ticks.png"
    ).exists()


def test_exp15_exact_poincare_map_sanity_rows_are_near_exact() -> None:
    rows = run_experiment()
    transition_rows = [row for row in rows if row["kind"] == "transition"]
    loop_rows = [row for row in rows if row["kind"] == "loop"]

    assert len(transition_rows) == 3
    assert len(loop_rows) == 1
    assert max(abs(float(row["fitted_beta_error"])) for row in transition_rows) < 1e-3
    assert max(float(row["rmse"]) for row in transition_rows) < 1e-3
    assert abs(float(loop_rows[0]["beta_composition_error"])) < 1e-3
    assert float(loop_rows[0]["translation_composition_error_norm"]) < 1e-3
