from __future__ import annotations

import os
import subprocess
import sys

import pytest


def test_discrete_radar_experiment_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")

    env = os.environ.copy()
    src_path = os.path.abspath("src")
    env["PYTHONPATH"] = (
        src_path
        if not env.get("PYTHONPATH")
        else f"{src_path}{os.pathsep}{env['PYTHONPATH']}"
    )

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp11_discrete_observer_radar_reconstruction.py",
            "--n-values",
            "20",
            "--tick-values",
            "8",
            "--repetitions",
            "1",
            "--seed",
            "17",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=os.getcwd(),
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert "Wrote event results" in result.stdout
    assert (tmp_path / "data" / "discrete_radar_reconstruction_events.csv").exists()
    assert (tmp_path / "data" / "discrete_radar_reconstruction_summary.csv").exists()
    assert (tmp_path / "figures" / "discrete_radar_time_scatter.png").exists()

