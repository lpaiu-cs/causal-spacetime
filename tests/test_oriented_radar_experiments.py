from __future__ import annotations

import os
import subprocess
import sys

import pytest


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    src_path = os.path.abspath("src")
    env["PYTHONPATH"] = (
        src_path
        if not env.get("PYTHONPATH")
        else f"{src_path}{os.pathsep}{env['PYTHONPATH']}"
    )
    return env


def test_single_observer_reflection_degeneracy_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")

    result = subprocess.run(
        [sys.executable, "experiments/exp12_single_observer_reflection_degeneracy.py"],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert "Wrote summary" in result.stdout
    assert os.path.exists("outputs/data/single_observer_reflection_degeneracy.csv")
    assert os.path.exists("outputs/figures/single_observer_reflection_degeneracy.png")


def test_oriented_radar_lorentz_map_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp13_oriented_radar_lorentz_map_recovery.py",
            "--n-values",
            "20",
            "--tick-values",
            "8",
            "--beta-values",
            "0.3",
            "--repetitions",
            "1",
            "--seed",
            "19",
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

    assert "Wrote event results" in result.stdout
    assert (tmp_path / "data" / "oriented_radar_lorentz_events.csv").exists()
    assert (tmp_path / "data" / "oriented_radar_lorentz_summary.csv").exists()
    assert (
        tmp_path / "figures" / "oriented_radar_lab_position_scatter.png"
    ).exists()
