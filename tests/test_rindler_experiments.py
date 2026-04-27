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


def test_exp16_rindler_horizon_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp16_rindler_horizon_reconstruction.py",
            "--T",
            "4.0",
            "--accelerations",
            "2.0",
            "--n-values",
            "30",
            "--tick-values",
            "8",
            "--repetitions",
            "1",
            "--seed",
            "29",
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
    assert (tmp_path / "data" / "rindler_horizon_reconstruction_events.csv").exists()
    assert (tmp_path / "data" / "rindler_horizon_reconstruction_summary.csv").exists()
    assert (tmp_path / "figures" / "rindler_accessibility_map.png").exists()


def test_exp17_inertial_vs_rindler_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp17_inertial_vs_rindler_accessibility.py",
            "--T",
            "4.0",
            "--n-events",
            "30",
            "--tick-count",
            "16",
            "--acceleration",
            "2.0",
            "--seed",
            "31",
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

    assert "Wrote comparison data" in result.stdout
    assert (tmp_path / "data" / "inertial_vs_rindler_accessibility.csv").exists()
    assert (tmp_path / "figures" / "inertial_vs_rindler_accessibility.png").exists()
