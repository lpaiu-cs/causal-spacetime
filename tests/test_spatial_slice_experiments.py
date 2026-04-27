from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp45_spatial_slice_exact_sanity import run_experiment


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    src_path = os.path.abspath("src")
    env["PYTHONPATH"] = (
        src_path
        if not env.get("PYTHONPATH")
        else f"{src_path}{os.pathsep}{env['PYTHONPATH']}"
    )
    return env


def _skip_without_matplotlib() -> None:
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")


def test_exp41_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp41_radar_time_order_from_brackets.py",
            "--n-values",
            "80",
            "--tick-values",
            "16",
            "--repetitions",
            "1",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=90,
    )

    assert "Wrote radar-time order data" in result.stdout
    assert (tmp_path / "data" / "radar_time_order_from_brackets.csv").exists()


def test_exp42_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp42_same_slice_distance_order_preservation.py",
            "--n-values",
            "80",
            "--tick-values",
            "16",
            "--bin-width-values",
            "2",
            "--repetitions",
            "1",
            "--max-pairs-per-slice",
            "20",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=90,
    )

    assert "Wrote same-slice distance-order data" in result.stdout
    assert (tmp_path / "data" / "same_slice_distance_order_preservation.csv").exists()


def test_exp43_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp43_sliced_observer_order_null_baseline.py",
            "--N",
            "90",
            "--tick-values",
            "16",
            "--bin-width-values",
            "2",
            "--constraint-counts",
            "40",
            "--repetitions",
            "1",
            "--steps",
            "20",
            "--restarts",
            "1",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=90,
    )

    assert "Wrote sliced observer-order null baseline data" in result.stdout
    assert (tmp_path / "data" / "sliced_observer_order_null_baseline.csv").exists()


def test_exp44_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp44_slice_width_sensitivity.py",
            "--N",
            "90",
            "--tick-count",
            "16",
            "--bin-width-values",
            "1",
            "2",
            "--repetitions",
            "1",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=90,
    )

    assert "Wrote slice-width sensitivity data" in result.stdout
    assert (tmp_path / "data" / "slice_width_sensitivity.csv").exists()


def test_exp45_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
