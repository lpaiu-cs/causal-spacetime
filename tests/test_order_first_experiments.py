from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp30_ordinal_exact_sanity import run_experiment


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


def test_exp25_radar_return_order_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp25_radar_return_distance_order.py",
            "--tick-values",
            "16",
            "32",
            "--target-counts",
            "20",
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
        timeout=60,
    )

    assert "Wrote radar-return distance order data" in result.stdout
    assert (tmp_path / "data" / "radar_return_distance_order.csv").exists()
    assert (
        tmp_path / "figures" / "radar_return_order_inversion_vs_ticks.png"
    ).exists()


def test_exp26_metric_representation_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    env = _env_with_src()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp26_metric_representation_scale_invariance.py",
        ],
        check=True,
        cwd=os.getcwd(),
        env={**env, "MPLBACKEND": "Agg"},
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert "Wrote metric representation scale-invariance data" in result.stdout


def test_exp27_ratio_stability_cli_smoke() -> None:
    _skip_without_matplotlib()

    result = subprocess.run(
        [sys.executable, "experiments/exp27_ratio_stability_from_calibration.py"],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert "Wrote ratio stability data" in result.stdout


def test_exp28_oriented_chart_order_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp28_oriented_chart_distance_order_preservation.py",
            "--n-values",
            "60",
            "--tick-values",
            "16",
            "--repetitions",
            "1",
            "--pair-count",
            "20",
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

    assert "Wrote oriented chart distance-order data" in result.stdout
    assert (
        tmp_path / "data" / "oriented_chart_distance_order_preservation.csv"
    ).exists()


def test_exp29_metric_representability_cli_smoke() -> None:
    result = subprocess.run(
        [sys.executable, "experiments/exp29_metric_representability_diagnostics.py"],
        check=True,
        cwd=os.getcwd(),
        env=_env_with_src(),
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert "Wrote metric representability diagnostics" in result.stdout


def test_exp30_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
