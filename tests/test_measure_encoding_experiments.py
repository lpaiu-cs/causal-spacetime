from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp24_measure_sprinkling_exact_sanity import run_experiment


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


def test_exp21_physical_measure_sprinkling_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp21_physical_measure_sprinkling.py",
            "--n-values",
            "40",
            "--repetitions",
            "1",
            "--pairs-per-repetition",
            "5",
            "--seed",
            "43",
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

    assert "Wrote pair results" in result.stdout
    assert (tmp_path / "data" / "physical_measure_sprinkling_pairs.csv").exists()
    assert (tmp_path / "data" / "physical_measure_sprinkling_summary.csv").exists()
    assert (tmp_path / "figures" / "physical_measure_volume_scatter.png").exists()


def test_exp22_local_measure_profile_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp22_local_measure_profile_estimation.py",
            "--n-values",
            "40",
            "--repetitions",
            "1",
            "--num-bins",
            "6",
            "--seed",
            "47",
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

    assert "Wrote bin results" in result.stdout
    assert (tmp_path / "data" / "local_measure_profile_bins.csv").exists()
    assert (tmp_path / "data" / "local_measure_profile_summary.csv").exists()
    assert (tmp_path / "figures" / "local_measure_profile_shape.png").exists()


def test_exp23_thinning_stability_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp23_thinning_coarse_graining_stability.py",
            "--N",
            "120",
            "--repetitions",
            "1",
            "--keep-probabilities",
            "1.0",
            "0.5",
            "--pairs-per-repetition",
            "5",
            "--seed",
            "53",
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

    assert "Wrote pair results" in result.stdout
    assert (tmp_path / "data" / "thinning_coarse_graining_pairs.csv").exists()
    assert (tmp_path / "data" / "thinning_coarse_graining_summary.csv").exists()
    assert (
        tmp_path / "figures" / "thinning_volume_error_vs_keep_probability.png"
    ).exists()


def test_exp24_exact_measure_sanity_rows_are_accurate() -> None:
    rows = run_experiment(T=2.0, scale=1.5)

    assert len(rows) == 5
    assert max(float(row["absolute_error"]) for row in rows) < 1e-3
