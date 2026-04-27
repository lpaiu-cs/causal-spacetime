from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp20_conformal_volume_exact_sanity import run_experiment


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    src_path = os.path.abspath("src")
    env["PYTHONPATH"] = (
        src_path
        if not env.get("PYTHONPATH")
        else f"{src_path}{os.pathsep}{env['PYTHONPATH']}"
    )
    return env


def test_exp18_conformal_order_ambiguity_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp18_conformal_order_ambiguity.py",
            "--N",
            "40",
            "--seed",
            "37",
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

    assert "Wrote summary" in result.stdout
    assert (tmp_path / "data" / "conformal_order_ambiguity_summary.csv").exists()
    assert (tmp_path / "figures" / "conformal_order_ambiguity_scales.png").exists()


def test_exp19_weighted_conformal_volume_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp19_weighted_conformal_volume_reconstruction.py",
            "--n-values",
            "40",
            "--repetitions",
            "1",
            "--pairs-per-repetition",
            "5",
            "--seed",
            "41",
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

    assert "Wrote pair results" in result.stdout
    assert (tmp_path / "data" / "weighted_conformal_volume_pairs.csv").exists()
    assert (tmp_path / "data" / "weighted_conformal_volume_summary.csv").exists()
    assert (tmp_path / "figures" / "weighted_conformal_volume_scatter.png").exists()


def test_exp20_exact_conformal_sanity_rows_are_accurate() -> None:
    rows = run_experiment(T=2.0, scale=1.5)

    assert len(rows) == 3
    assert max(float(row["absolute_error"]) for row in rows) < 1e-3
