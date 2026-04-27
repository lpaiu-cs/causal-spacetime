from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import pytest

from causal_spacetime_lab.dimension import (
    estimate_myrheim_meyer_dimension,
    myrheim_meyer_relation_fraction,
    relation_fraction,
)


def test_relation_fraction_convention() -> None:
    C = np.array(
        [
            [False, True, True],
            [False, False, True],
            [False, False, False],
        ]
    )

    assert relation_fraction(C) == pytest.approx(3 / 6)


def test_myrheim_meyer_curve_known_d2_value() -> None:
    assert myrheim_meyer_relation_fraction(2.0) == pytest.approx(0.25)


@pytest.mark.parametrize("dimension", [2.0, 3.0, 4.0])
def test_dimension_inversion_recovers_exact_theoretical_values(
    dimension: float,
) -> None:
    fraction = myrheim_meyer_relation_fraction(dimension)

    estimate = estimate_myrheim_meyer_dimension(fraction)

    assert estimate == pytest.approx(dimension, abs=1e-5)


def test_dimension_reconstruction_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
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
            "experiments/exp10_dimension_reconstruction.py",
            "--dims",
            "2",
            "--n-values",
            "50",
            "--repetitions",
            "1",
            "--seed",
            "13",
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

    assert "Wrote results" in result.stdout
    assert (tmp_path / "data" / "dimension_reconstruction_results.csv").exists()
    assert (tmp_path / "data" / "dimension_reconstruction_summary.csv").exists()
    assert (tmp_path / "figures" / "dimension_estimate_vs_N.png").exists()

