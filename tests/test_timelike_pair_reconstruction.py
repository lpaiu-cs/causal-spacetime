from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import pytest

from causal_spacetime_lab.causal import causal_matrix_1p1
from causal_spacetime_lab.estimators import estimate_tau_from_interval_count_1p1
from causal_spacetime_lab.intervals import alexandrov_interval_indices
from causal_spacetime_lab.metrics import minkowski_tau_1p1


def test_internal_pair_interval_reconstruction_on_deterministic_example() -> None:
    events = np.array(
        [
            [-0.6, 0.0],
            [0.0, 0.0],
            [0.6, 0.0],
        ]
    )
    C = causal_matrix_1p1(events)
    interval_count = alexandrov_interval_indices(C, 0, 2).size
    true_tau = minkowski_tau_1p1(events[0], events[2])
    rho = interval_count / (0.5 * true_tau**2)

    estimate = estimate_tau_from_interval_count_1p1(interval_count, rho)

    assert estimate == pytest.approx(true_tau)


def test_timelike_pair_experiment_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
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
            "experiments/exp07_timelike_pair_reconstruction_convergence.py",
            "--n-values",
            "40",
            "--repetitions",
            "1",
            "--pairs-per-repetition",
            "5",
            "--max-pairs-for-chain",
            "2",
            "--seed",
            "11",
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

    assert "Wrote pair results" in result.stdout
    assert (tmp_path / "data" / "timelike_pair_reconstruction_pairs.csv").exists()
    assert (tmp_path / "data" / "timelike_pair_reconstruction_summary.csv").exists()
    assert (
        tmp_path / "figures" / "timelike_pair_reconstruction_scatter.png"
    ).exists()
    assert (
        tmp_path / "figures" / "timelike_pair_reconstruction_error_vs_N.png"
    ).exists()
