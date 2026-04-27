from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp53_cross_slice_transport_exact_sanity import run_experiment


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


def test_exp46_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp46_cross_slice_predicate_undefined.py",
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

    assert "Wrote cross-slice undefined predicate data" in result.stdout
    assert (tmp_path / "data" / "cross_slice_predicate_undefined.csv").exists()


def test_exp47_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp47_sliced_constraint_graph_decomposition.py",
            "--N",
            "80",
            "--tick-count",
            "16",
            "--bin-width-values",
            "2",
            "--constraint-count",
            "40",
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

    assert "Wrote sliced constraint graph data" in result.stdout
    assert (tmp_path / "data" / "sliced_constraint_graph_decomposition.csv").exists()


def test_exp48_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp48_slice_local_embedding_validation.py",
            "--N-values",
            "80",
            "--tick-values",
            "16",
            "--bin-width-values",
            "2",
            "--constraint-count",
            "40",
            "--repetitions",
            "1",
            "--steps",
            "10",
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

    assert "Wrote slice-local embedding validation data" in result.stdout
    assert (tmp_path / "data" / "slice_local_embedding_validation.csv").exists()


def test_exp49_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp49_slice_gauge_dependence.py",
            "--N",
            "80",
            "--tick-count",
            "16",
            "--bin-width",
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

    assert "Wrote slice gauge-dependence data" in result.stdout
    assert (tmp_path / "data" / "slice_gauge_dependence.csv").exists()


def test_exp50_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp50_anchor_constrained_transport.py",
            "--N",
            "80",
            "--tick-count",
            "16",
            "--bin-width",
            "2",
            "--constraint-count",
            "40",
            "--repetitions",
            "1",
            "--steps",
            "10",
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

    assert "Wrote anchor-constrained transport data" in result.stdout
    assert (tmp_path / "data" / "anchor_constrained_transport.csv").exists()


def test_exp51_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp51_persistence_dependent_velocity.py",
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

    assert "Wrote persistence-dependent velocity data" in result.stdout
    assert (tmp_path / "data" / "persistence_dependent_velocity.csv").exists()


def test_exp52_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp52_noisy_transport_sensitivity.py",
            "--N",
            "80",
            "--tick-count",
            "16",
            "--bin-width",
            "2",
            "--noise-levels",
            "0.0",
            "0.05",
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

    assert "Wrote noisy transport sensitivity data" in result.stdout
    assert (tmp_path / "data" / "noisy_transport_sensitivity.csv").exists()


def test_exp53_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
