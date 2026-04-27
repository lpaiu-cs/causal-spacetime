from __future__ import annotations

import os
import subprocess
import sys

import pytest

from experiments.exp66_persistence_matching_exact_sanity import run_experiment


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


def test_exp60_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp60_persistence_predicate_undefined.py",
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

    assert "Wrote persistence predicate data" in result.stdout
    assert (tmp_path / "data" / "persistence_predicate_undefined.csv").exists()


def test_exp61_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp61_symmetric_persistence_ambiguity.py",
            "--slice-count",
            "3",
            "--object-counts",
            "4",
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

    assert "Wrote symmetric persistence ambiguity data" in result.stdout
    assert (tmp_path / "data" / "symmetric_persistence_ambiguity.csv").exists()


def test_exp62_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp62_relational_persistence_matching_recovery.py",
            "--slice-count",
            "5",
            "--object-counts",
            "4",
            "--motion-scales",
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
        timeout=60,
    )

    assert "Wrote relational persistence matching data" in result.stdout
    output = tmp_path / "data" / "relational_persistence_matching_recovery.csv"
    assert output.exists()


def test_exp63_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp63_partial_label_constrained_persistence.py",
            "--slice-count",
            "5",
            "--object-count",
            "5",
            "--known-fractions",
            "0.0",
            "0.5",
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

    assert "Wrote partial-label persistence data" in result.stdout
    assert (tmp_path / "data" / "partial_label_constrained_persistence.csv").exists()


def test_exp64_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp64_crossing_persistence_failure.py",
            "--slice-count",
            "6",
            "--object-count",
            "5",
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

    assert "Wrote crossing persistence failure data" in result.stdout
    assert (tmp_path / "data" / "crossing_persistence_failure.csv").exists()


def test_exp65_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()

    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp65_persistence_hypothesis_dependence.py",
            "--slice-count",
            "5",
            "--object-count",
            "5",
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

    assert "Wrote persistence-hypothesis dependence data" in result.stdout
    assert (tmp_path / "data" / "persistence_hypothesis_dependence.csv").exists()


def test_exp66_exact_sanity_rows_pass() -> None:
    rows = run_experiment()

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
