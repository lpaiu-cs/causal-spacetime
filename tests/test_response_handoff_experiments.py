from __future__ import annotations

import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.abspath("experiments"))

from experiments.exp133_response_handoff_exact_sanity import (
    run_experiment as run_exp133,
)
from experiments.exp138_response_handoff_preregistration_rules import (
    run_experiment as run_exp138,
)
from experiments.exp139_response_handoff_manifest_read_exact_sanity import (
    run_experiment as run_exp139,
)


def _env_with_src() -> dict[str, str]:
    env = os.environ.copy()
    paths = [os.path.abspath("src"), os.path.abspath("experiments")]
    if env.get("PYTHONPATH"):
        paths.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


def _skip_without_matplotlib() -> None:
    try:
        import matplotlib.pyplot  # noqa: F401
    except Exception as exc:
        pytest.skip(f"matplotlib is not importable in this interpreter: {exc}")


def test_exp133_exact_sanity(tmp_path) -> None:  # type: ignore[no-untyped-def]
    rows = run_exp133(tmp_path)

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)


def test_exp134_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp134_response_handoff_threshold_sensitivity.py",
            "--reference-length",
            "48",
            "--emission-positions",
            "6",
            "12",
            "18",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--min-agreement-values",
            "0.7",
            "--min-null-z-values",
            "0.0",
            "--min-bootstrap-values",
            "0.6",
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
    assert "Wrote response handoff threshold sensitivity" in result.stdout


def test_exp135_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp135_response_handoff_protocol_selection.py",
            "--reference-length",
            "48",
            "--emission-positions",
            "6",
            "12",
            "18",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--min-margins",
            "0.0",
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
    assert "Wrote response handoff protocol selection" in result.stdout


def test_exp136_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp136_response_handoff_manifest_export.py",
            "--reference-length",
            "48",
            "--emission-positions",
            "6",
            "12",
            "18",
            "--layer-delay-ranks",
            "3",
            "5",
            "--targets-per-layer",
            "3",
            "--max-manifests",
            "2",
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
    assert "Wrote response handoff manifest export" in result.stdout


def test_exp137_cli_smoke(tmp_path) -> None:  # type: ignore[no-untyped-def]
    _skip_without_matplotlib()
    result = subprocess.run(
        [
            sys.executable,
            "experiments/exp137_response_handoff_failure_catalog.py",
            "--target-count",
            "20",
            "--protocol-count-values",
            "2",
            "--reachable-probabilities",
            "0.3",
            "--unique-rank-counts",
            "3",
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
    assert "Wrote response handoff failure catalog" in result.stdout


def test_exp138_table_experiment() -> None:
    rows = run_exp138()

    assert rows
    assert any(row["rule_name"] == "heldout_not_used_for_fitting" for row in rows)


def test_exp139_exact_sanity(tmp_path) -> None:  # type: ignore[no-untyped-def]
    rows = run_exp139(tmp_path)

    assert rows
    assert all(float(row["passed"]) == 1.0 for row in rows)
