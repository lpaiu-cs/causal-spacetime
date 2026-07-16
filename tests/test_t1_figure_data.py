"""Regression tests for the tracked theory-figure data table.

The committed docs/theory/t1_figure_data.json must regenerate from the
script at its pinned configurations (numerically, to floating tolerance,
so the check is platform-independent; the byte-identical statement of
manuscript Section 12 is per-machine), and the two claims the
manuscript's Figure 6 visualizes -- zero band violations and the ~1/K
RMSE slope -- must hold on the regenerated data. Everything here is
deterministic; a failure means drift, not noise.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "theory"
sys.path.insert(0, str(EXPERIMENT_DIR))

from t1_figure_data import (  # noqa: E402
    RESULTS_PATH,
    band_table,
    resolution_table,
)


def _assert_close(regenerated, committed, path=""):
    """Structural equality with float tolerance (platform-independent)."""

    if isinstance(committed, dict):
        assert isinstance(regenerated, dict), path
        assert regenerated.keys() == committed.keys(), path
        for key in committed:
            _assert_close(regenerated[key], committed[key], f"{path}.{key}")
    elif isinstance(committed, list):
        assert len(regenerated) == len(committed), path
        for i, (a, b) in enumerate(zip(regenerated, committed, strict=True)):
            _assert_close(a, b, f"{path}[{i}]")
    elif isinstance(committed, float) and not isinstance(committed, bool):
        assert math.isclose(regenerated, committed, rel_tol=1e-9, abs_tol=1e-9), (
            path, regenerated, committed,
        )
    else:
        assert regenerated == committed, (path, regenerated, committed)


def test_band_table_matches_the_committed_artifact():
    regenerated = band_table()
    committed = json.loads(RESULTS_PATH.read_text(encoding="utf-8"))["band"]
    _assert_close(regenerated, committed)
    assert regenerated["violations"] == 0
    assert regenerated["n_reachable"] == 400
    residuals = np.asarray(regenerated["residual"], dtype=float)
    assert residuals.min() >= -1.0 and residuals.max() < 1.0


def test_resolution_table_matches_the_committed_artifact():
    regenerated = resolution_table()
    committed = json.loads(RESULTS_PATH.read_text(encoding="utf-8"))[
        "resolution"
    ]
    _assert_close(regenerated, committed)
    assert -1.3 < regenerated["rmse_loglog_slope"] < -0.7
    ladder = regenerated["ladder"]
    assert [row["ticks"] for row in ladder] == [12, 24, 48, 96, 192, 384]
    assert all(row["within_bound"] for row in ladder)
    rmse = [row["rmse"] for row in ladder]
    assert rmse == sorted(rmse, reverse=True)
