"""Deterministic data table behind the manuscript's theory figure.

Paper B's Section 8 figure shows three verified T1 results: the Model-D
quantization band (Section 8.1), the 1/K resolution law (Section 8.3), and
the protocol-dependent density scaling (Section 8.4). The density panel
reads the already-tracked docs/theory/t1_g2_density_scaling_results.json;
this script materializes the other two panels' data by re-running the
corresponding t1_verification.py checks' computations at their pinned
configurations (band: seed 0, 400 targets, K = 96; resolution ladder:
seed 2, K = 12..384) and writes the tracked table
docs/theory/t1_figure_data.json (audit-trail convention of PR #8: the
committed artifact regenerates byte-identically from this script).

The script asserts the two claims the figure visualizes -- zero band
violations and a log-log RMSE slope within the pinned CI window -- so the
tracked table can never silently disagree with the manuscript.

Usage:
    python t1_figure_data.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from t1_verification import (
    _sample_reachable_targets,
    bracket_widths_ranks,
    check_resolution_scaling,
    predicted_center,
    residuals_in_band,
)

ROOT = Path(__file__).resolve().parents[2]
RESULTS_PATH = ROOT / "docs" / "theory" / "t1_figure_data.json"

BAND_CONFIG = {"seed": 0, "n_targets": 400, "span": 1.4, "ticks": 96, "x0": 0.0}
RESOLUTION_CONFIG = {
    "seed": 2,
    "n_targets": 200,
    "tick_ladder": [12, 24, 48, 96, 192, 384],
    "span": 1.4,
}


def band_table() -> dict:
    """Per-target band data at the pinned controlled configuration
    (identical sampling and measurement to check_quantization_band)."""

    cfg = BAND_CONFIG
    rng = np.random.default_rng(cfg["seed"])
    delta = cfg["span"] / (cfg["ticks"] - 1)
    targets = _sample_reachable_targets(
        rng, cfg["n_targets"], cfg["x0"], cfg["span"]
    )
    widths = bracket_widths_ranks(
        targets, cfg["x0"], cfg["span"], cfg["ticks"]
    )
    reachable = ~np.isnan(widths)
    abs_dx = np.abs(targets[reachable, 1] - cfg["x0"])
    measured = widths[reachable]
    center = predicted_center(abs_dx, delta)
    residuals = measured - center
    violations = int((~residuals_in_band(residuals)).sum())
    assert violations == 0, f"band violated: {violations}"
    return {
        "config": cfg,
        "delta": delta,
        "n_reachable": int(reachable.sum()),
        "violations": violations,
        "abs_dx_over_delta": (abs_dx / delta).round(6).tolist(),
        "measured_w": measured.round(6).tolist(),
        "residual": residuals.round(6).tolist(),
    }


def resolution_table() -> dict:
    """The K-ladder RMSE rows and slope (check_resolution_scaling verbatim)."""

    cfg = RESOLUTION_CONFIG
    result = check_resolution_scaling(
        seed=cfg["seed"],
        n_targets=cfg["n_targets"],
        tick_ladder=tuple(cfg["tick_ladder"]),
        span=cfg["span"],
    )
    assert result["passed"], result
    return {
        "config": cfg,
        "ladder": result["ladder"],
        "rmse_loglog_slope": result["rmse_loglog_slope"],
    }


def main() -> None:
    results = {
        "band": band_table(),
        "resolution": resolution_table(),
        "density_panel_source": "docs/theory/t1_g2_density_scaling_results.json",
        "status": (
            "deterministic regeneration of the CI-pinned T1 checks at their "
            "pinned configurations; figure data only, no new measurement"
        ),
    }
    RESULTS_PATH.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(
        f"band: {results['band']['n_reachable']} targets, "
        f"{results['band']['violations']} violations; "
        f"resolution slope {results['resolution']['rmse_loglog_slope']:.4f}"
    )
    print(f"wrote {RESULTS_PATH}")


if __name__ == "__main__":
    main()
