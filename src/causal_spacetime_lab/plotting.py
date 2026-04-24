"""Plotting helpers for experiment outputs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def save_timelike_reconstruction_plot(
    rows: Sequence[Mapping[str, float]],
    output_path: str | Path,
) -> Path:
    """Save an error plot for the timelike reconstruction experiment."""

    if not rows:
        raise ValueError("rows must not be empty")

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    n_events = np.asarray([row["n_events"] for row in rows], dtype=float)
    interval_errors = np.asarray(
        [row["interval_tau_abs_error"] for row in rows],
        dtype=float,
    )
    chain_errors = np.asarray(
        [row["chain_tau_abs_error"] for row in rows],
        dtype=float,
    )

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(n_events, interval_errors, marker="o", label="Interval cardinality")
    ax.plot(n_events, chain_errors, marker="s", label="Longest chain estimate")
    ax.set_xscale("log")
    ax.set_xlabel("Sprinkled interior events")
    ax.set_ylabel("Absolute proper-time error")
    ax.set_title("Timelike proper-time reconstruction in a 1+1D causal diamond")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)

    return path

