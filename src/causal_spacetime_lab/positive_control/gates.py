"""Representability analysis, metrics, and gate evaluation for PC-V1.

Metric formulas follow docs/prereg/pc_v1_positive_control.md Section 8.
The positional-argsort disagreement metric is deliberately absent
(forbidden; see Section 7 of the preregistration).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from causal_spacetime_lab.embedding_stability import pairwise_order_stability
from causal_spacetime_lab.ordinal_embedding import (
    embedding_distance_order_error,
    fit_ordinal_embedding_gradient_descent,
    quadruplet_violation_rate,
)
from causal_spacetime_lab.positive_control.dissimilarity import (
    build_constraint_split,
    flip_constraint_sides,
    margin_from_probe_quantile,
    profile_dissimilarity_matrix,
)
from causal_spacetime_lab.positive_control.echo_profiles import EchoProfileMatrix


class InstrumentIntegrityError(RuntimeError):
    """Raised when the executed fit budget differs from the requested one."""


@dataclass(frozen=True)
class RepresentabilityFitPolicy:
    """Preregistered fit and analysis policy (PC-V1 Sections 6-8)."""

    embedding_dims: tuple[int, ...] = (1, 2, 3)
    gate_dim: int = 1
    steps: int = 1500
    restarts: int = 5
    learning_rate: float = 0.05
    stability_fits: int = 4
    train_constraints: int = 4000
    heldout_constraints: int = 1000
    margin_quantile: float = 0.25
    min_common_columns: int = 4
    pair_train_fraction: float = 0.8
    stability_comparisons: int = 5000
    truth_comparisons: int = 10_000


def _checked_fit(
    n_points: int,
    dim: int,
    constraints: NDArray[np.int_],
    policy: RepresentabilityFitPolicy,
    seed: int,
    restarts: int,
) -> NDArray[np.float64]:
    coords, diagnostics = fit_ordinal_embedding_gradient_descent(
        n_points,
        dim,
        constraints,
        steps=policy.steps,
        learning_rate=policy.learning_rate,
        seed=seed,
        restarts=restarts,
    )
    if diagnostics["steps"] != float(policy.steps) or diagnostics[
        "restarts"
    ] != float(restarts):
        raise InstrumentIntegrityError(
            "executed fit budget differs from requested budget"
        )
    return coords


def analyze_profiles(
    profiles: EchoProfileMatrix,
    truth_x: NDArray[np.float64] | None,
    policy: RepresentabilityFitPolicy,
    seed: int,
    include_flip_control: bool = False,
) -> list[dict[str, float]]:
    """Run the representability analysis on one profile matrix.

    Returns one row per embedding dimension. ``truth_x`` is the true spatial
    coordinate column for geometric scenes, or None for controls.
    """

    n = profiles.target_count
    dissimilarity = profile_dissimilarity_matrix(
        profiles,
        policy.min_common_columns,
    )
    margin = margin_from_probe_quantile(
        dissimilarity,
        quantile=policy.margin_quantile,
        seed=seed + 3,
    )
    split = build_constraint_split(
        dissimilarity,
        policy.train_constraints,
        policy.heldout_constraints,
        margin,
        train_fraction=policy.pair_train_fraction,
        seed=seed + 5,
    )

    rows: list[dict[str, float]] = []
    for dim in policy.embedding_dims:
        fit_seed = seed + 100 * dim
        coords = _checked_fit(
            n, dim, split.train, policy, fit_seed, policy.restarts
        )
        train_violation = quadruplet_violation_rate(coords, split.train)
        heldout_violation = quadruplet_violation_rate(coords, split.heldout)

        stability_embeddings = [coords]
        for replica in range(policy.stability_fits):
            stability_embeddings.append(
                _checked_fit(
                    n,
                    dim,
                    split.train,
                    policy,
                    fit_seed + 1000 * (replica + 1),
                    restarts=1,
                )
            )
        stability = pairwise_order_stability(
            stability_embeddings,
            num_pair_comparisons=policy.stability_comparisons,
            seed=seed + 7,
        )

        truth_error = float("nan")
        if truth_x is not None:
            truth_error = embedding_distance_order_error(
                coords,
                np.asarray(truth_x, dtype=float).reshape(n, 1),
                num_pair_comparisons=policy.truth_comparisons,
                seed=seed + 9,
            )

        flip_heldout = float("nan")
        if include_flip_control:
            flipped = flip_constraint_sides(split.train, seed=seed + 13)
            flip_coords = _checked_fit(
                n, dim, flipped, policy, fit_seed + 77, policy.restarts
            )
            flip_heldout = quadruplet_violation_rate(flip_coords, split.heldout)

        rows.append(
            {
                "embedding_dim": float(dim),
                "margin": float(margin),
                "train_constraint_count": float(split.train.shape[0]),
                "heldout_constraint_count": float(split.heldout.shape[0]),
                "train_violation": float(train_violation),
                "heldout_violation": float(heldout_violation),
                "generalization_gap": float(heldout_violation - train_violation),
                "restart_order_disagreement": float(
                    stability["mean_order_disagreement"]
                ),
                "truth_distance_order_error": float(truth_error),
                "flip_control_heldout": float(flip_heldout),
            }
        )
    return rows


@dataclass(frozen=True)
class GateThresholds:
    """Frozen gate thresholds (PC-V1 Section 10)."""

    gate_dim: int
    heldout_max: float
    null_gap_min: float
    stability_max: float
    truth_error_max: float
    frozen_commit: str = "UNFROZEN"


HARD_HELDOUT_CEILING = 0.35
HARD_NULL_GAP_FLOOR = 0.10
HARD_STABILITY_CEILING = 0.30


def propose_thresholds(
    structured_rows: list[dict[str, float]],
    null_gaps: list[float],
    gate_dim: int,
) -> GateThresholds:
    """Mechanical Stage A threshold construction (PC-V1 Section 10)."""

    at_dim = [
        row for row in structured_rows if int(row["embedding_dim"]) == gate_dim
    ]
    if not at_dim or not null_gaps:
        raise ValueError("no Stage A rows at the gate dimension")
    heldout = np.asarray([row["heldout_violation"] for row in at_dim])
    stability = np.asarray([row["restart_order_disagreement"] for row in at_dim])
    truth = np.asarray([row["truth_distance_order_error"] for row in at_dim])
    gaps = np.asarray(null_gaps, dtype=float)

    def round_up(value: float) -> float:
        return float(np.ceil(value / 0.05) * 0.05)

    def round_down(value: float) -> float:
        return float(np.floor(value / 0.05) * 0.05)

    return GateThresholds(
        gate_dim=gate_dim,
        heldout_max=min(
            HARD_HELDOUT_CEILING, round_up(float(np.quantile(heldout, 0.9)))
        ),
        null_gap_min=max(
            HARD_NULL_GAP_FLOOR, round_down(float(np.quantile(gaps, 0.1)))
        ),
        stability_max=min(
            HARD_STABILITY_CEILING, round_up(float(np.quantile(stability, 0.9)))
        ),
        truth_error_max=round_up(float(np.quantile(truth, 0.9))),
    )


def evaluate_gates(
    structured_row: dict[str, float],
    null_gap: float,
    thresholds: GateThresholds,
    apply_truth_gate: bool,
) -> dict[str, float | str]:
    """Evaluate G1-G4 for one seed at the gate dimension."""

    g1 = structured_row["heldout_violation"] <= thresholds.heldout_max
    g2 = null_gap >= thresholds.null_gap_min
    g3 = (
        structured_row["restart_order_disagreement"] <= thresholds.stability_max
    )
    g4 = True
    if apply_truth_gate:
        g4 = (
            structured_row["truth_distance_order_error"]
            <= thresholds.truth_error_max
        )
    return {
        "g1_heldout_pass": float(g1),
        "g2_null_gap_pass": float(g2),
        "g3_stability_pass": float(g3),
        "g4_truth_pass": float(g4),
        "all_gates_pass": float(g1 and g2 and g3 and g4),
        "null_gap": float(null_gap),
    }


def save_thresholds(thresholds: GateThresholds, path: Path) -> None:
    """Write a thresholds JSON (used by the freeze procedure)."""

    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "gate_dim": thresholds.gate_dim,
        "heldout_max": thresholds.heldout_max,
        "null_gap_min": thresholds.null_gap_min,
        "stability_max": thresholds.stability_max,
        "truth_error_max": thresholds.truth_error_max,
        "frozen_commit": thresholds.frozen_commit,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def load_frozen_thresholds(path: Path) -> GateThresholds:
    """Load frozen thresholds; refuse to proceed if absent."""

    if not path.exists():
        raise FileNotFoundError(
            f"frozen thresholds not found at {path}; Stage B/C may only run "
            "after the PC-V1 freeze procedure (preregistration Section 12)"
        )
    payload = json.loads(path.read_text(encoding="utf-8"))
    return GateThresholds(
        gate_dim=int(payload["gate_dim"]),
        heldout_max=float(payload["heldout_max"]),
        null_gap_min=float(payload["null_gap_min"]),
        stability_max=float(payload["stability_max"]),
        truth_error_max=float(payload["truth_error_max"]),
        frozen_commit=str(payload["frozen_commit"]),
    )
