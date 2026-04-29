"""Preregistered configuration for manifest-family comparison diagnostics."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FamilyComparisonConfig:
    """Fixed diagnostic thresholds for manifest-family comparison."""

    dims: list[int]
    steps: int
    restarts: int
    learning_rate: float
    null_repetitions: int
    restart_count: int
    heldout_violation_threshold: float
    generalization_gap_threshold: float
    min_null_separation: float
    max_restart_std: float
    max_latent_order_disagreement: float
    seed: int


def default_family_comparison_config() -> FamilyComparisonConfig:
    """Return preregistered default family-comparison settings."""

    return FamilyComparisonConfig(
        dims=[1, 2, 3],
        steps=800,
        restarts=3,
        learning_rate=0.05,
        null_repetitions=5,
        restart_count=8,
        heldout_violation_threshold=0.20,
        generalization_gap_threshold=0.10,
        min_null_separation=0.10,
        max_restart_std=0.10,
        max_latent_order_disagreement=0.30,
        seed=0,
    )
