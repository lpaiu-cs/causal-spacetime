"""End-to-end smoke and gate-logic tests for the PC-V1 pipeline."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from causal_spacetime_lab.positive_control.echo_profiles import (
    measure_bracket_echo_profiles,
    shuffle_profile_columns,
)
from causal_spacetime_lab.positive_control.gates import (
    GateThresholds,
    RepresentabilityFitPolicy,
    analyze_profiles,
    evaluate_gates,
    load_frozen_thresholds,
    propose_thresholds,
    save_thresholds,
)
from causal_spacetime_lab.positive_control.rewire import (
    geometry_free_scene,
    transitive_closure,
)
from causal_spacetime_lab.positive_control.scene import (
    PositiveControlSceneConfig,
    build_positive_control_scene,
)

SMOKE_CONFIG = PositiveControlSceneConfig(
    n_events=400,
    chain_positions=(-0.2, -0.07, 0.07, 0.2),
    ticks_per_chain=32,
    target_band_t=0.10,
    target_band_x=0.20,
    max_targets=16,
    min_targets=8,
    min_bracketing_chains=4,
    seed=3,
)

# Steps/restarts must stay at or above the measured adequacy floor: with
# 200 steps the structured fit reads as a null fit (heldout ~0.33), while
# 800 steps recovers the planted geometry (heldout ~0.04) — the exact
# underpowered-optimizer failure mode documented in the program review.
SMOKE_POLICY = RepresentabilityFitPolicy(
    embedding_dims=(1,),
    steps=800,
    restarts=3,
    stability_fits=2,
    train_constraints=400,
    heldout_constraints=120,
    stability_comparisons=2000,
    truth_comparisons=4000,
)


def test_structured_beats_destructive_null_and_recovers_truth_order() -> None:
    scene = build_positive_control_scene(SMOKE_CONFIG)
    profiles = measure_bracket_echo_profiles(scene)
    truth_x = scene.target_coords[:, 1]

    structured = analyze_profiles(profiles, truth_x, SMOKE_POLICY, seed=3)[0]
    shuffled = analyze_profiles(
        shuffle_profile_columns(profiles, seed=14), None, SMOKE_POLICY, seed=3
    )[0]

    # Smoke scale (4 chains): centering shrinks the absolute separation; the
    # full 6-chain separation is validated in Stage A. Assert the ordering.
    assert structured["heldout_violation"] < 0.10
    assert shuffled["heldout_violation"] > structured["heldout_violation"] + 0.05
    assert structured["truth_distance_order_error"] < 0.25
    assert structured["restart_order_disagreement"] < 0.35
    assert np.isnan(shuffled["truth_distance_order_error"])


def test_geometry_free_control_order_is_transitive_with_chains_intact() -> None:
    from causal_spacetime_lab.positive_control.rewire import (
        geometric_post_closure_density,
    )

    scene = build_positive_control_scene(SMOKE_CONFIG)
    control, density = geometry_free_scene(scene, seed=5)
    closed = control.causal
    assert np.array_equal(closed, transitive_closure(closed))
    for chain in control.chain_index_arrays:
        assert bool(np.all(closed[chain[:-1], chain[1:]]))
    assert control.causal_digest != scene.causal_digest

    # The repaired control matches the geometric post-closure density rather
    # than percolating to a near-complete (degenerate) order.
    geometric_density, _, _ = geometric_post_closure_density(scene)
    assert abs(density - geometric_density) <= 0.05
    assert density < 0.9


def test_threshold_proposal_respects_hard_floors_and_round_trips(
    tmp_path,
) -> None:
    structured_rows = [
        {
            "embedding_dim": 1.0,
            "heldout_violation": 0.9,
            "restart_order_disagreement": 0.9,
            "truth_distance_order_error": 0.12,
        }
    ]
    thresholds = propose_thresholds(structured_rows, [0.5], gate_dim=1)
    assert thresholds.heldout_max == pytest.approx(0.35)
    assert thresholds.stability_max == pytest.approx(0.30)
    assert thresholds.null_gap_min == pytest.approx(0.50)
    assert thresholds.truth_error_max == pytest.approx(0.15)

    path = tmp_path / "thresholds.json"
    save_thresholds(thresholds, path)
    loaded = load_frozen_thresholds(path)
    assert loaded == dataclasses.replace(thresholds, frozen_commit="UNFROZEN")


def test_stage_bc_refuse_to_run_without_frozen_thresholds(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        load_frozen_thresholds(tmp_path / "missing.json")


def test_gate_evaluation_directions() -> None:
    thresholds = GateThresholds(
        gate_dim=1,
        heldout_max=0.2,
        null_gap_min=0.1,
        stability_max=0.3,
        truth_error_max=0.15,
    )
    passing = {
        "heldout_violation": 0.1,
        "restart_order_disagreement": 0.2,
        "truth_distance_order_error": 0.1,
    }
    result = evaluate_gates(passing, 0.3, thresholds, apply_truth_gate=True)
    assert result["all_gates_pass"] == 1.0

    blocked = evaluate_gates(
        {**passing, "heldout_violation": 0.5}, 0.0, thresholds, apply_truth_gate=False
    )
    assert blocked["g1_heldout_pass"] == 0.0
    assert blocked["g2_null_gap_pass"] == 0.0
    assert blocked["all_gates_pass"] == 0.0
