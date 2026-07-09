"""Compute the shared-scalar confound before/after data for Paper B Figure 3.

On the PC-V1 Stage C specificity seeds (200-219 range; the frozen run used
200-209), take the same density-matched geometry-free random-order scenes and
run the identical d=1 pipeline twice: once on RAW bracket-width dissimilarity
(center_references=False) and once on PARALLAX dissimilarity
(center_references=True). The parallax column reproduces the frozen Stage C
random_order result; the raw column is the counterfactual that would have
false-passed the gate. Deterministic; writes confound_data.csv.

Usage: python docs/paper/paper_b/figures/compute_confound_data.py
"""

from __future__ import annotations

import csv
from pathlib import Path

from causal_spacetime_lab.ordinal_embedding import (
    fit_ordinal_embedding_gradient_descent,
    quadruplet_violation_rate,
)
from causal_spacetime_lab.positive_control.dissimilarity import (
    build_constraint_split,
    margin_from_probe_quantile,
    profile_dissimilarity_matrix,
)
from causal_spacetime_lab.positive_control.echo_profiles import (
    measure_bracket_echo_profiles,
)
from causal_spacetime_lab.positive_control.gates import RepresentabilityFitPolicy
from causal_spacetime_lab.positive_control.rewire import (
    geometric_post_closure_density,
    geometry_free_scene,
)
from causal_spacetime_lab.positive_control.scene import (
    PositiveControlSceneConfig,
    SceneValidityError,
    build_positive_control_scene,
)

SEEDS = range(200, 210)  # PC-V1 Stage C specificity seeds
OUT = Path("docs/paper/paper_b/figures/confound_data.csv")


def _heldout_at_d1(profiles, policy, seed, center: bool) -> float:
    dissimilarity = profile_dissimilarity_matrix(
        profiles, policy.min_common_columns, center_references=center
    )
    margin = margin_from_probe_quantile(
        dissimilarity, quantile=policy.margin_quantile, seed=seed + 3
    )
    split = build_constraint_split(
        dissimilarity,
        policy.train_constraints,
        policy.heldout_constraints,
        margin,
        train_fraction=policy.pair_train_fraction,
        seed=seed + 5,
    )
    coords, _ = fit_ordinal_embedding_gradient_descent(
        profiles.target_count,
        1,
        split.train,
        steps=policy.steps,
        learning_rate=policy.learning_rate,
        seed=seed + 100,
        restarts=policy.restarts,
    )
    return quadruplet_violation_rate(coords, split.heldout)


def main() -> None:
    policy = RepresentabilityFitPolicy()
    rows: list[dict] = []
    for seed in SEEDS:
        scene = build_positive_control_scene(PositiveControlSceneConfig(seed=seed))
        target_density, _, _ = geometric_post_closure_density(scene)
        try:
            control, _ = geometry_free_scene(scene, seed + 41)
        except SceneValidityError as error:
            print(f"seed {seed}: skipped ({error})")
            continue
        profiles = measure_bracket_echo_profiles(control)
        raw = _heldout_at_d1(profiles, policy, seed, center=False)
        parallax = _heldout_at_d1(profiles, policy, seed, center=True)
        rows.append({"seed": seed, "raw": raw, "parallax": parallax})
        print(f"seed {seed}: raw(false-pass)={raw:.3f}  parallax(block)={parallax:.3f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["seed", "raw", "parallax"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {OUT} ({len(rows)} seeds)")


if __name__ == "__main__":
    main()
