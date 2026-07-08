"""Shared runner logic for the PC-V1 positive-control experiments."""

from __future__ import annotations

import csv
import math
import subprocess
from pathlib import Path

from causal_spacetime_lab.positive_control.echo_profiles import (
    measure_bracket_echo_profiles,
    relabel_targets,
    shuffle_profile_columns,
)
from causal_spacetime_lab.positive_control.gates import (
    RepresentabilityFitPolicy,
    analyze_profiles,
)
from causal_spacetime_lab.positive_control.rewire import geometry_free_scene
from causal_spacetime_lab.positive_control.scene import (
    PositiveControlSceneConfig,
    SceneValidityError,
    build_positive_control_scene,
)

FROZEN_THRESHOLDS_PATH = Path("docs/prereg/frozen/pc_v1_thresholds.json")
DEFAULT_OUTPUT_DIR = Path("outputs/positive_control")


def primary_scene_config(seed: int) -> PositiveControlSceneConfig:
    """Preregistered primary configuration (PC-V1 Section 4)."""

    return PositiveControlSceneConfig(seed=seed)


def smoke_scene_config(seed: int) -> PositiveControlSceneConfig:
    """Small engineering-smoke configuration (not preregistered; no gates)."""

    return PositiveControlSceneConfig(
        n_events=400,
        chain_positions=(-0.2, -0.07, 0.07, 0.2),
        ticks_per_chain=32,
        target_band_t=0.10,
        target_band_x=0.20,
        max_targets=16,
        min_targets=8,
        min_bracketing_chains=4,
        seed=seed,
    )


def smoke_fit_policy() -> RepresentabilityFitPolicy:
    """Reduced fit policy for engineering smoke runs only.

    Steps/restarts stay at the adequacy floor measured during skeleton
    validation (200 steps reads as a null fit; 800 steps recovers the
    signal) — only constraint counts and comparisons are reduced.
    """

    return RepresentabilityFitPolicy(
        embedding_dims=(1, 2),
        steps=800,
        restarts=3,
        stability_fits=2,
        train_constraints=500,
        heldout_constraints=150,
        stability_comparisons=2000,
        truth_comparisons=4000,
    )


def git_describe() -> str:
    """Return a short code-version string for provenance rows."""

    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def parse_seed_spec(spec: str) -> list[int]:
    """Parse '0-9' or '1,4,7' style seed specifications."""

    seeds: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part[1:]:
            low, high = part.split("-", 1)
            seeds.extend(range(int(low), int(high) + 1))
        else:
            seeds.append(int(part))
    if not seeds:
        raise ValueError("seed specification is empty")
    return seeds


def run_seed_conditions(
    config: PositiveControlSceneConfig,
    policy: RepresentabilityFitPolicy,
    conditions: list[str],
    stage: str,
    code_version: str,
) -> list[dict[str, float | str]]:
    """Run the requested conditions for one seed, returning flat rows.

    Scene-validity and constraint-pool failures are recorded as
    status != ok rows (preregistration Section 3), never silently skipped.
    """

    seed = config.seed
    base: dict[str, float | str] = {
        "stage": stage,
        "seed": float(seed),
        "code_version": code_version,
    }
    try:
        scene = build_positive_control_scene(config)
    except SceneValidityError as error:
        return [
            {
                **base,
                "condition": condition,
                "status": f"scene_invalid: {error}",
            }
            for condition in conditions
        ]

    provenance = scene.provenance_row()
    measured = measure_bracket_echo_profiles(scene)
    truth_x = scene.target_coords[:, 1]

    rows: list[dict[str, float | str]] = []
    for condition in conditions:
        try:
            if condition == "structured":
                metric_rows = analyze_profiles(
                    measured, truth_x, policy, seed, include_flip_control=True
                )
                extra: dict[str, float | str] = {}
            elif condition == "column_shuffled":
                metric_rows = analyze_profiles(
                    shuffle_profile_columns(measured, seed + 11), None, policy, seed
                )
                extra = {}
            elif condition == "column_shuffled_repeat":
                metric_rows = analyze_profiles(
                    shuffle_profile_columns(measured, seed + 211), None, policy, seed
                )
                extra = {}
            elif condition == "relabel_symmetry":
                relabeled, truth_permuted = relabel_targets(
                    measured, truth_x, seed + 31
                )
                metric_rows = analyze_profiles(
                    relabeled, truth_permuted, policy, seed
                )
                extra = {}
            elif condition == "random_order":
                control_scene, density = geometry_free_scene(scene, seed + 41)
                control_profiles = measure_bracket_echo_profiles(control_scene)
                metric_rows = analyze_profiles(
                    control_profiles, None, policy, seed
                )
                extra = {"achieved_density": density}
            elif condition == "random_order_shuffled":
                control_scene, density = geometry_free_scene(scene, seed + 41)
                control_profiles = measure_bracket_echo_profiles(control_scene)
                metric_rows = analyze_profiles(
                    shuffle_profile_columns(control_profiles, seed + 61),
                    None,
                    policy,
                    seed,
                )
                extra = {"achieved_density": density}
            else:
                raise ValueError(f"unknown condition: {condition}")
        except (SceneValidityError, ValueError) as error:
            rows.append(
                {
                    **base,
                    **provenance,
                    "condition": condition,
                    "status": f"analysis_invalid: {error}",
                }
            )
            continue

        for metric_row in metric_rows:
            rows.append(
                {
                    **base,
                    **provenance,
                    **extra,
                    "condition": condition,
                    "status": "ok",
                    **metric_row,
                }
            )
    return rows


def heldout_at_gate_dim(
    rows: list[dict[str, float | str]],
    condition: str,
    seed: int,
    gate_dim: int,
) -> float | None:
    """Extract one condition's held-out violation at the gate dimension."""

    for row in rows:
        if (
            row.get("condition") == condition
            and row.get("status") == "ok"
            and float(row.get("seed", math.nan)) == float(seed)
            and int(float(row.get("embedding_dim", -1))) == gate_dim
        ):
            return float(row["heldout_violation"])
    return None


def structured_row_at_gate_dim(
    rows: list[dict[str, float | str]],
    condition: str,
    seed: int,
    gate_dim: int,
) -> dict[str, float] | None:
    """Extract one condition's full metric row at the gate dimension."""

    for row in rows:
        if (
            row.get("condition") == condition
            and row.get("status") == "ok"
            and float(row.get("seed", math.nan)) == float(seed)
            and int(float(row.get("embedding_dim", -1))) == gate_dim
        ):
            return {
                key: float(value)
                for key, value in row.items()
                if isinstance(value, (int, float))
            }
    return None


def write_rows_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    """Write rows with a stable union header."""

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def cohens_d(sample_a: list[float], sample_b: list[float]) -> float:
    """Cohen's d between two samples (a minus b)."""

    if len(sample_a) < 2 or len(sample_b) < 2:
        return float("nan")
    mean_a = sum(sample_a) / len(sample_a)
    mean_b = sum(sample_b) / len(sample_b)
    var_a = sum((v - mean_a) ** 2 for v in sample_a) / (len(sample_a) - 1)
    var_b = sum((v - mean_b) ** 2 for v in sample_b) / (len(sample_b) - 1)
    pooled = math.sqrt(
        ((len(sample_a) - 1) * var_a + (len(sample_b) - 1) * var_b)
        / (len(sample_a) + len(sample_b) - 2)
    )
    if pooled == 0.0:
        return float("inf") if mean_a != mean_b else 0.0
    return (mean_a - mean_b) / pooled
