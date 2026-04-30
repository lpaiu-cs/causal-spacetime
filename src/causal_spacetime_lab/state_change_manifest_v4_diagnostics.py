"""Diagnostics for preregistered protocol-invariant v4 manifests."""

from __future__ import annotations

from pathlib import Path

from causal_spacetime_lab.state_change_manifest_v3_protocol_diagnostics import (
    compute_v3_protocol_coverage_rows,
    compute_v3_protocol_failed_accounting_rows,
    compute_v3_protocol_latent_order_stability_rows,
    compute_v3_protocol_manifest_family_fit_rows,
    compute_v3_protocol_null_taxonomy_rows,
    compute_v3_protocol_restart_stability_rows,
    compute_v3_protocol_stricter_criteria_rows,
    manifest_metadata_audit_rows,
    summarize_v3_protocol_fit_rows,
)


def compute_v4_protocol_manifest_family_fit_rows(
    manifest_dir: Path,
    *,
    dims: list[int],
    steps: int,
    restarts: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Fit protocol-bearing v4 manifests and return per-manifest rows."""

    return compute_v3_protocol_manifest_family_fit_rows(
        manifest_dir,
        dims=dims,
        steps=steps,
        restarts=restarts,
        seed=seed,
    )


def summarize_v4_protocol_fit_rows(
    fit_rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    """Summarize protocol-bearing v4 fit rows by family and dimension."""

    return summarize_v3_protocol_fit_rows(fit_rows)


def compute_v4_protocol_null_taxonomy_rows(
    manifest_dir: Path,
    *,
    embedding_dim: int,
    null_repetitions: int,
    steps: int,
    restarts: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute null taxonomy rows for protocol-bearing v4 manifests."""

    return compute_v3_protocol_null_taxonomy_rows(
        manifest_dir,
        embedding_dim=embedding_dim,
        null_repetitions=null_repetitions,
        steps=steps,
        restarts=restarts,
        seed=seed,
    )


def compute_v4_protocol_stricter_criteria_rows(
    manifest_dir: Path,
    *,
    dims: list[int],
    steps: int,
    restarts: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute stricter diagnostic rows without carry-forward decisions."""

    return compute_v3_protocol_stricter_criteria_rows(
        manifest_dir,
        dims=dims,
        steps=steps,
        restarts=restarts,
        seed=seed,
    )


def compute_v4_protocol_failed_accounting_rows(
    manifest_dir: Path,
) -> list[dict[str, float | str]]:
    """Report failed/report-only accounting for v4 manifests."""

    return compute_v3_protocol_failed_accounting_rows(manifest_dir)


def compute_v4_protocol_coverage_rows(
    manifest_dir: Path,
) -> list[dict[str, float | str]]:
    """Extract v4 target and pair-node coverage rows."""

    return compute_v3_protocol_coverage_rows(manifest_dir)


def compute_v4_protocol_restart_stability_rows(
    manifest_dir: Path,
    *,
    embedding_dim: int,
    restart_count: int,
    steps: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute v4 restart-stability rows."""

    return compute_v3_protocol_restart_stability_rows(
        manifest_dir,
        embedding_dim=embedding_dim,
        restart_count=restart_count,
        steps=steps,
        seed=seed,
    )


def compute_v4_protocol_latent_order_stability_rows(
    manifest_dir: Path,
    *,
    embedding_dim: int,
    restart_count: int,
    steps: int,
    seed: int,
) -> list[dict[str, float | str]]:
    """Compute v4 latent-order stability rows."""

    return compute_v3_protocol_latent_order_stability_rows(
        manifest_dir,
        embedding_dim=embedding_dim,
        restart_count=restart_count,
        steps=steps,
        seed=seed,
    )


def manifest_v4_metadata_audit_rows(
    manifest_dir: Path,
) -> list[dict[str, float | str]]:
    """Audit protocol/profile/provenance metadata on v4 manifests."""

    return manifest_metadata_audit_rows(manifest_dir)
