from __future__ import annotations

from causal_spacetime_lab.state_change_manifest_v4_diagnostics import (
    compute_v4_protocol_coverage_rows,
    compute_v4_protocol_latent_order_stability_rows,
    compute_v4_protocol_restart_stability_rows,
    manifest_v4_metadata_audit_rows,
)


def test_manifest_v4_metadata_audit_rows_detect_metadata(m44_manifest_dir) -> None:
    rows = manifest_v4_metadata_audit_rows(m44_manifest_dir)

    assert rows
    assert all(float(row["has_measurement_protocol_id"]) == 1.0 for row in rows)
    assert all(float(row["has_measurement_protocol_hash"]) == 1.0 for row in rows)
    assert all(float(row["has_profile_metadata"]) == 1.0 for row in rows)
    assert all(float(row["has_handoff_provenance"]) == 1.0 for row in rows)


def test_compute_v4_protocol_coverage_rows_returns_finite_rows(
    m44_manifest_dir,
) -> None:
    rows = compute_v4_protocol_coverage_rows(m44_manifest_dir)

    assert rows
    assert all(float(row["target_coverage_fraction"]) >= 0.0 for row in rows)
    assert all(float(row["pair_node_coverage_fraction"]) >= 0.0 for row in rows)


def test_v4_restart_and_latent_order_stability_rows(m44_manifest_dir) -> None:
    restart_rows = compute_v4_protocol_restart_stability_rows(
        m44_manifest_dir,
        embedding_dim=1,
        restart_count=2,
        steps=20,
        seed=0,
    )
    latent_rows = compute_v4_protocol_latent_order_stability_rows(
        m44_manifest_dir,
        embedding_dim=1,
        restart_count=2,
        steps=20,
        seed=0,
    )

    assert restart_rows
    assert latent_rows
    assert "restart_std" in restart_rows[0]
    assert "latent_order_disagreement" in latent_rows[0]
