from __future__ import annotations

import numpy as np

from causal_spacetime_lab.state_change_manifest_v3_protocol_diagnostics import (
    compute_v3_protocol_coverage_rows,
    compute_v3_protocol_failed_accounting_rows,
    compute_v3_protocol_latent_order_stability_rows,
    compute_v3_protocol_restart_stability_rows,
)


def test_compute_v3_protocol_coverage_rows_returns_finite_coverage(
    m41_manifest_dir,
) -> None:
    rows = compute_v3_protocol_coverage_rows(m41_manifest_dir)

    assert rows
    assert np.isfinite(float(rows[0]["target_coverage_fraction"]))
    assert np.isfinite(float(rows[0]["pair_node_coverage_fraction"]))


def test_compute_v3_protocol_restart_stability_rows_return_restart_std(
    m41_manifest_dir,
) -> None:
    rows = compute_v3_protocol_restart_stability_rows(
        m41_manifest_dir,
        embedding_dim=2,
        restart_count=2,
        steps=50,
        seed=0,
    )

    assert rows
    assert "restart_std" in rows[0]


def test_compute_v3_protocol_latent_order_stability_rows_return_disagreement(
    m41_manifest_dir,
) -> None:
    rows = compute_v3_protocol_latent_order_stability_rows(
        m41_manifest_dir,
        embedding_dim=2,
        restart_count=2,
        steps=50,
        seed=0,
    )

    assert rows
    assert "latent_order_disagreement" in rows[0]


def test_compute_v3_protocol_failed_accounting_rows_include_provenance(
    m41_manifest_dir,
) -> None:
    rows = compute_v3_protocol_failed_accounting_rows(m41_manifest_dir)

    assert rows
    assert rows[0]["handoff_provenance_type"]

