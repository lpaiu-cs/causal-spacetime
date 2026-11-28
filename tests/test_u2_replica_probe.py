"""Finite geometry/algebra and interface checks, not asymptotic evidence."""

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import linalg

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments/exploratory"))
from subspacing_observer import local_pair_correlator, observe  # noqa: E402
from subspacing_wave_sim import OrderOracle, PanelPrimitive, WaveGeometry  # noqa: E402
from u2_spectral_and_replica_probe import (  # noqa: E402
    aggregate_replicas,
    complex_list,
    conformal_basis,
    constraint_jacobian,
    diffeomorphism_basis,
    orthogonal_compression,
    overlap,
    pair_margin_gradient,
    spectral_cell,
    transported_scores,
)


def test_seven_conformal_modes_fix_tips_and_exact_null_pairs():
    rng = np.random.default_rng(17)
    left = rng.uniform(-0.2, 0.2, (40, 4))
    right = left.copy()
    right[:, 0] += rng.uniform(0.1, 0.5, 40)
    right[:, 2:] += rng.uniform(-0.2, 0.2, (40, 2))
    right[:, 1] += np.sum((right[:, 2:] - left[:, 2:]) ** 2, axis=1) / (
        2 * (right[:, 0] - left[:, 0])
    )
    margin, gradient = pair_margin_gradient(left, right)
    fields = conformal_basis(np.concatenate((left, right))).reshape(80, 4, 7)
    variation = np.einsum("ni,nik->nk", gradient, fields[40:] - fields[:40])
    assert np.max(abs(margin)) < 1e-15
    assert np.max(abs(variation)) < 1e-12
    assert np.linalg.matrix_rank(fields.reshape(-1, 7)) == 7
    assert (
        np.max(abs(conformal_basis(np.array([[-0.5, 0, 0, 0], [0.5, 1, 0, 0]]))))
        < 1e-15
    )


def test_finite_margin_conformal_leakage_is_not_zero():
    left = np.array([[-0.2, 0.2, 0.02, -0.03]])
    null = np.array([[0.2, 0.2 + 0.1**2 / 0.8, 0.12, -0.03]])
    leakage = []
    for eta in (1e-2, 1e-3, 1e-4):
        right = null.copy()
        right[:, 1] += eta
        _, gradient = pair_margin_gradient(left, right)
        fields = conformal_basis(np.concatenate((left, right))).reshape(2, 4, 7)
        leakage.append(np.linalg.norm(gradient @ (fields[1] - fields[0])))
    assert leakage[0] > 1e-6
    assert 8 < leakage[0] / leakage[1] < 12
    assert 8 < leakage[1] / leakage[2] < 12


def test_margin_gradient_by_finite_difference():
    left, right = np.array([-0.3, 0.1, 0.03, 0.02]), np.array([0.2, 0.8, -0.1, 0.05])
    _, analytic = pair_margin_gradient(left, right)
    numerical = []
    for axis in np.eye(4):
        numerical.append(
            (
                pair_margin_gradient(left, right + 1e-6 * axis)[0]
                - pair_margin_gradient(left, right - 1e-6 * axis)[0]
            )
            / 2e-6
        )
    np.testing.assert_allclose(analytic, numerical, atol=1e-10)


def test_two_dimensional_null_faces_and_diffeomorphisms():
    left = np.array([[-0.2, 0.2], [-0.3, 0.4]])
    right = np.array([[-0.2, 0.8], [0.3, 0.4]])
    fields = diffeomorphism_basis(np.concatenate((left, right)), 8).reshape(4, 2, 16)
    np.testing.assert_allclose(fields[2, 0] - fields[0, 0], 0, atol=1e-14)
    np.testing.assert_allclose(fields[3, 1] - fields[1, 1], 0, atol=1e-14)
    tips = np.array([[-0.5, 0], [0.5, 1]])
    assert np.max(abs(diffeomorphism_basis(tips))) < 1e-14
    points = np.array([[-0.2, 0.2], [-0.199, 0.8], [0.2, 0.201]])
    jac, counts = constraint_jacobian(points, 0.01, tips=False)
    h = (jac.T @ jac).toarray()
    assert h[0, 0] > 90  # near the u face
    assert h[1, 1] > 90  # near the v face
    assert counts["near_pair_count"] == 2
    reverse_face, reverse_counts = constraint_jacobian(
        np.array([[-0.2, 0.8], [-0.199, 0.2]]), 0.01, tips=False
    )
    assert reverse_counts["near_pair_count"] == 1
    assert (reverse_face.T @ reverse_face)[0, 0] > 90
    # Smooth diffeomorphisms are only approximately null at finite bandwidth.
    lattice = np.column_stack(
        (np.linspace(-0.45, 0.45, 30), np.linspace(0.05, 0.95, 30))
    )
    modes = diffeomorphism_basis(lattice)
    wide, _ = constraint_jacobian(lattice, 0.08)
    narrow, _ = constraint_jacobian(lattice, 0.008)
    assert np.linalg.norm(narrow @ modes) < np.linalg.norm(wide @ modes) / 10


@pytest.mark.parametrize("dim", [2, 4])
def test_psd_symmetry_and_strict_projection(dim):
    rng = np.random.default_rng(24)
    points = (
        WaveGeometry(0).sprinkle(16, rng)
        if dim == 4
        else rng.random((16, 2)) - [0.5, 0]
    )
    jac, _ = constraint_jacobian(points, 0.05)
    h = (jac.T @ jac).toarray()
    np.testing.assert_allclose(h, h.T, atol=1e-12)
    assert linalg.eigvalsh(h)[0] > -1e-10
    basis = conformal_basis(points)
    q = linalg.null_space(basis.T)
    compressed = orthogonal_compression(h, basis)
    assert compressed.shape == (len(h) - basis.shape[1],) * 2
    np.testing.assert_allclose(
        linalg.eigvalsh(compressed), linalg.eigvalsh(q.T @ h @ q), atol=1e-10
    )
    # Omitting tips leaves exact translations; fixed-tip rows lift these.
    free, _ = constraint_jacobian(points, 0.05, tips=False)
    translation = np.tile(np.eye(dim), (len(points), 1))
    np.testing.assert_allclose(free @ translation, 0, atol=1e-14)
    assert np.linalg.norm(jac @ translation) > 0


def test_spectral_cpu_and_gpu_agree():
    cp = pytest.importorskip("cupy")
    if cp.cuda.runtime.getDeviceCount() == 0:
        pytest.skip("CUDA device required")
    cpu = spectral_cell(24, 4, 0.05, 3, "cpu")
    gpu = spectral_cell(24, 4, 0.05, 3, "gpu")
    np.testing.assert_allclose(cpu["lambda_1_to_16"], gpu["lambda_1_to_16"], atol=1e-9)
    assert cpu["lambda_min_perp"] == pytest.approx(gpu["lambda_min_perp"], abs=1e-9)


def test_overlap_and_pair_uncertainty():
    score = np.array([1 + 2j, 3 - 1j])
    assert overlap(complex_list(score), complex_list(score)) == pytest.approx(7.5)
    assert abs(overlap(complex_list(score), complex_list(1j * score))) < 1e-14
    rows = []
    for value in (1.0, 3.0, 5.0, 7.0):
        rows.append(
            {
                "status": "completed",
                "N": 1024,
                "r": 32,
                "metrics": {"blind_energy_excess": value},
                "arms": {
                    arm: dict(
                        phase_cos_mean=0.0,
                        phase_coherence_abs=0.2,
                        u_error_rms_over_lambda=2.0,
                    )
                    for arm in ("wave1", "wave2")
                },
            }
        )
    rows.append({"status": "failed", "N": 1024, "r": 32})
    result = aggregate_replicas(rows)[0]
    assert result["completed_pairs"] == 4
    assert result["alpha"] == 2
    stats = result["metrics"]["blind_energy_excess"]
    assert stats["mean"] == 4
    assert stats["se_of_mean"] == pytest.approx(np.std([1, 3, 5, 7], ddof=1) / 2)


def test_diagnostic_output_uses_only_blind_inputs():
    """Diagnostic opt-in preserves scores and aligns returned probe IDs/moments."""
    cp = pytest.importorskip("cupy")
    if cp.cuda.runtime.getDeviceCount() == 0:
        pytest.skip("CUDA device required")
    geometry = WaveGeometry(4)
    points = geometry.sprinkle(131072, np.random.default_rng(53), "gpu")
    oracle = OrderOracle(geometry, points, "gpu")
    ids = np.arange(len(points))
    arguments = (oracle.query, ids[:8192], ids[106496:], ids[8192:106496], [64])
    plain = observe(*arguments, q=4, count_query=oracle.counts, max_probes=64)
    baseline = observe(
        *arguments, q=4, count_query=oracle.counts, max_probes=64, orthogonalize=False
    )
    assert baseline["orthogonalize"] is False and baseline["dtheta_norm"] == 0
    assert plain["orthogonalize"] is True
    detailed = observe(
        *arguments, q=4, count_query=oracle.counts, max_probes=64, diagnostics=True
    )
    assert "diagnostics" not in plain
    np.testing.assert_array_equal(plain["score"], detailed["score"])
    diagnostic = detailed["diagnostics"]
    np.testing.assert_array_equal(
        diagnostic["moments"],
        oracle.counts(
            ids[106496:],
            diagnostic["probe_ids"],
            detailed["early_ids"],
            detailed["late_ids"],
        ),
    )
    assert diagnostic["points"].shape == (detailed["accepted"], 4)
    assert not np.shares_memory(diagnostic["points"], points)


def test_truth_transport_identity_and_phase_rotation():
    from subspacing_observer import BackgroundModel, Calibration, spectrum

    early = np.array(
        [
            [-0.25, 0.25, 0.12, 0],
            [-0.18, 0.22, -0.06, 0.12],
            [-0.29, 0.34, -0.1, -0.07],
            [-0.12, 0.31, 0.08, -0.09],
        ]
    )
    late = early.copy()
    late[:, :2] = -early[:, :2] + [0, 1]
    calibration = Calibration(4, early, late, 0, 0, 0)
    rng = np.random.default_rng(19)
    points = np.column_stack(
        (
            rng.uniform(-0.05, 0.05, 32),
            rng.uniform(0.45, 0.55, 32),
            rng.uniform(-0.02, 0.02, (32, 2)),
        )
    )
    model = BackgroundModel(4)
    moments, jac = model.moments(points, early, late, True)
    residual = rng.normal(0, 1e-5, moments.shape)
    score = transported_scores(
        points, moments + residual, calibration, 64, points[:, 0]
    )
    reference = spectrum(points, residual, jac, calibration, [64])
    np.testing.assert_allclose(
        score, [reference["score"][0], reference["score_neg"][0]], atol=1e-16
    )
    shifted = transported_scores(
        points, moments + residual, calibration, 64, points[:, 0] + np.pi / 128
    )
    np.testing.assert_allclose(shifted, score * [-1j, 1j], atol=1e-16)


@pytest.mark.parametrize("eta", [0, -1, np.nan, np.inf])
def test_invalid_eta(eta):
    with pytest.raises(ValueError, match="eta"):
        constraint_jacobian(np.array([[0.0, 0.5], [0.1, 0.6]]), eta)


@pytest.mark.parametrize("dc", [0.1, np.inf])
def test_local_pair_correlator_matches_explicit_pairs(dc):
    rng = np.random.default_rng(590)
    u = rng.uniform(-0.25, 0.25, 27)
    u[1] = u[0]  # Different events in one bin must remain in the sum.
    z = rng.normal(size=len(u)) + 1j * rng.normal(size=len(u))
    frequencies = np.arange(-128, 129, 16)
    delta = u[:, None] - u[None, :]
    kernel = np.exp(-0.5 * (delta / dc) ** 2)
    np.fill_diagonal(kernel, 0)
    expected = np.array(
        [
            np.sum(
                kernel * (z[:, None] * z[None, :].conj() * np.exp(-1j * f * delta)).real
            )
            for f in frequencies
        ]
    )
    actual = local_pair_correlator(u, z, frequencies, dc)
    assert np.max(abs(actual - expected)) < (1e-12 if np.isinf(dc) else 1e-4)
    if np.isinf(dc):
        np.testing.assert_allclose(
            actual,
            abs(z @ np.exp(-1j * np.outer(u, frequencies))) ** 2 - np.sum(abs(z) ** 2),
            atol=1e-12,
            rtol=0,
        )
    single = local_pair_correlator(u[:1], z[:1], frequencies, dc)
    np.testing.assert_allclose(single, 0, atol=1e-12)


def test_ensemble_intervals_and_threshold():
    from run_whitened_ensemble import aggregate, error_statistics, separation, wilson

    np.testing.assert_allclose(wilson(1, 8), [0.02241749, 0.47088818], atol=1e-8)
    stats = separation([2, 3, 4, 5], [0, 1, 2, 3])
    assert stats["difference"] == 2
    assert stats["ci95"][0] < 2 < stats["ci95"][1]
    assert separation([1], [2])["status"] == "insufficient_samples"
    assert error_statistics([1, 2, 3])["median"] == 2
    rows = []
    for mode in ("background", "circular_wave"):
        for i in range(8):
            rows.append(
                dict(
                    part=1,
                    N=2**22,
                    r=64,
                    mode=mode,
                    status="completed",
                    peak_hit=False,
                    q_err=i / 100,
                    whitened_target_snr=1,
                    whitened_peak_z=float(i + (mode == "circular_wave")),
                )
            )
    result = aggregate(rows)["part1"][0]["empirical_threshold"]
    assert result["whitened_peak_z"] == 6
    assert result["false_positives"] == 1
    assert result["wave_detections"] == 2


def test_batched_panel_conversion_matches_individual_legendre_conversion():
    from numpy.polynomial import legendre, polynomial

    def function(u):
        return np.exp(5j * u) / np.cos(2 * u) ** 2

    primitive = PanelPrimitive(function, 37)
    nodes, _ = legendre.leggauss(8)
    starts = np.arange(37) / 37 - 0.5
    values = function(starts[:, None] + (nodes + 1) / 74)
    integral = legendre.legint(legendre.legfit(nodes, values.T, 7), axis=0) / 74
    reference = np.stack([legendre.leg2poly(integral[:, j]) for j in range(37)])
    reference[:, 0] -= polynomial.polyval(-1, reference.T)
    np.testing.assert_allclose(primitive.coeff, reference, rtol=0, atol=1e-15)
