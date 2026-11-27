"""Independent small checks for the exploratory detector (not detection power)."""

import ast
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments" / "exploratory"))
from subspacing_observer import (  # noqa: E402
    FUTURE,
    PAST,
    BackgroundModel,
    Calibration,
    anchor_theta_jacobian,
    bump_window,
    calibrate,
    flat_gram_initialization,
    neyman_orthogonalize,
    observe,
    order_counts,
    select_anchors,
    spectrum,
    tangent_project,
    triangulate,
)
from subspacing_wave_sim import OrderOracle, WaveGeometry  # noqa: E402


def anchors():
    early = np.array(
        [
            [-0.25, 0.25, 0.12, 0],
            [-0.18, 0.22, -0.06, 0.12],
            [-0.29, 0.34, -0.10, -0.07],
            [-0.12, 0.31, 0.08, -0.09],
            [-0.32, 0.20, -0.02, 0.06],
        ]
    )
    late = np.array(
        [
            [0.24, 0.76, 0.10, 0.03],
            [0.17, 0.79, -0.12, 0.06],
            [0.31, 0.68, 0.04, -0.13],
            [0.12, 0.72, -0.07, -0.05],
            [0.28, 0.84, 0.03, 0.11],
        ]
    )
    return early, late


def test_blindness():
    source = (ROOT / "experiments/exploratory/subspacing_observer.py").read_text()
    tree = ast.parse(source)
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(a.name for a in node.names)
        if isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
    assert all(
        "wave_sim" not in name and "positive_control" not in name for name in imports
    )
    assert "importlib" not in source and "__import__" not in source
    assert "a_X" not in [
        node.id for node in ast.walk(tree) if isinstance(node, ast.Name)
    ]


@pytest.mark.parametrize("q", [0, 0.25, 1, 4])
def test_background_population_and_jacobian(q):
    early, late = anchors()
    points = np.array([[-0.05, 0.44, 0.03, -0.04], [0.06, 0.55, -0.06, 0.03]])
    model, geometry = BackgroundModel(q), WaveGeometry(q)
    actual, jac = model.moments(points, early, late, True)
    truth = geometry.population(points, early, late)
    assert np.max(abs(actual - truth)) < 1e-12
    for j in range(4):
        delta = np.zeros_like(points)
        delta[:, j] = 1e-6
        fd = (
            model.moments(points + delta, early, late)
            - model.moments(points - delta, early, late)
        ) / 2e-6
        assert np.max(abs(fd - jac[:, :, j])) < 1e-7
    calibration = Calibration(q, early, late, 0, 0, 0)
    recovered, residual, _ = triangulate(actual, calibration, iterations=10)
    assert np.max(abs(recovered - points)) < 1e-9
    assert np.max(abs(residual)) < 1e-12


def test_wave_fixture_agreement():
    path = str(ROOT / "experiments/positive_control")
    sys.path.insert(0, path)
    try:
        fixture_class = importlib.import_module(
            "order_wave_null_coordinate"
        ).WaveFixture
        for q in (0.05, 1, 4):
            geometry = WaveGeometry(q, "circular_wave", 1 / 128)
            fixture = fixture_class(q, 128, tighter=True)
            u = np.linspace(-0.49, 0.45, 33)
            assert (
                np.max(abs(geometry.matrix(-0.5, u) - fixture.matrix(-0.5, u))) < 1e-11
            )
            expected = fixture.sol.sol(u)[0]
            assert np.max(abs(geometry.primitive(u) - expected)) < 1e-12
            assert (
                abs(geometry.volume - np.pi * fixture.segment(-0.5, 0.5, 2)[0]) < 1e-11
            )
    finally:
        sys.path.remove(path)


@pytest.mark.parametrize("mode", ["background", "circular_wave"])
def test_physical_volume_sprinkling(mode):
    geometry = WaveGeometry(4, mode, 1 / 32)
    points = geometry.sprinkle(16000, np.random.default_rng(14))
    u, v, y = points[:, 0], points[:, 1], points[:, 2:]
    li, ri = (
        np.linalg.inv(geometry.matrix(-0.5, u)),
        np.linalg.inv(geometry.matrix(u, 0.5)),
    )
    lower = 0.5 * np.einsum("ni,nij,nj->n", y, li, y)
    upper = 1 - 0.5 * np.einsum("ni,nij,nj->n", y, ri, y)
    assert np.all((v > lower) & (v < upper))
    # Rosenblatt probability transforms must be independent uniforms.
    s = 0.5 * np.einsum("ni,nij,nj->n", y, li + ri, y)
    uniforms = np.column_stack(
        (
            geometry.cdf(u) / geometry.cdf.prefix[-1],
            2 * s - s * s,
            (v - lower) / (upper - lower),
        )
    )
    assert np.max(abs(uniforms.mean(axis=0) - 0.5)) < 0.009
    assert np.max(abs((uniforms**2).mean(axis=0) - 1 / 3)) < 0.009
    assert np.max(abs(np.corrcoef(uniforms.T) - np.eye(3))) < 0.025
    assert abs(np.corrcoef(u[:-1], u[1:])[0, 1]) < 0.025
    # Independent Cartesian midpoint cubature of a(u)^2*dv*dy integrates
    # the conditional radial measure: E[s]=1/3, not the ellipse-uniform 1/2.
    assert abs(s.mean() - 1 / 3) < 0.009
    if mode == "background":
        assert abs(geometry.volume - np.pi * BackgroundModel(4).gamma) < 1e-12


def test_flat_gram_and_joint_population_recovery():
    early, late = anchors()
    points = np.vstack((early, late))
    for q in (0, 1, 4):
        model = BackgroundModel(q)
        pf = np.column_stack(
            (model.interval(PAST, points), model.interval(points, FUTURE))
        )
        cross = model.interval(early[:, None], late[None])
        if q == 0:
            t, spatial = flat_gram_initialization(pf, cross)
            expected_t = (points[:, 0] + points[:, 1] + 0.5) / 2
            true_spatial = np.column_stack(
                ((points[:, 0] - points[:, 1] + 0.5) / 2, points[:, 2:] / np.sqrt(2))
            )
            assert np.max(abs(t - expected_t)) < 1e-12
            assert (
                np.max(abs(spatial @ spatial.T - true_spatial @ true_spatial.T)) < 1e-9
            )
        fit = calibrate(pf, cross, q=None, max_nfev=250)
        assert fit.rms < 1e-10
        assert abs(fit.q - q) < 2e-5
        assert fit.early[0, 2] > 0 and abs(fit.early[0, 3]) < 1e-14
        assert fit.early[1, 3] > 0


def test_nonzero_wave_signal_and_tangent_projection():
    early, late = anchors()
    q, lam = 1, 1 / 64
    geometry = WaveGeometry(q, "circular_wave", lam)
    rng = np.random.default_rng(28)
    points = np.column_stack(
        (
            rng.uniform(-0.08, 0.08, 192),
            rng.uniform(0.42, 0.58, 192),
            rng.uniform(-0.04, 0.04, (192, 2)),
        )
    )
    model = BackgroundModel(q)
    mask = np.all(model.moments(points, early, late) > 1e-5, axis=1)
    points = points[mask]
    mu = geometry.population(points, early, late, 128)
    calibration = Calibration(q, early, late, 0, 0, 0)
    recovered, residual, jac = triangulate(mu, calibration, initial=points)
    projected = tangent_project(residual, jac)
    assert np.max(abs(np.einsum("nki,nk->ni", jac, projected))) < 1e-14
    result = spectrum(recovered, residual, jac, calibration, np.arange(32, 129, 4))
    assert abs(result["score"][8]) / lam**2 > 1e-3
    orth, _ = neyman_orthogonalize(recovered, residual, jac, calibration, free_q=True)
    corrected = spectrum(recovered, orth, jac, calibration, np.arange(32, 129, 4))
    assert corrected["amplitude"][8] / lam**2 > 1e-3
    assert corrected["amplitude"][8] > 0.9 * result["amplitude"][8]
    assert corrected["peak_frequency"] == 64


@pytest.mark.parametrize("free_q", [False, True])
@pytest.mark.parametrize("q", [0, 1e-6, 1, 4])
def test_anchor_theta_differential(q, free_q):
    early, late = anchors()
    points = np.array([[-0.05, 0.44, 0.03, -0.04], [0.06, 0.55, -0.06, 0.03]])
    model = BackgroundModel(q)
    derivative = anchor_theta_jacobian(model, points, early, late, free_q)
    stacked = np.vstack((early, late))
    for column, index in enumerate(np.flatnonzero(np.arange(stacked.size) != 3)):
        delta = np.zeros(stacked.size)
        delta[index] = 1e-6
        plus, minus = stacked + delta.reshape(-1, 4), stacked - delta.reshape(-1, 4)
        fd = (
            model.moments(points, plus[: len(early)], plus[len(early) :])
            - model.moments(points, minus[: len(early)], minus[len(early) :])
        ) / 2e-6
        assert np.max(abs(fd - derivative[..., column])) < 1e-7
    if free_q:
        low, high = max(0, q - 1e-6), q + 1e-6
        fd = (
            BackgroundModel(high).moments(points, early, late)
            - BackgroundModel(low).moments(points, early, late)
        ) / (high - low)
        assert np.max(abs(fd - derivative[..., -1])) < 1e-7


@pytest.mark.parametrize("free_q", [False, True])
def test_neyman_zero_and_linear_anchor_removal(free_q):
    early, late = anchors()
    rng = np.random.default_rng(162)
    points = rng.uniform(
        [-0.08, 0.42, -0.04, -0.04], [0.08, 0.58, 0.04, 0.04], (256, 4)
    )
    model = BackgroundModel(4)
    moments, jac = model.moments(points, early, late, True)
    calibration = Calibration(4, early, late, 0, 0, 0)
    zero, norm = neyman_orthogonalize(
        points, np.zeros_like(moments), jac, calibration, free_q
    )
    assert np.max(abs(zero)) < 1e-15 and norm == 0
    nuisance = anchor_theta_jacobian(model, points, early, late, free_q)
    planted = nuisance @ rng.normal(0, 1e-5, nuisance.shape[-1])
    corrected, _ = neyman_orthogonalize(points, planted, jac, calibration, free_q)
    assert np.max(abs(corrected)) < 1e-11
    result = spectrum(points, zero, jac, calibration, np.arange(24, 225, 4))
    assert result["shot_sd"] == result["c_ref_estimate"] == 0
    assert result["whitened_snr"] == 0
    assert np.all(np.isfinite(result["whitened_amplitude"]))
    assert np.all(np.isfinite(result["whitened_sigma_bg"]))
    result = spectrum(points, planted, jac, calibration, np.arange(24, 225, 4))
    np.testing.assert_allclose(
        result["whitened_amplitude"], result["amplitude"] / result["whitened_sigma_bg"]
    )
    assert result["shot_sd"] > 0


def test_shear_differential_and_compact_window():
    early, late = anchors()
    points = np.array([[-0.03, 0.46, 0.04, -0.02], [0.04, 0.54, -0.03, 0.02]])
    model = BackgroundModel(4)
    symbol = model.shear_response(points, early, late)
    endpoints = [(PAST, 1), (FUTURE, -1)]
    endpoints += [(a, 1) for a in early] + [(b, -1) for b in late]
    for channel, (anchor, sign) in enumerate(endpoints):
        first, second = (anchor, points) if sign == 1 else (points, anchor)
        dm = model.m(second[..., 0]) - model.m(first[..., 0])
        dy = second[..., 2:] - first[..., 2:]
        for beta in (1, 1j):
            tensor = np.array(
                [
                    [complex(beta).real, complex(beta).imag],
                    [complex(beta).imag, -complex(beta).real],
                ]
            )
            values = []
            for epsilon in (-1e-6, 1e-6):
                matrix = (
                    dm[:, None, None] * np.eye(2)
                    - sign
                    * epsilon
                    * tensor[None]
                    / np.cos(model.k * points[:, 0])[:, None, None] ** 2
                )
                kap = (
                    second[..., 1]
                    - first[..., 1]
                    - 0.5 * np.einsum("ni,nij,nj->n", dy, np.linalg.inv(matrix), dy)
                )
                values.append(
                    model.alpha(second[..., 0] - first[..., 0])
                    / model.gamma
                    * np.maximum(kap, 0) ** 2
                )
            np.testing.assert_allclose(
                (values[1] - values[0]) / 2e-6,
                (symbol[:, channel] * beta).real,
                atol=1e-8,
            )
    rng = np.random.default_rng(91)
    box = rng.uniform([-0.5, 0, -0.3, -0.3], [0.5, 1, 0.3, 0.3], (10000, 4))
    active = bump_window(box) > 0
    assert active.any()
    for q in (0, 0.25, 1, 4):
        model = BackgroundModel(q)
        assert np.all(model.interval(PAST, box[active]) > 0)
        assert np.all(model.interval(box[active], FUTURE) > 0)


def test_blind_split_protocol():
    pytest.importorskip("cupy")
    geometry = WaveGeometry(4)
    points = geometry.sprinkle(131072, np.random.default_rng(53), "gpu")
    oracle = OrderOracle(geometry, points, "gpu")
    ids = np.arange(len(points))
    result = observe(
        oracle.query,
        ids[:8192],
        ids[106496:],
        ids[8192:106496],
        np.arange(32, 193, 8),
        q=None,
        count_query=oracle.counts,
        max_probes=64,
    )
    assert 0 <= result["q_hat"] < np.pi**2
    assert result["accepted"] >= 8
    assert result["n_ref"] == len(ids[106496:])
    assert result["tangent_leakage"] < 1e-12
    with pytest.raises(ValueError, match="disjoint"):
        observe(oracle.query, ids[:50], ids[:20], ids[100:200], [64], q=4)


def test_cpu_counts_and_strict_boundary():
    geometry = WaveGeometry(1)
    points = geometry.sprinkle(128, np.random.default_rng(2))
    oracle = OrderOracle(geometry, points)
    order, reverse = oracle.query(np.arange(128), np.arange(128))
    assert np.array_equal(order.T, reverse)
    assert not np.diag(order).any()
    refs, probes = np.arange(20, 100), np.arange(100, 128)
    early, late = np.arange(5), np.arange(5, 10)
    counts = order_counts(oracle.query, refs, probes, early, late, tile=17)
    for j, x in enumerate(probes):
        expected = [np.mean(order[refs, x]), np.mean(order[x, refs])]
        expected.extend(np.mean(order[a, refs] & order[refs, x]) for a in early)
        expected.extend(np.mean(order[x, refs] & order[refs, b]) for b in late)
        np.testing.assert_allclose(counts[j], expected, atol=0, rtol=0)
    special = np.array([[0, 0.2, 0, 0], [0, 0.3, 0, 0], [0, 0.3, 0.1, 0]])
    fw, _ = OrderOracle(geometry, special).query(np.arange(3), np.arange(3))
    assert fw[0, 1] and not fw[0, 2]
    with pytest.raises(ValueError):
        order_counts(oracle.query, refs, refs)


def test_gpu_cpu_parity():
    cp = pytest.importorskip("cupy")
    if cp.cuda.runtime.getDeviceCount() == 0:
        pytest.skip("no CUDA device")
    geometry = WaveGeometry(4, "circular_wave", 1 / 96)
    points = geometry.sprinkle(512, np.random.default_rng(5))
    cpu, gpu = OrderOracle(geometry, points), OrderOracle(geometry, points, "gpu")
    ids = np.arange(len(points))
    expected, actual = cpu.query(ids, ids), gpu.query(ids, ids)
    assert all(np.array_equal(a, b) for a, b in zip(expected, actual, strict=True))
    refs, probes = np.arange(20, 400), np.arange(400, 512)
    a, b = np.arange(5), np.arange(5, 10)
    np.testing.assert_array_equal(
        gpu.counts(refs, probes, a, b), order_counts(cpu.query, refs, probes, a, b)
    )
    # Cross both the 8192-reference and 4096-probe tile boundaries.
    tiled_points = geometry.sprinkle(12320, np.random.default_rng(51))
    tiled_cpu = OrderOracle(geometry, tiled_points)
    tiled_gpu = OrderOracle(geometry, tiled_points, "gpu")
    refs, probes = np.arange(10, 8203), np.arange(8203, 12320)
    np.testing.assert_array_equal(
        tiled_gpu.counts(refs, probes, a, b),
        order_counts(tiled_cpu.query, refs, probes, a, b, tile=1024),
    )
    u = np.linspace(-0.5, 0.5, 1000)
    assert (
        np.max(abs(cp.asnumpy(geometry.primitive.gpu(u)) - geometry.primitive(u)))
        < 1e-12
    )
    generated = geometry.sprinkle(512, np.random.default_rng(53), "gpu")
    rng = np.random.default_rng(53)
    draws = rng.random((512, 4))
    draws = draws[rng.permutation(512)]
    u, v, y = generated[:, 0], generated[:, 1], generated[:, 2:]
    li = np.linalg.inv(geometry.matrix(-0.5, u))
    ri = np.linalg.inv(geometry.matrix(u, 0.5))
    bottom = 0.5 * np.einsum("ni,nij,nj->n", y, li, y)
    top = 1 - 0.5 * np.einsum("ni,nij,nj->n", y, ri, y)
    s = 0.5 * np.einsum("ni,nij,nj->n", y, li + ri, y)
    np.testing.assert_allclose(
        geometry.cdf(u) / geometry.cdf.prefix[-1], draws[:, 0], atol=1e-12, rtol=0
    )
    np.testing.assert_allclose(s, 1 - np.sqrt(draws[:, 1]), atol=1e-12, rtol=0)
    np.testing.assert_allclose(
        (v - bottom) / (top - bottom), draws[:, 3], atol=1e-11, rtol=0
    )


@pytest.mark.parametrize("q", [-0.1, np.pi**2, np.nan])
def test_invalid_geometry(q):
    with pytest.raises(ValueError):
        WaveGeometry(q)
    with pytest.raises(ValueError):
        BackgroundModel(q)


@pytest.mark.parametrize("q", [0, 0.25, 1, 4])
def test_global_probe_initialization(q):
    early, late = anchors()
    model = BackgroundModel(q)
    points = WaveGeometry(q).sprinkle(65536, np.random.default_rng(909))
    mu = model.moments(points, early, late)
    active = np.all(mu > 1e-5, axis=1)
    assert active.sum() >= 50
    points, mu = points[active], mu[active]
    recovered, residual, _ = triangulate(
        mu, Calibration(q, early, late, 0, 0, 0), initial=None
    )
    assert np.max(np.abs(recovered - points)) < 1e-11
    assert np.max(np.abs(residual)) < 1e-12


def test_spectrum_parity_invariance():
    early, late = anchors()
    rng = np.random.default_rng(513)
    points = rng.uniform(
        [-0.08, 0.42, -0.04, -0.04], [0.08, 0.58, 0.04, 0.04], (128, 4)
    )
    model = BackgroundModel(4)
    _, jac = model.moments(points, early, late, True)
    residual = rng.normal(0, 1e-4, jac.shape[:2])
    frequencies = np.arange(24, 257, 4)
    original = spectrum(
        points, residual, jac, Calibration(4, early, late, 0, 0, 0), frequencies
    )
    parity = np.array([1, 1, 1, -1])
    reflected = spectrum(
        points * parity,
        residual,
        jac * parity,
        Calibration(4, early * parity, late * parity, 0, 0, 0),
        frequencies,
    )
    for a, b in (
        ("amplitude_pos", "amplitude_neg"),
        ("amplitude_neg", "amplitude_pos"),
        ("amplitude", "amplitude"),
        ("whitened_amplitude", "whitened_amplitude"),
    ):
        np.testing.assert_allclose(original[a], reflected[b], atol=1e-12, rtol=0)
    assert original["peak_signed_frequency"] == -reflected["peak_signed_frequency"]
    assert (
        original["whitened_peak_signed_frequency"]
        == -reflected["whitened_peak_signed_frequency"]
    )


def test_reflected_gram_calibration(monkeypatch):
    import subspacing_observer as observer

    original = observer.flat_gram_initialization
    reflection = np.diag([1, 1, -1])
    assert np.linalg.det(reflection) == -1

    def reflected(pf, cross):
        t, spatial = original(pf, cross)
        return t, spatial @ reflection

    monkeypatch.setattr(observer, "flat_gram_initialization", reflected)
    early, late = anchors()
    points = np.vstack((early, late))
    model = BackgroundModel(4)
    pf = np.column_stack((model.interval(PAST, points), model.interval(points, FUTURE)))
    cross = model.interval(early[:, None], late[None])
    fit = calibrate(pf, cross, q=None, max_nfev=250)
    assert fit.rms < 1e-10
    assert abs(fit.q - 4) < 2e-5
    recovered = np.vstack((fit.early, fit.late))
    # Compare the inferred geometry modulo its unobservable transverse O(2) gauge.
    np.testing.assert_allclose(recovered[:, :2], points[:, :2], atol=2e-5, rtol=0)
    np.testing.assert_allclose(
        recovered[:, 2:] @ recovered[:, 2:].T,
        points[:, 2:] @ points[:, 2:].T,
        atol=2e-5,
        rtol=0,
    )


def test_spectrum_zero_symbol(monkeypatch):
    early, late = anchors()
    points = np.tile([0.01, 0.5, 0.02, 0.01], (8, 1))
    _, jac = BackgroundModel(1).moments(points, early, late, True)
    rng = np.random.default_rng(13)
    shear = rng.normal(size=jac.shape[:2]) + 1j * rng.normal(size=jac.shape[:2])
    shear[0] = 0
    monkeypatch.setattr(BackgroundModel, "shear_response", lambda *args: shear)
    residual = rng.normal(size=jac.shape[:2])
    residual[0] = 1e20  # Zero sensitivity must contribute exactly zero, even here.
    calibration = Calibration(1, early, late, 0, 0, 0)
    frequencies = np.arange(24, 257, 4)
    result = spectrum(points, residual, jac, calibration, frequencies)
    residual[0] = 0
    expected = spectrum(points, residual, jac, calibration, frequencies)
    np.testing.assert_array_equal(result["amplitude"], expected["amplitude"])
    assert np.all(np.isfinite(result["amplitude"]))
    assert result["peak_amplitude"] < 10
    shear[:] = 0
    result = spectrum(points, residual, jac, calibration, frequencies)
    assert np.all(result["amplitude"] == 0)
    assert result["s_perp_norm_sq_sum"] == 0


def test_two_tier_anchor_selection():
    geometry = WaveGeometry(4)
    points = geometry.sprinkle(8192, np.random.default_rng(53))
    oracle = OrderOracle(geometry, points)
    pilot = np.arange(len(points))
    core_a, core_b = select_anchors(oracle.query, pilot)
    early, late = select_anchors(oracle.query, pilot, calibration_count=10)
    np.testing.assert_array_equal(early[:4], core_a)
    np.testing.assert_array_equal(late[:4], core_b)
    assert len(early) == len(late) == 10
    assert oracle.query(early, late)[0].all()


@pytest.mark.parametrize("q", [0, 1e-320, 1e-12, 0.25, 1, 4, 9.8])
def test_fast_alpha_matches_hypergeometric(q):
    from scipy.special import hyp0f1

    model = BackgroundModel(q)
    length = np.r_[0, np.geomspace(1e-10, 1, 1000)]
    z = -q * length**2 / 4
    a, b = hyp0f1(2.5, z), hyp0f1(1.5, z)
    expected = length**2 * a / (6 * b)
    derivative = (
        2 * length * a / b
        - length**3
        * q
        / 2
        * (hyp0f1(3.5, z) / (2.5 * b) - a * hyp0f1(2.5, z) / (1.5 * b**2))
    ) / 6
    np.testing.assert_allclose(model.alpha(length), expected, atol=1e-14, rtol=2e-12)
    np.testing.assert_allclose(
        model.alpha_prime(length), derivative, atol=1e-14, rtol=2e-12
    )


def test_calibration_reports_unweighted_rms():
    early, late = anchors()
    points = np.vstack((early, late))
    model = BackgroundModel(4)
    rng = np.random.default_rng(741)
    pf = np.column_stack((model.interval(PAST, points), model.interval(points, FUTURE)))
    cross = model.interval(early[:, None], late[None])
    pf += rng.normal(0, 1e-5, pf.shape)
    cross += rng.normal(0, 1e-5, cross.shape)
    fit = calibrate(pf, cross, q=4, starts=1)
    recovered = np.vstack((fit.early, fit.late))
    pred = np.column_stack(
        (model.interval(PAST, recovered), model.interval(recovered, FUTURE))
    )
    raw = np.r_[
        (pred - pf).ravel(),
        (model.interval(fit.early[:, None], fit.late[None]) - cross).ravel(),
    ]
    assert fit.rms > 1e-7
    assert fit.rms == pytest.approx(np.sqrt(np.mean(raw**2)), abs=1e-15)


def test_calibration_seeds_each_q_bucket(monkeypatch):
    import subspacing_observer as observer

    solve = observer.least_squares
    guesses = []

    def record(fun, initial, **kwargs):
        if len(initial) > 6:
            guesses.append(initial[-1])
        return solve(fun, initial, **kwargs)

    monkeypatch.setattr(observer, "least_squares", record)
    early, late = anchors()
    points = np.vstack((early, late))
    model = BackgroundModel(1)
    pf = np.column_stack((model.interval(PAST, points), model.interval(points, FUTURE)))
    calibrate(
        pf,
        model.interval(early[:, None], late[None]),
        q=None,
        orientations=4,
        starts=1,
        max_nfev=1,
    )
    np.testing.assert_array_equal(guesses, [0.05, 0.5, 1.5, 4.0, 7.0])


def test_campaign_collects_both_modes_and_preserves_history(monkeypatch):
    import json
    from concurrent.futures import ThreadPoolExecutor
    from tempfile import TemporaryDirectory

    import run_subspacing_pilot as runner

    monkeypatch.setattr(runner, "ProcessPoolExecutor", ThreadPoolExecutor)

    with TemporaryDirectory(prefix="subspacing-test-", dir=ROOT) as directory:
        assert Path(directory).resolve().is_relative_to(ROOT.resolve())
        output = Path(directory) / "campaign.json"
        monkeypatch.setattr(runner, "primitive_benchmark", lambda backend: {})
        monkeypatch.setattr(runner, "stage_one", lambda: [])
        monkeypatch.setattr(runner, "campaign_analysis", lambda summary: {})

        def dataset(n, q, mode, r, seed, backend, max_probes):
            return [
                dict(
                    N=n,
                    q_true=q,
                    mode=mode,
                    r=r,
                    seed=seed,
                    known_q=known,
                    status="completed",
                )
                for known in (True, False)
            ]

        monkeypatch.setattr(runner, "finite_dataset", dataset)
        monkeypatch.setattr(sys, "argv", ["pilot", "--output", str(output)])
        runner.main()
        completed = json.loads(output.read_text())
        assert len(completed["stage2"]) == 36
        assert len(completed["stage3"]) == 70
        assert len({r["seed"] for r in completed["stage2"] + completed["stage3"]}) == 53
        assert sum(r["N"] == 2**23 for r in completed["stage3"]) == 10
        monkeypatch.setattr(
            sys, "argv", ["pilot", "--output", str(output), "--budget-seconds", "0"]
        )
        runner.main()
        skipped = json.loads(output.read_text())
        assert skipped["previous_campaigns"] == [completed]
        assert len(skipped["stage2"]) == 18 and len(skipped["stage3"]) == 35
        assert all(
            r["status"] == "not_run_budget"
            for r in skipped["stage2"] + skipped["stage3"]
        )
