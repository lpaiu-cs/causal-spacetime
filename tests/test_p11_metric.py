"""Regressions for the P11 metric instrument (prereg v1.6, Section 8).

The load-bearing pins: the chain estimator must equal the shared
definition's (L - 2)/sqrt(N) under the frozen rho = N/2 convention,
the verdict table must be single-valued by precedence, the power
literals must implement the frozen selection rule, and every seed
window must be private and fresh against the full documented list.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p11_metric import (  # noqa: E402
    DELTA_EQ,
    K_PAIRS,
    P11_LADDER,
    PILOT_BLOCKS,
    SCORE_OFFSET,
    SKIP_CAP,
    STAGE_A_BLOCKS,
    STRIDE,
    TAU_BAND,
    VERIFY_BASE,
    VERIFY_COUNT,
    continuum_uv,
    draw_disjoint_pairs,
    eligible_pool,
    power_requirements,
    run_sample,
    score_pair,
    supports_disjoint,
    verdict,
)


def test_continuum_uv_is_the_rank_grid():
    pi = np.array([2, 0, 3, 1])
    u, v = continuum_uv(pi)
    assert np.allclose(u, np.arange(4) / 4.0)
    assert np.allclose(v, pi / 4.0)


def test_tau_chain_equals_shared_definition_with_endpoints_subtracted():
    """Identity permutation: the interval (0, k) contains k - 1
    interior events forming a full chain, so the closed chain length is
    k + 1 and tau_hat must be exactly (k + 1 - 2)/sqrt(N)."""

    n, k = 100, 40
    pi = np.arange(n)
    u, v = continuum_uv(pi)
    scored = score_pair(u, v, 0, k)
    assert scored["chain_length"] == k + 1
    assert scored["tau_chain"] == pytest.approx((k - 1) / np.sqrt(n))
    # and the volume estimator through the shared metrics definition
    assert scored["m_open"] == k - 1
    assert scored["tau_vol"] == pytest.approx(2.0 * np.sqrt((k - 1) / n))
    # truth via the frozen convention
    assert scored["tau_true"] == pytest.approx(2.0 * k / n)


def test_conditioned_interior_mean_uses_rank_gaps():
    n, k = 100, 40
    pi = np.arange(n)
    u, v = continuum_uv(pi)
    scored = score_pair(u, v, 0, k)
    assert scored["m_conditioned"] == pytest.approx(
        (k - 1) * (k - 1) / (n - 2)
    )


def test_eligible_pool_matches_band_by_construction():
    rng = np.random.default_rng(7)
    pi = rng.permutation(300)
    u, v = continuum_uv(pi)
    pool = eligible_pool(u, v)
    assert pool.shape[0] > 0
    taus = 2.0 * np.sqrt(
        (u[pool[:, 1]] - u[pool[:, 0]]) * (v[pool[:, 1]] - v[pool[:, 0]])
    )
    assert np.all(taus >= TAU_BAND[0]) and np.all(taus <= TAU_BAND[1])
    assert np.all(u[pool[:, 1]] > u[pool[:, 0]])
    assert np.all(v[pool[:, 1]] > v[pool[:, 0]])


def test_supports_disjoint_is_box_disjointness():
    a = (0.0, 0.2, 0.0, 0.2)
    assert supports_disjoint(a, (0.3, 0.5, 0.0, 0.2))   # u-disjoint
    assert supports_disjoint(a, (0.0, 0.2, 0.3, 0.5))   # v-disjoint
    assert not supports_disjoint(a, (0.1, 0.3, 0.1, 0.3))  # overlap


def test_draw_rejection_cap_marks_incomplete():
    """A pool of mutually overlapping candidates can never yield six
    disjoint supports; the frozen cap must fire and mark the sample
    incomplete rather than loop."""

    u = np.linspace(0.0, 1.0, 40)
    v = np.linspace(0.0, 1.0, 40)
    pool = np.array([[0, 39 - k] for k in range(10)])  # nested boxes
    pairs, complete, rejections = draw_disjoint_pairs(
        pool, u, v, np.random.default_rng(0)
    )
    assert not complete
    assert len(pairs) == 1
    assert rejections == 9  # pool exhausted before the cap


def test_run_sample_completes_and_scores_at_bottom_rung():
    record, complete = run_sample(600, 190000, single_stream=True)
    assert complete
    assert record["pool_size"] > 0
    assert np.isfinite(record["y"])
    assert record["min_m_open"] > 0


def test_verification_mode_smoke_completeness():
    """20 single-stream samples at the bottom rung: the Section 4 pin
    demands 1998/2000; this smoke demands 20/20 at test speed (a
    failure here means the full pin is hopeless)."""

    complete_count = sum(
        run_sample(600, 190000 + k, single_stream=True)[1]
        for k in range(20)
    )
    assert complete_count == 20


def test_spacelike_pool_and_dual_box_normalization():
    """Reversed permutation: every pair is spacelike. For the extreme
    pair (0, N-1) the dual box is the whole open square and its
    unanchored LIS is the estimator with NO endpoint correction:
    d_hat must equal L/sqrt(N) exactly."""

    from p11_metric import (
        eligible_pool_spacelike,
        run_sample_spacelike,
        score_pair_spacelike,
    )

    n = 60
    pi = np.arange(n)[::-1]
    u, v = continuum_uv(pi)
    pool = eligible_pool_spacelike(u, v)
    dists = 2.0 * np.sqrt(
        (u[pool[:, 1]] - u[pool[:, 0]]) * (v[pool[:, 0]] - v[pool[:, 1]])
    )
    assert pool.shape[0] > 0
    assert np.all(dists >= TAU_BAND[0]) and np.all(dists <= TAU_BAND[1])
    assert np.all(v[pool[:, 0]] > v[pool[:, 1]])   # spacelike ordering

    scored = score_pair_spacelike(u, v, 0, n - 1)
    assert scored["tau_chain"] == pytest.approx(
        scored["chain_length"] / np.sqrt(n)
    )
    assert scored["tau_true"] == pytest.approx(
        2.0 * np.sqrt((u[n - 1] - u[0]) * (v[0] - v[n - 1]))
    )

    record, complete = run_sample_spacelike(600, 600000, single_stream=True)
    assert complete and np.isfinite(record["y"])


def test_dual_box_is_reproduced_by_pure_order_forcing():
    """Section 10 amendment: box membership must be order data. Start
    from the single oriented incomparable pair x ⊏ y and propagate the
    spatial orientation by Gallai's Gamma-forcing, using ONLY the
    causal relation (soundness derivable case by case in 1+1D):

        a ⊏ b  and  a || c  and  b comparable to c   =>   a ⊏ c
        a ⊏ b  and  c || b  and  a comparable to c   =>   c ⊏ b

    The forced set {z : x ⊏ z ⊏ y, z || x, z || y} must equal the
    rank-window box the implementation uses. No coordinates enter the
    forcing — only the causal matrix."""

    from p11_metric import eligible_pool_spacelike, score_pair_spacelike

    rng = np.random.default_rng(41)
    checked = 0
    for _trial in range(6):
        n = 200   # dense enough for the forcing witnesses; at n = 80
        #           one of six trials fails to percolate (soundness
        #           always holds) -- the Section 10 caveat, measured
        pi = rng.permutation(n)
        u, v = continuum_uv(pi)
        pool = eligible_pool_spacelike(u, v)
        if pool.shape[0] == 0:
            continue
        i, j = (int(x) for x in pool[rng.integers(0, pool.shape[0])])
        # causal matrix from ranks -- the order relation itself
        causal = (u[:, None] < u[None, :]) & (v[:, None] < v[None, :])
        incomp = ~causal & ~causal.T & ~np.eye(n, dtype=bool)
        comp = causal | causal.T

        # Vectorized fixed point over the two Gamma rules:
        #   rule 1: (O @ comp) masked to incomparable pairs
        #   rule 2: (comp @ O) masked to incomparable pairs
        oriented = np.zeros((n, n), dtype=bool)
        oriented[i, j] = True
        while True:
            grown = (oriented
                     | ((oriented @ comp) & incomp)
                     | ((comp @ oriented) & incomp))
            if (grown == oriented).all():
                break
            oriented = grown
        forced_box = {
            z for z in range(n)
            if oriented[i, z] and oriented[z, j]
            and incomp[i, z] and incomp[j, z]
        }
        rank_box = {
            z for z in range(n)
            if u[i] < u[z] < u[j] and v[j] < v[z] < v[i]
        }
        assert forced_box <= rank_box      # soundness, unconditionally
        assert forced_box == rank_box      # completeness at this density
        scored = score_pair_spacelike(u, v, i, j)
        assert scored["m_open"] == len(rank_box)
        checked += 1
    assert checked >= 5


def test_dual_box_certificate_certifies_scored_samples():
    """The production certificate: every scored pair's box membership
    must be order-forced via the two-step Gamma witnesses. At the
    experiment's densities the certificate should hold; a certified
    member is a box member by soundness, so certifying all members
    verifies exactly what the scorer consumes."""

    from p11_metric import (
        dual_box_order_certificate,
        eligible_pool_spacelike,
        run_sample_spacelike,
    )

    rng = np.random.default_rng(97)
    certified = 0
    for _trial in range(4):
        n = 300
        pi = rng.permutation(n)
        u, v = continuum_uv(pi)
        pool = eligible_pool_spacelike(u, v)
        i, j = (int(x) for x in pool[rng.integers(0, pool.shape[0])])
        certified += int(dual_box_order_certificate(u, v, i, j))
    assert certified == 4

    record, complete = run_sample_spacelike(600, 600100,
                                            single_stream=True)
    assert complete and record["box_order_certified"]


def test_certificate_tiers_pinned_on_the_two_production_cases():
    """The two Stage B samples the anchored tier flagged, pinned:
    (600, 425200) has a pair the anchored closure misses but the full
    Gamma fixed point certifies (the truncation artifact — anchors
    accumulate across generations); (1200, 442200) has one pair of
    six that BOTH tiers refuse — a genuine forcing failure, the
    Section 10 caveat in production, carried as the sample flag."""

    from p11_metric import (
        SCORE_OFFSET,
        draw_disjoint_pairs,
        dual_box_order_certificate,
        eligible_pool_spacelike,
        full_gamma_certificate,
        run_sample_spacelike,
    )

    def pairs_for(n, seed):
        rng = np.random.default_rng(seed)
        pi = rng.permutation(n)
        u, v = continuum_uv(pi)
        pool = eligible_pool_spacelike(u, v)
        draw_rng = np.random.default_rng(seed + SCORE_OFFSET)
        pairs, complete, _ = draw_disjoint_pairs(pool, u, v, draw_rng)
        assert complete
        return u, v, pairs

    u, v, pairs = pairs_for(600, 425200)
    anchored = [dual_box_order_certificate(u, v, i, j) for i, j in pairs]
    assert not all(anchored)
    assert all(
        a or full_gamma_certificate(u, v, i, j)
        for a, (i, j) in zip(anchored, pairs, strict=True)
    )
    record, complete = run_sample_spacelike(600, 425200)
    assert complete and record["box_order_certified"]

    u, v, pairs = pairs_for(1200, 442200)
    tiered = [
        dual_box_order_certificate(u, v, i, j)
        or full_gamma_certificate(u, v, i, j)
        for i, j in pairs
    ]
    assert tiered.count(False) == 1
    record, complete = run_sample_spacelike(1200, 442200)
    assert complete and not record["box_order_certified"]


def test_stage_b_windows_are_private_and_fresh():
    from p11_metric import (
        PILOT_B_BLOCKS,
        STAGE_B_BLOCKS,
        VERIFY_B_BASE,
    )

    assert PILOT_B_BLOCKS == {600: (336000, 220), 2400: (380000, 220)}
    assert STAGE_B_BLOCKS == {600: (424000, 80), 1200: (440000, 80),
                              2400: (456000, 80)}
    spans = []
    for base, slots in {**PILOT_BLOCKS, **STAGE_A_BLOCKS,
                        **PILOT_B_BLOCKS, **STAGE_B_BLOCKS}.values():
        spans.append(set(range(base, base + STRIDE * slots)))
    for x in spans:
        for w in spans:
            assert x is w or not (x & w)
    experimental_top = max(max(s) for s in spans)
    for n in P11_LADDER:
        vspan = set(range(VERIFY_B_BASE[n], VERIFY_B_BASE[n] + VERIFY_COUNT))
        assert min(vspan) > experimental_top
        for s in spans:
            assert not (vspan & s)


def test_verdict_precedence_is_single_valued():
    assert verdict(-0.05, -0.01, True) == "IMPROVES"      # also in margin
    assert verdict(0.01, 0.05, True) == "DEGRADES"        # also in margin
    assert verdict(-0.05, 0.05, True) == "FLAT-WITHIN-MARGIN"
    assert verdict(-0.10, 0.05, True) == "INCONCLUSIVE"
    assert verdict(-0.05, 0.05, False) == "UNRESOLVED"
    assert verdict(-0.05, -0.01, False) == "IMPROVES"
    assert abs(DELTA_EQ - 0.067) < 1e-12


def _gaussian_pilot(sigma: float, seed: int, n: int = 200) -> np.ndarray:
    y = np.random.default_rng(seed).standard_normal(n)
    y = (y - y.mean()) / y.std(ddof=1)          # exact sample sigma
    return sigma * y


def test_bonett_bound_is_above_the_variance_and_kurtosis_adaptive():
    from p11_metric import bonett_variance_bound

    gaussian = bonett_variance_bound(_gaussian_pilot(0.10, 1))
    assert gaussian["s2"] == pytest.approx(0.01)
    # Gaussian-kurtosis inflation sits near 1.18, never below 1
    assert 1.05 < gaussian["bound"] / gaussian["s2"] < 1.35

    # heavy tails must widen the bound: same variance, Laplace-ish y
    rng = np.random.default_rng(2)
    heavy = rng.laplace(size=200)
    heavy = 0.10 * (heavy - heavy.mean()) / heavy.std(ddof=1)
    heavy_bound = bonett_variance_bound(heavy)
    assert heavy_bound["g4"] > gaussian["g4"]
    assert (heavy_bound["bound"] / heavy_bound["s2"]
            > gaussian["bound"] / gaussian["s2"])


def test_power_selection_rule_branches():
    """v1.8: sizing from summed per-rung Bonett bounds; branch
    structure exercised with exact-sample-sigma Gaussian pilots."""

    # small sigmas -> n_eq drives, flat available
    p = power_requirements(_gaussian_pilot(0.07, 3), _gaussian_pilot(0.07, 4))
    assert p["flat_available"] and p["n_per_rung"] == p["n_eq"]
    assert p["n_eq"] <= 60 and not p["infeasible"]
    # moderate sigmas -> flat unaffordable, superiority sizing
    p = power_requirements(_gaussian_pilot(0.15, 5), _gaussian_pilot(0.15, 6))
    assert not p["flat_available"] and not p["infeasible"]
    assert p["n_per_rung"] == max(p["n_sup"], 12)
    # large sigmas -> superiority breaches the cap -> infeasible
    p = power_requirements(_gaussian_pilot(0.40, 7), _gaussian_pilot(0.40, 8))
    assert p["infeasible"] and p["n_per_rung"] is None
    # variances enter as sums; the bound exceeds the point estimate
    p = power_requirements(_gaussian_pilot(0.05, 9), _gaussian_pilot(0.13, 10))
    assert p["s2"] == pytest.approx(0.05 ** 2 + 0.13 ** 2)
    assert p["s2_90"] > p["s2"]


def test_stage_gate_rejects_stale_prerequisite(tmp_path, monkeypatch):
    """A prerequisite artifact stamped by an older implementation (or
    left behind by a crashed rerun) must be refused, not consumed."""

    import p11_metric as mod

    (tmp_path / "p11_verification_summary.json").write_text(
        '{"code_version": "0ld5tamp", "pin_passed": true}\n',
        encoding="utf-8",
    )
    monkeypatch.setattr(mod, "_diag_code_version", lambda: "abc1234")
    with pytest.raises(SystemExit, match="implementation changed"):
        mod.run_pilot(tmp_path)


def test_preflight_refuses_a_dirty_worktree(monkeypatch, tmp_path):
    import p11_metric as mod

    monkeypatch.setattr(mod, "_diag_code_version",
                        lambda: "abc1234-dirty")
    with pytest.raises(SystemExit, match="dirty"):
        mod.run_verify(tmp_path)


def test_stage_pass_gate_requires_improves_and_reachable_stamp(
        monkeypatch, tmp_path):
    """Section 6: pilot-B may not consume its quarantined windows
    unless Stage A's FROZEN record reads IMPROVES with a stamp
    reachable from HEAD."""

    import subprocess

    import p11_metric as mod

    monkeypatch.setattr(mod, "FROZEN_P11_DIR", tmp_path)
    # missing record refuses
    with pytest.raises(SystemExit, match="has not run"):
        mod._require_stage_pass("p11_stage_a_summary.json", "Stage A")
    # failing verdict refuses
    (tmp_path / "p11_stage_a_summary.json").write_text(
        '{"verdict": "DEGRADES", "code_version": "whatever"}\n',
        encoding="utf-8",
    )
    with pytest.raises(SystemExit, match="not\\s+IMPROVES"):
        mod._require_stage_pass("p11_stage_a_summary.json", "Stage A")
    # a stamp that is not in the object store refuses too -- but for its
    # OWN reason. The message must not say "not an ancestor", because a
    # missing object says nothing about ancestry: that conflation is what
    # three findings on PRs #35 and #36 were built on, and what this
    # repository's own CI run then reproduced against a shallow clone.
    (tmp_path / "p11_stage_a_summary.json").write_text(
        '{"verdict": "IMPROVES", "code_version": "0000000"}\n',
        encoding="utf-8",
    )
    with pytest.raises(SystemExit) as absent:
        mod._require_stage_pass("p11_stage_a_summary.json", "Stage A")
    assert "not present in this checkout" in str(absent.value)
    assert "NOT evidence that the history was flattened" in str(absent.value)
    assert "not an ancestor" not in str(absent.value)
    # and a stamp that IS present but off the ancestry gets the other
    # message, so both branches are covered rather than just the one the
    # environment happens to produce
    monkeypatch.setattr(mod, "stamp_reachability",
                        lambda _stamp: "not-ancestor")
    with pytest.raises(SystemExit, match="present but is not"):
        mod._require_stage_pass("p11_stage_a_summary.json", "Stage A")
    monkeypatch.undo()
    monkeypatch.setattr(mod, "FROZEN_P11_DIR", tmp_path)
    # passing verdict with a reachable stamp is accepted
    head = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    (tmp_path / "p11_stage_a_summary.json").write_text(
        f'{{"verdict": "IMPROVES", "code_version": "{head}"}}\n',
        encoding="utf-8",
    )
    result = mod._require_stage_pass("p11_stage_a_summary.json", "Stage A")
    assert result["verdict"] == "IMPROVES"


def test_preflight_refuses_an_unknown_stamp(monkeypatch, tmp_path):
    """git-less environments stamp 'unknown', which carries no
    provenance and would even satisfy the prerequisite equality check
    against another unknown — refused like dirty."""

    import p11_metric as mod

    monkeypatch.setattr(mod, "_diag_code_version", lambda: "unknown")
    with pytest.raises(SystemExit, match="unknown"):
        mod.run_verify(tmp_path)


def test_calibrated_bound_records_coverage_and_never_shrinks():
    from p11_metric import (
        BONETT_Z,
        bonett_variance_bound,
        calibrated_variance_bound,
    )

    gaussian = calibrated_variance_bound(_gaussian_pilot(0.10, 11), "t1")
    assert 0.0 <= gaussian["coverage_at_nominal"] <= 1.0
    assert gaussian["z_used"] >= BONETT_Z
    assert (gaussian["bound"]
            >= bonett_variance_bound(_gaussian_pilot(0.10, 11))["bound"]
            - 1e-12)
    # calibration flag consistent with z
    assert gaussian["calibrated"] == (gaussian["z_used"] > BONETT_Z)

    # a heavy-tailed pilot: if nominal coverage misses the target, the
    # calibrated z must exceed nominal and coverage at z_used meets it
    rng = np.random.default_rng(12)
    heavy = rng.standard_t(df=3, size=200)
    result = calibrated_variance_bound(0.1 * heavy, "t2")
    if result["coverage_at_nominal"] < 0.95:
        assert result["z_used"] > BONETT_Z


def test_windows_are_private_fresh_and_reserve_sized():
    """Stride-200 privacy, the documented spans, and freshness against
    every seed range the programme has ever used."""

    used = (set(range(0, 10)) | set(range(100, 120)) | set(range(400, 420))
            | set(range(500, 520)) | set(range(820, 960))
            | set(range(1000, 1020)) | set(range(1100, 1160))
            | set(range(9001, 9060)) | set(range(30000, 30380))
            | set(range(40000, 40169)) | set(range(41000, 41060))
            | set(range(43000, 55000)))
    blocks = []
    for n, (base, slots) in {**PILOT_BLOCKS, **STAGE_A_BLOCKS}.items():
        span = set(range(base, base + STRIDE * slots))
        blocks.append((n, base, span))
        assert not (span & used), (n, base)
        assert max(0, 3, 5, 9, 61, 100, SCORE_OFFSET) < STRIDE
    for _, _, a in blocks:
        for _, _, b in blocks:
            assert a is b or not (a & b)
    assert PILOT_BLOCKS[600] == (200000, 220)
    assert PILOT_BLOCKS[2400] == (244000, 220)
    assert STAGE_A_BLOCKS == {600: (288000, 80), 1200: (304000, 80),
                              2400: (320000, 80)}
    experimental_bottom = min(
        base for base, _ in {**PILOT_BLOCKS, **STAGE_A_BLOCKS}.values()
    )
    for n in P11_LADDER:
        verify_span = set(range(VERIFY_BASE[n], VERIFY_BASE[n] + VERIFY_COUNT))
        assert max(verify_span) < experimental_bottom
        assert not (verify_span & used)
    assert SKIP_CAP == 20 and K_PAIRS == 6


def test_stage_c_fiducials_are_off_diagonal_and_frozen():
    """Section 11 v2: every target fiducial must sit off u = v, where
    the two-post discriminant vanishes identically -- the defect that
    killed v1's unenumerated endpoint grid."""

    from p11_metric import ANCHOR_FIDUCIALS, TARGET_FIDUCIALS

    assert ANCHOR_FIDUCIALS == ((0.15, 0.15), (0.85, 0.85), (0.15, 0.85))
    assert TARGET_FIDUCIALS == ((0.38, 0.44), (0.38, 0.56), (0.50, 0.41),
                                (0.50, 0.59), (0.62, 0.44), (0.62, 0.56))
    for fu, fv in TARGET_FIDUCIALS:
        assert abs(fu - fv) >= 0.05           # off the degenerate axis
        assert 0.35 <= fu <= 0.65 and 0.35 <= fv <= 0.65
        # strictly inside the anchor cone, so eligibility is generic
        assert 0.15 < fu < 0.85 and 0.15 < fv < 0.85


def test_two_post_roots_clamp_the_negative_discriminant():
    """The frozen D := max(D, 0) path: a measurement pair that admits
    no real intersection must still yield the tangency root, never an
    exception and never a completeness rejection."""

    from p11_metric import two_post_roots

    # exact-geometry case: a target at (0.5, 0.45) between the posts
    u1, v1, u2, v2 = 0.15, 0.15, 0.85, 0.85
    A = (0.5 - u1) * (0.45 - v1)
    B = (u2 - 0.5) * (v2 - 0.45)
    roots = two_post_roots(A, B, u1, v1, u2, v2)
    assert roots
    assert min(abs(r[0] - 0.5) + abs(r[1] - 0.45) for r in roots) < 1e-9

    # inflate both measurements so no real intersection exists
    bad = two_post_roots(A * 4.0, B * 4.0, u1, v1, u2, v2)
    assert bad and all(np.isfinite(r).all() for r in bad)


def test_gauss_newton_sign_convention_and_exactness():
    """The sign multiplies the PRODUCT, not each difference (the bug
    that inverted the third constraint). With exact measurements the
    fit must return the true point."""

    from p11_metric import gauss_newton_fit

    true_u, true_v = 0.5, 0.42
    a1, a2, a3 = (0.15, 0.15), (0.85, 0.85), (0.15, 0.85)
    A = (true_u - a1[0]) * (true_v - a1[1])
    B = (a2[0] - true_u) * (a2[1] - true_v)
    # a3 is spacelike to the target: the product is negative
    prod3 = (true_u - a3[0]) * (true_v - a3[1])
    assert prod3 < 0
    C = -prod3
    fit = gauss_newton_fit((0.45, 0.5), (A, B, C),
                           ((a1[0], a1[1], 1.0), (a2[0], a2[1], 1.0),
                            (a3[0], a3[1], -1.0)))
    assert abs(fit[0] - true_u) < 1e-9 and abs(fit[1] - true_v) < 1e-9


def test_stage_c_sample_is_metric_only_not_rank_readout():
    """The anti-circularity property, checked on data: a rank-readout
    estimator would have error identically zero. Stage C's error must
    be strictly positive and finite, and its record must carry the
    order certificate."""

    from p11_metric import run_sample_coordinates

    record, complete = run_sample_coordinates(600, 700000)
    assert complete
    assert 0.0 < record["median_coord_error"] < 1.0
    assert np.isfinite(record["y"])
    assert "box_order_certified" in record


def test_stage_c_windows_are_private_and_fresh():
    from p11_metric import (
        PILOT_C_BLOCKS,
        STAGE_C_BLOCKS,
        VERIFY_C_BASE,
    )

    earlier = set()
    for base, slots in (list(PILOT_BLOCKS.values())
                        + list(STAGE_A_BLOCKS.values())):
        earlier |= set(range(base, base + STRIDE * slots))
    c_spans = []
    for base, slots in (list(PILOT_C_BLOCKS.values())
                        + list(STAGE_C_BLOCKS.values())):
        span = set(range(base, base + STRIDE * slots))
        assert not (span & earlier)
        c_spans.append(span)
    for a in c_spans:
        for b in c_spans:
            assert a is b or not (a & b)
    # verification-C is quarantined above every experimental window
    top = max(base + STRIDE * slots
              for base, slots in (list(PILOT_C_BLOCKS.values())
                                  + list(STAGE_C_BLOCKS.values())))
    for n in P11_LADDER:
        vspan = set(range(VERIFY_C_BASE[n], VERIFY_C_BASE[n] + 2000))
        assert min(vspan) >= top


def test_published_calibrated_bound_matches_the_power_input():
    """Deferred P2 (1): the artifact must publish the bound power
    actually consumed. The per-rung calibrated bounds must sum to the
    power step's s2_90 exactly -- the mismatch that recurred in all
    three pilots."""

    from p11_metric import calibrated_variance_bound, power_requirements

    y_b = _gaussian_pilot(0.12, 21)
    y_t = _gaussian_pilot(0.14, 22)
    published = (calibrated_variance_bound(y_b, "bottom")["bound"]
                 + calibrated_variance_bound(y_t, "top")["bound"])
    assert power_requirements(y_b, y_t)["s2_90"] == pytest.approx(published)
