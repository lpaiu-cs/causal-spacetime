"""Regressions for P12 Stage B (prereg Section 10.11).

Section 10.11 lists eight tests and says each is "written against the
case it is meant to catch", so every test below names its case. Five of
them are ported regressions for defects that already bit once — P13's
C5 (both arms gated) and C2 (two-sided quantile), and the addendum
review's C11/C21 (a gate that recorded its failure and continued), C17
(seed blocks checked by eye), C10 (a statistic that dropped the twin's
variance) and C25 (a stale expected value inside the staleness
regression).

Nothing here runs a campaign. The expensive paths are exercised on a
handful of seeds only where a cheap structural assertion cannot reach
the behaviour.
"""

from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (
    Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
)
sys.path.insert(0, str(EXPERIMENT_DIR))

from p12_stage_b import (  # noqa: E402
    CROSS_ARM_TOLERANCE,
    DELTA_EQ_B,
    DELTA_STAR_B,
    ETA_LO,
    FROZEN_P13V3_DIR,
    M_TARGET,
    M_TOLERANCE,
    N_CAP,
    N_EQ_COEFF,
    N_FLOOR,
    P12B_LADDER,
    R_TAU2_TRUE,
    R_TRUE,
    RECOVERY_THRESHOLD,
    REL_SD_A,
    REL_SD_B,
    SPENT_RANGES,
    TAU_ELL,
    WINDOW_CEIL,
    WINDOW_FLOOR,
    XHALF,
    _abort_on_gate_failure,
    assert_windows_disjoint_and_fresh,
    cross_arm_m_gate,
    g_hat,
    g_of_s,
    invert_g,
    m_gate,
    paired_recovery,
    power_requirements_b,
    seed_blocks,
    stage_b_record,
    volume_closed_form,
)

# --------------------------------------------------------------------
# helpers: synthetic per-sample records, so the gate and aggregation
# tests do not need a campaign
# --------------------------------------------------------------------


def _records(n, mean_g, mean_m, mean_tau=1.0, spread=0.0, seed=0):
    rng = np.random.default_rng(seed)
    jitter = rng.normal(0.0, spread, n) if spread else np.zeros(n)
    return [
        {"mean_g": float(mean_g + jitter[i]), "mean_m": float(mean_m),
         "mean_tau_hat": float(mean_tau)}
        for i in range(n)
    ]


# ====================================================================
# 1. the closed-form volume against quadrature, and position
#    independence at fixed tau
# ====================================================================

def _volume_quadrature(eta_p, x_p, eta_q, x_q, ell=1.0, n=128):
    """Independent integration of Omega^2 over the diamond, using
    neither the closed form nor tau."""

    u_p, v_p = eta_p + x_p, eta_p - x_p
    u_q, v_q = eta_q + x_q, eta_q - x_q
    nodes, weights = np.polynomial.legendre.leggauss(n)
    u = 0.5 * (u_q - u_p) * nodes + 0.5 * (u_q + u_p)
    v = 0.5 * (v_q - v_p) * nodes + 0.5 * (v_q + v_p)
    wu = 0.5 * (u_q - u_p) * weights
    wv = 0.5 * (v_q - v_p) * weights
    eta = 0.5 * (u[:, None] + v[None, :])
    return float(np.einsum("i,j,ij->", wu, wv, ell ** 2 / eta ** 2) * 0.5)


def _tau_curved(eta_p, x_p, eta_q, x_q, ell=1.0):
    z = ((eta_p ** 2 + eta_q ** 2 - (x_p - x_q) ** 2)
         / (2.0 * eta_p * eta_q))
    return ell * math.acosh(z)


@pytest.mark.parametrize(
    ("eta_p", "x_p", "deta", "dx"),
    [(-1.9, 0.0, 0.6, 0.1), (-1.5, 0.7, 0.4, 0.05),
     (-1.2, -0.4, 0.15, 0.02), (-6.0, 2.0, 3.0, 0.5),
     (-3.0, -1.0, 1.2, 0.3)],
)
def test_volume_closed_form_matches_quadrature(eta_p, x_p, deta, dx):
    """The case: Section 5 refused the small-diamond EXPANSION because
    its coefficient is "the kind of unverified constant that killed
    Stage C v1". If the closed form were wrong, every recovery would be
    biased by a factor nothing else in the pipeline could reveal.

    Review C12 is why the tolerance is asserted at 1e-9 rather than
    reported: the design check's first version measured 4.1e-9 against
    a pinned 1e-9, recorded it, and continued.
    """

    eta_q, x_q = eta_p + deta, x_p + dx
    tau = _tau_curved(eta_p, x_p, eta_q, x_q)
    closed = volume_closed_form(tau)
    quad = _volume_quadrature(eta_p, x_p, eta_q, x_q)
    assert abs(closed - quad) / closed < 1e-9


def test_volume_is_position_independent_at_fixed_tau():
    """The case: maximal symmetry REQUIRES the diamond volume to depend
    on tau alone. A position-dependent volume would mean the derivation
    dropped a term, and the ladder would read different curvature at
    different depths in the patch.
    """

    target, volumes = 0.45, []
    for eta_p in (-1.9, -1.6, -1.3):
        lo, hi = 1e-6, abs(eta_p) - 1e-6
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if _tau_curved(eta_p, 0.0, eta_p + mid, 0.0) < target:
                lo = mid
            else:
                hi = mid
        deta = 0.5 * (lo + hi)
        volumes.append(_volume_quadrature(eta_p, 0.0, eta_p + deta, 0.0))
    spread = (max(volumes) - min(volumes)) / float(np.mean(volumes))
    assert spread < 1e-9
    assert abs(np.mean(volumes) - volume_closed_form(target)) < 1e-9


# ====================================================================
# 2. G strictly decreasing, and G^-1(G(s)) = s to 1e-10
# ====================================================================

def test_g_is_strictly_decreasing_and_starts_at_one_half():
    """The case: the inversion is single-valued only because G is
    monotone. If it were not, `invert_g`'s bisection would return one of
    several branches, silently and seed-dependently.
    """

    ss = np.unique(np.concatenate([np.linspace(1e-6, 1.0, 400),
                                   np.linspace(1.0, 6.0, 400)]))
    gs = np.array([g_of_s(float(s)) for s in ss])
    assert np.all(np.diff(gs) < 0.0)
    assert g_of_s(1e-12) == pytest.approx(0.5, abs=1e-12)
    assert g_of_s(0.0) == pytest.approx(0.5, abs=1e-12)


@pytest.mark.parametrize("s", [0.05, 0.15, 0.375, 0.75, 1.0, 2.0, 4.0])
def test_g_inverse_round_trips(s):
    """The case: Section 10.11 item 2 pins the round trip at 1e-10.
    The bisection has a fixed iteration count, so a widened bracket or
    a changed floor could quietly cost precision the recovery needs --
    at the operating point the amplification is 12.9x, so inversion
    slop is multiplied.
    """

    recovered = invert_g(g_of_s(s))
    assert recovered is not None
    assert abs(recovered - s) < 1e-10


def test_operating_point_inverts_to_the_true_curvature():
    """The case: the whole chain g -> s -> R must be exact on exact
    inputs, or no amount of statistics helps. tau/ell = 1.5 gives
    R tau^2 = 4.5 and R = 2.
    """

    g_true = volume_closed_form(TAU_ELL) / TAU_ELL ** 2
    s = invert_g(g_true)
    assert s == pytest.approx(TAU_ELL / 2.0, abs=1e-10)
    assert 8.0 * s ** 2 == pytest.approx(R_TAU2_TRUE, rel=1e-10)
    assert 8.0 * s ** 2 / TAU_ELL ** 2 == pytest.approx(R_TRUE, rel=1e-10)


# ====================================================================
# 3. the boundary convention
# ====================================================================

@pytest.mark.parametrize("g", [0.5, 0.5 + 1e-12, 0.6, 1.0, 2.0])
def test_boundary_returns_none_at_or_above_the_flat_value(g):
    """The case: noise pushes `g` above the flat value 1/2, where no
    positive curvature is representable. Section 10.2 freezes this as
    the EDGE OF THE PARAMETER SPACE -- R_hat = 0, relative error
    exactly 1 -- and not as a clamp, so it must be reported rather than
    smoothed.
    """

    assert invert_g(g) is None


def test_boundary_samples_score_exactly_one_and_are_counted():
    """The case: review C13 found "no single sample recovers R" false
    against the artifact, because the medians EXCLUDED undefined
    inversions while the power block included them at error 1. So the
    boundary population must appear in the count AND in the error array
    with the value 1 exactly.
    """

    # twin g below curved g => Q > 1 => Q/2 > 1/2 => boundary
    curved = _records(4, mean_g=0.90, mean_m=19.0)
    twin = _records(4, mean_g=0.80, mean_m=19.0)
    out = paired_recovery(curved, twin)
    assert out["boundary_hits"] == 4
    assert out["boundary_fraction"] == pytest.approx(1.0)
    assert np.all(out["rel_error"] == 1.0)
    assert out["count_below_one"] == 0

    # A sample INSIDE the domain is not counted as boundary, and error
    # exactly 1 is then reserved for the boundary: this fixture sits in
    # the domain but far from the truth (its tau_hat is wrong), so the
    # error is large WITHOUT being the boundary's 1 -- which is the
    # distinction C13 collapsed.
    inside = paired_recovery(curved, _records(4, mean_g=2.0, mean_m=19.0))
    assert inside["boundary_hits"] == 0
    assert np.all(np.isfinite(inside["rel_error"]))
    assert not np.any(inside["rel_error"] == 1.0)

    # and a genuinely recovering sample lands below 1 and is counted
    g_true = g_of_s(TAU_ELL / 2.0)
    good = paired_recovery(
        _records(4, mean_g=g_true, mean_m=76.0, mean_tau=TAU_ELL),
        _records(4, mean_g=0.5, mean_m=76.0, mean_tau=TAU_ELL))
    assert good["boundary_hits"] == 0
    assert good["count_below_one"] == 4
    assert good["fraction_below_one"] == pytest.approx(1.0)


def test_recovery_is_exact_when_the_ratio_is_the_truth():
    """The case: the aggregation must not introduce its own bias. Fed a
    ratio equal to the continuum truth and the true tau, the paired
    statistic has to return relative error zero.
    """

    g_true = g_of_s(TAU_ELL / 2.0)
    curved = _records(3, mean_g=g_true, mean_m=76.0, mean_tau=TAU_ELL)
    twin = _records(3, mean_g=0.5, mean_m=76.0, mean_tau=TAU_ELL)
    out = paired_recovery(curved, twin)
    assert np.all(out["rel_error"] < 1e-9)


# ====================================================================
# 4. the m-gate covering BOTH arms (P13 review C5, ported)
# ====================================================================

def test_m_gate_fails_on_a_drifting_twin_beside_a_clean_curved_arm():
    """The case, verbatim from P13's review C5: Section 5 gives the twin
    "the same eligibility conditions, K, rejection cap, fill rule and
    m-gate", and P13's first implementation gated only the curved rows.
    A twin rung could then drift in discreteness and still have its
    contrast used as a control. In Stage B the twin is INSIDE the
    estimator as the normalizer, so an ungated twin is a silent bias
    rather than a decorative one.
    """

    curved = {r: _records(5, 1.0, M_TARGET[r]) for r in P12B_LADDER}
    twin = {r: _records(5, 1.0, M_TARGET[r]) for r in P12B_LADDER}
    assert m_gate(curved, "curved")["passed"]
    assert m_gate(twin, "twin")["passed"]

    # drift ONE twin rung well past the tolerance, leave curved clean
    drifted = dict(twin)
    drifted[1200] = _records(5, 1.0, M_TARGET[1200] * 1.12)
    curved_gate = m_gate(curved, "curved")
    twin_gate = m_gate(drifted, "twin")
    assert curved_gate["passed"], "the curved arm is clean by construction"
    assert not twin_gate["passed"], "the twin drift must be caught"
    assert not twin_gate["per_rung"]["1200"]["passed"]
    # and it aborts rather than being recorded
    with pytest.raises(SystemExit, match="ABORTS on its frozen m gates"):
        _abort_on_gate_failure([curved_gate, twin_gate])


def test_m_gate_targets_are_per_rung_because_the_ladder_moves_m():
    """The case: P13's gate form does NOT transfer literally, and a
    thoughtless port would gate against a grand mean. Stage B sweeps
    density, so m runs 19 -> 76 BY DESIGN and a grand-mean test would
    fail on the design itself. This pins the per-rung form so a future
    "harmonization" with P13 cannot reintroduce the wrong gate.
    """

    exact = {r: _records(3, 1.0, M_TARGET[r]) for r in P12B_LADDER}
    assert m_gate(exact, "curved")["passed"]
    grand = float(np.mean([M_TARGET[r] for r in P12B_LADDER]))
    # every rung AT the grand mean would pass a P13-style gate and must
    # fail this one, since none of them is at its own target
    at_grand = {r: _records(3, 1.0, grand) for r in P12B_LADDER}
    gate = m_gate(at_grand, "curved")
    assert not gate["passed"]
    assert not gate["per_rung"][str(P12B_LADDER[0])]["passed"]
    assert not gate["per_rung"][str(P12B_LADDER[-1])]["passed"]


# ====================================================================
# 4b. the CROSS-ARM m gate at +/- 1%, and it must ABORT
# ====================================================================

def test_cross_arm_gate_catches_arms_offset_from_one_another():
    """The case, from review C11: two arms each internally flat but
    offset from each other. Q_hat divides one by the other, so both
    within-arm gates pass while numerator and denominator carry
    different discreteness biases. C11 measured the cost: a 2-3.4%
    offset gave a ~1.5% bias mismatch, and with Q/(1-Q) ~ 11 that is
    ~17% injected into the recovery.
    """

    curved = {r: _records(4, 1.0, M_TARGET[r]) for r in P12B_LADDER}
    matched = {r: _records(4, 1.0, M_TARGET[r]) for r in P12B_LADDER}
    assert cross_arm_m_gate(curved, matched)["passed"]

    offset = {r: _records(4, 1.0, M_TARGET[r] * 0.97)
              for r in P12B_LADDER}
    # each arm is internally perfect against its own targets...
    assert m_gate(curved, "curved")["passed"]
    inner = m_gate(offset, "twin")
    assert inner["passed"], (
        "a 3% offset is inside the 5% within-arm tolerance, which is "
        "exactly why the cross-arm gate has to exist")
    # ...and the cross-arm gate is what sees it
    gate = cross_arm_m_gate(curved, offset)
    assert not gate["passed"]
    for rung in P12B_LADDER:
        assert gate["per_rung"][str(rung)]["cross_arm_offset"] == (
            pytest.approx(-0.03, abs=1e-9))


def test_cross_arm_gate_failure_aborts_and_writes_nothing():
    """The case, from review C21: the first fix for C11 RECORDED
    `gate_passed: false` and continued, so nothing stopped the artifact
    from being written and copied into frozen/. That is review C12's
    shape reappearing inside the repair for C11. Section 10.11 item 4b
    requires asserting the abort, not the flag.
    """

    curved = {r: _records(4, 1.0, M_TARGET[r]) for r in P12B_LADDER}
    offset = {r: _records(4, 1.0, M_TARGET[r] * 1.05)
              for r in P12B_LADDER}
    gate = cross_arm_m_gate(curved, offset)
    assert not gate["passed"]
    with pytest.raises(SystemExit) as excinfo:
        _abort_on_gate_failure([gate])
    assert "nothing written" in str(excinfo.value)


def test_cross_arm_tolerance_is_tighter_than_the_within_arm_one():
    """The case: if the two tolerances were ever equalized, the
    cross-arm gate would stop being able to see the offsets the
    within-arm gate lets through, which is the only reason it exists.
    """

    assert CROSS_ARM_TOLERANCE < M_TOLERANCE


# ====================================================================
# 4c. the paired statistic must carry the twin's variation
# ====================================================================

def test_paired_statistic_has_larger_variance_than_a_fixed_denominator():
    """The case, from review C10: dividing every curved sample by the
    twin's RUNG MEAN gives a statistic whose stored variance contains
    none of the twin arm's sampling variation, while the power
    calculation treats the values as independent -- so the projected n
    comes out too optimistic. Section 10.11 item 4c says to assert the
    inequality rather than trust the docstring, so this constructs the
    substituted version on the SAME data and compares.
    """

    # TWO fixture properties are load-bearing, and both were got wrong
    # on the way to this test, so they are written down rather than left
    # to look arbitrary.
    #
    # (a) No sample may reach the boundary. Boundary hits are pinned at
    #     error exactly 1, which compresses dispersion, so a fixture
    #     straddling the edge would let the comparison turn on the clamp.
    #
    # (b) The fixture must sit AWAY from exact recovery. `y` is
    #     log10 of the error, so at zero error it has a logarithmic
    #     singularity, and there `log|Q - Q_true| = log(sigma) +
    #     log|Z|` -- the variance of `y` becomes INDEPENDENT of sigma
    #     and the comparison measures distribution shape instead of the
    #     twin's contribution. Centring on the truth therefore makes
    #     this test pass or fail for the wrong reason. The campaign's
    #     actual errors are 26% / 26% / 15%, so a 2% bias mismatch --
    #     which is also what Section 10.6 measures -- puts the fixture
    #     where the statistic really lives.
    g_true = g_of_s(TAU_ELL / 2.0)
    curved = _records(400, g_true * 0.98, 76.0, mean_tau=TAU_ELL,
                      spread=0.002, seed=1)
    twin = _records(400, 0.5, 76.0, mean_tau=TAU_ELL,
                    spread=0.002, seed=2)

    paired = paired_recovery(curved, twin)
    assert paired["boundary_hits"] == 0, "fixture must stay off the edge"
    assert 0.1 < paired["median_rel_error"] < 0.9, (
        "fixture must sit away from exact recovery -- see (b)")

    # the C10 construction: one fixed estimated denominator
    twin_mean = float(np.mean([r["mean_g"] for r in twin]))
    substituted = paired_recovery(
        curved, [{"mean_g": twin_mean, "mean_m": 76.0,
                  "mean_tau_hat": TAU_ELL} for _ in curved])

    var_paired = float(np.var(paired["y"], ddof=1))
    var_fixed = float(np.var(substituted["y"], ddof=1))
    assert var_fixed < var_paired, (
        f"substituting the twin's rung mean understates the variance: "
        f"{var_fixed:.6f} vs paired {var_paired:.6f}")
    # and by roughly the predicted amount: for a ratio of two arms with
    # comparable relative dispersion the variances add in quadrature, so
    # dropping one arm should remove about half. A ratio near 1 would
    # mean the twin's variation was not really being carried.
    assert 0.3 < var_fixed / var_paired < 0.8, (
        f"expected the fixed denominator to lose about half the "
        f"variance, got {var_fixed / var_paired:.3f}")


def test_paired_statistic_uses_each_sample_s_own_tau_hat():
    """The case: if the aggregation reached for a rung-level tau_hat,
    the statistic would smuggle a shared random quantity back in by a
    different door. Two samples differing ONLY in tau_hat must give
    different y.
    """

    g_true = g_of_s(TAU_ELL / 2.0)
    curved = [
        {"mean_g": g_true, "mean_m": 76.0, "mean_tau_hat": TAU_ELL},
        {"mean_g": g_true, "mean_m": 76.0, "mean_tau_hat": 1.2 * TAU_ELL},
    ]
    twin = _records(2, 0.5, 76.0, mean_tau=TAU_ELL)
    out = paired_recovery(curved, twin)
    assert out["rel_error"][0] != out["rel_error"][1]
    assert out["rel_error"][0] < 1e-9


# ====================================================================
# 4d. seed-block disjointness, enumerated
# ====================================================================

def test_seed_blocks_are_pairwise_disjoint_and_enumerated():
    """The case, from review C17: the design check's previous layout put
    rung 600's calibration iteration 1 on exactly rung 1200's iteration
    0, and 14 cross-rung collisions followed behind a comment that
    claimed separation. So the blocks are enumerated and checked, both
    arms and every stage.
    """

    blocks = seed_blocks()
    # every stage and arm is represented, so the check cannot pass by
    # enumerating too little
    labels = {label for label, _, _ in blocks}
    for rung in P12B_LADDER:
        assert f"verify-{rung}" in labels
        assert f"stage-curved-{rung}" in labels
        assert f"stage-twin-{rung}" in labels
    assert len(blocks) == 3 + 2 + 2 + 3 + 3

    for i, (la, loa, hia) in enumerate(blocks):
        assert loa <= hia
        for lb, lob, hib in blocks[i + 1:]:
            assert hia < lob or hib < loa, (
                f"blocks overlap: {la} [{loa}..{hia}] vs "
                f"{lb} [{lob}..{hib}]")


def test_disjointness_check_aborts_rather_than_reporting():
    """The case: C17's collisions were behind a comment. A checker that
    returns a flag would be the same defect one layer up, so this pins
    that the helper raises.
    """

    import p12_stage_b as B

    original = dict(B.STAGE_B_TWIN_BLOCKS)
    try:
        # land the twin's bottom block exactly on the curved one
        B.STAGE_B_TWIN_BLOCKS[600] = B.STAGE_B_BLOCKS[600]
        with pytest.raises(SystemExit, match="overlap"):
            B.assert_windows_disjoint_and_fresh()
    finally:
        B.STAGE_B_TWIN_BLOCKS.clear()
        B.STAGE_B_TWIN_BLOCKS.update(original)
    # and the real layout passes
    assert B.assert_windows_disjoint_and_fresh()


# ====================================================================
# 5. Delta*_B recomputed from the frozen a, b and m ladder (C25)
# ====================================================================

def test_delta_star_b_is_recomputable_from_its_derivation():
    """The case, from review C25: Section 10.11 named `-0.2516` here,
    the value from the SUPERSEDED 400-sample check, while 10.4 and the
    frozen artifact give `-0.2508`. As specified the test would have
    rejected a correct implementation -- a stale expected value inside a
    regression whose whole job is to catch staleness. So the expected
    value is DERIVED from the frozen coefficients rather than typed
    twice.
    """

    def rel_sd(m):
        return math.sqrt(REL_SD_A / m + REL_SD_B * m ** (-2.0 / 3.0))

    m_bot, m_top = M_TARGET[P12B_LADDER[0]], M_TARGET[P12B_LADDER[-1]]
    derived = math.log10(rel_sd(m_top) / rel_sd(m_bot))
    assert derived == pytest.approx(DELTA_STAR_B, abs=5e-5)
    assert DELTA_STAR_B == pytest.approx(-0.2508, abs=1e-9)
    # and it is STEEPER than the chain-only rate, which is the evidence
    # Section 10.4 offers that it was derived and not borrowed
    chain_only = -math.log10(4.0) / 3.0
    assert derived < chain_only
    # P11's rate must not have been reused wholesale
    assert abs(DELTA_STAR_B - (-0.2007)) > 0.04


def test_delta_eq_b_is_one_third_of_the_target():
    """The case: P11's convention is delta_eq = |Delta*|/3. If the two
    constants ever drift apart the equivalence row silently changes
    meaning.
    """

    assert DELTA_EQ_B == pytest.approx(abs(DELTA_STAR_B) / 3.0, abs=5e-5)


# ====================================================================
# 6. the equivalence coefficient is TWO-SIDED (P13 review C2, ported)
# ====================================================================

def test_equivalence_coefficient_is_two_sided():
    """The case, from P13's review C2: a within-margin verdict needs
    BOTH interval bounds inside the margin, i.e. |Delta_hat| <
    delta_eq - 1.960 se, so with delta_eq / se = 1.960 + z the power at
    Delta = 0 is P(|Z| < z) = 2 Phi(z) - 1, NOT Phi(z). P13 had this
    slot "upgraded" with a one-sided quantile, which broke a convention
    that was already right. This derives the delivered power back out
    of the constant so the slip cannot recur here.
    """

    z = math.sqrt(N_EQ_COEFF) * DELTA_EQ_B - 1.960
    assert z == pytest.approx(1.645, abs=1e-3)
    delivered = 2.0 * 0.5 * (1.0 + math.erf(z / math.sqrt(2.0))) - 1.0
    assert delivered == pytest.approx(0.90, abs=2e-3), (
        "1.645 must deliver a CENTRAL 90%; if this reads 0.95 the "
        "quantile has been re-interpreted as one-sided")


def test_power_selection_rule_and_cap():
    """The case: Section 10.10 freezes the selection rule and a cap of
    1300, with FLAT-WITHIN-MARGIN available only if n_eq fits. A cap
    quietly widened would make an unaffordable verdict look purchasable.
    """

    assert (N_FLOOR, N_CAP) == (12, 1300)

    tiny = np.random.default_rng(11).normal(0.0, 0.05, 200)
    power = power_requirements_b(tiny, tiny.copy())
    assert power["n_eq"] <= N_CAP
    assert power["flat_available"] is True
    assert power["infeasible"] is False
    assert N_FLOOR <= power["n_per_rung"] <= N_CAP

    huge = np.random.default_rng(12).normal(0.0, 3.0, 200)
    wide = power_requirements_b(huge, huge.copy())
    assert wide["n_eq"] > N_CAP
    assert wide["flat_available"] is False
    if wide["n_sup"] > N_CAP:
        assert wide["infeasible"] is True
        assert wide["n_per_rung"] is None


# ====================================================================
# 7. window privacy and freshness, against the RANGES themselves
# ====================================================================

def test_windows_are_fresh_against_every_documented_spent_range():
    """The case: Section 10.11 item 7 asks this be asserted against the
    ranges rather than against a floor, so that a window landing BELOW
    a spent range is caught too. A floor test would pass a block placed
    at 500000, which P11 and P12 already spent.
    """

    blocks = assert_windows_disjoint_and_fresh()
    assert blocks
    for label, lo, hi in blocks:
        assert WINDOW_FLOOR <= lo <= hi <= WINDOW_CEIL, label
        for name, slo, shi in SPENT_RANGES:
            assert hi < slo or shi < lo, f"{label} overlaps {name}"


def test_freshness_check_catches_a_block_below_a_spent_range():
    """The case: this is the defect a floor comparison cannot see, so it
    is constructed explicitly rather than assumed impossible.
    """

    import p12_stage_b as B

    original = dict(B.STAGE_B_BLOCKS)
    try:
        # inside P12's design-check space, i.e. BELOW Stage B's window
        B.STAGE_B_BLOCKS[600] = (2_500_000, 1320)
        with pytest.raises(SystemExit) as excinfo:
            B.assert_windows_disjoint_and_fresh()
        message = str(excinfo.value)
        assert "leaves Stage B's window" in message or "spent" in message
    finally:
        B.STAGE_B_BLOCKS.clear()
        B.STAGE_B_BLOCKS.update(original)


def test_spent_ranges_cover_the_campaigns_the_prereg_names():
    """The case: the freshness test is only as good as its list. If a
    campaign's range were dropped from SPENT_RANGES the check would
    still pass while reusing seeds, so the list's contents are pinned.
    """

    named = {name for name, _, _ in SPENT_RANGES}
    assert any("design check" in n for n in named)
    assert any("v1" in n for n in named)
    assert any("v2" in n for n in named)
    assert any("v3" in n for n in named)
    for _name, lo, hi in SPENT_RANGES:
        assert lo <= hi
        assert hi < WINDOW_FLOOR, "a spent range must end below Stage B's"


# ====================================================================
# 8. the co-requirement returns RATE-ONLY
# ====================================================================

def _recovery(rel_error):
    return {"rel_error": rel_error,
            "rel_error_ci": [rel_error - 0.05, rel_error + 0.05]}


def test_rate_only_when_the_rate_passes_and_recovery_fails():
    """The case, frozen in Section 10.7: "If (i) passes and (ii) fails,
    the record reads RATE-ONLY and states in the same sentence that no
    curvature was recovered -- the gate measured a decaying bias and
    nothing else." Section 5 predicted this exact failure mode ("a gate
    that passes while recovering nothing"), and P13's Section 11 record
    is what it cost to learn that a verdict word gets read as recovery.
    """

    out = stage_b_record(-0.25, [-0.32, -0.18], True, _recovery(0.62))
    assert out["rate_verdict"] == "IMPROVES"
    assert out["recovery_passed"] is False
    assert out["outcome"] == "RATE-ONLY"
    assert "NO CURVATURE WAS RECOVERED" in out["statement"]
    # the size travels with the word (Section 10.9)
    assert "0.620" in out["statement"]


def test_pass_requires_both_halves():
    """The case: the co-requirement is a conjunction. A design that
    returned the good word on either half alone would be the thing
    Section 10.7 was added to prevent.
    """

    both = stage_b_record(-0.25, [-0.32, -0.18], True, _recovery(0.15))
    assert both["outcome"] == "RECOVERS-CURVATURE"
    assert both["recovery_passed"] is True
    # the scope sentence Section 10.9 requires beside any pass
    assert "conformally flat" in both["statement"]
    assert "d >= 3" in both["statement"]

    # recovery alone, with the rate gate inconclusive, is not a pass
    weak = stage_b_record(-0.01, [-0.09, 0.07], True, _recovery(0.10))
    assert weak["rate_verdict"] != "IMPROVES"
    assert weak["outcome"] == weak["rate_verdict"]
    assert weak["outcome"] != "RECOVERS-CURVATURE"


def test_recovery_threshold_is_evaluated_at_exactly_the_frozen_value():
    """The case: 0.25 is a design choice Section 10.7 states as one, and
    the headroom is 1.92 sigma. A threshold nudged after seeing data is
    what preregistration forbids, so the boundary is pinned on both
    sides.
    """

    assert RECOVERY_THRESHOLD == 0.25
    at = stage_b_record(-0.25, [-0.32, -0.18], True,
                        _recovery(RECOVERY_THRESHOLD))
    assert at["recovery_passed"] is True, "the threshold is inclusive"
    just_over = stage_b_record(-0.25, [-0.32, -0.18], True,
                               _recovery(RECOVERY_THRESHOLD + 1e-9))
    assert just_over["recovery_passed"] is False
    assert just_over["outcome"] == "RATE-ONLY"


# ====================================================================
# the cross-stage gate, and the constants it must not confuse
# ====================================================================

def test_cross_stage_gate_reads_p13_campaign_v3_not_v1():
    """The case: `frozen/p13/` holds campaign v1's CONFOUNDED record and
    `frozen/p13v3/` holds the CURVATURE-ROBUST one Section 10.10 names
    ("Stage A-13C"). A gate pointed at the obvious directory would abort
    on a real record, and a gate that accepted any verdict would wave
    the CONFOUNDED one through.
    """

    import json

    assert FROZEN_P13V3_DIR.name == "p13v3"
    v3 = json.loads((FROZEN_P13V3_DIR / "p13_stage_a_summary.json")
                    .read_text(encoding="utf-8"))
    assert v3["verdict"] == "CURVATURE-ROBUST"

    v1_dir = FROZEN_P13V3_DIR.parent / "p13"
    v1 = json.loads((v1_dir / "p13_stage_a_summary.json")
                    .read_text(encoding="utf-8"))
    assert v1["verdict"] == "CONFOUNDED", (
        "if v1's record ever reads CURVATURE-ROBUST this test is stale, "
        "not passing")


def test_prerequisite_stamps_are_reachable_from_head():
    """The case: a stamp that is not an ancestor of HEAD cannot have its
    provenance audited here, which is what the merge policy protects.
    Section 10.10 requires BOTH prerequisites reachable.
    """

    import json

    stamps = []
    for directory, name in (
        (FROZEN_P13V3_DIR.parent / "p12", "p12_stage_a_summary.json"),
        (FROZEN_P13V3_DIR, "p13_stage_a_summary.json"),
    ):
        artifact = json.loads((directory / name).read_text(encoding="utf-8"))
        stamps.append(artifact["code_version"])
    for stamp in stamps:
        assert subprocess.run(
            ["git", "merge-base", "--is-ancestor", stamp, "HEAD"],
            capture_output=True,
        ).returncode == 0, f"{stamp} is not an ancestor of HEAD"


def test_operating_point_constants_come_from_p13s_own_tables():
    """The case: Section 10.8 says "P13's 1.50 rung VERBATIM". Copying
    the numbers instead of importing them is how two modules drift apart
    silently, so this pins the identity.
    """

    import p13_tau_ell as P13

    assert TAU_ELL == 1.50
    assert (ETA_LO, XHALF) == P13.PATCH[1.50][:2]
    assert ETA_LO == -7.0 and XHALF == 17.0
    # the ladder is the top rung's intensity scaled 4x end to end
    import p12_stage_b as B

    assert B.LADDER_RHO[2400] == P13.PATCH[1.50][2]
    assert B.LADDER_RHO[600] == pytest.approx(19.025)
    assert B.LADDER_RHO[1200] == pytest.approx(38.05)
    # the twin intensities are the addendum's CALIBRATED literals, not
    # P13's value scaled -- scaling left the arms 2-3% apart in m and
    # that offset dominated the recovery error (review C11)
    assert B.TWIN_RHO[600] == pytest.approx(2.1948)
    assert B.TWIN_RHO[1200] == pytest.approx(4.3370)
    assert B.TWIN_RHO[2400] == pytest.approx(8.6990)
    assert B.TWIN_RHO[2400] != pytest.approx(P13.TWIN_RHO[1.50])


def test_g_hat_is_order_plus_count_only():
    """The case: the claim being tested by the whole stage is "order +
    number = geometry". If g_hat took a length it would be circular, so
    its signature and behaviour are pinned: doubling rho at fixed
    (m, tau_hat) must halve g.
    """

    assert g_hat(38, 2.0, 3.0) == pytest.approx(38 / (2.0 * 9.0))
    assert g_hat(38, 4.0, 3.0) == pytest.approx(
        0.5 * g_hat(38, 2.0, 3.0))


# ====================================================================
# the amended Verification-B pin, re-derived rather than trusted
# ====================================================================

def _binom_upper_tail(k, n, p):
    """P(X >= k) for X ~ Bin(n, p), exact finite sum (no scipy here)."""

    if p >= 1.0:
        return 1.0 if k <= n else 0.0
    if p <= 0.0:
        return 1.0 if k <= 0 else 0.0
    return min(1.0, sum(math.comb(n, j) * p ** j * (1.0 - p) ** (n - j)
                        for j in range(k, n + 1)))


def _clopper_pearson_lower(k, n, alpha=0.05):
    """One-sided lower bound: the p solving P(X >= k | n, p) = alpha."""

    if k >= n:
        return alpha ** (1.0 / n)
    lo, hi = 0.0, k / n
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _binom_upper_tail(k, n, mid) > alpha:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _negbinom_failures_upper_tail(c, successes, p):
    """P(failures > c) before the successes-th success, exact sum."""

    return max(0.0, 1.0 - sum(
        math.comb(f + successes - 1, f) * p ** successes * (1.0 - p) ** f
        for f in range(0, c + 1)))


def _frozen_curved_completions():
    import json

    path = (FROZEN_P13V3_DIR.parent / "p12"
            / "p12_stage_b_ensemble_check.json")
    ens = json.loads(path.read_text(encoding="utf-8"))
    return {int(rung): (ens["curved"][rung]["n_complete"],
                        ens["curved"][rung]["n_attempted"])
            for rung in ("600", "1200", "2400")}


def test_verification_pins_are_derived_from_frozen_data():
    """The case, and it is review C25's shape caught before it bit: a
    preregistered threshold typed as a literal drifts from its own
    derivation. So the pins are RECOMPUTED here from the frozen
    design-check completions by the stated rule, and the literals in the
    module are asserted equal to the result.

    The rule:  pin = min(1998, largest k with
                         P(X >= k | Bin(2000, p_lo)) >= 0.999)
    with p_lo the one-sided 95% Clopper-Pearson lower bound on the frozen
    curved completion. Inputs are frozen data only, so this test would
    have produced the same pins before any Stage B sample was drawn.
    """

    import p12_stage_b as B

    for rung, (k, n) in _frozen_curved_completions().items():
        p_lo = _clopper_pearson_lower(k, n)
        pin = B.VERIFY_COUNT
        while pin > 0 and _binom_upper_tail(
                pin, B.VERIFY_COUNT, p_lo) < 0.999:
            pin -= 1
        expected = min(1998, pin)
        assert B.VERIFY_PIN_B[rung] == expected, (
            f"rung {rung}: frozen {k}/{n} gives p_lo={p_lo:.6f} and "
            f"pin {expected}, module says {B.VERIFY_PIN_B[rung]}")


def test_each_pin_is_clearable_at_its_own_frozen_rate():
    """The case, which is the defect being repaired: 1998 of 2000 demands
    0.999 and SITS ON THE MEDIAN of Bin(2000, 0.999), so at exactly the
    rate it demands it fails half the time. Every derived pin must be
    clearable with probability >= 0.999 at the pessimistic end of its own
    frozen rate, and the inherited value must be shown NOT to be -- if
    that second half ever passed, the amendment would have lost its
    reason.
    """

    import p12_stage_b as B

    for rung, (k, n) in _frozen_curved_completions().items():
        p_lo = _clopper_pearson_lower(k, n)
        clears = _binom_upper_tail(B.VERIFY_PIN_B[rung], B.VERIFY_COUNT,
                                  p_lo)
        assert clears >= 0.999, (
            f"rung {rung} pin {B.VERIFY_PIN_B[rung]} is cleared only "
            f"{clears:.4f} of the time at p_lo={p_lo:.6f}")

    # and the inherited single value is not clearable at the bottom rung
    k, n = _frozen_curved_completions()[600]
    p_lo = _clopper_pearson_lower(k, n)
    inherited = _binom_upper_tail(1998, B.VERIFY_COUNT, p_lo)
    assert inherited < 0.5, (
        f"the inherited pin 1998 clears {inherited:.4f} of the time at "
        "the bottom rung's frozen rate -- if this is no longer far below "
        "0.999 the amendment is unjustified and should be reverted")
    # the median property, stated directly
    at_demanded_rate = _binom_upper_tail(1998, B.VERIFY_COUNT, 0.999)
    assert 0.5 < at_demanded_rate < 0.8, (
        "1998 should sit near the median of the very rate it demands")


def test_pins_are_per_rung_and_only_ever_loosen():
    """The case: a single pin across a 4x density sweep is the same
    category of error as a grand-mean m gate, and a derivation that could
    TIGHTEN a pin could be used to make a later stage easier to reach by
    lowering an earlier one. Neither is allowed.
    """

    import p12_stage_b as B

    assert set(B.VERIFY_PIN_B) == set(P12B_LADDER)
    assert all(pin <= 1998 for pin in B.VERIFY_PIN_B.values())
    # the rungs genuinely differ, which is why one value cannot serve
    assert len(set(B.VERIFY_PIN_B.values())) > 1


def test_the_pin_is_a_smoke_test_and_the_fill_rule_is_the_real_gate():
    """The case, and the reason the loosening is defensible rather than
    convenient. The pin cannot be the feasibility gate: at n = 2000 the
    healthy rate and the rate at which the campaign starts aborting are
    only about 2.7 sigma apart, so no threshold can be both clearable and
    a reliable feasibility test. This pins the two facts that make the
    division of labour correct -- the campaign is comfortable at the
    frozen rate, and it does abort once the rate collapses.
    """

    need, slots, skipcap = 1300, 1320, 20
    assert slots - need == skipcap, (
        "the slot allowance and the skip cap must agree, since together "
        "they are the feasibility gate")

    # comfortable at the frozen bottom-rung rate
    assert _negbinom_failures_upper_tail(skipcap, need, 0.997) < 1e-4
    # and it really does abort when completion collapses
    assert _negbinom_failures_upper_tail(skipcap, need, 0.98) > 0.5

    # the derived bottom-rung pin catches the collapse that matters...
    import p12_stage_b as B

    caught = 1.0 - _binom_upper_tail(B.VERIFY_PIN_B[600],
                                     B.VERIFY_COUNT, 0.98)
    assert caught > 0.99, f"collapse to 0.98 caught only {caught:.4f}"
    # ...and is NOT expected to catch a drift the campaign survives, so
    # that limitation is asserted rather than left as a surprise
    survivable = 1.0 - _binom_upper_tail(B.VERIFY_PIN_B[600],
                                         B.VERIFY_COUNT, 0.995)
    assert survivable < 0.1
    assert _negbinom_failures_upper_tail(skipcap, need, 0.995) < 1e-3
