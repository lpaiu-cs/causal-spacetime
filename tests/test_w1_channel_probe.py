"""W1 channel probe: the load-bearing seals are (1) frozen-equivalence
-- the generalized channel predicate must reproduce the FROZEN vacuum
predicate bit-identically on the vacuum signature -- and (2) the
curvature decomposition pin: each channel's Ricci/Weyl split in the
Brinkmann conventions (R_uiuj = -1/2 d_i d_j H), so the channel
labels cannot drift from what the legs integrate (PR review R1 caught
a factor-2 profile error in the first cut's docstring)."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (Path(__file__).resolve().parents[1]
                  / "experiments" / "positive_control")
sys.path.insert(0, str(EXPERIMENT_DIR))

import p14_plane_wave as pw  # noqa: E402
import p14_probe_p3c as p3c  # noqa: E402
import probe_seed_ledger as ledger  # noqa: E402
import w1_channel_probe as w1  # noqa: E402
from p14_plane_wave import Slab, arms, sprinkle  # noqa: E402

_ARTIFACT = (Path(__file__).resolve().parents[1] / "docs" / "prereg"
             / "p14_w1_channel_results.json")


def _operating_point():
    label, w, du, dv, dx, dy = p3c.POINT
    slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
    curved, flat = arms(slab, w)
    rho = p3c.E_N / slab.coordinate_volume
    return w, rho, curved, flat


def test_w1_stream_is_observed_and_reruns_are_replays():
    """W1's results exist, so its seed is observed-spent: the ledger
    hands it out only through the replay path, and a fresh allocation
    of it must abort."""

    assert ledger.replay_scalar("w1_exploration") == w1.SEED
    assert w1.SEED in ledger.spent_scalars()
    assert ledger.S3_PILOT_SEED in ledger.spent_scalars()
    assert ledger.S3_SEED in ledger.spent_scalars()


def test_channel_curvature_decomposition_is_pinned():
    """vacuum: pure Weyl at eigenvalue magnitude w^2 with R_uu = 0.
    ricci: pure trace with R_uu = +2 w^2 (NEC) and zero Weyl.
    mixed: exactly HALF of each pure arm -- the license for the
    midpoint contrast."""

    w = 1.3
    vac = w1.channel_curvature(w1.VACUUM, w)
    ric = w1.channel_curvature(w1.RICCI, w)
    mix = w1.channel_curvature(w1.MIXED, w)

    assert vac["ricci_uu"] == 0.0
    assert vac["weyl_eigenvalues"] == (-w * w, w * w)
    assert vac["h_coefficients"] == (w * w, -w * w)

    assert ric["ricci_uu"] == pytest.approx(2 * w * w)
    assert ric["ricci_uu"] > 0.0
    assert ric["weyl_eigenvalues"] == (0.0, 0.0)

    assert mix["ricci_uu"] == pytest.approx(0.5 * ric["ricci_uu"])
    assert mix["weyl_eigenvalues"][0] == pytest.approx(
        0.5 * vac["weyl_eigenvalues"][0])
    assert mix["weyl_eigenvalues"][1] == pytest.approx(
        0.5 * vac["weyl_eigenvalues"][1])


def test_vacuum_signature_reproduces_the_frozen_predicate():
    """Every pair of a sprinkle at the frozen operating point: same
    related/None verdict, same margin interval, same escalation flag
    as the frozen `causal_relation` -- bit-identical."""

    w, rho, curved, flat = _operating_point()
    rng = np.random.default_rng(97)
    pts = sprinkle(flat, rho, rng)
    pts = pts[np.argsort(pts[:, 0], kind="stable")][:90]
    checked = 0
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            frozen = pw.causal_relation(curved, pts[i], pts[j])
            ours = w1.channel_relation(w1.VACUUM, pts[i], pts[j], w)
            assert ours.related == frozen.related
            assert ours.margin_lower == frozen.margin_lower
            assert ours.margin_upper == frozen.margin_upper
            assert ours.escalated == frozen.escalated
            checked += 1
    assert checked > 3000


def test_ricci_signature_is_xy_symmetric():
    """H = -w^2 (x^2 + y^2) is symmetric under x <-> y, so the
    predicate must be too."""

    w, rho, curved, flat = _operating_point()
    rng = np.random.default_rng(98)
    pts = sprinkle(flat, rho, rng)[:60]
    swapped = pts[:, [0, 1, 3, 2]]
    for i in range(len(pts)):
        for j in range(len(pts)):
            if i == j:
                continue
            a = w1.channel_relation(w1.RICCI, pts[i], pts[j], w)
            b = w1.channel_relation(w1.RICCI, swapped[i], swapped[j], w)
            assert a.related == b.related


def test_mixed_signature_reduces_to_flat_on_the_y_axis():
    """With both events at y = 0 the focusing leg costs nothing, so
    the mixed channel must agree with the frozen FLAT arm."""

    w, rho, curved, flat = _operating_point()
    rng = np.random.default_rng(99)
    pts = sprinkle(flat, rho, rng)[:60].copy()
    pts[:, 3] = 0.0
    disagreements = ambiguous = 0
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            frozen = pw.causal_relation(flat, pts[i], pts[j])
            ours = w1.channel_relation(w1.MIXED, pts[i], pts[j], w)
            if frozen.related is None or ours.related is None:
                ambiguous += 1
            elif ours.related != frozen.related:
                disagreements += 1
    assert disagreements == 0
    assert ambiguous == 0


def test_conjugate_straddle_raises_for_focusing_channels():
    w = 1.0
    p = np.array([0.0, 0.0, 0.1, 0.1])
    q = np.array([np.pi, 0.0, 0.1, 0.1])
    with pytest.raises(ValueError):
        w1.channel_relation(w1.RICCI, p, q, w)


def test_artifact_midpoint_residual_matches_its_own_readings():
    """The committed artifact's headline contrast recomputes from its
    stored per-reading deltas: midpoint_residual =
    mixed - (vacuum + ricci)/2, elementwise and in summary."""

    if not _ARTIFACT.exists():
        pytest.skip("W1 artifact not present in this checkout")
    r = json.loads(_ARTIFACT.read_text(encoding="utf-8"))
    # The committed W1 artifact is a deterministic replay of the
    # observed stream and must say so at the artifact boundary.
    assert r["run_kind"] == "replay"
    assert "40000211" in r["replay_of"]
    vac = r["delta"]["vacuum"]["per_reading"]
    ric = r["delta"]["ricci"]["per_reading"]
    mix = r["delta"]["mixed"]["per_reading"]
    stored = r["contrasts"]["midpoint_residual"]
    recomputed = [m - 0.5 * (v + c)
                  for v, c, m in zip(vac, ric, mix, strict=True)]
    assert len(stored["per_reading"]) == len(recomputed)
    for got, want in zip(stored["per_reading"], recomputed, strict=True):
        assert math.isclose(got, want, rel_tol=0, abs_tol=1e-15)
    assert math.isclose(stored["mean"],
                        sum(recomputed) / len(recomputed),
                        rel_tol=1e-12, abs_tol=1e-18)
