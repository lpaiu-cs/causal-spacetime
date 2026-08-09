"""W1 channel probe: the load-bearing seal is frozen-equivalence --
the generalized channel predicate must reproduce the FROZEN vacuum
predicate bit-identically on the vacuum signature, so that the ricci
and mixed censuses differ from the frozen chain only in the channel
signature, never in the float path, error model, or escalation."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
sys.path.insert(0, str(EXPERIMENT_DIR))

import p14_plane_wave as pw  # noqa: E402
import p14_probe_p3c as p3c  # noqa: E402
import w1_channel_probe as w1  # noqa: E402
from p14_plane_wave import Slab, arms, sprinkle  # noqa: E402
from seed_windows import (  # noqa: E402
    P11_P13_SPENT_RANGES,
    P12_ALLOCATION_DECADE,
    assert_point_seeds_fresh,
)

_RANGES = P11_P13_SPENT_RANGES + (P12_ALLOCATION_DECADE,)


def _operating_point():
    label, w, du, dv, dx, dy = p3c.POINT
    slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
    curved, flat = arms(slab, w)
    rho = p3c.E_N / slab.coordinate_volume
    return w, rho, curved, flat


def test_w1_seed_is_fresh_against_the_full_ledger():
    assert_point_seeds_fresh({"w1_exploration": w1.SEED},
                             w1.SPENT_SCALARS, _RANGES, "W1")


def test_s3_seed_is_now_spent():
    with pytest.raises(SystemExit):
        assert_point_seeds_fresh({"w1_exploration": 40_000_201},
                                 w1.SPENT_SCALARS, _RANGES, "W1")


def test_vacuum_signature_reproduces_the_frozen_predicate_bit_identically():
    """Every pair of a sprinkle at the frozen operating point: same
    related/None verdict, same margin interval, same escalation flag
    as the frozen `causal_relation`."""

    w, rho, curved, flat = _operating_point()
    rng = np.random.default_rng(97)
    pts = sprinkle(flat, rho, rng)
    pts = pts[np.argsort(pts[:, 0], kind="stable")][:90]
    checked = 0
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            frozen = pw.causal_relation(curved, pts[i], pts[j])
            ours = w1.channel_relation(w1.VACUUM, pts[i], pts[j], w)
            assert ours.related is frozen.related or ours.related == frozen.related
            assert ours.margin_lower == frozen.margin_lower
            assert ours.margin_upper == frozen.margin_upper
            assert ours.escalated == frozen.escalated
            checked += 1
    assert checked > 3000


def test_ricci_signature_is_xy_symmetric():
    """H = -(w^2/2)(x^2 + y^2) is symmetric under x <-> y, so the
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
