"""S3 exploration probe: seed-ledger freshness, both predicates'
geometry, ambiguity bounds, and reading determinism. Exploration-grade
-- the probe has no frozen gate, so these pin the parts a wrong
implementation would silently corrupt: the aggregated seed discipline,
the flat chord, the Shapiro sign, the census interval, and
reproducibility."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (Path(__file__).resolve().parents[1]
                  / "experiments" / "positive_control")
sys.path.insert(0, str(EXPERIMENT_DIR))

import probe_seed_ledger as ledger  # noqa: E402
import s1_schwarzschild_cost as s1  # noqa: E402
import s3_schwarzschild_probe as s3  # noqa: E402
from seed_windows import assert_point_seeds_fresh  # noqa: E402


def test_s3_seed_is_fresh_against_the_aggregated_ledger():
    ledger.assert_probe_seed_fresh("s3_exploration")


def test_aggregated_ledger_contains_the_seeds_the_hand_copy_missed():
    """PR review regression: the first cut hand-copied the spent set
    and missed these; the aggregation must carry them."""

    spent = ledger.spent_scalars("s3_exploration")
    for missed in (20_260_831, 20_260_832, 20_260_841, 20_260_842):
        assert missed in spent
    assert ledger.W1_SEED in spent


def test_spent_scalars_would_be_rejected():
    """The S1 bench seed and a C1 execution seed must both abort,
    proving the aggregated set actually guards."""

    spent = ledger.spent_scalars("s3_exploration")
    for taken in (40_000_101, 40_000_061, 20_260_841):
        with pytest.raises(SystemExit):
            assert_point_seeds_fresh({"s3_exploration": taken}, spent,
                                     ledger.SPENT_RANGES, "S3")


def test_p12_decade_would_be_rejected():
    with pytest.raises(SystemExit):
        assert_point_seeds_fresh({"s3_exploration": 34_000_061},
                                 ledger.spent_scalars("s3_exploration"),
                                 ledger.SPENT_RANGES, "S3")


def test_flat_predicate_is_the_exact_chord_cone():
    """Radial pair r=10 -> r=12 at identical angles: chord 2, so
    dt=2.1 is related and dt=1.9 is not; a transverse pair follows
    the spherical law of cosines."""

    p = np.array([0.0, 10.0, 0.5, 0.3])
    assert s3.flat_relation(p, np.array([2.1, 12.0, 0.5, 0.3]))
    assert not s3.flat_relation(p, np.array([1.9, 12.0, 0.5, 0.3]))
    q = np.array([0.0, 10.0, 0.5, 0.3 + np.pi / 2])
    chord = s3.flat_chord(p, q)
    assert 0.0 < chord < 20.0
    later = np.array([chord + 1e-6, 10.0, 0.5, p[3] + np.pi / 2])
    assert s3.flat_relation(p, later)


def test_schwarzschild_radial_flight_exceeds_flat():
    """The Shapiro sign, deterministically: the radial coordinate
    flight time (exact tortoise difference) exceeds the flat chord
    time, so the curved census can only lose radial-pair relations
    near the cone -- the mechanism behind the expected delta < 0.
    Radial pairs only; the full-domain sign is what the probe
    measures, in the frozen coordinates."""

    t_min, err = s1.flight_time(10.0, 20.0, 0.0)
    assert err == 0.0
    assert t_min > 10.0


def test_reading_reports_coherent_bounds_and_is_deterministic():
    """Same seed, same reading -- identical censuses (the
    restart-on-interruption contract); and the interval census obeys
    f_lower <= f_upper with the gap equal to ambiguous/pairs."""

    a = s3.reading(np.random.default_rng(12345), n_events=25)
    b = s3.reading(np.random.default_rng(12345), n_events=25)
    assert a == b
    fm_lo, fm_hi, f0, ambiguous, escalated = a
    pairs = 25 * 24 // 2
    assert 0.0 <= fm_lo <= fm_hi <= 1.0
    assert fm_hi - fm_lo == pytest.approx(ambiguous / pairs)
    assert 0.0 <= f0 <= 1.0
