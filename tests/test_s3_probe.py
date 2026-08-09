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


def test_s3_official_seed_is_a_fresh_allocation():
    assert ledger.assert_fresh_scalar("s3_exploration") == s3.SEED
    assert s3.SEED == 40_000_231


def test_smoke_stream_is_separate_and_cannot_burn_the_official_seed():
    """PR review R3: `--smoke` draws only the dedicated observed
    stream; the official allocation stays untouched and unobserved,
    and re-allocating the smoke stream as fresh must abort."""

    assert s3.SMOKE_SEED != s3.SEED
    assert s3.SMOKE_SEED == 40_000_221
    assert s3.SMOKE_SEED in ledger.spent_scalars()
    assert s3.SEED not in ledger.spent_scalars()
    with pytest.raises(SystemExit):
        assert_point_seeds_fresh({"s3_smoke": s3.SMOKE_SEED},
                                 ledger.spent_scalars(),
                                 ledger.SPENT_RANGES, "S3")


def test_aggregated_ledger_carries_every_review_caught_omission():
    """PR review regressions: the hand-copy missed the 2026083x/4x
    seeds, the first union missed S1's BENCH_SEED + 1, and the pilot
    and W1 streams are observed-spent."""

    spent = ledger.spent_scalars()
    for missed in (20_260_831, 20_260_832, 20_260_841, 20_260_842,
                   40_000_102):
        assert missed in spent
    assert ledger.S3_PILOT_SEED in spent
    assert ledger.W1_SEED in spent


def test_observed_scalars_would_be_rejected_as_fresh():
    """The pilot's own seed, the S1 bench seed, and a C1 execution
    seed must all abort a fresh allocation -- reuse of an observed
    stream can never pass as fresh (PR review R2)."""

    spent = ledger.spent_scalars()
    for taken in (40_000_201, 40_000_101, 40_000_061, 20_260_841):
        with pytest.raises(SystemExit):
            assert_point_seeds_fresh({"s3_exploration": taken}, spent,
                                     ledger.SPENT_RANGES, "S3")


def test_p12_decade_would_be_rejected():
    with pytest.raises(SystemExit):
        assert_point_seeds_fresh({"s3_exploration": 34_000_061},
                                 ledger.spent_scalars(),
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
