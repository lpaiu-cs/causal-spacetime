"""S3 exploration probe: seed-ledger freshness, both predicates'
geometry, and reading determinism. Exploration-grade -- the probe has
no frozen gate, so these pin the parts a wrong implementation would
silently corrupt: the seed discipline, the flat chord, the Shapiro
sign, and reproducibility."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = Path(__file__).resolve().parents[1] / "experiments" / "positive_control"
sys.path.insert(0, str(EXPERIMENT_DIR))

import s1_schwarzschild_cost as s1  # noqa: E402
import s3_schwarzschild_probe as s3  # noqa: E402
from seed_windows import (  # noqa: E402
    P11_P13_SPENT_RANGES,
    P12_ALLOCATION_DECADE,
    assert_point_seeds_fresh,
)

_RANGES = P11_P13_SPENT_RANGES + (P12_ALLOCATION_DECADE,)


def test_s3_seed_is_fresh_against_the_full_ledger():
    assert_point_seeds_fresh({"s3_exploration": s3.SEED},
                             s3.SPENT_SCALARS, _RANGES, "S3")


def test_spent_scalars_would_be_rejected():
    """The ledger regression: the S1 bench seed and a C1 execution
    seed must both abort, proving the spent set actually guards."""

    for spent in (40_000_101, 40_000_061):
        with pytest.raises(SystemExit):
            assert_point_seeds_fresh({"s3_exploration": spent},
                                     s3.SPENT_SCALARS, _RANGES, "S3")


def test_p12_decade_would_be_rejected():
    with pytest.raises(SystemExit):
        assert_point_seeds_fresh({"s3_exploration": 34_000_061},
                                 s3.SPENT_SCALARS, _RANGES, "S3")


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
    assert s3.flat_relation(p, np.array([chord + 1e-6, 10.0, 0.5, p[3] + np.pi / 2]))


def test_schwarzschild_radial_flight_exceeds_flat():
    """The Shapiro sign, deterministically: the radial coordinate
    flight time (exact tortoise difference) exceeds the flat chord
    time, so the curved census can only lose radial-pair relations
    near the cone -- the mechanism behind the expected delta < 0."""

    t_min, err = s1.flight_time(10.0, 20.0, 0.0)
    assert err == 0.0
    assert t_min > 10.0


def test_reading_is_deterministic():
    """Same seed, same reading -- byte-identical censuses (the
    restart-on-interruption contract for the exploration)."""

    a = s3.reading(np.random.default_rng(12345), n_events=25)
    b = s3.reading(np.random.default_rng(12345), n_events=25)
    assert a == b
