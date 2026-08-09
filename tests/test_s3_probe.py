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


def test_s3_stream_is_observed_and_reruns_are_replays():
    """The official artifact exists, so 40_000_231 is observed-spent:
    the ledger hands it out only through the replay path, and a fresh
    allocation of it must abort."""

    assert ledger.replay_scalar("s3_exploration") == s3.SEED
    assert s3.SEED == 40_000_231
    assert s3.SEED in ledger.spent_scalars()
    with pytest.raises(SystemExit):
        assert_point_seeds_fresh({"s3_exploration": s3.SEED},
                                 ledger.spent_scalars(),
                                 ledger.SPENT_RANGES, "S3")


def test_smoke_stream_is_separate_from_the_official_stream():
    """PR review R3: `--smoke` draws only the dedicated observed
    stream, and re-allocating it as fresh must abort."""

    assert s3.SMOKE_SEED != s3.SEED
    assert s3.SMOKE_SEED == 40_000_221
    assert s3.SMOKE_SEED in ledger.spent_scalars()
    with pytest.raises(SystemExit):
        assert_point_seeds_fresh({"s3_smoke": s3.SMOKE_SEED},
                                 ledger.spent_scalars(),
                                 ledger.SPENT_RANGES, "S3")


def test_official_artifact_pins_seed_and_clean_lineage():
    """The committed artifact: seed is the observed stream, entry and
    exit git states match and are clean, the delta bounds recompute
    from the raw f arrays, and identified_ci95 is exactly the outer
    endpoints of the two bound summaries."""

    artifact = (Path(__file__).resolve().parents[1] / "docs" / "prereg"
                / "p14_s3_probe_results.json")
    if not artifact.exists():
        pytest.skip("S3 artifact not present in this checkout")
    import json
    r = json.loads(artifact.read_text(encoding="utf-8"))
    assert r["params"]["seed"] == ledger.replay_scalar("s3_exploration")
    assert r["code"]["start"] == r["code"]["end"]
    assert r["code"]["start"]["dirty"] is False
    f0 = r["f_flat"]["per_reading"]
    fl = r["f_schwarzschild_lower"]["per_reading"]
    fh = r["f_schwarzschild_upper"]["per_reading"]
    lo = r["delta_lower"]["per_reading"]
    hi = r["delta_upper"]["per_reading"]
    assert len(f0) == len(fl) == len(fh) == len(lo) == len(hi)
    for a, b, c, d, e in zip(f0, fl, fh, lo, hi, strict=True):
        assert abs(d - (b - a)) < 1e-15
        assert abs(e - (c - a)) < 1e-15
    assert r["identified_ci95"] == [
        r["delta_lower"]["ci95_student_t"][0],
        r["delta_upper"]["ci95_student_t"][1]]


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
