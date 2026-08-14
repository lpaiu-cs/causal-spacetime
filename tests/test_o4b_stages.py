"""The O4b execution order, the unbiased prefix, and what survives a
failure.

O4 ran `G3a -> G1/G2 -> G3b`, spent twelve hours, discovered the
instrumentation was broken and kept nothing. These tests hold the
three corrections: G3b first, every drawn point in the estimator, and
a failure that preserves what was already paid for.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o4b_checkpoint as ck  # noqa: E402
import o4b_stages as st  # noqa: E402


class _Acc:
    """Stands in for the empirical-Bernstein accumulator: these tests
    are about WHICH points reach it, not about the interval."""

    def __init__(self) -> None:
        self.values: list[float] = []

    def add(self, z: float) -> None:
        self.values.append(z)

    @property
    def n(self) -> int:
        return len(self.values)


def test_g3b_runs_before_the_rest_of_g1():
    """The order is part of the freeze, not a convenience."""

    assert st.STAGES == ("g3a", "g3b", "g1", "g2")
    assert st.STAGES.index("g3b") < st.STAGES.index("g1")


def test_a_g3a_failure_has_spent_no_fresh_seed():
    """A property of the order, not a promise about it: G3a is a
    frozen case table and a deterministic solver."""

    assert st.fresh_seed_touched("g3a") is False
    assert st.fresh_seed_touched("g3b") is True
    assert st.fresh_seed_touched("g1") is True


def test_every_drawn_point_enters_the_estimator():
    """Including the ones that are not candidates and the ones that
    are not eligible. Filtering either way estimates something else."""

    prefix = st.Prefix(_Acc(), n_avail=4)
    prefix.observe(0.0)                                  # Z = 0
    prefix.observe(0.5, lambda: {"eligible": False,
                                 "fully_testable": False,
                                 "reason": "W_robust <= 0"})
    prefix.observe(0.25, lambda: {"eligible": True,
                                  "fully_testable": True})
    assert prefix.acc.values == [0.0, 0.5, 0.25]


def test_the_estimator_takes_the_point_before_any_predicate_is_asked():
    """Judging first and accumulating second is the same code with
    the bias in it, so the order is checked directly."""

    seen: list[str] = []

    class Watching(_Acc):
        def add(self, z: float) -> None:
            seen.append("accumulate")
            super().add(z)

    def judge() -> dict:
        seen.append("judge")
        return {"eligible": True, "fully_testable": True}

    st.Prefix(Watching(), n_avail=1).observe(0.5, judge)
    assert seen == ["accumulate", "judge"]


def test_zero_window_points_are_accumulated_but_are_not_candidates():
    """The two tallies over one stream. `L_S1 > 0` is the availability
    report's denominator; it is not the estimator's sample."""

    prefix = st.Prefix(_Acc(), n_avail=2)
    for z in (0.0, 0.0, 0.4, 0.6):
        prefix.observe(z, lambda: {"eligible": True,
                                   "fully_testable": True})
    assert prefix.acc.n == 4
    assert prefix.zero_window == 2
    assert prefix.candidates == 2
    report = prefix.report()
    assert report["accumulated_points"] == 4
    assert report["candidates"] == 2
    assert report["accumulates_every_drawn_point"] is True


def test_a_non_candidate_is_never_judged():
    """It is not part of the availability question, and asking would
    charge the budget for an answer nothing reads."""

    def judge() -> dict:                             # pragma: no cover
        raise AssertionError("a Z = 0 point was judged")

    st.Prefix(_Acc(), n_avail=1).observe(0.0, judge)


def test_rates_are_withheld_from_an_incomplete_prefix():
    """A rate from a prefix that stopped early has a denominator
    chosen by the stopping -- the exact defect the fixed prefix
    exists to avoid."""

    prefix = st.Prefix(_Acc(), n_avail=100)
    for _ in range(10):
        prefix.observe(0.5, lambda: {"eligible": True,
                                     "fully_testable": True})
    report = prefix.report()
    assert report["complete"] is False
    assert "rates_withheld" in report
    assert "eligible_rate" not in report
    assert "fully_testable_rate" not in report


def test_a_complete_prefix_reports_on_the_fixed_denominator():
    prefix = st.Prefix(_Acc(), n_avail=4)
    for i in range(4):
        prefix.observe(0.5, lambda i=i: {
            "eligible": True,
            "fully_testable": i < 3,
            "reason": "construction-unavailable"})
    report = prefix.report()
    assert report["complete"] is True
    assert report["eligible_rate"] == 1.0
    assert report["fully_testable_rate"] == 0.75
    assert report["not_fully_testable_by_reason"] == {
        "construction-unavailable": 1}
    assert "fixed prefix, not the scan" in report["rate_denominator"]


def test_a_g3b_failure_preserves_the_g1_prefix():
    """The statistic is already paid for and is valid independently of
    the instrumentation check. Discarding it is what O4 did."""

    failure = st.StageFailure(
        "g3b", "INVALID", "contract mismatch at cluster 12,041",
        preserved={"g1_prefix": {"n": 8_192, "mean_z": 0.31}})
    record = st.incident(failure, {"freeze_sha": "0" * 40})
    assert record["preserved"]["g1_prefix"]["n"] == 8_192
    assert record["verdict"] is None
    assert record["outcome"] == "INVALID"
    assert record["termination_reason"] == (
        "contract mismatch at cluster 12,041")


def test_an_incident_never_carries_a_verdict():
    """Fail-closed means no scientific sentence is published, not that
    nothing is recorded."""

    with pytest.raises(ValueError, match="not a failure outcome"):
        st.StageFailure("g3b", "VERDICT", "whatever")
    with pytest.raises(ValueError, match="not a frozen stage"):
        st.StageFailure("g4", "INVALID", "whatever")


def test_invalid_and_inconclusive_are_not_interchangeable():
    """The first indicts the instrument; the second indicts nothing."""

    assert set(st.OUTCOMES) == {"VERDICT", "INVALID", "INCONCLUSIVE",
                                "ABORT"}
    source = " ".join(Path(st.__file__).read_text(
        encoding="utf-8").split())
    assert "they are not interchangeable" in source
    assert "the first indicts the instrument" in source
    assert "indicts nothing" in source
    for outcome in ("INVALID", "INCONCLUSIVE", "ABORT"):
        assert st.StageFailure("g3b", outcome, "x").outcome == outcome


def test_a_checkpoint_carries_the_rng_position_not_a_draw_count(
        tmp_path):
    """Resuming has to continue the same stream. A count reproduces it
    only if every consumer drew exactly as before."""

    np = pytest.importorskip("numpy")
    rng = np.random.default_rng(40_000_401)
    rng.random(1_000)

    class _B:
        def state(self) -> dict:
            return {"calls": 2_000, "wall_s": 3.5}

    path = st.write_checkpoint(
        tmp_path / "ck.json", "g1_chunk",
        freeze_sha="a" * 40, digest="b" * 64, seed=40_000_401,
        rng=rng, samples=1_000,
        statistics={"mean_z": 0.25}, budget=_B())
    record = ck.read(path)
    assert record["rng_position"]["bit_generator"] == "PCG64"
    assert record["rng_position"] != np.random.default_rng(
        40_000_401).bit_generator.state
    assert record["partial"] is True
    assert record["budget"]["calls"] == 2_000


def test_the_new_campaign_scalars_are_fresh_and_301_is_untouched():
    """A retired stream is never re-entered, and `40,000,301` stays
    withdrawn -- reusing it would make the withdrawal look like a
    deferral."""

    import probe_seed_ledger as ledger

    assert ledger.assert_fresh_scalar("o4b_g1_audit") == 40_000_401
    assert ledger.assert_fresh_scalar("o4b_g2_leakage") == 40_000_411
    assert 40_000_301 not in ledger.spent_scalars()
    assert 40_000_301 not in ledger.FRESH_PROBE_SCALARS.values()
    assert 40_000_281 in ledger.spent_scalars()      # O4, retired
    assert 40_000_291 in ledger.spent_scalars()
