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


def test_the_availability_counters_freeze_at_the_n_availth_candidate():
    """The R6 defect, and the whole reason the fixed prefix exists.

    The scan has to continue past `N_avail` -- G3b needs `K_G3B`
    fully-testable clusters and a prefix with even one unavailable
    candidate will not have supplied them. Letting the counters run on
    puts results-chosen points in the numerator over a fixed
    denominator; a rate above 1 is the visible form, a slightly
    inflated rate the invisible one."""

    prefix = st.Prefix(_Acc(), n_avail=3)
    for _ in range(5):
        prefix.observe(0.5, lambda: {"eligible": True,
                                     "fully_testable": True})
    assert prefix.acc.n == 5                 # the estimator took all 5
    assert prefix.candidates == 3            # the report froze at 3
    assert prefix.eligible == 3
    assert prefix.fully_testable == 3
    assert prefix.scan_candidates == 2
    assert prefix.scan_fully_testable == 2

    report = prefix.report()
    assert report["eligible_rate"] == 1.0    # not 5/3
    assert report["fully_testable_rate"] == 1.0
    assert report["beyond_the_prefix"]["candidates"] == 2


def test_the_k_g3b_count_spans_the_prefix_and_the_scan():
    """Only the AVAILABILITY report is confined to the prefix. The
    contract sample is every fully-testable cluster."""

    prefix = st.Prefix(_Acc(), n_avail=2)
    for i in range(6):
        prefix.observe(0.5, lambda i=i: {"eligible": True,
                                         "fully_testable": i % 2 == 0,
                                         "reason": "unavailable"})
    assert prefix.fully_testable == 1        # of the first 2
    assert prefix.scan_fully_testable == 2   # of the next 4
    assert prefix.total_fully_testable == 3
    assert prefix.report()[
        "total_fully_testable_for_k_g3b"] == 3
    assert prefix.report()["fully_testable_rate"] == 0.5


def test_zero_window_points_past_the_prefix_are_also_split():
    """They still enter the estimator; they just stop entering the
    fixed report."""

    prefix = st.Prefix(_Acc(), n_avail=1)
    prefix.observe(0.0)
    prefix.observe(0.5, lambda: {"eligible": True,
                                 "fully_testable": True})
    prefix.observe(0.0)
    assert prefix.acc.n == 3
    assert prefix.zero_window == 1
    assert prefix.scan_zero_window == 1


def test_a_reused_run_context_cannot_overrule_the_failure():
    """The R5 defect, and the same boundary as the checkpoint's. The
    incident is write-once provenance: a context carrying `stage` or
    `verdict` would file a G3b failure as a G1 result."""

    failure = st.StageFailure("g3b", "INVALID", "contract mismatch")
    for key in st.INCIDENT_KEYS:
        with pytest.raises(ValueError, match="the incident itself owns"):
            st.incident(failure, {key: "g1"})
    record = st.incident(failure, {"freeze_sha": "0" * 40})
    assert record["stage"] == "g3b"
    assert record["verdict"] is None


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
        rng=rng, rng_stream="o4b_g1_audit", samples=1_000,
        statistics={"mean_z": 0.25}, budget=_B())
    record = ck.read(path)
    assert record["rng_position"]["bit_generator"] == "PCG64"
    assert record["rng_position"] != np.random.default_rng(
        40_000_401).bit_generator.state
    assert record["rng_stream"] == "o4b_g1_audit"
    assert record["partial"] is True
    assert record["budget"]["calls"] == 2_000


def test_the_document_freezes_the_counter_freeze():
    doc = (_REPO / "docs" / "prereg"
           / "p14_o4_g3_prereg_reopen.md").read_text(encoding="utf-8")
    assert "### 4.2 availability 카운터는 `N_avail` 번째 후보에서" in doc
    assert "eligible_rate > 1" in doc
    assert "beyond_the_prefix" in doc
    assert "availability 보고뿐" in doc


def test_the_campaign_scalars_are_retired_and_301_is_untouched():
    """After the 2026-08-15 run the O4b scalars are OBSERVED, not
    fresh: the run drew both streams to their frozen sample sizes, so
    they are spent whether or not a verdict was published. A retired
    stream is never re-entered, and `40,000,301` stays withdrawn --
    reusing it would make the withdrawal look like a deferral."""

    import probe_seed_ledger as ledger
    import pytest

    # spent, and now in OBSERVED under their functional names
    assert ledger.OBSERVED_PROBE_SCALARS["o4b_g1_audit"] == 40_000_401
    assert ledger.OBSERVED_PROBE_SCALARS["o4b_g2_leakage"] == 40_000_411
    assert 40_000_401 in ledger.spent_scalars()
    assert 40_000_411 in ledger.spent_scalars()
    assert ledger.FRESH_PROBE_SCALARS == {}

    # a retired name cannot be re-entered as a fresh allocation
    for retired in ("o4b_g1_audit", "o4b_g2_leakage"):
        with pytest.raises(KeyError):
            ledger.assert_fresh_scalar(retired)
        # but a labelled reproduction can read the observed seed
        assert ledger.replay_scalar(retired) in (40_000_401, 40_000_411)

    assert 40_000_301 not in ledger.spent_scalars()
    assert 40_000_301 not in ledger.FRESH_PROBE_SCALARS.values()
    assert 40_000_281 in ledger.spent_scalars()      # O4, retired
    assert 40_000_291 in ledger.spent_scalars()
