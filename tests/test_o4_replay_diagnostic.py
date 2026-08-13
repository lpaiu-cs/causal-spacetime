"""The O4 G3 replay diagnostic: contract, boundaries, and the algebra.

Two kinds of test live here. The first kind enforces the PI's boundary
list mechanically -- a replay may not observe, may not touch a fresh or
withdrawn scalar, may not reach the reservation authority, and may not
emit anything shaped like a result. The second kind pins the failure
identity the diagnostic exists to expose, so a later edit cannot quietly
change what "the outside probe's second leg" means.

Nothing here asserts a value the campaign was supposed to produce.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import o4_replay_diagnostic as rd  # noqa: E402
import o4_sizing as sz  # noqa: E402
import o4_volume_audit as o4  # noqa: E402
import probe_seed_ledger as ledger  # noqa: E402

_SOURCE = (_REPO / "experiments" / "oracle"
           / "o4_replay_diagnostic.py").read_text(encoding="utf-8")

#: Key substrings that would mean a scientific quantity had leaked into
#: a replay artifact. The same guard the incident record carries, plus
#: the words a verdict would arrive under.
_FORBIDDEN = ("mean", "var", "half_width", "v_s1", "identified",
              "leak", "cp_upper", "concordant", "discordant", "power",
              "band", "verdict", "gate", "status", "estimate")


@pytest.fixture
def covering(monkeypatch):
    """Shrink the frozen stress set so a one-cluster walk counts as
    covering it. The suppression rule is `probed == frozen`, and this
    is the only way to exercise the covering branch without paying for
    100,000 clusters."""

    monkeypatch.setitem(o4.FROZEN, "g3_clusters", 1)


def _walk_keys(node, path=""):
    if isinstance(node, dict):
        for key, value in node.items():
            yield f"{path}.{key}", key
            yield from _walk_keys(value, f"{path}.{key}")
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from _walk_keys(value, f"{path}[{i}]")


# ------------------------------------------------------- it is a replay

def test_the_artifact_declares_itself_a_replay():
    art = rd.assemble(1)
    assert art["run_kind"] == "replay"
    assert art["replay_of"] == "PR #70 incident"
    assert art["replayed_stream"]["seed"] == 40_000_281
    assert art["replayed_stream"]["via"].endswith("replay_scalar")


def test_no_key_anywhere_could_be_read_as_a_result():
    art = rd.assemble(1)
    offenders = [(where, key) for where, key in _walk_keys(art)
                 if any(bad in key.lower() for bad in _FORBIDDEN)]
    assert offenders == [], (
        f"a replay artifact must not carry result-shaped keys: "
        f"{offenders}")


def test_the_scope_says_what_was_not_computed():
    art = rd.assemble(1)
    absent = art["scope"]["not_computed"]
    assert "G1 statistics" in absent
    assert "G2 leakage" in absent
    assert "any verdict" in absent
    assert art["scope"]["does_not_stop_at_first_undecided"] is True


def test_counts_are_labelled_diagnostic_not_estimates(covering):
    art = rd.assemble(1)
    says = art["findings"]["counts_are"]
    assert "not estimates" in says
    assert "diagnostic frequencies" in says


# --------------------------------------------------- seeds and streams

def test_the_diagnostic_never_allocates_a_fresh_scalar():
    assert "assert_fresh_scalar" not in _SOURCE
    assert ledger.FRESH_PROBE_SCALARS == {}


def test_the_withdrawn_scalar_stays_unspent():
    assert rd._WITHDRAWN_UNSPENT == 40_000_301
    assert rd._WITHDRAWN_UNSPENT not in ledger.spent_scalars()


def test_the_ledger_is_verified_rather_than_trusted(monkeypatch):
    monkeypatch.setattr(rd, "replay_scalar", lambda name: 12345)
    with pytest.raises(SystemExit, match="does not match the incident"):
        rd.verify_replay_surface()


def test_a_live_fresh_allocation_stops_the_replay(monkeypatch):
    monkeypatch.setattr(rd, "FRESH_PROBE_SCALARS", {"future": 1})
    with pytest.raises(SystemExit, match="fresh pool is not empty"):
        rd.verify_replay_surface()


def test_spending_the_withdrawn_scalar_stops_the_replay(monkeypatch):
    monkeypatch.setattr(
        rd, "spent_scalars",
        lambda: frozenset({rd._WITHDRAWN_UNSPENT}))
    with pytest.raises(SystemExit, match="withdrawn unspent"):
        rd.verify_replay_surface()


def test_g2_is_not_replayed():
    assert "o4_aborted_g2" not in _SOURCE.replace(
        "o4_aborted_g2 is not replayed", "")
    assert "run_g2" not in _SOURCE


def _identifiers() -> set[str]:
    """Every name the module actually uses as code.

    Read from the AST, not from the text: the module DESCRIBES the
    reservation in prose, and a substring scan cannot tell a sentence
    about a ref from a call that writes one."""

    import ast
    tree = ast.parse(_SOURCE)
    return {node.attr for node in ast.walk(tree)
            if isinstance(node, ast.Attribute)} | {
        node.id for node in ast.walk(tree)
        if isinstance(node, ast.Name)} | {
        alias.name.split(".")[0] for node in ast.walk(tree)
        if isinstance(node, ast.Import) for alias in node.names}


def test_the_reservation_is_neither_claimed_nor_deleted():
    used = _identifiers()
    for forbidden in ("reserve_remote", "probe_reservation_namespace",
                      "remote_reservation", "reservation_authority",
                      "_make_commit", "subprocess"):
        assert forbidden not in used, (
            f"a replay must not reach the reservation authority, but "
            f"the module calls {forbidden!r}")


# ------------------------------------------------- the replayed surface

def test_the_replay_surface_matches_the_executed_freeze():
    surface = rd.verify_replay_surface()
    snapshot = json.loads(rd._EXECUTED.read_text(encoding="utf-8"))
    assert set(surface["files_verified"]) == (
        set(snapshot["files"]) - {rd._LEDGER_REL})
    assert surface["declared_exception"]["path"] == rd._LEDGER_REL


def test_every_numeric_file_is_checked_not_just_some(monkeypatch,
                                                     tmp_path):
    """The exemption is one named file, and it is the only one."""

    snapshot = json.loads(rd._EXECUTED.read_text(encoding="utf-8"))
    for rel in snapshot["files"]:
        if rel == rd._LEDGER_REL:
            continue
        broken = dict(snapshot)
        broken["files"] = dict(snapshot["files"]) | {rel: "00" * 32}
        path = tmp_path / f"snapshot-{rel.replace('/', '_')}.json"
        path.write_text(json.dumps(broken), encoding="utf-8")
        monkeypatch.setattr(rd, "_EXECUTED", path)
        with pytest.raises(SystemExit, match="differ from the executed"):
            rd.verify_replay_surface()


def test_environment_drift_stops_the_replay(monkeypatch, tmp_path):
    snapshot = json.loads(rd._EXECUTED.read_text(encoding="utf-8"))
    broken = dict(snapshot)
    broken["environment"] = dict(snapshot["environment"]) | {
        "numpy": "0.0.0-not-this-one"}
    path = tmp_path / "snapshot.json"
    path.write_text(json.dumps(broken), encoding="utf-8")
    monkeypatch.setattr(rd, "_EXECUTED", path)
    with pytest.raises(SystemExit, match="environment drift"):
        rd.verify_replay_surface()


def test_a_missing_snapshot_is_fatal_not_skipped(monkeypatch,
                                                tmp_path):
    monkeypatch.setattr(rd, "_EXECUTED", tmp_path / "absent.json")
    with pytest.raises(SystemExit, match="nothing to reproduce"):
        rd.verify_replay_surface()


# ------------------------------------------------ the prefix is bounded

def test_the_prefix_stops_at_the_clusters_asked_for():
    got = list(rd.stress_clusters(rd._STREAM_SEED, 3,
                                  o4.FROZEN["tol"]))
    assert [c[0] for c in got] == [0, 1, 2]


def test_the_prefix_is_the_frozen_g1_draw():
    """Same generator, same chunking, same acceptance test -- the
    clusters must be G1's own, not a re-derivation."""

    import numpy as np
    rng = np.random.default_rng(rd._STREAM_SEED)
    rs, ths = o4._draw(rng, o4.CHUNK, sz.R_LO, sz.R_HI, sz.PSI_MAX)
    expected = []
    for r, th in zip(rs, ths, strict=True):
        ell, t1, t2 = o4._ell(float(r), float(th), o4.FROZEN["tol"])
        if ell > 0.0:
            expected.append((float(r), float(th), t1, t2))
        if len(expected) == 2:
            break
    got = list(rd.stress_clusters(rd._STREAM_SEED, 2,
                                  o4.FROZEN["tol"]))
    assert [c[1:] for c in got] == expected


def test_probing_past_the_frozen_stress_set_is_refused(monkeypatch):
    monkeypatch.setattr(
        sys, "argv",
        ["o4_replay_diagnostic.py", "--clusters",
         str(o4.FROZEN["g3_clusters"] + 1)])
    with pytest.raises(SystemExit, match="stress set is frozen"):
        rd.main()


# ------------------------------------------------------- the algebra

def test_the_outside_probe_second_leg_carries_the_fixed_offset():
    """The identity the incident record derived without executing:
    at the outside probe, the `x->q` leg's decision distance is the
    fixed offset itself, with no geometry in it at all. That is why
    `err >= 1e-6` decides the outcome on its own."""

    index, r, th, t1, t2 = next(iter(
        rd.stress_clusters(rd._STREAM_SEED, 1, o4.FROZEN["tol"])))
    rec = rd.probe_cluster(index, r, th, t1, t2, o4.FROZEN["tol"])
    outside = next(p for p in rec["probes"]
                   if p["probe"] == rd.OUTSIDE)
    leg = next(lg for lg in outside["legs"] if lg["leg"] == rd.LEG_XQ)
    off = o4.FROZEN["g3_stress_offset"]
    assert leg["abs_dt_minus_t_min"] == pytest.approx(off, abs=1e-12)
    assert leg["undecided"] is (leg["err"] >= leg[
        "abs_dt_minus_t_min"])


def test_the_midpoint_probe_carries_half_the_window():
    index, r, th, t1, t2 = next(iter(
        rd.stress_clusters(rd._STREAM_SEED, 1, o4.FROZEN["tol"])))
    rec = rd.probe_cluster(index, r, th, t1, t2, o4.FROZEN["tol"])
    half = 0.5 * rec["window"]["L"]
    mid = next(p for p in rec["probes"] if p["probe"] == rd.MIDPOINT)
    for leg in mid["legs"]:
        assert leg["abs_dt_minus_t_min"] == pytest.approx(half,
                                                          rel=1e-9)


def test_the_readout_agrees_with_the_predicate():
    """Every leg re-derives what `causal_relation` decided on; if the
    two ever parted the diagnostic would be reading the wrong call."""

    for index, r, th, t1, t2 in rd.stress_clusters(
            rd._STREAM_SEED, 5, o4.FROZEN["tol"]):
        rec = rd.probe_cluster(index, r, th, t1, t2,
                               o4.FROZEN["tol"])
        for probe in rec["probes"]:
            for leg in probe["legs"]:
                assert leg["rederivation_agrees"], (probe, leg)


def test_dpsi_is_the_predicates_angle_not_theta():
    """The predicate recovers the angle through `acos(cos theta)`, so
    the diagnostic must use that value and not the `theta` handed to
    `_ell` -- otherwise `t_min` could differ in the last ulp from the
    one the decision actually used."""

    import numpy as np
    p = np.array([0.0, sz.R_IN, 0.0, 0.0])
    x = np.array([1.0, 15.0, 0.4, 0.0])
    assert rd._dpsi(p, x) == math.acos(max(-1.0, min(1.0,
                                                     math.cos(0.4))))


def test_a_negative_dt_is_recorded_as_the_short_circuit():
    leg = rd._leg(sz.R_IN, 15.0, -1.0, 0.4, o4.FROZEN["tol"],
                  rd.LEG_PX, False)
    assert leg["short_circuited_negative_dt"] is True
    assert leg["undecided"] is False
    assert leg["rederivation_agrees"] is True
    assert "t_min" not in leg


# ------------------------------------------- the walk does not give up

def _fake_leg(leg: str, undecided: bool) -> dict:
    """A leg carrying every field the real one does, so a consumer that
    reads the record (`main`'s report, the artifact) is exercised on
    the same shape."""

    gap = 1e-7 if undecided else 1.0
    err = 1e-6 if undecided else 1e-9
    return {"leg": leg, "dt": 1.0, "dpsi": 0.1, "t_min": 1.0 - gap,
            "err": err, "abs_dt_minus_t_min": gap,
            "decision_margin": gap - err,
            "short_circuited_negative_dt": False,
            "undecided": undecided,
            "returned": "undecided" if undecided else "true",
            "rederivation_agrees": True}


def _fake_cluster(index: int, *, undecided_probe=None,
                  mismatch_probe=None) -> dict:
    probes = []
    for kind in (rd.MIDPOINT, rd.OUTSIDE):
        undecided = kind == undecided_probe
        probes.append({
            "probe": kind, "t_x": 1.0, "want": kind == rd.MIDPOINT,
            "undecided": undecided,
            "boolean_mismatch": (not undecided
                                 and kind == mismatch_probe),
            "legs": [_fake_leg(rd.LEG_PX, undecided),
                     _fake_leg(rd.LEG_XQ, False)],
        })
    return {"cluster_index": index, "r": 15.0, "theta": 0.1,
            "T1": 1.0, "T2": 1.0,
            "window": {"lo": 1.0, "hi": 7.5, "L": 6.5},
            "probes": probes}


@pytest.fixture
def synthetic(monkeypatch):
    """Five clusters, three of them defective, with the solver and the
    draw stubbed out so the walk's own logic is what is under test."""

    plan = {
        0: {},
        1: {"undecided_probe": rd.OUTSIDE},
        2: {},
        3: {"undecided_probe": rd.MIDPOINT},
        4: {"mismatch_probe": rd.OUTSIDE},
    }
    monkeypatch.setattr(
        rd, "stress_clusters",
        lambda seed, clusters, tol, progress=None: (
            (i, 15.0, 0.1, 1.0, 1.0) for i in range(clusters)))
    monkeypatch.setattr(
        rd, "probe_cluster",
        lambda index, r, th, t1, t2, tol: _fake_cluster(
            index, **plan[index]))
    return plan


def test_the_walk_continues_past_every_undecided(synthetic):
    """Clusters 3 and 4 lie past the first undecided at cluster 1, and
    both are still reached -- G3 could not have seen either."""

    found = rd.walk(5, o4.FROZEN["tol"])
    assert found["clusters_probed"] == 5
    assert {s["cluster_index"] for s in found["first_site_per_cause"]
            } == {1, 3, 4}


def test_the_first_undecided_is_the_first_one_not_the_worst(synthetic):
    found = rd.walk(5, o4.FROZEN["tol"])
    assert found["first_undecided"]["cluster_index"] == 1
    assert found["first_undecided"]["cause"] == rd.CAUSE_OUTSIDE


def test_every_cause_is_counted_separately(synthetic, monkeypatch):
    monkeypatch.setitem(o4.FROZEN, "g3_clusters", 5)
    counts = {row["cause"]: row["clusters"]
              for row in rd.walk(5, o4.FROZEN["tol"])["counts"]}
    assert counts == {rd.CAUSE_MIDPOINT: 1, rd.CAUSE_OUTSIDE: 1,
                      rd.CAUSE_MISMATCH: 1}


def test_a_partial_walk_reports_no_frequency_at_all(synthetic):
    """A count over the leading prefix is not a count over the frozen
    stress set, and a published partial artifact would read as though
    it were."""

    found = rd.walk(5, o4.FROZEN["tol"])
    assert found["covers_frozen_stress_set"] is False
    for withheld in ("counts", "counts_are", "clusters_with_no_cause"):
        assert withheld not in found
    assert "not a count over the frozen stress set" in (
        found["counts_withheld"])


def test_a_partial_walk_still_reports_its_sites(synthetic):
    """Suppression is of frequencies, not of facts: "cluster 1 was
    undecided here, with these numbers" is exact about a prefix."""

    found = rd.walk(5, o4.FROZEN["tol"])
    assert found["first_undecided"]["cluster_index"] == 1
    assert len(found["first_site_per_cause"]) == 3


def test_the_console_does_not_print_a_partial_frequency(
        monkeypatch, capsys, synthetic):
    monkeypatch.setattr(sys, "argv",
                        ["o4_replay_diagnostic.py", "--clusters", "5"])
    rd.main()
    out = capsys.readouterr().out
    assert "clean clusters" not in out
    assert "not a count over the frozen stress set" in out


def test_each_cause_records_its_first_site(synthetic):
    firsts = {s["cause"]: s["cluster_index"]
              for s in rd.walk(5, o4.FROZEN["tol"])
              ["first_site_per_cause"]}
    assert firsts == {rd.CAUSE_MIDPOINT: 3, rd.CAUSE_OUTSIDE: 1,
                      rd.CAUSE_MISMATCH: 4}


def test_a_full_walk_does_report_its_frequencies(synthetic,
                                                 monkeypatch):
    monkeypatch.setitem(o4.FROZEN, "g3_clusters", 5)
    found = rd.walk(5, o4.FROZEN["tol"])
    assert found["covers_frozen_stress_set"] is True
    assert found["clusters_with_no_cause"] == 2
    assert "counts_withheld" not in found


# ------------------------------------------------------- publication

def test_nothing_is_written_without_an_out_path(monkeypatch, tmp_path):
    monkeypatch.setattr(sys, "argv",
                        ["o4_replay_diagnostic.py", "--clusters", "1"])
    rd.main()
    assert list(tmp_path.iterdir()) == []


def test_publication_is_no_clobber(monkeypatch, tmp_path, synthetic):
    out = tmp_path / "diag.json"
    monkeypatch.setattr(sys, "argv",
                        ["o4_replay_diagnostic.py", "--clusters", "5",
                         "--out", str(out)])
    rd.main()
    assert json.loads(out.read_text(encoding="utf-8"))["run_kind"] == (
        "replay")
    with pytest.raises(SystemExit, match="write-once"):
        rd.main()


def test_the_published_artifact_names_the_incident_it_replays(
        monkeypatch, tmp_path, synthetic):
    out = tmp_path / "diag.json"
    monkeypatch.setattr(sys, "argv",
                        ["o4_replay_diagnostic.py", "--clusters", "5",
                         "--out", str(out)])
    rd.main()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert art["incident_record"] == "docs/prereg/p14_o4_incident.json"
    assert (_REPO / art["incident_record"]).exists()
    assert "never an observation" in art["not_a_result"]
