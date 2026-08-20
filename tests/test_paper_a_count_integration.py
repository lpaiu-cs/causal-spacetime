"""Paper A <-> prediction-anchored count-stage integration contracts.

Section 6.8 and the claim boundary now CLAIM the per-rung Poisson-count
verdicts, so the numbers they print must be the artifacts' numbers and
the ladder they describe must be the ladder that ran. These tests hold
three lines: the cited figures are derived from the artifacts (not
typed), every rung actually returned CONCORDANT, and the promotion is
bounded -- per-rung only, no joint verdict, no interpolation.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
_PREREG = REPO / "docs" / "prereg"
_PAPER = REPO / "docs" / "paper" / "paper_a"

sys.path.insert(0, str(REPO / "experiments" / "oracle"))

import s6_rungs as s6  # noqa: E402

MANUSCRIPT = _PAPER / "manuscript.md"
CLAIM_BOUNDARY = _PAPER / "claim_boundary.md"

#: (mu as printed, artifact) for the three executed rungs
RUNGS = (
    ("0.1333", "p14_o5_count.json"),
    ("0.1867", "p14_s6_m14_count.json"),
    ("0.2400", "p14_s6_m18_count.json"),
    ("0.4000", "p14_s6_m30_count.json"),
)


def _art(name: str) -> dict:
    return json.loads((_PREREG / name).read_text(encoding="utf-8"))


def _flat(path: Path) -> str:
    return " ".join(path.read_text(encoding="utf-8").split())


def test_every_cited_rung_actually_returned_concordant():
    for _, name in RUNGS:
        a = _art(name)
        assert a["decision"]["verdict"] == "CONCORDANT", name


def _sec68_rows() -> dict:
    """The Section 6.8 table, parsed into {mu: [cells]} so a value can
    be checked against the rung it is printed ON -- a whole-document
    substring search would still pass if two rungs' rows were swapped."""

    rows = {}
    for line in MANUSCRIPT.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) == 11 and re.fullmatch(r"0\.\d{4}", cells[0]):
            rows[cells[0]] = cells
    return rows


def test_the_manuscript_table_matches_each_rung_row_by_row():
    """Every printed cell is re-derived from THAT rung's artifact:
    mass, certified volume, N, K_certain, U_amb, C, D, band and the
    verdict must all line up on one row."""

    rows = _sec68_rows()
    assert set(rows) == {mu for mu, _ in RUNGS}, sorted(rows)

    for mu, name in RUNGS:
        a = _art(name)
        d, scan, fz = a["decision"], a["scan"], a["frozen_config"]
        cells = rows[mu]
        # the mass and the compactness on this row must be the same
        # rung under the ladder's one exact path -- neither is trusted
        assert f"{s6.mu(float(cells[1])):.4f}" == mu, (mu, "mu<->M")
        assert cells[2] == f"[{fz['v_lo']:.6f}, {fz['v_hi']:.6f}]", (mu, "V")
        assert cells[4] == f"{scan['n_total']:,}", (mu, "N")
        assert cells[5] == f"{scan['k_certain']:,}", (mu, "K")
        assert cells[6] == str(scan["u_ambiguous"]), (mu, "U")
        assert cells[7] == (f"[{d['c_lo_volume']:.6f}, "
                            f"{d['c_hi_volume']:.6f}]"), (mu, "C")
        assert cells[8] == f"[{d['d_lo']:.6f}, {d['d_hi']:+.6f}]", (mu, "D")
        assert cells[9] == f"{d['band']:.6f}", (mu, "band")
        assert d["verdict"] in cells[10], (mu, "verdict")


def test_the_row_check_would_catch_a_swap():
    """Guard the guard: the deeper rungs differ in every checked
    column, so swapping their rows really would fail the row-by-row
    test above -- the check is not vacuous."""

    rows = _sec68_rows()
    mus = [mu for mu, _ in RUNGS]
    by_mu = dict(RUNGS)
    a14, a18 = _art(by_mu[mus[1]]), _art(by_mu[mus[2]])
    assert a14["scan"]["k_certain"] != a18["scan"]["k_certain"]
    assert a14["scan"]["u_ambiguous"] != a18["scan"]["u_ambiguous"]
    assert a14["decision"]["band"] != a18["decision"]["band"]
    assert rows[mus[1]][5] != rows[mus[2]][5]
    assert rows[mus[1]][6] != rows[mus[2]][6]


def test_the_claim_boundary_prints_each_rung_as_one_statement():
    """The claim boundary states mu, K, U, D and the band together, so
    the correspondence is checked as a single string per rung."""

    flat = _flat(CLAIM_BOUNDARY)
    for mu, name in RUNGS:
        a = _art(name)
        d, scan = a["decision"], a["scan"]
        stmt = (f"`mu = {mu}` K = {scan['k_certain']:,}, "
                f"U = {scan['u_ambiguous']}, "
                f"D = [{d['d_lo']:.6f}, {d['d_hi']:+.6f}] "
                f"in +-{d['band']:.6f}")
        assert stmt in flat, (mu, stmt)
        assert name in flat, name


def test_the_containment_the_paper_claims_actually_holds():
    """CONCORDANT is not taken on trust: the identified discrepancy is
    re-checked against the band from the serialized values."""

    for _, name in RUNGS:
        d = _art(name)["decision"]
        assert d["d_lo"] >= -d["band"], name
        assert d["d_hi"] <= d["band"], name


def test_the_ladder_is_the_three_executed_rungs_in_order():
    """The paper's mu ladder must be the executed masses, ascending,
    and genuinely distinct -- no rung is an isometric copy."""

    mus = [float(mu) for mu, _ in RUNGS]
    assert mus == sorted(mus) and len(set(mus)) == 4
    flat = _flat(MANUSCRIPT)
    assert "no rung is an isometric" in flat
    assert "mu = 2M/r_c" in flat


def test_the_stated_rung_COUNT_tracks_the_executed_ladder():
    """Prose that COUNTS rungs must move when the ladder does. Promoting
    a rung used to update the claim row while leaving "the three rung
    verdicts" and "carries the three per-rung count verdicts" behind, so
    the claim boundary and the archived-state authority disagreed about
    whether the new artifact is inside the claim."""

    word = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}[len(RUNGS)]
    stale = {w for n, w in
             {2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}.items()
             if n != len(RUNGS)}

    cb = _flat(CLAIM_BOUNDARY)
    assert f"the {word} rung verdicts are" in cb
    # the claim boundary states the range bound without a numeral
    # ("between or beyond the tested rungs"), which is count-agnostic
    # by construction -- pin that it stays that way
    assert "between or beyond the tested rungs" in cb
    for w in stale:
        assert f"the {w} rung verdicts" not in cb, w

    man = _flat(_PAPER / "artifact_manifest.md")
    assert f"carries the {word} per-rung count" in man
    for w in stale:
        assert f"carries the {w} per-rung count" not in man, w

    ms = _flat(MANUSCRIPT)
    assert f"{word.capitalize()} rungs were executed" in ms
    assert f"on each of the {word} preregistered rungs" in ms
    assert f"beyond the {word} tested values" in ms
    for w in stale:
        assert f"{w.capitalize()} rungs were executed" not in ms, w
        assert f"beyond the {w} tested values" not in ms, w


def test_the_promotion_is_bounded_in_both_documents():
    """Executed must not grow into more than was run."""

    cb = _flat(CLAIM_BOUNDARY)
    assert "NO joint or composite cross-rung verdict" in cb
    assert "no interpolation, no extrapolation" in cb
    assert "NOT a confirmation of causal-set theory" in cb
    ms = _flat(MANUSCRIPT)
    assert "no joint primary verdict" in ms
    assert "there is no joint primary verdict" in ms.lower()
