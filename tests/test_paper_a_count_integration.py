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
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
_PREREG = REPO / "docs" / "prereg"
_PAPER = REPO / "docs" / "paper" / "paper_a"

MANUSCRIPT = _PAPER / "manuscript.md"
CLAIM_BOUNDARY = _PAPER / "claim_boundary.md"

#: (mu as printed, artifact) for the three executed rungs
RUNGS = (
    ("0.1333", "p14_o5_count.json"),
    ("0.1867", "p14_s6_m14_count.json"),
    ("0.2400", "p14_s6_m18_count.json"),
)


def _art(name: str) -> dict:
    return json.loads((_PREREG / name).read_text(encoding="utf-8"))


def _flat(path: Path) -> str:
    return " ".join(path.read_text(encoding="utf-8").split())


def test_every_cited_rung_actually_returned_concordant():
    for _, name in RUNGS:
        a = _art(name)
        assert a["decision"]["verdict"] == "CONCORDANT", name


def test_the_manuscript_prints_the_artifact_numbers_not_typed_ones():
    """K, U, D and the band in the Section 6.8 table are re-derived from
    each artifact and must appear verbatim."""

    flat = _flat(MANUSCRIPT)
    for mu, name in RUNGS:
        a = _art(name)
        d = a["decision"]
        assert mu in flat, mu
        assert f"{a['scan']['k_certain']:,}" in flat, name
        assert f"[{d['d_lo']:.6f}, {d['d_hi']:+.6f}]" in flat, name
        assert f"{d['band']:.6f}" in flat, name


def test_the_claim_boundary_prints_the_artifact_numbers():
    flat = _flat(CLAIM_BOUNDARY)
    for mu, name in RUNGS:
        a = _art(name)
        d = a["decision"]
        assert f"`mu = {mu}` K = {a['scan']['k_certain']:,}" in flat, name
        assert (f"D = [{d['d_lo']:.6f}, {d['d_hi']:+.6f}] "
                f"in +-{d['band']:.6f}") in flat, name
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
    assert mus == sorted(mus) and len(set(mus)) == 3
    flat = _flat(MANUSCRIPT)
    assert "no rung is an isometric" in flat
    assert "mu = 2M/r_c" in flat


def test_the_promotion_is_bounded_in_both_documents():
    """Executed must not grow into more than was run."""

    cb = _flat(CLAIM_BOUNDARY)
    assert "NO joint or composite cross-rung verdict" in cb
    assert "no interpolation, no extrapolation" in cb
    assert "NOT a confirmation of causal-set theory" in cb
    ms = _flat(MANUSCRIPT)
    assert "no joint primary verdict" in ms
    assert "there is no joint primary verdict" in ms.lower()
