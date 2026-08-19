"""O4b downstream-integration contract tests.

The recovered O4b verdict is now cited by live consumers — the Paper A
manuscript, the claim boundary, the artifact manifest, and the two oracle
theory notes. These tests hold three lines at once: the FORMAL VERDICT
(the frozen stage-verdict table applied to the published result yields
CONCORDANT, and the cited numbers are the artifact's numbers verbatim),
the CONSUMER TEXT (the pre-certification phrasings — "partially derived",
"remains open there", "not-yet-designed" — must not return, while the
still-open Poisson-count boundary must stay stated), and the HISTORY
(the O4 abort and the O4b incident keep their original no-verdict
wording; the recovery is linked, never retro-graded).
"""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
_PREREG = REPO / "docs" / "prereg"
_PAPER = REPO / "docs" / "paper" / "paper_a"
_THEORY = REPO / "docs" / "theory"

RESULTS = _PREREG / "p14_o4b_results.json"
MANUSCRIPT = _PAPER / "manuscript.md"
CLAIM_BOUNDARY = _PAPER / "claim_boundary.md"
ARTIFACT_MANIFEST = _PAPER / "artifact_manifest.md"
ORACLE_NOTE = _THEORY / "schwarzschild_volume_oracle_note.md"
ORACLE_CERT = _THEORY / "schwarzschild_volume_oracle_certification.md"

#: The frozen figures the ruling records. Compared as exact reprs, so a
#: consumer citing a rounded or drifted number fails here.
G1_LO = "56.448806185841875"
G1_HI = "56.9822829864225"
DISC_LO = "-0.8992123405928396"
DISC_HI = "0.7695461186219887"
BAND = "1.7034113309135284"
G2_UPPER = "0.14195058753928652"
G2_BUDGET = "0.14195094424279403"


def _results() -> dict:
    return json.loads(RESULTS.read_text(encoding="utf-8"))


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


# ------------------------------------------------- the formal verdict

def test_the_frozen_verdict_table_yields_concordant():
    """The frozen O4 stage-verdict table (p14_o4_volume_audit.md §5):
    G1 and G2 both concordant, G3 valid, frozen n reached -> CONCORDANT.
    The published artifact sits exactly on that row."""

    r = _results()
    assert r["g1"]["status"] == "concordant"
    assert r["g2"]["status"] == "concordant"
    assert r["g1"]["n"] == 26_200_000          # the frozen n, reached
    assert r["g2"]["n"] == 1_072_696
    # G3a passed and the availability scan completed -- G3 not invalid
    assert r["availability"]["complete"] is True
    table = _read(_PREREG / "p14_o4_volume_audit.md")
    assert "| G1 and G2 both concordant | `CONCORDANT` |" in table


def test_the_result_is_the_canonical_artifact_with_the_frozen_numbers():
    r = _results()
    assert r["kind"] == "results"
    assert r["run_kind"] == "recovered_completed_campaign"
    assert repr(r["g1"]["v_s1_lo"]) == G1_LO
    assert repr(r["g1"]["v_s1_hi"]) == G1_HI
    assert [repr(x) for x in r["g1"]["identified_discrepancy"]] == [
        DISC_LO, DISC_HI]
    assert repr(r["g1"]["band_abs"]) == BAND
    assert repr(r["g2"]["leak_upper_abs"]) == G2_UPPER
    assert repr(r["g2"]["budget_abs"]) == G2_BUDGET
    assert r["g2"]["leaking_points"] == 0


def test_the_recovery_flags_are_explicit():
    """No new seed / no solver / no resampling / no gate change -- the
    ruling's recovery contract, carried in the artifact itself."""

    rc = _results()["recovery"]
    assert rc["no_new_seed"] is True
    assert rc["no_solver_call"] is True
    assert rc["no_resampling"] is True
    assert rc["no_gate_change"] is True
    assert rc["reproduced_preserved_verdicts_bit_exact"] is True


# --------------------------------------- the consumers cite the artifact

def test_the_manuscript_cites_the_artifact_numbers_verbatim():
    text = _read(MANUSCRIPT)
    for token in (G1_LO, G1_HI, DISC_LO, DISC_HI, BAND,
                  G2_UPPER, G2_BUDGET):
        assert token in text, token
    assert "CONCORDANT" in text
    assert "recovered_completed_campaign" in text
    assert "3.25%" in text


def test_the_claim_boundary_cites_the_artifact_numbers_verbatim():
    text = _read(CLAIM_BOUNDARY)
    for token in (G1_LO, G1_HI, DISC_LO, DISC_HI, BAND,
                  G2_UPPER, G2_BUDGET):
        assert token in text, token
    assert "CONCORDANT" in text
    assert "no new seed, no solver call" in text


def test_the_claim_boundary_states_the_forbidden_promotions():
    """The audit must not silently grow: each forbidden reading is
    named as a non-claim. Whitespace-normalized, so markdown line
    wrapping cannot hide a needle."""

    flat = " ".join(_read(CLAIM_BOUNDARY).split())
    assert "NOT a prediction-anchored promotion of S4/S5" in flat
    assert "NOT a C1/C2 joint verdict" in flat
    assert "NOT a Poisson causal-set count verification" in flat
    assert "NOT mass- or domain-generality" in flat
    assert "NOT a confirmation of causal-set theory" in flat
    assert "NOT complete separation or general volume accuracy" in flat
    assert "Nothing in Sections 6-6.7 is upgraded" in flat


# ------------------------------------------- stale phrasing stays out

def test_the_pre_certification_phrasings_do_not_return():
    """The corrected consumers may not regress to the pre-O3 state.
    'partially derived' survives only as the note's historical mention
    of the phrasing it replaced."""

    for path in (MANUSCRIPT, CLAIM_BOUNDARY):
        text = _read(path)
        assert "partially derived" not in text, path.name
        assert "remains open there" not in text, path.name
    note = " ".join(_read(ORACLE_NOTE).split())
    assert "not-yet-designed" not in note
    # the one permitted survival: the note explains the OLD phrasing
    assert '"partially derived theoretical path" phrasing described' in note
    cert = " ".join(_read(ORACLE_CERT).split())
    assert "Paper A is unchanged by this arc until" not in cert


def test_the_poisson_count_stage_is_stated_as_executed_per_rung():
    """The count stage is no longer open -- it ran across the mass
    ladder and is claimed per rung. This guard therefore flips
    direction: it now holds the promotion DOWN, so that "executed"
    never grows into a joint verdict, an interpolated statement, or a
    claim about causal-set theory."""

    for path, needle in (
            (MANUSCRIPT, "Poisson-count stage"),
            (CLAIM_BOUNDARY, "Prediction-anchored Poisson count, per rung"),
            (ORACLE_NOTE, "a separate, since-executed"),
            (ORACLE_CERT, "executed across the S6 mass ladder")):
        assert needle in " ".join(_read(path).split()), path.name

    flat = " ".join(_read(CLAIM_BOUNDARY).split())
    assert "NO joint or composite cross-rung verdict" in flat
    assert "no interpolation, no extrapolation" in flat
    assert "NOT a confirmation of causal-set theory" in flat
    assert "one fixed intensity per rung" in flat
    # and the auxiliary audit is still only an instrument statement
    assert "NOT a Poisson causal-set count verification" in flat


def test_s4_s5_verdicts_do_not_use_the_oracle():
    """The sentence the ruling keeps: the Section 6.7 verdicts are
    operationally anchored and consumed no oracle."""

    assert "uses no diamond-volume oracle" \
        in " ".join(_read(MANUSCRIPT).split())
    assert "no diamond-volume oracle is used in the S4/S5 verdicts" \
        in " ".join(_read(CLAIM_BOUNDARY).split())


# --------------------------------------------------- history preserved

def test_the_o4_abort_record_is_not_retro_graded():
    """The O4 campaign was never recovered -- its statistics were not
    persisted -- and its record keeps saying so."""

    rec = json.loads((_PREREG / "p14_o4_incident.json")
                     .read_text(encoding="utf-8"))
    assert rec["verdict"] is None
    assert rec["termination_reason"] == "g3-undecided"
    doc = _read(_PREREG / "p14_o4_incident.md")
    assert "no scientific verdict" in doc.lower() \
        or "no verdict" in doc.lower()
    assert "recovered in `p14_o4b_results.json`" not in doc


def test_the_o4b_incident_keeps_its_verdict_and_links_the_recovery():
    """The historical record stands -- verdict null, ABORT -- and the
    recovery is linked as a later event, not a regrade."""

    rec = json.loads((_PREREG / "p14_o4b_incident.json")
                     .read_text(encoding="utf-8"))
    assert rec["verdict"] is None
    assert rec["outcome"] == "ABORT"
    doc = _read(_PREREG / "p14_o4b_incident.md")
    flat = " ".join(doc.split())
    assert "This record grades nothing." in doc
    assert "subsequently recovered in `p14_o4b_results.json`" in flat
    assert "Nothing below > is retroactively regraded" in flat


def test_the_manifest_carries_the_o4b_bundle_and_the_inventory_rows():
    text = _read(ARTIFACT_MANIFEST)
    assert "Auxiliary O4b instrument-audit bundle" in text
    for rel in ("p14_o3_volume.json",
                "p14_o4b_executed_freeze_manifest.json",
                "p14_o4b_incident.json", "p14_o4b_checkpoint.json",
                "p14_o4b_results.json", "o4b_recover.py",
                "test_o4b_recover.py", "test_o4b_wiring_fixes.py"):
        assert rel in text, rel
    assert "No campaign rerun stands behind any row" in text
