"""S5 prereg contract tests: frozen constants, DeLong-vs-brute
equality, exact CP bounds, the four-way outcome partition with strict
boundaries, the BA train/test rules, the falsifiability battery, the
power certification table, and the seed discipline."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (Path(__file__).resolve().parents[1]
                  / "experiments" / "positive_control")
sys.path.insert(0, str(EXPERIMENT_DIR))

import probe_seed_ledger as ledger  # noqa: E402
import s5_schwarzschild_c2 as s5  # noqa: E402

_S3 = (Path(__file__).resolve().parents[1] / "docs" / "prereg"
       / "p14_s3_probe_results.json")


def _s3_arms():
    r = json.loads(_S3.read_text(encoding="utf-8"))
    return (np.array(r["f_schwarzschild_lower"]["per_reading"]),
            np.array(r["f_flat"]["per_reading"]))


def test_s5_seeds_are_fresh_distinct_allocations():
    assert ledger.assert_fresh_scalar("s5_curved") == s5.SEED_CURVED
    assert ledger.assert_fresh_scalar("s5_flat") == s5.SEED_FLAT
    assert s5.SEED_CURVED == 40_000_251
    assert s5.SEED_FLAT == 40_000_261
    assert s5.SEED_CURVED != s5.SEED_FLAT
    assert s5.SEED_CURVED not in ledger.spent_scalars()
    assert s5.SEED_FLAT not in ledger.spent_scalars()


def test_frozen_constants_match_the_document():
    assert s5.EPS_AUC == 0.10
    assert s5.EPS_BA == 0.10
    assert s5.N_ARM == 300
    assert s5.E_N == 300
    assert s5.CANDIDATE_N == (300, 600, 1000)
    assert s5.POWER_TARGET == 0.90
    assert s5.POWER_B == 20_000
    assert s5.POWER_SEED_BASE == (781, 4860)


def test_delong_equals_brute_force_auc():
    c0, f0 = _s3_arms()
    rng = np.random.default_rng(21)
    for _ in range(5):
        c = c0[rng.integers(0, 300, 120)]
        f = f0[rng.integers(0, 300, 150)]
        auc, lo, hi = s5.delong_auc_ci(c, f)
        assert auc == pytest.approx(s5.auc_brute(c, f), abs=1e-12)
        assert lo < auc < hi


def test_cp_table_exact_endpoints_and_monotonicity():
    t = s5._cp_table(10, 0.0125)
    assert t[0][0] == 0.0
    assert t[10][1] == 1.0
    assert t[10][0] == pytest.approx(0.0125 ** (1 / 10), abs=1e-9)
    for k in range(10):
        assert t[k][0] <= t[k + 1][0]
        assert t[k][1] <= t[k + 1][1]
        assert t[k][0] < t[k][1]


def test_outcome_partition_and_strict_boundaries():
    assert s5.auc_outcome((0.61, 0.99)) == "DETECTED"
    assert s5.auc_outcome((0.45, 0.55)) == "EQUIVALENT-AT-MARGIN"
    assert s5.auc_outcome((0.10, 0.35)) == "DIRECTION-REVERSED"
    assert s5.auc_outcome((0.55, 0.65)) == "INCONCLUSIVE"
    # exact boundaries never pass a branch
    assert s5.auc_outcome((0.60, 0.99)) == "INCONCLUSIVE"
    assert s5.auc_outcome((0.40, 0.55)) == "INCONCLUSIVE"
    assert s5.auc_outcome((0.45, 0.60)) == "INCONCLUSIVE"
    assert s5.auc_outcome((0.10, 0.40)) == "INCONCLUSIVE"


def test_identified_agreement_rule():
    """Disagreeing bound series resolve to INCONCLUSIVE, agreeing
    series pass through."""

    c0, f0 = _s3_arms()
    out, _ = s5.stage_outcome(c0, c0, f0)
    assert out == "DETECTED"
    shifted = c0 + 0.05  # upper series pushed toward the flat arm
    out2, det = s5.stage_outcome(c0, shifted, f0)
    assert det["outcomes"][0] == "DETECTED"
    assert det["outcomes"][1] != "DETECTED"
    assert out2 == "INCONCLUSIVE"


def test_ba_threshold_tie_goes_to_the_flat_side():
    thr = 0.5
    curved_test = np.array([0.5, 0.49, 0.51])
    flat_test = np.array([0.5, 0.51, 0.49])
    ba, lo, hi = s5.ba_cp_bonferroni(curved_test, flat_test, thr)
    # curved: only 0.49 is < thr (tie NOT curved); flat: 0.5 and 0.51
    assert ba == pytest.approx(0.5 * (1 / 3 + 2 / 3), abs=1e-12)
    assert lo < ba < hi


def test_negative_controls_prove_the_gates_can_fail():
    nc = s5.negative_controls()
    assert nc["detected_on_null"] is False
    assert nc["equivalent_on_effect"] is False
    assert nc["effect_outcome"] == "DETECTED"
    assert nc["reversed_outcome"] == "DIRECTION-REVERSED"
    assert nc["ba_gate_on_null"] is False


def test_power_certification_fast_sanity():
    c0, f0 = _s3_arms()
    for row in ("auc_effect", "auc_null", "ba_effect"):
        k = s5.power_row(row, 300, c0, f0, b=400)
        assert k >= 360, (row, k)


@pytest.mark.slow
def test_power_certification_full_table_and_selection():
    """The frozen certification: n = 300 rows exactly as pinned in
    the prereg doc, and the minimum-n rule selects 300."""

    c0, f0 = _s3_arms()
    chosen, table = s5.select_n(c0, f0)
    assert chosen == 300
    rows = table[300]
    assert rows["auc_effect"]["successes"] == 20_000
    assert rows["auc_null"]["successes"] == 19_564
    assert rows["ba_effect"]["successes"] == 20_000
    assert rows["auc_null"]["cp95_lower"] == pytest.approx(0.976080,
                                                           abs=5e-6)
    for r in rows.values():
        assert r["cp95_lower"] >= s5.POWER_TARGET


def test_sd_upper_factor_is_the_exact_chi2_value():
    """PR #60 execution review: the frozen doc says the EXACT
    chi-square one-sided 95% upper factor; the Wilson-Hilferty
    approximation (relative error -5.65e-7 at df=299, optimistic)
    is rejected. Pinned against the review's independent constant,
    and the gamma implementation is verified through two closed-form
    identities."""

    assert s5._sd_upper_factor(299) == pytest.approx(
        1.0724931824796697, abs=1e-12)
    assert s5._sd_upper_factor(299) > 1.0
    for x in (0.1, 0.5, 1.7, 4.0):
        assert s5._gammp(0.5, x) == pytest.approx(
            math.erf(math.sqrt(x)), abs=1e-13)
        assert s5._gammp(1.0, x) == pytest.approx(
            1.0 - math.exp(-x), abs=1e-13)


def test_frozen_sentences_are_verbatim_in_the_document():
    """PR #60 review: the runner SENTENCES dict is the canonical
    source and every sentence appears byte-identically (unwrapped) in
    the prereg document's section 8."""

    doc = (Path(__file__).resolve().parents[1] / "docs" / "prereg"
           / "p14_s5_schwarzschild_c2.md").read_text(encoding="utf-8")
    for key, sentence in s5.SENTENCES.items():
        assert sentence in doc, key


def test_cli_is_fail_closed(monkeypatch):
    """PR #60 review: an unknown argument or --help must exit at the
    CLI boundary and can never fall through into the fresh-seed
    campaign (argparse rejects unknowns with exit code 2 and handles
    --help with exit code 0, both before any seed is touched)."""

    monkeypatch.setattr(sys, "argv", ["s5_schwarzschild_c2", "--smok"])
    with pytest.raises(SystemExit) as bad:
        s5.main()
    assert bad.value.code == 2
    monkeypatch.setattr(sys, "argv", ["s5_schwarzschild_c2", "--help"])
    with pytest.raises(SystemExit) as help_exit:
        s5.main()
    assert help_exit.value.code == 0


def test_freeze_manifest_matches_the_working_tree():
    m = json.loads(s5._FREEZE_MANIFEST.read_text(encoding="utf-8"))
    files = m["files"]
    assert len(files) == 8
    assert "experiments/positive_control/s5_schwarzschild_c2.py" in files
    assert "docs/prereg/p14_s3_probe_results.json" in files
    repo = Path(__file__).resolve().parents[1]
    for rel, want in files.items():
        assert s5._sha256(repo / rel) == want, rel
    s5.verify_freeze("test")
