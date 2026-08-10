"""S4 prereg contract tests: the frozen constants, the gate logic
(including the falsifiability battery and strict boundary rule), the
seed discipline, and the power-certification table. The doc is
docs/prereg/p14_s4_schwarzschild_c1.md; anything asserted here is
frozen there."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

EXPERIMENT_DIR = (Path(__file__).resolve().parents[1]
                  / "experiments" / "positive_control")
sys.path.insert(0, str(EXPERIMENT_DIR))

import probe_seed_ledger as ledger  # noqa: E402
import s4_schwarzschild_c1 as s4  # noqa: E402


def test_s4_stream_is_observed_and_reruns_are_replays():
    """The campaign artifact exists, so 40_000_241 is observed-spent:
    the ledger hands it out only through the replay path, and a fresh
    allocation of it must abort."""

    assert ledger.replay_scalar("s4_campaign") == s4.SEED
    assert s4.SEED == 40_000_241
    assert s4.SEED in ledger.spent_scalars()
    assert s4.SMOKE_SEED == 40_000_221
    assert s4.SMOKE_SEED in ledger.spent_scalars()
    from seed_windows import assert_point_seeds_fresh
    with pytest.raises(SystemExit):
        assert_point_seeds_fresh({"s4_campaign": s4.SEED},
                                 ledger.spent_scalars(),
                                 ledger.SPENT_RANGES, "S4")


def test_replay_write_ownership_cannot_touch_the_fresh_artifact():
    """PR #58 review: the write boundary. In the current ledger state
    every non-smoke run is a replay, replay output is owned by the
    replay path regardless of what exists, and a 'fresh' run while
    the fresh artifact exists refuses instead of overwriting the
    original lineage."""

    assert "s4_campaign" not in ledger.FRESH_PROBE_SCALARS
    assert s4.artifact_policy("replay", True) == s4._REPLAY_ARTIFACT
    assert s4.artifact_policy("replay", False) == s4._REPLAY_ARTIFACT
    assert s4.artifact_policy("fresh_observation", False) == s4._ARTIFACT
    with pytest.raises(SystemExit):
        s4.artifact_policy("fresh_observation", True)


def test_executed_freeze_snapshot_matches_the_historical_blobs():
    """Paper A integration review: the CURRENT manifest hashes the
    post-result replay surface, so the manifest that governed the
    EXECUTED campaign is preserved as an immutable snapshot and
    verified here directly against the historical blobs at the freeze
    commit -- every digest it froze must equal the blob at ceed85d,
    and the executed artifact's lineage must be exactly that commit."""

    import hashlib
    import json
    import subprocess
    repo = Path(__file__).resolve().parents[1]
    snap = repo / "docs" / "prereg" / "p14_s4_executed_freeze_manifest.json"
    if not snap.exists():
        pytest.skip("executed-freeze snapshot not present")
    freeze = "ceed85d8bc0eb8e3cdddb071a2381defe1d3fe15"

    def blob(rel):
        return subprocess.run(
            ["git", "cat-file", "blob", f"{freeze}:{rel}"], cwd=repo,
            capture_output=True, check=True).stdout

    hist = blob("docs/prereg/p14_s4_freeze_manifest.json")
    assert snap.read_bytes().replace(b"\r\n", b"\n") == hist
    m = json.loads(hist.decode("utf-8"))
    assert len(m["files"]) == 8
    for rel, want in m["files"].items():
        assert hashlib.sha256(blob(rel)).hexdigest() == want, rel
    r = json.loads((repo / "docs" / "prereg" / "p14_s4_results.json")
                   .read_text(encoding="utf-8"))
    assert r["code"]["start"]["rev"] == freeze
    # the current replay-surface manifest is intentionally different
    cur = (repo / "docs" / "prereg" / "p14_s4_freeze_manifest.json")
    assert cur.read_bytes().replace(b"\r\n", b"\n") != hist


def test_campaign_artifact_pins_the_frozen_outcome():
    """The committed artifact: the observed seed, the exact freeze
    checkout at entry AND exit, run_kind = fresh_observation (a replay
    must not silently replace it), and the CONFIRMED verdict
    recomputes from the stored raw per-reading arrays through the
    frozen gate functions."""

    import json
    artifact = (Path(__file__).resolve().parents[1] / "docs" / "prereg"
                / "p14_s4_results.json")
    if not artifact.exists():
        pytest.skip("S4 results not present in this checkout")
    r = json.loads(artifact.read_text(encoding="utf-8"))
    assert r["params"]["seed"] == ledger.replay_scalar("s4_campaign")
    assert r["code"]["start"] == r["code"]["end"]
    assert r["code"]["start"]["dirty"] is False
    assert r["code"]["start"]["rev"] == \
        "ceed85d8bc0eb8e3cdddb071a2381defe1d3fe15"
    assert r["run_kind"] == "fresh_observation"
    lower = np.array(r["delta_lower"]["per_reading"])
    upper = np.array(r["delta_upper"]["per_reading"])
    ci = s4.identified_ci(lower, upper)
    wci = s4.welch_identified_ci(lower, upper, s4.s3_block())
    assert list(ci) == r["identified_ci95"]
    assert list(wci) == r["welch_identified_ci95"]
    verdict, b_label = s4.stage_verdict(ci, wci)
    assert verdict == r["verdict"] == "CONFIRMED"
    assert b_label == r["gate_b"] == "REPLICATED"
    assert r["gate_a"] is True and s4.gate_a(ci)
    assert r["ambiguity"]["ambiguous"] == 0


def test_frozen_constants_match_the_document():
    assert s4.N_READINGS == 300
    assert s4.E_N == 300
    assert s4.EPS_DET == 0.0036
    assert s4.EPS_REP == 0.0012
    assert s4.KAPPA == 3.0e-4
    assert s4.SD_CONS == 0.00129861497
    assert s4.POWER_SEED == (781, 4840)
    assert s4.POWER_B == 20_000


def test_sd_cons_is_the_chi2_upper_bound_of_the_s3_sd():
    """SD_cons = SD_S3 * sqrt(299 / chi2_{0.05,299}); the quantile is
    checked through the Wilson-Hilferty approximation, which is
    accurate to ~1e-3 relative at df=299."""

    d3 = s4.s3_block()
    sd3 = float(d3.std(ddof=1))
    assert sd3 == pytest.approx(0.0012108375, abs=1e-9)
    df = 299
    z05 = -1.6448536269514722
    wh = df * (1.0 - 2.0 / (9.0 * df) + z05 * math.sqrt(2.0 / (9.0 * df))) ** 3
    sd_cons_wh = sd3 * math.sqrt(df / wh)
    assert s4.SD_CONS == pytest.approx(sd_cons_wh, rel=2e-3)
    assert s4.SD_CONS > sd3


def test_identified_ci_reduces_to_student_t_when_unambiguous():
    rng = np.random.default_rng(7)
    x = rng.normal(-0.036, 0.0012, 300)
    lo, hi = s4.identified_ci(x, x)
    from p14_probe_p2 import student_t_crit
    t = student_t_crit(299)
    h = t * x.std(ddof=1) / math.sqrt(300)
    assert lo == pytest.approx(x.mean() - h, abs=1e-15)
    assert hi == pytest.approx(x.mean() + h, abs=1e-15)


def test_boundary_rule_is_strict():
    """The strict-comparison rule lives in the gate predicates, so it
    is pinned by feeding intervals that sit EXACTLY on each boundary:
    none may pass. (A mean() round-trip cannot construct an exact
    float boundary, which is precisely why the rule is strict.)"""

    assert not s4.gate_a((-0.037, -s4.EPS_DET))
    assert not s4.gate_null((-s4.EPS_DET, 0.0))
    assert not s4.gate_null((0.0, s4.EPS_DET))
    assert s4.gate_b((-0.0005, s4.EPS_REP)) != "REPLICATED"
    assert s4.gate_b((-s4.EPS_REP, 0.0005)) != "REPLICATED"
    assert s4.gate_b((s4.EPS_REP, 0.002)) != "DISCORDANT"
    assert s4.gate_b((s4.EPS_REP, 0.002)) == "B-INCONCLUSIVE"


def test_negative_controls_prove_the_gates_can_fail():
    """The frozen falsifiability battery (doc section 7): every
    negative construction fails its gate, and the A-only construction
    reaches DETECTED-NOT-REPLICATED."""

    nc = s4.negative_controls()
    assert nc["gate_a_at_threshold"] is False
    assert nc["gate_a_inside"] is False
    assert nc["gate_b_offband"] != "REPLICATED"
    assert nc["null_at_2eps"] is False
    assert nc["a_only_stage"] == "DETECTED-NOT-REPLICATED"
    assert nc["a_only_b_label"] != "REPLICATED"


def test_stage_partition_covers_every_case_exactly_once():
    """Constructed CI/WCI pairs hit each verdict; the partition of
    doc section 4 is total."""

    conf = s4.stage_verdict((-0.037, -0.035), (-0.0005, 0.0005))
    assert conf == ("CONFIRMED", "REPLICATED")
    dnr = s4.stage_verdict((-0.041, -0.039), (-0.005, -0.003))
    assert dnr[0] == "DETECTED-NOT-REPLICATED"
    nod = s4.stage_verdict((-0.001, 0.001), (0.033, 0.037))
    assert nod[0] == "NO-DETECTION"
    inc = s4.stage_verdict((-0.005, -0.001), (0.02, 0.04))
    assert inc[0] == "INCONCLUSIVE"


def test_replicated_forces_gate_a_at_these_margins():
    """The frozen inclusion: on real-scale blocks, whenever Gate B
    says REPLICATED, Gate A holds (doc section 4 derivation)."""

    s3d = s4.s3_block()
    src, m3 = s4._power_source()
    rng = np.random.default_rng(11)
    checked = 0
    for _ in range(200):
        x = src[rng.integers(0, len(src), 300)] + m3
        wci = s4.welch_identified_ci(x, x, s3d)
        if s4.gate_b(wci) == "REPLICATED":
            assert s4.gate_a(s4.identified_ci(x, x))
            checked += 1
    assert checked > 150


def test_freeze_manifest_matches_the_working_tree():
    """The content-addressed freeze identity: every protocol-surface
    file matches its frozen digest, the runner pins ITSELF, and the
    S3 comparison block is pinned -- any change to any of them must
    consciously regenerate the manifest under review."""

    import json
    m = json.loads(s4._FREEZE_MANIFEST.read_text(encoding="utf-8"))
    files = m["files"]
    assert len(files) == 8
    for rel in ("experiments/positive_control/s4_schwarzschild_c1.py",
                "docs/prereg/p14_s3_probe_results.json",
                "docs/prereg/p14_s4_schwarzschild_c1.md"):
        assert rel in files
    repo = Path(__file__).resolve().parents[1]
    for rel, want in files.items():
        assert s4._sha256(repo / rel) == want, rel
    s4.verify_freeze("test")


def test_freeze_verification_refuses_drift(monkeypatch, tmp_path):
    """A doctored digest must abort -- the clean-but-later-tree
    scenario from the PR #57 review."""

    import json
    m = json.loads(s4._FREEZE_MANIFEST.read_text(encoding="utf-8"))
    m["files"]["experiments/positive_control/s4_schwarzschild_c1.py"] = \
        "0" * 64
    doctored = tmp_path / "manifest.json"
    doctored.write_text(json.dumps(m), encoding="utf-8")
    monkeypatch.setattr(s4, "_FREEZE_MANIFEST", doctored)
    with pytest.raises(SystemExit):
        s4.verify_freeze("test")


def test_power_certification_fast_sanity():
    """B=500 per row must already sit above the 90% target."""

    for j in range(5):
        assert s4.power_row(j, b=500) >= 450


@pytest.mark.slow
def test_power_certification_full_table():
    """The frozen table: 20000/20000 on every row, CP95 lower bound
    0.9998156 >= 0.90."""

    for j in range(5):
        assert s4.power_row(j) == s4.POWER_B
    assert s4.cp_lower_all_success(s4.POWER_B) == pytest.approx(
        0.9998156, abs=5e-8)
    assert s4.cp_lower_all_success(s4.POWER_B) >= 0.90
