"""O4 freeze contract tests.

The freeze-commit discipline, as for O3: the frozen gate, the
content-addressed manifest, the write-once/no-results boundary and the
fail-closed CLI are pinned BEFORE execution. What is new here is that
the STATISTICAL machinery carries contract tests of its own -- the
review ruling is explicit that coverage is bought by Maurer-Pontil
Theorem 4 itself, so what must be pinned is the identity of this
implementation with the theorem's statement, never a simulation.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "experiments" / "oracle"))
sys.path.insert(0, str(_REPO / "experiments" / "positive_control"))

import empirical_bernstein as eb  # noqa: E402
import exact_binomial as xb  # noqa: E402
import o4_sizing as sz  # noqa: E402
import o4_volume_audit as o4  # noqa: E402

# ------------------------------------------------------- freeze shell

def test_manifest_digests_match_the_working_tree():
    """The pinned set covers the certified surface, the statistics
    modules, the sizing certification, the runner, the preregistration,
    the INSTRUMENT UNDER AUDIT, and the O3 result the audit is measured
    against."""

    manifest = o4.verify_digests("test")
    assert set(manifest["files"]) == {
        "experiments/oracle/certified_interval.py",
        "experiments/oracle/certified_flight_time.py",
        "experiments/oracle/empirical_bernstein.py",
        "experiments/oracle/exact_binomial.py",
        "experiments/oracle/o4_sizing.py",
        "experiments/oracle/o4_volume_audit.py",
        "experiments/positive_control/s1_schwarzschild_cost.py",
        "experiments/positive_control/probe_seed_ledger.py",
        "docs/prereg/p14_o4_volume_audit.md",
        "docs/prereg/p14_o4_sizing.json",
        "docs/prereg/p14_o3_volume.json",
        "pyproject.toml",
    }


def test_environment_lock_covers_the_whole_apparatus(monkeypatch):
    manifest = json.loads(o4._MANIFEST.read_text(encoding="utf-8"))
    env = manifest["environment"]
    assert set(env) == {"python", "gmpy2", "mpfr", "gmp", "numpy"}
    monkeypatch.setattr(o4, "_environment", lambda: dict(env))
    o4.verify_environment("test", manifest)
    monkeypatch.setattr(o4, "_environment",
                        lambda: dict(env, numpy="0.0.0"))
    with pytest.raises(SystemExit, match="environment drift"):
        o4.verify_environment("test", manifest)


_REV = "0" * 40

#: The freeze merge the approved campaign ran from.
_FREEZE_REV = "1eb9461f1af739968403c633a7a44cbba9a9f948"


@pytest.fixture
def clean_outputs(monkeypatch, tmp_path):
    """The campaign HAS run, so `p14_o4_reservation.json` now exists in
    the tree and preflight's write-once check fires before anything
    else. Tests of the other checks redirect the two write-once paths
    to an empty directory; the live files are asserted on separately,
    in the incident record's own tests."""

    fakes = (tmp_path / "p14_o4_reservation.json",
             tmp_path / "p14_o4_results.json")
    monkeypatch.setattr(o4, "_RESERVATION", fakes[0])
    monkeypatch.setattr(o4, "_ARTIFACT", fakes[1])
    monkeypatch.setattr(o4, "_WRITE_ONCE", fakes)
    return fakes


@pytest.fixture
def offline(monkeypatch, clean_outputs):
    """Preflight's later steps talk to the real reservation authority
    -- including a push. Tests that are about the EARLIER checks must
    not reach it: CI caught `test_preflight_refuses_any_commit...`
    pushing a probe ref to the live repository."""

    def forbidden(*a, **k):                          # pragma: no cover
        raise AssertionError("this test must not touch the remote")

    monkeypatch.setattr(o4, "remote_reservation", lambda: None)
    monkeypatch.setattr(o4, "probe_reservation_namespace",
                        lambda: None)
    monkeypatch.setattr(o4, "reserve_remote", forbidden)


def test_the_freeze_commit_carried_neither_result_nor_reservation():
    """Historical, as for O3: the campaign has since run and left a
    reservation in the tree, so only the frozen commit can attest that
    the rules were sealed with nothing observed. The result artifact
    never appeared at all -- the run aborted in G3."""

    for rel in ("p14_o4_results.json", "p14_o4_reservation.json"):
        probe = subprocess.run(
            ["git", "cat-file", "-e", f"{_FREEZE_REV}:docs/prereg/{rel}"],
            cwd=_REPO, capture_output=True)
        assert probe.returncode != 0, rel
    assert not o4._ARTIFACT.exists(), "no verdict was ever published"


@pytest.mark.parametrize("which", [0, 1])
def test_preflight_refuses_when_a_write_once_output_exists(
        monkeypatch, tmp_path, which, offline):
    """Both the reservation and the result are terminal: the audit is
    write-once, so a second attempt on the same streams is refused
    whether or not the first one produced a verdict."""

    names = ("p14_o4_reservation.json", "p14_o4_results.json")
    fakes = [tmp_path / n for n in names]
    fakes[which].write_text("{}", encoding="utf-8")
    monkeypatch.setattr(o4, "_RESERVATION", fakes[0])
    monkeypatch.setattr(o4, "_ARTIFACT", fakes[1])
    monkeypatch.setattr(o4, "_WRITE_ONCE", tuple(fakes))
    monkeypatch.setattr(o4, "verify_freeze", lambda stage: {})
    monkeypatch.setattr(o4, "_git_state",
                        lambda: {"rev": _REV, "dirty": False})
    with pytest.raises(SystemExit, match="write-once"):
        o4.preflight(_REV)


def test_preflight_refuses_a_dirty_tree(monkeypatch, offline):
    monkeypatch.setattr(o4, "verify_freeze", lambda stage: {})
    monkeypatch.setattr(o4, "_git_state",
                        lambda: {"rev": _REV, "dirty": True})
    with pytest.raises(SystemExit, match="clean exact checkout"):
        o4.preflight(_REV)


def test_preflight_refuses_any_commit_the_approval_does_not_name(
        monkeypatch, offline):
    """Digest verification cannot certify the manifest itself: a later
    commit that edits a protocol file and re-pins the manifest in the
    SAME commit passes every digest check. The approved SHA is the only
    thing that distinguishes it (review R1)."""

    monkeypatch.setattr(o4, "verify_freeze", lambda stage: {})
    monkeypatch.setattr(o4, "_git_state",
                        lambda: {"rev": "a" * 40, "dirty": False})
    with pytest.raises(SystemExit, match="the approval does not name"):
        o4.preflight(_REV)
    for junk in ("", "abc", "a" * 39, "z" * 40, "A" * 41):
        with pytest.raises(SystemExit, match="full 40-hex"):
            o4.preflight(junk)
    # case-insensitive and whitespace-tolerant, but otherwise exact: a
    # matching SHA now gets past the rev check and dies later, on the
    # retired streams (see the un-runnable test below)
    with pytest.raises(KeyError, match="never re-entered"):
        o4.preflight("  " + "A" * 40 + "\n")


def test_the_executed_freeze_can_no_longer_start_a_campaign():
    """The streams this freeze names were drawn and retired, so its own
    preflight must refuse -- there is no configuration of the tree in
    which `1eb9461` may run again. A re-run needs a new freeze with new
    scalars."""

    out = subprocess.run(
        [sys.executable, str(_REPO / "experiments" / "oracle"
                             / "o4_volume_audit.py"),
         "--preflight", "--freeze-rev", _FREEZE_REV],
        cwd=_REPO, capture_output=True, text=True)
    assert out.returncode != 0
    combined = out.stdout + out.stderr
    assert "never re-entered" in combined or "dirty" in combined


def test_preflight_does_consult_the_authority(monkeypatch,
                                             clean_outputs):
    """Guard against the `offline` fixture hiding a regression: a
    preflight that clears every local check must still have asked the
    authority and proved it can write there."""

    called = []
    monkeypatch.setattr(o4, "verify_freeze", lambda stage: {})
    monkeypatch.setattr(o4, "_git_state",
                        lambda: {"rev": _REV, "dirty": False})
    monkeypatch.setattr(o4, "remote_reservation",
                        lambda: called.append("read"))
    monkeypatch.setattr(o4, "probe_reservation_namespace",
                        lambda: called.append("probe"))
    # the seed check is last and now refuses (streams retired), so the
    # call still proves the authority was consulted on the way there
    with pytest.raises(KeyError, match="never re-entered"):
        o4.preflight(_REV)
    assert called == ["read", "probe"]


def test_campaign_refuses_to_run_without_the_approved_sha(monkeypatch):
    script = str(_REPO / "experiments" / "oracle"
                 / "o4_volume_audit.py")
    out = subprocess.run([sys.executable, script, "--preflight"],
                         capture_output=True, text=True)
    assert out.returncode != 0
    assert "--freeze-rev is required" in out.stderr


def test_the_campaigns_own_outputs_do_not_read_as_a_dirty_tree(
        monkeypatch):
    """The reservation is written BEFORE the run, so the exit lineage
    check must not mistake it for protocol drift -- while any other
    modified path still must."""

    lines = ["?? docs/prereg/p14_o4_reservation.json",
             "?? docs/prereg/p14_o4_results.json"]

    def fake_run(cmd, **kw):
        out = ("\n".join(lines) if "status" in cmd else "f" * 40)
        return subprocess.CompletedProcess(cmd, 0, out + "\n", "")

    monkeypatch.setattr(o4.subprocess, "run", fake_run)
    assert o4._git_state() == {"rev": "f" * 40, "dirty": False}
    lines.append(" M experiments/oracle/o4_sizing.py")
    assert o4._git_state()["dirty"] is True


def test_cli_is_fail_closed():
    script = str(_REPO / "experiments" / "oracle"
                 / "o4_volume_audit.py")
    bad = subprocess.run([sys.executable, script, "--preflght"],
                         capture_output=True, text=True)
    assert bad.returncode == 2
    helped = subprocess.run([sys.executable, script, "--help"],
                            capture_output=True, text=True)
    assert helped.returncode == 0
    assert "--preflight" in helped.stdout


def test_smoke_is_capped_against_pre_observation():
    """A large smoke run on the FROZEN anchors would pre-observe V_S1
    at a precision comparable to the equivalence band -- the O3 lesson
    about collapsing the freeze boundary."""

    script = str(_REPO / "experiments" / "oracle"
                 / "o4_volume_audit.py")
    out = subprocess.run(
        [sys.executable, script, "--smoke",
         str(o4.SMOKE_MAX + 1)], capture_output=True, text=True)
    assert out.returncode != 0
    assert "collapse the freeze boundary" in out.stderr


# -------------------------------------- Maurer-Pontil Theorem 4

def test_variance_is_the_theorems_pairwise_form():
    """V_n = sum_{i<j}(Z_i-Z_j)^2 / (n(n-1)) is what Theorem 4 names;
    the streaming accumulator must compute exactly that."""

    xs = [0.0, 0.25, 0.5, 0.125, 1.0, 0.75, 0.03125]
    acc = eb.Accumulator()
    for x in xs:
        acc.add(x)
    assert acc.var == pytest.approx(eb.pairwise_variance(xs),
                                    rel=1e-12)


def test_half_width_is_the_literal_theorem_formula():
    n, var, delta = 26_200_000, 3.5e-4, 0.02
    ln = math.log(2.0 / delta)
    want = math.sqrt(2.0 * var * ln / n) + 7.0 * ln / (3.0 * (n - 1))
    assert eb.half_width(n, var, delta) == pytest.approx(want,
                                                         rel=1e-15)


def test_out_of_range_aborts_rather_than_clipping():
    """Clipping would bias the mean and void the [0,1] hypothesis, so
    the precondition is fail-closed rather than repaired."""

    acc = eb.Accumulator()
    with pytest.raises(eb.EBError, match=r"outside \[0, 1\]"):
        acc.add(1.0000001)
    with pytest.raises(eb.EBError, match="non-finite"):
        acc.add(float("nan"))


def test_interval_needs_two_samples():
    acc = eb.Accumulator()
    acc.add(0.5)
    with pytest.raises(eb.EBError, match="n >= 2"):
        eb.interval(acc, 0.02)


def test_statistic_is_invariant_to_stream_partition():
    """The i.i.d. unit is the sample POINT: a partitioned run must
    produce the same mean and V_n as a single-stream run over the same
    points, so the number of generating streams can never enter the
    inference."""

    xs = [(i * 37 % 101) / 200.0 for i in range(500)]
    whole = eb.Accumulator()
    for x in xs:
        whole.add(x)
    parts = []
    for lo, hi in ((0, 111), (111, 260), (260, 500)):
        a = eb.Accumulator()
        for x in xs[lo:hi]:
            a.add(x)
        parts.append(a)
    merged = eb.Accumulator()
    for a in parts:
        merged.merge(a)
    assert merged.n == whole.n
    assert merged.mean == pytest.approx(whole.mean, rel=1e-12)
    assert merged.var == pytest.approx(whole.var, rel=1e-10)


def test_rescaling_preserves_the_interval_monotonically():
    acc = eb.Accumulator()
    for x in (0.1, 0.2, 0.3, 0.4):
        acc.add(x)
    iv = eb.interval(acc, 0.02)
    lo, hi = iv.rescaled(sz.SCALE)
    assert lo == pytest.approx(sz.SCALE * iv.lo, rel=1e-15)
    assert hi == pytest.approx(sz.SCALE * iv.hi, rel=1e-15)
    with pytest.raises(eb.EBError):
        iv.rescaled(0.0)


# ------------------------------------------- Clopper-Pearson

def test_zero_count_matches_the_closed_form():
    for n in (100, 50_000, 1_072_696):
        for alpha in (0.05, 0.005):
            assert xb.cp_upper(0, n, alpha) == pytest.approx(
                xb.cp_upper_zero(n, alpha), rel=1e-9)


def test_cp_upper_is_exact_at_the_defining_equation():
    p = xb.cp_upper(3, 1000, 0.005)
    assert xb.binom_cdf(3, 1000, p) == pytest.approx(0.005, rel=1e-6)


def test_g3_cluster_bound_is_the_frozen_number():
    assert sz.g3_upper() == pytest.approx(2.995687e-05, rel=1e-6)


# ---------------------------------------------- sizing certification

def test_l_max_is_an_upward_bound_never_a_nearest_rounding():
    """A downward L_max would understate the variance ceiling and
    inflate the power claim, so the frozen constant must dominate the
    certified margin -- 1.5599927 (nearest at 7 places) would be
    BELOW it and is rejected here by construction."""

    cert = sz._CERT["margins"]["nonempty"]
    assert sz.L_MAX_UB >= float(cert.hi)
    assert sz.L_MAX_UB >= 1.5599927415085288
    assert 1.5599927 < float(cert.lo)      # the trap this test closes


def test_geometry_constants_sit_inside_their_enclosures():
    for value, iv in ((sz.B_BOX, sz.B_BOX_IV),
                      (sz.B_OUT, sz.B_OUT_IV),
                      (sz.SCALE, sz.SCALE_IV)):
        assert iv.lo_float() <= value <= iv.hi_float()


def test_tau_clears_the_full_width_floor():
    """The floor is the oracle's FULL certified width, not its
    half-width: the identified set absorbs the whole interval because
    V_true's position inside it is unknown."""

    floor = (sz.V_HI - sz.V_LO) / sz.V_REF
    assert floor == pytest.approx(0.01999426, rel=1e-5)
    assert sz.TAU > floor
    assert 0.015 < floor                   # tau = 1.5% is impossible


def test_power_is_certified_on_the_continuum_not_a_grid():
    """Interval-arithmetic branch and bound over the whole segment
    [V_lo, V_hi]; a finite grid is not a certification."""

    cert = sz.certify_power(sz.N_G1)
    assert cert["power_lower_bound"] >= sz.POWER_TARGET
    assert cert["worst_endpoint"] == "V_hi"
    assert cert["boxes"] >= 1


def test_frozen_size_dominates_the_minimum_certified_size():
    assert sz.N_G1 >= 26_012_722
    assert sz.certify_power(26_012_722)["power_lower_bound"] \
        < sz.certify_power(sz.N_G1)["power_lower_bound"]


def test_exp_bound_is_an_upper_bound_on_the_true_exponential():
    """Validity everywhere; tightness where the branch and bound
    actually decides. Truncating exp(x)'s series always UNDER-states
    it, so the reciprocal always OVER-states exp(-x) -- the direction
    the proof needs. The residual looseness grows with x, but by then
    exp(-x) is far below the 0.0967 the bound must beat."""

    from certified_interval import Iv
    for x in (0.0, 0.5, 2.30258509, 12.0):
        got = sz.exp_neg_upper(Iv(x)).hi_float()
        assert got >= math.exp(-x)
        assert got <= math.exp(-x) * (1.0 + 1e-9)
    far = sz.exp_neg_upper(Iv(40.0)).hi_float()
    assert math.exp(-40.0) <= far < 1e-17


def test_sizing_artifact_matches_the_module():
    art = json.loads(
        (_REPO / "docs" / "prereg" / "p14_o4_sizing.json")
        .read_text(encoding="utf-8"))
    assert art["g1"]["n"] == sz.N_G1
    assert art["g2"]["points"] == sz.g2_points()
    assert art["gate"]["tau"] == sz.TAU
    assert art["gate"]["familywise"] == pytest.approx(0.05)
    assert art["gate"]["total_sentence"] == pytest.approx(0.0325)
    assert art["geometry"]["L_max_upper_bound"] == sz.L_MAX_UB


def test_cap_headroom_is_stated_as_two_distinct_ratios():
    """Review item 4: 80M calls is 1.46x the expectation while 24 h is
    2.05x the wall clock -- the wording must not merge them."""

    b = json.loads(
        (_REPO / "docs" / "prereg" / "p14_o4_sizing.json")
        .read_text(encoding="utf-8"))["budget"]
    assert b["cap_calls"] == o4.FROZEN["max_calls"]
    assert b["cap_wall_s"] == o4.FROZEN["max_wall_s"]
    assert 1.4 < b["cap_calls_ratio"] < 1.5
    assert 2.0 < b["cap_wall_ratio"] < 2.1
    assert abs(b["cap_calls_ratio"] - b["cap_wall_ratio"]) > 0.5


# ------------------------------------------------ verdict composition

def _g1(status, powered=True):
    return {"status": status, "power_certified": powered}


def test_verdict_composition_covers_every_branch():
    ok3 = {"valid": True}
    assert o4.compose_verdict(_g1("concordant"),
                              {"status": "concordant"},
                              ok3) == "CONCORDANT"
    assert o4.compose_verdict(_g1("discordant"),
                              {"status": "concordant"},
                              ok3) == "DISCORDANT"
    assert o4.compose_verdict(_g1("concordant"),
                              {"status": "discordant"},
                              ok3) == "DISCORDANT"
    assert o4.compose_verdict(_g1("inconclusive"),
                              {"status": "concordant"},
                              ok3) == "INCONCLUSIVE"
    assert o4.compose_verdict(_g1("concordant", powered=False),
                              {"status": "concordant"},
                              ok3) == "INCONCLUSIVE"


def test_g3_failure_is_invalid_never_scientific_discordance():
    """The wrapper contract is instrumentation: a broken one voids the
    stage rather than reporting a physics disagreement."""

    bad3 = {"valid": False}
    for g1 in ("concordant", "discordant", "inconclusive"):
        for g2 in ("concordant", "discordant"):
            assert o4.compose_verdict(_g1(g1), {"status": g2},
                                      bad3) == "INVALID"


def test_discordance_survives_an_underpowered_run():
    """Disjointness is a coverage statement, so it holds at any n;
    concordance is not, so it needs the frozen n."""

    assert o4.compose_verdict(_g1("discordant", powered=False),
                              {"status": "concordant"},
                              {"valid": True}) == "DISCORDANT"


# -------------------------------------------------- frozen agreement

def test_frozen_config_matches_the_ruled_design():
    assert o4.FROZEN["tau"] == 0.030
    assert o4.FROZEN["delta_g1_per_side"] == 0.02
    assert o4.FROZEN["alpha_g2_per_end"] == 0.005
    assert o4.FROZEN["alpha_g3"] == 0.05
    assert o4.FROZEN["leak_budget_frac"] == 0.0025
    assert o4.FROZEN["n_g1"] == 26_200_000
    assert o4.FROZEN["g3_clusters"] == 100_000
    assert o4.FROZEN["max_calls"] == 80_000_000
    assert o4.FROZEN["max_wall_s"] == 86_400.0


def test_the_freeze_commit_carried_two_fresh_campaign_scalars():
    """A historical assertion, as for O3's no-result check: the freeze
    commit is what must have held these as FRESH. They were drawn on
    2026-08-12 and are retired now, so the live ledger cannot state
    this -- only the frozen blob can."""

    blob = subprocess.run(
        ["git", "show", f"{_FREEZE_REV}:experiments/positive_control/"
         f"probe_seed_ledger.py"],
        cwd=_REPO, capture_output=True, text=True, check=True).stdout
    frozen = blob.split("FRESH_PROBE_SCALARS")[1].split("}")[0]
    assert '"g1_audit": 40_000_281' in frozen
    assert '"g2_leakage": 40_000_291' in frozen


def test_the_smoke_stream_and_its_successors_stay_retired():
    from probe_seed_ledger import OBSERVED_PROBE_SCALARS
    for k in ("o4_smoke", "o4_smoke_g2", "o4_smoke_g3"):
        assert k in OBSERVED_PROBE_SCALARS
    assert (OBSERVED_PROBE_SCALARS["o4_smoke"] + 1
            == OBSERVED_PROBE_SCALARS["o4_smoke_g2"])
    assert (OBSERVED_PROBE_SCALARS["o4_smoke"] + 2
            == OBSERVED_PROBE_SCALARS["o4_smoke_g3"])


def test_g3_declares_the_lineage_it_actually_has():
    """G3 reuses G1's accepted points, so it must not also hold a seed
    it never reads: that would record a spent scalar as provenance for
    samples G1 produced (review R1)."""

    import inspect
    assert list(inspect.signature(o4.run_g3).parameters) == [
        "cfg", "stress", "budget"]
    src = inspect.getsource(o4)
    assert 'seeds["g3"]' not in src and "'g3':" not in src
    assert '"g3": ' not in src.split("def run_g3")[0]

    budget = o4._Budget(o4.FROZEN)
    out = o4.run_g3(o4.FROZEN | {"g3_clusters": 0}, [], budget)
    assert "inherited from G1" in out["seed_lineage"]


def _fake_remote(tmp_path, monkeypatch):
    """A real bare repo standing in for `origin`, so the reservation's
    create-only semantics are exercised against git itself rather than
    against a mock of what git is assumed to do."""

    bare, work = tmp_path / "remote.git", tmp_path / "work"
    subprocess.run(["git", "init", "-q", "--bare", str(bare)],
                   check=True)
    subprocess.run(["git", "init", "-q", str(work)], check=True)
    for cmd in (["config", "user.email", "o4@test"],
                ["config", "user.name", "o4"],
                ["commit", "-q", "--allow-empty", "-m", "freeze"],
                ["remote", "add", "origin", str(bare)]):
        subprocess.run(["git", "-C", str(work), *cmd], check=True)
    monkeypatch.setattr(o4, "_REPO", work)
    url = subprocess.run(["git", "-C", str(work), "remote", "get-url",
                          "origin"], check=True, capture_output=True,
                         text=True).stdout.strip()
    monkeypatch.setattr(o4, "_CANONICAL_AUTHORITY",
                        o4._normalise_remote(url))
    return work


@pytest.mark.parametrize("url", [
    "https://github.com/lpaiu-cs/causal-spacetime.git",
    "https://github.com/lpaiu-cs/causal-spacetime",
    "https://lpaiu-cs@github.com/lpaiu-cs/causal-spacetime.git",
    "git@github.com:lpaiu-cs/causal-spacetime.git",
    "ssh://git@github.com/lpaiu-cs/causal-spacetime.git",
    "git://github.com/lpaiu-cs/causal-spacetime.git",
    "https://github.com/lpaiu-cs/causal-spacetime/",
])
def test_every_spelling_of_the_canonical_repository_is_accepted(url):
    assert o4._normalise_remote(url) == o4._CANONICAL_AUTHORITY


@pytest.mark.parametrize("url", [
    "https://github.com/someone-else/causal-spacetime.git",
    "git@github.com:lpaiu-cs/causal-spacetime-mirror.git",
    "https://gitlab.com/lpaiu-cs/causal-spacetime.git",
    "/srv/mirrors/causal-spacetime.git",
])
def test_a_fork_or_mirror_is_not_the_reservation_authority(url):
    assert o4._normalise_remote(url) != o4._CANONICAL_AUTHORITY


def test_the_authority_is_an_identity_not_a_local_alias(monkeypatch,
                                                        tmp_path):
    """`origin` is a per-checkout alias: a fork or mirror holding the
    same approved SHA would find ITS OWN reservation ref empty and
    claim the same streams elsewhere, which the SHA check cannot see
    (review R3)."""

    work = _fake_remote(tmp_path, monkeypatch)
    o4.reservation_authority()          # the patched local authority
    monkeypatch.setattr(o4, "_CANONICAL_AUTHORITY",
                        "github.com/lpaiu-cs/causal-spacetime")
    with pytest.raises(SystemExit, match="not github.com/lpaiu-cs"):
        o4.reservation_authority()
    # and no path to the ref may skip the identity check
    with pytest.raises(SystemExit, match="not github.com/lpaiu-cs"):
        o4.remote_reservation()
    with pytest.raises(SystemExit, match="not github.com/lpaiu-cs"):
        o4.reserve_remote({"attempt": 1})
    with pytest.raises(SystemExit, match="not github.com/lpaiu-cs"):
        o4.probe_reservation_namespace()

    subprocess.run(["git", "-C", str(work), "remote", "remove",
                    "origin"], check=True)
    with pytest.raises(SystemExit, match="no `origin` remote"):
        o4.reservation_authority()


def test_reservation_is_global_not_local_to_a_checkout(monkeypatch,
                                                       tmp_path):
    """A working-tree file serialises only processes sharing that
    directory: two clones of the same approved SHA would each make
    their own reservation, and an aborted run whose worktree is
    discarded would take its reservation with it (review R2). The
    claim therefore lives on a remote ref."""

    _fake_remote(tmp_path, monkeypatch)
    assert o4.remote_reservation() is None
    first = o4.reserve_remote({"attempt": 1})
    assert o4.remote_reservation() == first

    # a second checkout -- fresh tree, same approved commit, no local
    # reservation file anywhere -- must still be refused
    with pytest.raises(SystemExit, match="already held"):
        o4.reserve_remote({"attempt": 2})
    # ... and refused even when it pushes the very same payload, which
    # a plain `git push` would have accepted as "everything up-to-date"
    with pytest.raises(SystemExit, match="already held"):
        o4.reserve_remote({"attempt": 1})
    assert o4.remote_reservation() == first


def test_an_unreachable_authority_refuses_rather_than_assumes_free(
        monkeypatch, tmp_path):
    _fake_remote(tmp_path, monkeypatch)
    (tmp_path / "remote.git").rename(tmp_path / "moved.git")
    with pytest.raises(SystemExit, match="reservation authority"):
        o4.remote_reservation()


def test_preflight_refuses_when_the_streams_are_already_claimed(
        monkeypatch, clean_outputs):
    monkeypatch.setattr(o4, "verify_freeze", lambda stage: {})
    monkeypatch.setattr(o4, "_git_state",
                        lambda: {"rev": _REV, "dirty": False})
    monkeypatch.setattr(o4, "remote_reservation", lambda: "d" * 40)
    with pytest.raises(SystemExit, match="already held"):
        o4.preflight(_REV)


def test_preflight_proves_the_claim_can_actually_be_made(monkeypatch,
                                                         tmp_path):
    """Push rights and the server's ref-namespace policy are only
    exercised by a real push; without the probe, a campaign could
    clear every check and then fail to claim its streams."""

    work = _fake_remote(tmp_path, monkeypatch)
    o4.probe_reservation_namespace()
    left = subprocess.run(["git", "-C", str(work), "ls-remote",
                           "origin"], capture_output=True, text=True)
    assert "probe" not in left.stdout, "the probe ref was not cleaned"

    (tmp_path / "remote.git").rename(tmp_path / "moved.git")
    with pytest.raises(SystemExit, match="must not start"):
        o4.probe_reservation_namespace()


def test_streams_are_reserved_before_the_first_draw(monkeypatch,
                                                    tmp_path):
    """An aborted campaign must still spend its seeds. Without a
    persistent reservation, a fail-closed abort leaves no artifact and
    the next attempt would read the already-observed streams as
    fresh."""

    reservation = tmp_path / "p14_o4_reservation.json"
    monkeypatch.setattr(o4, "_RESERVATION", reservation)
    monkeypatch.setattr(o4, "_ARTIFACT", tmp_path / "p14_o4_results.json")
    monkeypatch.setattr(o4, "preflight", lambda rev: {"git": {}})
    monkeypatch.setattr(o4, "_git_state",
                        lambda: {"rev": _REV, "dirty": False})
    monkeypatch.setattr(o4, "assert_fresh_scalar", lambda name: 7)
    claimed: list[dict] = []
    monkeypatch.setattr(o4, "reserve_remote",
                        lambda payload: claimed.append(payload) or "c" * 40)

    def abort(*a, **k):
        raise SystemExit("fail-closed: causal_relation returned "
                         "undecided at a G3 stress point")

    monkeypatch.setattr(o4, "assemble", abort)
    monkeypatch.setattr(sys, "argv",
                        ["o4", "--freeze-rev", _REV])
    with pytest.raises(SystemExit, match="fail-closed"):
        o4.main()

    assert reservation.exists(), "the abort left the streams unreserved"
    rec = json.loads(reservation.read_text(encoding="utf-8"))
    assert rec["seeds"] == {"g1": 7, "g2": 7}
    assert "g3" not in rec["seeds"]
    assert "even if the run aborts" in rec["note"]
    assert rec["reservation_object"] == "c" * 40
    # the global claim is made first: the local file only records it
    assert claimed and claimed[0]["seeds"] == {"g1": 7, "g2": 7}

    # and a second attempt on the same streams is now refused
    with pytest.raises(SystemExit, match="write-once"):
        o4._publish_write_once(reservation, "{}\n")


def test_prereg_pins_the_epistemic_boundary():
    doc = (_REPO / "docs" / "prereg" / "p14_o4_volume_audit.md") \
        .read_text(encoding="utf-8")
    plain = " ".join(doc.replace("**", "").split())
    assert "Nothing in Paper A is upgraded" in plain
    assert "not a Poisson-count audit" in plain
    assert "never called a certification" in plain
    assert "full certified width" in plain
    assert "three-way consistency" in plain
