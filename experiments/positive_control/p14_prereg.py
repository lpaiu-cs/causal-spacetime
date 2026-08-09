"""P14 preregistration -- stage runner (preparation commit, doc §8 step 1).

The PREREGISTERED stage over the completed P14 probe chain (PR #38-#48,
merge 1025c50). Design document: `docs/prereg/p14_prereg.md` (v0.4,
in-session review). Two claims at the frozen operating point
aniso-a1.0, never merged into one sentence:

- **C1 (paired ensemble mean).** theta_Delta = E[f_A - f_0] over
  paired readings of the SAME sprinkled points; confirmed when the
  95% t-CI of the mean paired difference sits entirely above
  +EPS_DELTA.
- **C2 (classifier ensemble discrimination).** The frozen s/AUC/BA
  rules separate curved from flat ensembles from single-poset global
  relation fractions -- an independent preregistered REPLICATION of
  P3-C's confirmation. Its non-confirmation never retroactively
  weakens the standing P3-C record (frozen sentences below).

Stage POSITIVE = C1 confirmed AND C2 confirmed, certified as a JOINT
event; joint equivalence is deliberately NOT defined (the joint null
rate is under the 0.90 floor at the frozen sizes). Non-positive
outcomes report each claim's own three-branch verdict with its own
frozen sentence.

**Execution order (frozen, doc §8).** Preparation commit P (this
code, no artifacts) -> preflight runs on a CLEAN checkout of P,
artifact records `code_version = P` -> final-freeze commit F carries
the certification artifact, the frozen n, and S1's cost sentence in
the freeze manifest -> the campaign runs on a clean checkout of F,
results record `code_version = F` -> results commit R. A test
asserts P is an ancestor of F is an ancestor of R. The campaign mode
REFUSES to run without the freeze manifest, and every stochastic
mode refuses to run on a dirty worktree.

Run:  python experiments/positive_control/p14_prereg.py preflight
      python experiments/positive_control/p14_prereg.py manifest <s1.txt>
      python experiments/positive_control/p14_prereg.py campaign
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import p12_stage_b
import p14_probe_p3c as p3c
from p14_plane_wave import Slab, arms, sprinkle
from p14_probe_p1 import clopper_pearson, relation_census
from p14_probe_p2 import student_t_crit
from seed_windows import assert_point_seeds_fresh

# ---------------------------------------------------------------------
# Frozen constants (prereg v0.4)
# ---------------------------------------------------------------------

#: Operating point and density BY IDENTITY with the probe chain --
#: one source, so a probe-side edit cannot silently diverge.
POINT = p3c.POINT
E_N = p3c.E_N

N_C1 = 3_000
N_C2 = p3c.N_ARM  # 4800/arm, identity with the certified P3-C size

#: C1 margin: the operational margin ANCHORED to the P3-C confirmed
#: block (review R-1) -- eps_s * pooled sd = 0.08061025 * 0.0044395871.
#: Frozen as a literal; a test cross-checks the product. Not a
#: physical margin, and never restandardized against stage data.
EPS_DELTA = 0.0003578762
EPS_DELTA_ANCHOR = (0.08061025, 0.0044395871)

#: Execution seeds -- exactly three, no spares (review R-4). Interrupted
#: runs restart the SAME seed; swapping seeds after an ambiguity
#: failure is forbidden (stage blocked instead). Placement (review
#: R-5): the v0.3-approved 20260861/71/72 fall inside P13 campaign
#: v3's LEDGER range [15000000..21503999], so house freshness cannot
#: be certified for them -- as it could not be for the date-styled
#: 2026xxxx probe seeds (a ledger-freshness violation on record; no
#: evidence of an exact integer stream reuse, and no retroactive
#: downgrade of the probe results). The stage sits at 40000000+,
#: past P12's full allocation decade [30000000..39999999] (P12
#: declared the space from 40000000 on as the next allocation).
CAMPAIGN_SEEDS = {"c1_paired": 40_000_061,
                  "c2_curved": 40_000_071,
                  "c2_flat": 40_000_072}

#: Bootstrap roots are ENTROPY VECTORS in the probe chain's [781, k]
#: convention (np.random.default_rng([781, k]) semantics), NOT spawn
#: keys of a shared root -- prereg v0.4 froze this meaning after the
#: review showed the two constructions give different streams.
BOOT_ROOTS = {"joint_effect": (781, 60),
              "c1_null": (781, 61),
              "c2_null": (781, 62)}
SPENT_BOOT_ROOTS = ((781, 41), (781, 50), (781, 51))
#: Children per replicate, frozen: joint spawns (c1, c2_curved,
#: c2_flat); the null branches spawn their own arms only.
BOOT_CHILDREN = {"joint_effect": 3, "c1_null": 1, "c2_null": 2}

B_JOINT = 4_000
B_C1_NULL = 20_000
B_C2_NULL = 20_000

#: Frozen reporting sentences. C2's three sentences never weaken the
#: standing P3-C confirmation (review v0.3): replication language
#: only, conflict recorded side by side, no retroactive cancellation.
C1_SENTENCES = {
    "confirmed": "paired ensemble 평균 이동이 ε_Δ를 넘는다",
    "equivalent": "이동은 ε_Δ 안이다",
    "inconclusive": "C1 inconclusive",
}
C2_SENTENCES = {
    "confirmed": "P3-C 분리를 독립적으로 재현했다.",
    "inconclusive": "이번 preregistered replication은 재확인하지 못했다; "
                    "기존 P3-C 확인 기록과 병기한다.",
    "equivalent": "이번 replication은 equivalence를 지지해 "
                  "기존 P3-C와 충돌한다.",
}

_REPO = Path(__file__).resolve().parents[2]
_P3E_ARTIFACT = _REPO / "docs" / "prereg" / "p14_probe_p3e_results.json"
_PREFLIGHT_ARTIFACT = _REPO / "docs" / "prereg" / "p14_prereg_preflight.json"
_FREEZE_MANIFEST = _REPO / "docs" / "prereg" / "p14_prereg_freeze.json"
_RESULTS_ARTIFACT = _REPO / "docs" / "prereg" / "p14_prereg_results.json"

#: The scalar spent set: every probe-chain seed, including the P3-C
#: confirmation streams 20260851/52 (consumed, hence spent).
_SPENT_SCALARS = frozenset(p3c.BURNED_SEEDS) \
    | frozenset(p3c.CAMPAIGN_SEEDS.values())
#: The house integer ledger: P11-P13 spent ranges plus P12's FULL
#: allocation decade -- the documented boundary, not just the used
#: ceiling 33263999 (review R-5: the check must be as strong as the
#: declared allocation).
_SPENT_RANGES = p12_stage_b.SPENT_RANGES + (
    ("P12 allocation decade", 30_000_000, 39_999_999),)


def assert_seed_layout() -> None:
    """Both stream kinds (prereg v0.4): scalar execution seeds against
    the spent scalars AND the house ranges; bootstrap roots against
    the spent entropy vectors. Structural spawn-key uniqueness is
    asserted exhaustively in the test suite."""

    assert_point_seeds_fresh(CAMPAIGN_SEEDS, _SPENT_SCALARS,
                             _SPENT_RANGES, "P14 prereg")
    roots = list(BOOT_ROOTS.values())
    assert len(set(roots)) == len(roots)
    clash = set(roots) & set(SPENT_BOOT_ROOTS)
    assert not clash, f"bootstrap roots reused: {sorted(clash)}"


def boot_layout(branch: str, b: int) -> list[tuple]:
    """The frozen spawn structure: root = SeedSequence(entropy vector),
    one child node per replicate, `BOOT_CHILDREN[branch]` grandchildren
    per replicate. Returns the per-replicate grandchild tuples."""

    root = np.random.SeedSequence(list(BOOT_ROOTS[branch]))
    return [tuple(rep.spawn(BOOT_CHILDREN[branch]))
            for rep in root.spawn(b)]


# ---------------------------------------------------------------------
# Provenance: execution commit and clean-checkout enforcement
# ---------------------------------------------------------------------


def code_version() -> str:
    """HEAD SHA of the checkout this run executes from. Aborts on a
    dirty worktree: constant-identity tests show the values match,
    only the recorded SHA proves WHICH commit executed (prereg
    v0.4), and a dirty tree makes the SHA a lie."""

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=_REPO,
                         capture_output=True, text=True,
                         check=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=_REPO,
                           capture_output=True, text=True,
                           check=True).stdout.strip()
    if dirty:
        raise SystemExit("stage computations run on a CLEAN checkout "
                         "only; commit or stash first")
    return sha


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------------
# C1: paired mean difference, t-interval, three branches
# ---------------------------------------------------------------------


def paired_samples(curved, flat, rho: float, seed: int, count: int):
    """Per-sprinkling PAIRED readings: sprinkle once, census the same
    points under both geometries (doc §4.1 -- identical intensity
    because det g = -1). Returns (delta, f_curved, f_flat,
    ambiguity dict) with per-reading ambiguous/escalated totals; the
    undecided-is-never-a-silent-False contract is inherited from the
    probe chain."""

    rng = np.random.default_rng(seed)
    delta = np.empty(count)
    fa = np.empty(count)
    f0 = np.empty(count)
    amb = {"curved": {"ambiguous": 0, "escalated": 0},
           "flat": {"ambiguous": 0, "escalated": 0}}
    for i in range(count):
        pts = sprinkle(flat, rho, rng)
        pairs = len(pts) * (len(pts) - 1) / 2
        ca = relation_census(curved, pts)
        c0 = relation_census(flat, pts)
        fa[i] = ca.related.sum() / pairs
        f0[i] = c0.related.sum() / pairs
        delta[i] = fa[i] - f0[i]
        amb["curved"]["ambiguous"] += ca.ambiguous
        amb["curved"]["escalated"] += ca.escalated
        amb["flat"]["ambiguous"] += c0.ambiguous
        amb["flat"]["escalated"] += c0.escalated
        if (i + 1) % 500 == 0:
            print(f"    c1 {seed}: {i + 1}/{count}", flush=True)
    return delta, fa, f0, amb


def c1_metrics(delta: np.ndarray, tcrit: float | None = None) -> dict:
    """Mean paired difference with the 95% two-sided Student-t CI on
    the between-sprinkling variance (doc §5.2's replication unit)."""

    n = len(delta)
    if tcrit is None:
        tcrit = student_t_crit(n - 1)
    mean = float(delta.mean())
    se = float(delta.std(ddof=1)) / math.sqrt(n)
    return {"mean": mean, "ci": (mean - tcrit * se, mean + tcrit * se)}


def c1_classify(ci: tuple[float, float], eps: float = EPS_DELTA) -> str:
    """Three branches, frozen positive direction (curved gains
    relations, theta_Delta > 0)."""

    lo, hi = ci
    if lo > eps:
        return "confirmed"
    if -eps < lo and hi < eps:
        return "equivalent"
    return "inconclusive"


# ---------------------------------------------------------------------
# Preflight: three branches, frozen B, CP-lower >= 0.90 each
# ---------------------------------------------------------------------


def _p3e_pairs(p3e_art: dict) -> tuple[np.ndarray, np.ndarray]:
    """The C1 effect source, frozen: P3-E's stored PAIRED per-sprinkling
    samples at the operating point (the paired pilot; its mean
    0.05018773 is below every unpaired block mean, which is why the
    unpaired minimum was rejected as a 'conservative' C1 effect)."""

    rec = next(r for r in p3e_art["ladder_g"] if r["label"] == POINT[0])
    return (np.asarray(rec["raw"]["f_curved"]),
            np.asarray(rec["raw"]["f_flat"]))


def preflight(p3e_art: dict) -> dict:
    """Certification at the frozen sizes: joint effect (B=4000),
    C1 centered null (B=20000), C2 null (B=20000); every branch must
    reach an exact CP 95% LOWER bound of 0.90. No certification
    failure may raise B or n automatically -- the stage blocks and
    the design review reopens (doc §4/§6)."""

    fa, f0 = _p3e_pairs(p3e_art)
    delta = fa - f0
    centered = delta - delta.mean()
    tcrit = student_t_crit(N_C1 - 1)

    t0 = time.perf_counter()
    joint = {"joint_confirmed": 0, "c1_confirmed": 0, "c2_confirmed": 0}
    for c1_ss, c2c_ss, c2f_ss in boot_layout("joint_effect", B_JOINT):
        r = np.random.default_rng(c1_ss)
        m1 = c1_metrics(delta[r.integers(0, len(delta), N_C1)], tcrit)
        ok1 = c1_classify(m1["ci"]) == "confirmed"
        ra = np.random.default_rng(c2c_ss)
        rb = np.random.default_rng(c2f_ss)
        m2 = p3c.metric_cis(fa[ra.integers(0, len(fa), N_C2)],
                            f0[rb.integers(0, len(f0), N_C2)])
        ok2 = p3c.classify(m2["ci_s"], m2["ci_auc"],
                           m2["ci_ba"]) == "confirmed"
        joint["c1_confirmed"] += ok1
        joint["c2_confirmed"] += ok2
        joint["joint_confirmed"] += ok1 and ok2

    c1_null = 0
    for (ss,) in boot_layout("c1_null", B_C1_NULL):
        r = np.random.default_rng(ss)
        m = c1_metrics(centered[r.integers(0, len(centered), N_C1)],
                       tcrit)
        c1_null += c1_classify(m["ci"]) == "equivalent"

    c2_null = 0
    for ca_ss, cf_ss in boot_layout("c2_null", B_C2_NULL):
        ra = np.random.default_rng(ca_ss)
        rb = np.random.default_rng(cf_ss)
        m = p3c.metric_cis(f0[ra.integers(0, len(f0), N_C2)],
                           f0[rb.integers(0, len(f0), N_C2)])
        c2_null += p3c.classify(m["ci_s"], m["ci_auc"],
                                m["ci_ba"]) == "equivalent"

    branches = {
        "joint_effect": {"reps": B_JOINT, "counts": joint,
                         "pass_key": "joint_confirmed"},
        "c1_null": {"reps": B_C1_NULL,
                    "counts": {"equivalent": c1_null},
                    "pass_key": "equivalent"},
        "c2_null": {"reps": B_C2_NULL,
                    "counts": {"equivalent": c2_null},
                    "pass_key": "equivalent"},
    }
    certified = True
    for name, block in branches.items():
        lo, hi = clopper_pearson(block["counts"][block["pass_key"]],
                                 block["reps"])
        block["pass_ci95_exact"] = [lo, hi]
        certified = certified and lo >= 0.90
    return {
        "code_version": code_version(),
        "point": POINT[0], "e_n": E_N,
        "n_c1": N_C1, "n_c2": N_C2,
        "eps_delta": EPS_DELTA,
        "eps_delta_anchor": list(EPS_DELTA_ANCHOR),
        "margins_c2": dict(p3c.P3C_MARGINS),
        "campaign_seeds": dict(CAMPAIGN_SEEDS),
        "boot_roots": {k: list(v) for k, v in BOOT_ROOTS.items()},
        "boot_children": dict(BOOT_CHILDREN),
        "effect_source": "p14_probe_p3e_results.json ladder_g "
                         "aniso-a1.0 raw paired samples",
        "branches": branches,
        "certified": bool(certified),
        "seconds": time.perf_counter() - t0,
    }


# ---------------------------------------------------------------------
# Freeze manifest (final-freeze commit F) and campaign (held)
# ---------------------------------------------------------------------


def write_freeze_manifest(s1_sentence: str) -> dict:
    """Step 4 of doc §8: records the certification artifact's digest
    and execution commit P, the S1 cost sentence, and the frozen
    sizes. Only the commit carrying this manifest is called the
    final freeze."""

    pre = json.loads(_PREFLIGHT_ARTIFACT.read_text(encoding="utf-8"))
    if not pre["certified"]:
        raise SystemExit("preflight is not certified; the stage is "
                         "blocked and the design review reopens")
    manifest = {
        "preflight_digest": _sha256(_PREFLIGHT_ARTIFACT),
        "preflight_code_version": pre["code_version"],
        "s1_sentence": s1_sentence.strip(),
        "n_c1": N_C1, "n_c2": N_C2,
        "eps_delta": EPS_DELTA,
        "campaign_seeds": dict(CAMPAIGN_SEEDS),
    }
    with open(_FREEZE_MANIFEST, "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=1)
        f.write("\n")
    return manifest


def run_campaign() -> dict:
    """The stage execution -- held until the final freeze exists. The
    manifest gate enforces doc §8's order mechanically: no campaign
    without a freeze manifest whose preflight digest matches the
    committed artifact."""

    if not _FREEZE_MANIFEST.exists():
        raise SystemExit("no freeze manifest: the campaign only runs "
                         "from a clean checkout of the final-freeze "
                         "commit (doc §8 step 5)")
    manifest = json.loads(_FREEZE_MANIFEST.read_text(encoding="utf-8"))
    if manifest["preflight_digest"] != _sha256(_PREFLIGHT_ARTIFACT):
        raise SystemExit("preflight artifact does not match the freeze "
                         "manifest digest -- stage blocked")
    assert_seed_layout()
    ver = code_version()

    label, w, du, dv, dx, dy = POINT
    slab = Slab(du=du, dv=dv, dx=dx, dy=dy)
    rho = E_N / slab.coordinate_volume
    curved, flat = arms(slab, w)

    delta, fa1, f01, amb1 = paired_samples(
        curved, flat, rho, CAMPAIGN_SEEDS["c1_paired"], N_C1)
    for arm_name, tot in amb1.items():
        assert tot["ambiguous"] == 0, \
            f"c1 {arm_name} ambiguity {tot['ambiguous']} != 0 -- " \
            "stage blocked (no seed swap; investigate first)"
    m1 = c1_metrics(delta)
    v1 = c1_classify(m1["ci"])

    a, amb_a, esc_a = p3c.arm_samples(
        curved, rho, CAMPAIGN_SEEDS["c2_curved"], N_C2)
    b, amb_b, esc_b = p3c.arm_samples(
        flat, rho, CAMPAIGN_SEEDS["c2_flat"], N_C2)
    assert amb_a == 0 and amb_b == 0, "c2 ambiguity != 0 -- stage blocked"
    m2 = p3c.metric_cis(a, b)
    v2 = p3c.classify(m2["ci_s"], m2["ci_auc"], m2["ci_ba"])

    art = {
        "code_version": ver,
        "point": label, "e_n": E_N,
        "n_c1": N_C1, "n_c2": N_C2,
        "eps_delta": EPS_DELTA,
        "margins_c2": dict(p3c.P3C_MARGINS),
        "seeds": dict(CAMPAIGN_SEEDS),
        "c1": {"metrics": m1, "verdict": v1,
               "sentence": C1_SENTENCES[v1],
               "ambiguity": amb1,
               "raw": {"delta": [float(v) for v in delta],
                       "f_curved": [float(v) for v in fa1],
                       "f_flat": [float(v) for v in f01]}},
        "c2": {"metrics": m2, "verdict": v2,
               "sentence": C2_SENTENCES[v2],
               "ambiguity": {"curved": {"ambiguous": amb_a,
                                        "escalated": esc_a},
                             "flat": {"ambiguous": amb_b,
                                      "escalated": esc_b}},
               "raw": {"f_curved": [float(v) for v in a],
                       "f_flat": [float(v) for v in b]}},
        "stage_positive": v1 == "confirmed" and v2 == "confirmed",
    }
    with open(_RESULTS_ARTIFACT, "w", encoding="utf-8",
              newline="\n") as f:
        json.dump(art, f, ensure_ascii=False, indent=1)
        f.write("\n")
    return art


def main(mode: str) -> None:
    if mode == "preflight":
        p3e_art = json.loads(_P3E_ARTIFACT.read_text(encoding="utf-8"))
        art = preflight(p3e_art)
        with open(_PREFLIGHT_ARTIFACT, "w", encoding="utf-8",
                  newline="\n") as f:
            json.dump(art, f, ensure_ascii=False, indent=1)
            f.write("\n")
        print(f"certified: {art['certified']} "
              f"({art['seconds']:.0f}s); artifact: "
              f"{_PREFLIGHT_ARTIFACT}")
    elif mode == "manifest":
        sentence = Path(sys.argv[2]).read_text(encoding="utf-8")
        m = write_freeze_manifest(sentence)
        print(f"freeze manifest written: {_FREEZE_MANIFEST}")
        print(f"  preflight_code_version: {m['preflight_code_version']}")
    elif mode == "campaign":
        art = run_campaign()
        print(f"stage_positive: {art['stage_positive']}; "
              f"c1: {art['c1']['verdict']}; c2: {art['c2']['verdict']}; "
              f"artifact: {_RESULTS_ARTIFACT}")
    else:
        raise SystemExit(f"unknown mode: {mode}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "preflight")
