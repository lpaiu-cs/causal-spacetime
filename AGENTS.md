# AGENTS.md

## Project identity

This repository is a research-oriented simulation lab for studying whether spacetime quantities can be operationally reconstructed from causal information structure.

The working research hypothesis is:

> Time and space are not derived from velocity as displacement divided by time. Instead, operational time and distance can be reconstructed from primitive causal or null information-accessibility relations, with metric scale requiring additional structure such as event density, clocks, or free-fall information.

This project should remain scientifically conservative. Do not overclaim that the simulations prove a new theory of physics. The simulations should test reconstruction procedures, failure cases, and known relativistic behavior.

## Main goals

Build a Python simulation codebase for:

1. Radar coordinate decomposition in special relativity.
2. Lorentz kinematics and length contraction as observer-dependent event selection.
3. Causal set sprinkling in Minkowski spacetime.
4. Timelike proper-time reconstruction from longest chains.
5. Interval-volume reconstruction from Alexandrov interval cardinality.
6. Spacelike-distance reconstruction experiments, initially as exploratory code.
7. Finite-speed lattice counterexamples showing that finite signal speed alone does not imply Lorentzian spacetime.

## Scientific constraints

Use natural units `c = 1` by default.

Distinguish carefully between:

- signal speed,
- causal order,
- null relation,
- observer protocol,
- proper time,
- coordinate time,
- spatial distance,
- metric scale.

Avoid language implying that speed is primitive in the usual `v = dx/dt` sense.

Prefer this language:

- "null information-accessibility structure"
- "causal order"
- "operational reconstruction"
- "radar decomposition"
- "event density"
- "causal interval cardinality"

Avoid this language unless explicitly discussed as a rejected formulation:

- "time is information-transfer speed"
- "space is information-transfer speed"
- "this proves spacetime is information"

### House rule: an accuracy gate is evaluated on the interval, not the point

**Adopted 2026-08-01. Derivation and price: `docs/prereg/p14_weyl_curvature.md` §6.**

When a preregistered gate is a claim about **accuracy** — "this instrument
recovers `X` to within `t`" — the verdict is decided by the confidence
interval, not by the point estimate. A point estimate inside `t` whose
interval crosses `t` does not support the sentence the record will quote.

The case that forced this. P12 §9.3's co-requirement (ii) was frozen as a
point comparison; the campaign's top rung read `0.1708` against a threshold
of `0.25` and passed, while its 95% interval reached `0.2644` — across the
line. §9.3 declined to re-score, correctly: converting a point rule to an
interval rule after seeing the interval is the same move as loosening a
threshold after seeing the data, in the opposite direction and equally
forbidden. So the rule is fixed here, in advance, for designs that do not
yet exist.

**What it costs, computed on the case that forced it rather than asserted.**
At P12 Stage B's realized variance (effective `SE = 0.0474`), clearing
`0.25` with the whole interval needs `SE <= 0.0404`, i.e. `1.38x` the
samples: `n = 1538` per rung per arm against the `1117` actually run — and
against P12's own frozen cap of `1300`. **Under this rule P12 would have
declared itself infeasible.** Adopt it knowing that, or do not adopt it.

Note where P12's `n` came from: `n_sup = 101` and `n_eq = 1117`, both from
the **rate** gate. The recovery gate's precision was never a sizing
constraint and simply came along for free. Under this rule it becomes one,
and on those numbers the binding one. **Every design from here sizes for
its accuracy gate explicitly, in the power section, before any data.**

Three things this rule is not:

- **Not retroactive.** P12's record stands on the rule it froze.
- **Not free.** See the paragraph above. A design adopting it and then
  discovering it cannot afford it has loosened a gate.
- **Not a licence to switch.** The direction is one-way: a design may not
  move from interval to point, ever, and a design that wants a point rule
  must argue for it in its own text before freezing.

Gates that are **not** accuracy claims — existence, sign, ordering, budget,
completion — keep whatever form their design freezes, and this rule says
nothing about them.

## Implementation standards

Use Python 3.11 or later.

Use:

- `numpy`
- `scipy` if needed
- `matplotlib`
- `networkx` only if useful and not performance-critical
- `pytest`
- `ruff`

Do not use heavyweight dependencies unless justified.

Keep the package modular.

Suggested structure:

```text
causal-spacetime-lab/
  README.md
  pyproject.toml
  AGENTS.md
  src/
    causal_spacetime_lab/
      __init__.py
      constants.py
      lorentz.py
      radar.py
      sprinkling.py
      causal.py
      chains.py
      intervals.py
      lattice.py
      metrics.py
      plotting.py
  experiments/
    exp01_radar_decomposition.py
    exp02_lorentz_length_contraction.py
    exp03_causalset_timelike_reconstruction.py
    exp04_interval_volume_reconstruction.py
    exp05_finite_speed_lattice_counterexample.py
  tests/
    test_lorentz.py
    test_radar.py
    test_sprinkling.py
    test_causal.py
    test_chains.py
  docs/
    research_notes.md
    references.md
```

## Merge policy: never squash, never rebase-merge

**Merge pull requests with a merge commit only.** Repository settings now
enforce this (`allow_squash_merge` and `allow_rebase_merge` are both off), but
the reason belongs here so it is not re-enabled by someone tidying settings.

Every frozen artifact in `docs/prereg/frozen/` records a `code_version`: the
short SHA of the commit that produced it. Those stamps are the programme's
provenance mechanism, and several records make claims about ORDER IN THE
HISTORY rather than about file contents — for example that a diagnostic's
predictions were committed before it ran, or that a design was frozen before
its implementation and before any data.

- **Squash** collapses a branch into one commit, so every stamp names an object
  no checkout can resolve and every ordering claim becomes an unverifiable
  assertion.
- **Rebase-merge** preserves the individual commits but REWRITES their SHAs,
  which breaks the stamps just as thoroughly and less visibly.

The cross-stage gate machinery also requires a *reachable* stamp, so both
strategies silently turn that precondition into a no-op rather than failing
loudly.

Corollaries for anyone working in this repository:

- Do not rewrite published history on a branch whose commits are already
  referenced by a `code_version` (no `rebase`, no `--amend`, no force-push).
- When a stage runs, commit first and run from a clean worktree, so the stamp
  the artifact records is a commit that exists. The runners enforce this with a
  preflight check that counts untracked files as dirty.
- A record may cite a stamp only if
  `git merge-base --is-ancestor <stamp> HEAD` succeeds, where **HEAD is the
  branch that carries the record, or any merge commit that has that branch as
  a parent.** The gate is a statement about the repository's real history, not
  about every possible export of it.

### What the gate does NOT mean, because this has been misread three times

Automated review environments commonly evaluate a pull request as a **synthetic
single-parent commit**: the branch's final tree grafted onto the fork point, or
onto the base tip, with the branch's individual commits absent. On such an
object `--is-ancestor <stamp> <synthetic>` exits 1 for every stamp on the
branch, because none of those commits are in its ancestry. That is a property
of the synthetic object, not evidence that the branch flattened its own history.

Three findings on PRs #35 and #36 (C1, C15, C23/C26) were built on exactly this
and reported as P1 provenance failures. The distinguishing check, which anyone
auditing should run before drawing the inference:

```
git cat-file -t <cited-sha>                  # does the object exist here at all?
git rev-list --parents -n 1 <cited-sha>      # one parent, or two?
git rev-list --parents -n 1 <branch-head>    # what the branch actually looks like
```

For PR #36 at head `5dd2537` the cited object did not exist in the repository,
while `git merge-base --is-ancestor 7ac5893 X` exited 0 for X in `HEAD`, the
pushed branch, `refs/pull/36/head`, **and GitHub's own `refs/pull/36/merge`**
(which has two parents, base tip and branch head, and therefore does preserve
the stamps).

The concern underneath the misreading is real, which is why it is answered by
configuration rather than by vigilance: squash and rebase merges are disabled at
the repository level, so no merge performed through GitHub can produce the
flattened history these findings describe. PR #35 is the worked example — merged
with a merge commit, and all six of its stamps verify reachable from
`origin/main` afterwards.
