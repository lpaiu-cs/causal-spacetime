# Pairwise Response-Profile Comparison

Milestone 29 defines admissible pairwise response-profile comparison
protocols. These protocols are pre-metric. They do not produce spatial
distance, physical distance, calibrated radar distance, or metric geometry.

## Why A New Protocol Is Needed

A single-reference response-order signature is a scalar preorder over targets.
It does not compare target-target pairs. A multi-reference response profile is
richer, but it still only records protocol-relative response ranks and
reachability.

For target `i`, write:

```text
P_i = (D_i^(1), ..., D_i^(m))
M_i = reachability mask
```

A pairwise comparison protocol `pi` defines a response-profile dissimilarity:

```text
delta_pi(i,j)
```

Then it can induce response-comparison constraints:

```text
(i,j) <_pi (k,l)
  iff delta_pi(i,j) < delta_pi(k,l)
```

Pairwise response-profile comparison is a chosen pre-metric protocol. It is
not target-target physical distance.

## Admissible Protocols

The implemented protocols include:

- profile-separation comparison,
- rank-gap mean comparison,
- rank-gap median comparison,
- reachability-mismatch comparison,
- combined rank-gap and mismatch comparison.

Each protocol must declare its missing-data policy:

- `common_reachable`: compare only columns where both targets are reachable,
- `require_all_reachable`: require both targets reachable in all protocol
  columns,
- `penalize_mismatch`: include a reachability mismatch penalty.

Different choices can produce different response-comparison orders. That is
protocol-choice dependence, not a failure of the finite diagnostic.

## Null Baselines

Before any future embedding attempt, null baselines are required. The
implemented baselines shuffle delay ranks, shuffle reachability masks, permute
target profiles, or rebuild random profiles with the same column marginals.

These baselines test whether profile structure is nontrivial under the chosen
comparison protocol. Null-baseline behavior does not prove geometry.

## Scope

Milestone 29 stops at response-comparison constraints. It does not fit
embeddings, recover target-target distance, or interpret profile
dissimilarity as spatial structure.

Milestone 30 validates constraint pools produced by these comparison
protocols. It checks held-out protocol agreement, bootstrap stability,
null-baseline separation, and coverage before any future representation test.
Validated pools remain pre-metric response-comparison outputs.

Milestone 31 exports handoff manifests for eligible pools. Those manifests are
frozen input specifications and not fitted representations.
