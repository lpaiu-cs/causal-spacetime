# P14 S1 — the cost of causality on Schwarzschild: a price, not a verdict

Per `p14_weyl_curvature.md` §8.1 (approved 2026-08-01) and the P14
preregistration review's scope reduction: **S1 records the cost of
the causal-predicate COMPONENT on this solver, patch, and N. Even an
affordable price leaves the Schwarzschild volume/campaign path
separately unresolved (no closed-form diamond volume); an
unaffordable one holds only this solver-domain-budget path.**

## The solver and the frozen domain

`p prec q` by direct-null-geodesic time-of-flight in the exterior
patch: shell `r/M ∈ [10, 20]`, polar cap `1.0 rad` (pairwise
`Δψ ≤ 2.0`), `t`-extent `40 M`. The impact parameter is bisected on
the no-turn/one-turn families in `u = 1/r`; quadratures use the
exact cubic factorization `R(u) = (u_t − u) Q(u)` (cancellation-free
at the turning point) with Richardson-style error bounds; radial
pairs are the exact tortoise difference; `|Δt − T_min|` inside the
error bound is UNDECIDED, never silently classified (P1's interval
discipline). Correctness is pinned against the flat limit, the
radial closed form, the weak-field Shapiro delay in SCHWARZSCHILD
coordinates (`T = d + 2M ln[(r1+r2+d)/(r1+r2−d)] − M(s1/r1 + s2/r2)`
— the bare isotropic log form disagrees by `~0.2 M` here and was the
first draft's wrong check), symmetry, monotonicity, and the
patch-safety bound (worst direct perihelion `> 5 M`, photon sphere
at `3 M`).

## The price (this host; artifact `p14_s1_cost.json`)

| tol | μs/pair |
|---|---|
| 1e-06 | 605 |
| 1e-08 (default) | **768** |
| 1e-10 | 949 |

Plane-wave predicate on the same host at aniso-a1.0: `2.11 μs/pair`
→ **price ratio ≈ 364×**. (An earlier draft quoted 1181 μs/574×;
~35% of that was the perihelion root-finder iterating past float
convergence — solver waste, not physics, and the review caught it.)
Projection at the default tolerance (single core, `O(n²)` relations
per sample, P11–P13-like sizes):

| n per sample | s/sample | hours per 4800 samples |
|---|---|---|
| 1000 | 384 | **512** |
| 2000 | 1536 | 2047 |
| 4000 | 6144 | 8192 |

1500 sampled pairs at every tolerance rung: `373` related, `0`
undecided. Wall-clock is host-dependent; the artifact records the
host and the measured microseconds. Benchmark seeds `40000101/102`
are cost-only streams (never statistics), placed above every ledger
range.

## The sentence this buys (for the preregistration's freeze manifest)

> S1은 해당 solver·patch·N에서 인과 술어 구성요소의 비용을 기록한다:
> shell r∈[10,20]M·cap 2 rad·tol 1e-8에서 **~0.77 ms/pair —
> plane-wave 술어의 ~360×, n=1000 기준 4800-표본 캠페인 환산
> ~510시간(단일 코어)**. 감당 가능하더라도 Schwarzschild
> 부피·캠페인 경로는 별도 미해결이며, 감당 불가능하면 해당
> solver·도메인·예산 경로만 보류한다.

What it does NOT buy: a verdict on the generalization path, a
volume oracle, or any claim beyond this solver-domain-budget
triple.
