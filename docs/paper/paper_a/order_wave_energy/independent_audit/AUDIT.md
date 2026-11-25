# Independent audit of the order-only resolution bracket

2026-11-24. Audit record for the three arguments in this folder that bound,
for the exact circular vacuum plane-wave family with `N = lambda^-alpha`, the
exponents at which the causal order alone can or cannot distinguish the wave
from its effective background.

This is an audit by AI referees. No human mathematician has reviewed these
arguments. "CORRECT" below means that an independent referee re-derived every
step and found no false step; it does not mean a published or peer-reviewed
theorem. Nothing here is a finite-sample result, and no claim is made beyond
this one family with known q.

## 1. What was audited and how

| # | Claim | Source documents |
|---|---|---|
| 1 | For every fixed `alpha < 4/3` the full order is asymptotically indistinguishable (`TV -> 0`) | `quantized_order_bound.md`, `full_order_resolution_probe.md` |
| 2 | With known `0 < q <= 1`, an explicit order-only test separates when `N lambda^3 / log N -> infinity` | `local_null_challenge.md`, `local_null_results.md` |
| 3 | Draft: for each fixed integer `K >= 1` an order-only test separates when `alpha > 2 + 2/K`, hence every fixed power above 2 | `gram_fourier_noise.md`, with `gram_fourier_challenge.md`, `gram_fourier_results.md` |

All three were derived by one AI author on 2026-09-24 (branch
`codex/order-wave-correlation`) and had not been audited.

Method.

- One referee per argument, each a sub-agent with a fresh context. Argument 3
  was split in two (background noise bound and scan; wave-side coupling and
  signal), each half taking the other half's displayed conclusions as
  hypotheses.
- Each referee received the paths of the source documents and the claim to
  referee. None received the orchestrator's assessment or another referee's
  work, and each was told to ignore status sentences in the sources.
- Each referee re-derived the formulas and wrote its own check code. Only one
  script (`proof2_shared_reference/s6_author_routine_crosscheck.py`) calls the
  author's code, deliberately, to feed it the referee's numbers.
- The orchestrator read the three arguments and re-derived their main chains
  before the reports arrived (Section 4), then re-ran every referee script
  from a clean worktree (Section 8).
- A separate sub-agent searched the literature (Section 6).

Grades: ERROR = a step or the conclusion is false as stated. GAP = a step is
not justified by the text. MINOR = the statement is true and the missing
piece is routine, or it is a wording, scope or practicality note.

## 2. Verdicts

| # | Verdict | Stated precisely |
|---|---|---|
| 1 | **CORRECT** | For every fixed `K >= 1` and `0 < q < pi^2` there is `C_K(q)` with `TV <= C_K [N lam^2 + (N lam^4)^(1/(2K+1)) + (N^(3K+1) lam^(4K))^(1/(2K+1))]`, uniformly in phase. Hence `TV -> 0` at every fixed `alpha < 4/3`. Nothing is proved at `alpha = 4/3`. |
| 2 | **CORRECT**, with `lambda -> 0` added to the statement | For fixed known `0 < q <= 1`, `lambda -> 0` and `N lambda^3 / log N -> infinity`, the stated test has false-positive probability at most `N^-2` for every `N >= 4` and power tending to one. |
| 3 | **CORRECT as an asymptotic existence argument**, both halves | The noise bound (5), the wave-side bound (10), the signal bound (11) and the scan logic are each true. The condition delivered is exactly `alpha > 2 + 2/K`, strictly. The argument is for a regularised statistic that the draft describes but does not write down, and one sentence of the draft is literally false without affecting any inequality (Section 3). |

No referee reported an ERROR or an unrepaired GAP. Every referee graded its
findings MINOR and supplied the short argument for each compressed step.

Combined, with `alpha_low` and `alpha_high` as defined in
[next_question.md](next_question.md):

    4/3 <= alpha_low <= alpha_high <= 2.

Reports: [proof 1](proof1_lower_bound/report.md),
[proof 2](proof2_shared_reference/report.md),
[draft 3, noise and scan](proof3a_null_noise/report.md),
[draft 3, coupling and signal](proof3b_wave_signal/report.md).

## 3. Findings to carry into the source documents

The source documents carry frozen stamps and this audit does not edit them.
These are the corrections to make when they are next revised.

### Argument 1 (lower bound)

1. Link the constants: the non-edge claim needs
   `(L - 1 - sqrt 3) h >= c0 delta`, with `c0 = 2 sqrt 2 C_E`,
   `C_E = 2 sqrt(q) (4 sec^2(sqrt(q)/2) - 2)`. With `L = 1` the claim fails
   (154 violations in the referee's control).
2. State the hypothesis of eq. (3) (crossing set only) and say why an edge has
   one orientation under both kernels (`T` is a common time function).
3. Write the smoothing lemma of Section 4 as a lemma with its constant,
   `L = [(m+1) rbar + max_j rbar'_j] / h^2`. The referee's slicing proof is in
   its report, Section 3.5.
4. State the Hellinger convention `H^2 = sum (sqrt p - sqrt q)^2` and that the
   number of components with an edge is at most the number of uncertain pairs.
5. `K = 1` is valid and returns the old `N lambda -> 0` bound. The chosen `h`
   is not the minimiser of (7) when `N h` dominates; this costs rate only.
6. The "strict separation" sentence also relies on `kappa(q) > 0` from
   `coordinate_oracle_barrier.md`, which was not audited here.

### Argument 2 (shared references)

1. Add `lambda -> 0` to the boxed statement. A sequence with `lambda` fixed
   and `N -> infinity` satisfies the hypothesis as written and is not covered.
2. Replace the imported bias expansion by the uniform bound
   `0 <= K_0 - Kbar_lambda <= C lambda^2 min(l, lambda)^2 / l` (referee report,
   Section 6.1). Only this upper bound is needed.
3. Say "the least integer `B` with `2^B >= 32 N^4`".
4. The rectangular domain `R` forces `[a_X/p_Y^2, a_Y/p_X^2]` into
   `(p0(q), 1/10)`, a window 4.1% wide at `q = 1` and 1.0% at `q = 0.25`. The
   test as specified therefore skips the frozen `q = 0.25` witnesses at
   `r = 1024, 2048` and the `q = 1` witnesses at `r <= 512` at any reference
   count. This does not affect the asymptotic claim. It does mean the "10 of
   12" population witnesses in `local_null_results.md` are not all certifiable
   even with unlimited data, and that file should say so.
5. Exactly one power of `log N` is needed, through the Bernstein term.

### Argument 3 (Gram-Fourier draft)

1. **Write the regularised statistic down.** The cutoff must act on all 18 raw
   moments, in particular on the anchor-interval moments `b`, not through the
   weight `w` alone. The statistic frozen in `gram_fourier_challenge.md`
   eq. (2) has only `w` and is not the one the noise bound covers. A probe
   that is not after anchor 3 has `H != 0` (the referee's example: 0.0029).
   Section 4 gives one explicit form.
2. With the cutoff in `b`, the identity (2) and the bound (5) hold for every
   anchor triple. "Good anchor" and the pilot-failure term are not needed
   under the background.
3. **Correct the false sentence in Section 6**: `|a_2x a_2y| = |G12|` holds
   only if `G22 = G11`. On the stated box the supremum is `2.0063e-5`, not
   `2e-5`. The certificate value `Im D_A >= 0.01069797` is overstated by
   `1.8e-7`; the corrected bound is `0.0106978`. Every target inequality
   survives: `Im D_A > 0.01` with margin `7e-4`, the grid coefficient
   `0.0101679 > 0.009`, and `|D_A| < 0.05` unchanged.
4. Write the Section 5 comparison as the identity
   `S_n - E F(mu_lam) = T1 + T2 + T3 + T4` (referee report 3b, Section 3.4)
   and add the remark that for probes outside an enlarged active set the
   terms vanish except on an event of probability `exp(-c n)`.
5. "Anchor pair slacks stay positive" means `b_ij > 0`; say so.
6. The null vote bound is `(K C_K / log n)^2`; eq. (5) lists `1/n` twice and
   needs `n_P >= c n`; the elementary check is exact only at
   `omega = 2 pi j`. The sentence "independent noise would lose (3)" is
   misleading: independent references give a smaller variance, and what (3)
   controls is the common-mode term of shared references.
7. The scan cap is the only binding requirement. A grid up to `sqrt(n)` keeps
   the null control and would give `alpha > 2 + 2/(K+1)`.

## 4. Orchestrator's cross-check

Re-derived by the orchestrator before reading the reports:

- The exact vacuum condition: for `gamma = a^2 R diag(e^(2 beta), e^(-2 beta)) R^T`
  with `theta' = 1/(2 lambda)`, `tr[(S^-1 S')^2] = 8 theta'^2 sinh^2(2 beta)`,
  so `a''/a = -q` exactly when `sinh(2 beta) = 2 lambda sqrt(q)`.
- Argument 1: the identity `dT^2 - |dz|^2 = 2 b dv - |dy|^2` behind the
  sandwich; the Hellinger step from the sandwich; the count
  `N^2 h (delta/h)^2`; the exponent algebra.
- Argument 3: the pathwise `1/omega` from bounded variation along a
  background null generator; the index count giving `n^(-k/2)`; the remainder
  device `A(mu_0 + t eps) = t J(t)`; the three replacements behind (10) and
  the orders `lambda/sqrt(n) + 1/n` of each.
- Exponent arithmetic of all three arguments and of the design draft:
  `orchestrator_checks.py`, 58 checks.

Two points of synthesis that no single referee could make.

**The interface between the two halves of argument 3.** Both referees named
the same least-trusted step: the regularised statistic must be one object for
which (5), (10) and (11) all hold. One explicit choice, supplied here and not
checked by a referee:

    A_reg(mu) = chi_anchor(mu_anchor) * prod_i chi_b(b_i) * w(u, v, R^2) * H,

where `chi_b` is smooth, zero for `b <= b_min/2` and one for `b >= b_min`,
`b_min` is the infimum of the population `b_{A_i X}` over active probes and
admissible anchors (positive, because the certified active-probe slack is at
least 0.100), and `chi_anchor` is a smooth cutoff in the 12 anchor moments,
one on the inner anchor neighbourhood and zero outside the outer one. The
support of `w` already keeps the probe moments `(p, f, a)` in a compact part
of the calibration's regular domain. Then `A_reg` is smooth with bounded
derivatives; it equals `w H` on the population vectors of both hypotheses for
active probes and accepted anchors once `lambda` is small, so the population
signal is unchanged; and on its support every `b` is positive, so every pair
is timelike, every Gram entry is a true inner product and `H = 0` for every
background probe. The phase `u(mu)` needs only a smooth extension off the
support. Existence is therefore not in doubt; what is missing is the text.

**Exponent 2 is set by the leading term, not by the remainder.** The referee
of the noise bound found that a grid capped at `sqrt(n)` would still be
controlled, and that in its 1+1 analogue the remainder bound is never
attained. The condition `alpha > 2` comes from comparing the signal
`lambda^2` with the leading shared-reference noise `1/(omega sqrt(n))` at
`omega = 1/lambda`, that is `lambda sqrt(n) -> infinity`. The referee of the
wave side found independently that the wave-specific error terms sit at the
same order `lambda/sqrt(n)` with no slack. So 2 is intrinsic to this
statistic. Increasing `K` only removes the artefact `2/K`.

## 5. These are statements about limits

Every referee measured how far the constants are from any reachable size.

| Argument | Measured | Consequence |
|---|---|---|
| 1 | `alpha = 1.3` needs `K >= 14` with exponent margin at most 0.0035; `c0` = 6.4, 18.1, 132, 13532 for `q` = 0.25, 1, 4, 9 | no finite-`lambda` content |
| 2 | certifying the `1/lambda = 2048` witness needs about `1e24` events with ideal constants, and a reference count between `1e25` and `1e26` with the implemented enclosures, against `lambda^-3 = 8.6e9` | no rejection on any sampled order is possible |
| 3 | `sd(u_hat)` about `74.5/sqrt(n)` at `q = 1`, so about `2e6` references before `u_hat` is inside the window; the anchor neighbourhood has probability about `7e-38` per candidate triple; `1/lambda` enters the scan near `lambda ~ 1e-40` for `alpha = 2 + 2/K + 0.1` | the complete instrument cannot be exercised at any size |

The bracket is therefore a set of asymptotic statements. It cannot be backed
by a finite-N campaign, and a reader has to accept the proofs. The source
documents already say that no finite-N power is claimed; these numbers say
how far away it is.

## 6. Prior art

Full report: [prior_art/report.md](prior_art/report.md). A search that fails
to retrieve a result is not a proof of novelty; the report lists its queries
and what it could not access.

- No item was found that states a threshold, minimax rate or
  total-variation indistinguishability law in N versus a wavelength for
  sprinkled causal orders, or that treats a high-frequency (Isaacson/Burnett)
  regime in a causal-set setting.
- Madsen, arXiv:2610.36695: a Riemann-tensor reconstruction rate
  `(L_rho/rho)^(1/[2(3d+4)])` for a fixed smooth geometry, with
  geometry-dependent constants. No lower bound, no wavelength, no wave.
- Ibrahim and Bombelli, arXiv:2611.09056: a heuristic element count
  (`1e8` to `1e9`) for one geodesic-deviation method on a linearised wave at
  wavelength comparable to the probed region. No lower bound; the authors say
  their estimate does not exclude other methods.
- Closest in form: Bombelli gr-qc/0002053 (a statistical distance between
  n-point order laws, evaluated for n = 2, 3) and Mueller arXiv:2503.01719
  Theorem 3 (an existential total-variation indistinguishability at a fixed
  number of points).
- Closest in technique: detection thresholds for random geometric graphs
  against Erdos-Renyi. The present problem is not a special case, because
  both hypotheses here are geometric and the relation is deterministic.
- The usual statement that structure below the mean spacing `rho^(-1/d)` is
  not encoded would put the threshold at `alpha = 4`. Every number in the
  bracket is below 4. No retrieved paper states or refutes this for a
  coherent wave filling the region.

## 7. Is the gap worth pursuing

**Recommendation: yes, starting with the one question that decides which
side to attack.** Detail, routes and stop rules are in
[next_question.md](next_question.md).

Reasons.

1. The bracket survived: three arguments, four independent referees, no
   error, 565 checks passing.
2. As far as a bounded search can tell, the question is open.
3. The gap is not a matter of polishing either proof. The lower-bound method
   cannot use larger cells (the uncertain-pair graph percolates) and the
   upper-bound statistic's exponent is set by its leading term. Both endpoints
   are properties of a method. The interval between them is one structural
   question: whether the order stores null-margin information only at the
   count level `N^(-1/2)` or below it.
4. Either answer is a law. If 2 is sharp, the order's resolution of a plane
   wave is exactly the count limit `N lambda^2 ~ 1`. If it is not, the order
   resolves below the count limit and the attaining observable is not a
   smooth functional of counts.

Risks, stated plainly.

- Asymptotic only (Section 5). This line produces proofs, not simulations.
- AI-only audit. The three places a human should look first are the general-K
  smoothing lemma in argument 1, the uniform bias bound in argument 2, and
  the regularised statistic together with the tip estimate (8) in argument 3.
- One family: exact circular plane wave, known `q <= 1`, fixed-N iid sampling.
  Nothing is claimed for unknown q, other waveforms, Poisson N or generic
  spacetimes.
- The heuristic accounting for the gap is inconclusive. Different assumptions
  give every value from 4/3 to 2.

What not to do next: do not run a finite-N power campaign of the Gram-Fourier
instrument (Section 5 shows it cannot be exercised), and do not add
population rows to either existing statistic (their exponents are fixed).

## 8. Provenance

- Scripts and three of the four referee reports were committed at `f05d673`
  before any stamped run.
- `run_all_checks.py` then ran from a clean worktree (it refuses otherwise,
  counting untracked files as dirty) and recorded `code_version = f05d673` in
  [run_logs/summary.json](run_logs/summary.json), with one stdout log per
  script. CPU only.

| Folder | Scripts | PASS | FAIL |
|---|---:|---:|---:|
| orchestrator | 1 | 58 | 0 |
| `proof1_lower_bound` | 8 | 301 | 0 |
| `proof2_shared_reference` | 6 | 89 | 0 |
| `proof3a_null_noise` | 6 | 64 | 0 |
| `proof3b_wave_signal` | 6 | 53 | 0 |
| total | 27 | 565 | 0 |

- The counts equal those each referee reported from its own runs.
- Before the stamped run the orchestrator applied `ruff format` and
  `ruff check --fix` to the referee scripts and added a file-level exemption
  for long lines and percent formatting, so that the repository lint stays
  clean. The stamped run is of the formatted scripts.
- The referees could not write report files; they returned their reports as
  text and the orchestrator saved them unchanged apart from one marked
  sentence each. The fourth report and this record were added after the run.
- To reproduce, from the repository root on a clean worktree:
  `python docs/paper/paper_a/order_wave_energy/independent_audit/run_all_checks.py --scipy-python <interpreter with NumPy and SciPy>`.

## 9. What this audit did not do

- It did not compute any total-variation distance between order laws.
- It did not exercise any of the three instruments on a sampled order.
- It did not audit `coordinate_oracle_barrier.md`, the earlier sufficient
  bounds at exponents 8 and 4, or the provenance sections (hashes, commit
  order) of the source documents.
- It did not re-fetch the prior-art sources.
- It did not edit any source document.
- It did not freeze or run any design. The two drafts in `next_question.md`
  are drafts.

## 10. Addendum, 2026-11-26: argument 3's statistic written out

Added after the audit, at the user's request ("논증 3 마무리"). Sections 1
to 9 above are unchanged.

- **What was done.** Section 3's items 1 to 7 for argument 3 were carried
  into a revision, [gram_fourier_noise_r2.md](../gram_fourier_noise_r2.md).
  The frozen original is not edited; its SHA-256 still matches the
  anchor-guard certificate.
  - The regularised statistic `A_reg` (r2, eq. (1a)) is written as one
    object.
  - Its cutoffs act on all 18 raw moments: the anchor numbers, the three
    anchor-pair `b_ij`, and the three anchor-to-probe `b_iX` from below.
  - The noise bound (5), the wave-side bound (10) and the signal bound (11)
    are applied to that same object.
  - The false sentence `|a_2x a_2y| = |G12|` is corrected, and `B` is
    defined.
- **Review.** A fresh-context AI referee:
  [argument3_revision/review.md](argument3_revision/review.md). Verdict:
  **CORRECT WITH MINOR FIXES.** It found:
  - one statement false as stated: compact support on `R^18`, false because
    large `b_iX` is not cut off. It is now stated on the bounded set
    `Q = [-1, 2]^18`, and nothing downstream changes;
  - seven minor items;
  - an overclaiming headline.

  All were corrected; r2 Section 9 lists them. Deterministic checks:
  [argument3_revision/checks.py](argument3_revision/checks.py) and the
  referee's [ref_checks.py](argument3_revision/ref_checks.py). Nothing was
  sampled.
- **What this does not change.**
  - The bracket `4/3 <= alpha_low <= alpha_high <= 2` is this audit's
    statement (Section 2). It is not strengthened.
  - Argument 3's remaining steps are still at outline level in its text,
    with the referee reports of Section 2 supplying them:
    - the index count of (5);
    - the tip estimates (6) to (9);
    - the uniform closeness of moments at the tips;
    - the expansion behind (11);
    - the U-statistic and Hoeffding bounds.
  - The certificate is the author's code. Its one wrong step is corrected
    by hand, not re-certified.
  - The source documents of arguments 1 and 2 have not been revised. Their
    items in Section 3 are still missing from their text, including the
    constant link of argument 1, without which its claim fails at `L = 1`.
  - The review is AI-level. No human has checked it.
