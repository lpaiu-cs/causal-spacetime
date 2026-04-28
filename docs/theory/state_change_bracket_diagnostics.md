# State-Change Bracket Diagnostics

Milestone 21 studies order-level brackets induced by a selected reference
chain in a finite state-change trigger network.

This is a reference-protocol diagnostic. It does not reconstruct metric radar
distance, calibrated clock time, finite-speed spatial geometry, or a physical
observer.

## Order-Level Brackets

Given a reference chain

```text
R = (r_0, r_1, ..., r_m)
```

and a target event `e`, define:

```text
p_R(e) = max { i : r_i ≺_T e }
q_R(e) = min { j : e ≺_T r_j }
```

`p_R(e)` is the predecessor reference position and `q_R(e)` is the successor
reference position. These are positions along the chosen reference chain, not
clock labels.

## Two-Sided Accessibility

An event is two-sided accessible relative to `R` when both bracket positions
exist:

```text
p_R(e) exists and q_R(e) exists
```

Reference-chain events themselves are excluded by default from target
bracketing. They are the reference backbone rather than bracketed targets.

## Rank-Level Quantities

For two-sided accessible events:

```text
T_rank(e) = p_R(e) + q_R(e)
W_rank(e) = q_R(e) - p_R(e)
```

`T_rank` is an order-level radar-time rank. `W_rank` is an order-level
bracket-width rank. The bracket-width rank is not metric distance and not
physical spatial distance.

Rank slices are coarse bins of `T_rank`. They are reference-chain-relative
order bins, not global time slices.

## Reference-Source Dependence

Local-system reference chains, greedy order-only chains, longest order chains,
and random baselines can produce different bracket coverage and rank orderings.
Milestone 21 records this as protocol-reference choice dependence.

That dependence is not a failure. It identifies which bracket-rank statements
are defined only after a reference protocol is chosen.

## Relation To Earlier Milestones

Milestone 19 provides the finite state-change trigger network. Milestone 20
scores candidate reference chains. Milestone 21 asks what a selected reference
chain can order-access through predecessor and successor brackets.

Milestone 22 extends this by fixing a reference emission position first. It
then asks which target events are after that emission and return to a later
reference position. This produces same-emission echo-return positions and
echo-delay ranks, still without metric distance or calibrated clock values.

Milestone 23 uses controlled echo-response motifs to validate those echo
diagnostics. It plants target events with known return positions and records
exact recovery, shortcut returns, and reference-protocol visibility.
