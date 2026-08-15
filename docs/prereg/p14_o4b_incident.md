# O4b incident — a publication-wiring defect, no scientific verdict

**This record grades nothing.** The O4b campaign ran to the end of both
gates and then stopped, one step before publication, on an
implementation defect in the runner's own exit path. No result file was
written and no verdict is published. The gate statistics that the run
computed are preserved here and in the checkpoint as *non-verdict*
partials, and whether they can become a result is the subject of a
separate recovery audit — not settled by this document.

## 1. What happened

The approved campaign ran once from the clean, exact checkout of the
freeze branch head `715865abc684224785e71a3130e17c50db35f947` (merged in
PR #75), against manifest digest
`cec650b9391af0fc11e4b6bb94455cdbc1a037c18ecec83149d0d4693e7d7be2`. In
order it:

- passed G3a — the 71-case wrapper contract, 6,537 metered solver calls,
  charged to the campaign budget before any fresh stream was opened;
- claimed `refs/o4b/reservation` = `46acee340bc247511546964b2925953721d5bb59`,
  the point of no return, spending `40,000,401` / `40,000,411`;
- ran G3b on the fixed prefix, then G1 to 26,200,000 points and G2 to
  1,072,696 points, each to its frozen sample size;
- reached the pre-publication provenance re-read and **raised
  `SystemExit` before writing any result.**

EXIT=1 after 43,531.6 s. Neither cap bound (55,351,929 of 80,000,000
calls; 43,531.6 of 86,400 s). This was **not** a fail-closed scientific
stop and **not** a cap: it was a defect in the code that verifies the
reservation before publishing.

## 2. The defect — a return value that was never returned

`Campaign.run` records whatever the post-G3a claim callback returns and
carries it to the exit check that must confirm the reservation ref still
holds *this* attempt's object before anything is published:

```python
self.reservation_object = on_g3a_passed()   # in Campaign.run
...
if on_before_publish is not None:
    self._staged("g2", lambda: on_before_publish(self.reservation_object))
```

But the real callback in `main` stores the claim object in an outer dict
and returns `None`:

```python
def claim() -> None:                 # returns None
    claimed["object"] = reservation.claim(payload)
```

So the ref was claimed correctly — the remote still holds
`46acee34…`, exactly the attempt object — yet `self.reservation_object`
was `None`. Eleven hours later the exit check received `None` and the
guard for that case fired:

```python
def verify_before_publish(obj):
    if obj is None:
        raise SystemExit("internal: no claim object to verify")
```

The reservation provenance the exit check exists to confirm was in fact
intact. What failed was not the check but the wiring that feeds it: the
object never reached it.

## 3. Why the suite did not catch it

Every success-path test drove `Campaign.run` with a *returning* stub for
the claim callback, e.g. `on_g3a_passed=lambda: "c" * 40`. The suite
therefore only ever exercised "the world where the callback returns the
object". The production path `main → claim → Campaign.run`, where the
callback returns `None`, was never run end-to-end in a test. Nineteen
review rounds hardened the failure boundaries; none ran the happy path
of the actual `main`. This was a preventable integration-test gap, and
the fix belongs with a regression test that binds the real `claim`
signature to `Campaign.run`.

## 4. What is preserved, and what it is not

The eleven hours of computation are not lost. `p14_o4b_incident.json`
and `p14_o4b_checkpoint.json` carry, as non-verdict partials:

- **G1** — n = 26,200,000; interval `[56.448806, 56.982283]`; identified
  discrepancy `[-0.899212, 0.769546]`; band 1.703411.
- **G2** — n = 1,072,696; leaking points 0; leak-upper 0.141951 within
  budget 0.141951.
- **availability** — fixed prefix 100,000 / 100,000 complete;
  eligible_rate 1.0; fully_testable_rate 1.0; zero-window points
  4,073,746.
- the PCG64 state, the freeze SHA, the manifest digest, the reservation
  object, the seeds, and the budget ledger.

**These are not a verdict.** The runner writes its single result artifact
only after the exit provenance check passes; that check never passed, so
no gate here has a *status* in the scientific sense, and this document
does not assign one. Reporting these as "the gates passed" would repeat
the O4 abort's retracted "G2 appears to have passed" — function
completion is not gate passage. Whether the preserved statistics support
a result is exactly what the separate recovery audit must decide.

## 5. Further runner-record provenance defects (recorded, not fixed)

The runner's machine-readable records carry three provenance labels that
are wrong or misleading. None of them changes a single computed number,
and the *ground truth* is present in the same files — but a careless
reader (especially an automated recovery audit) could misread them, so
each is named here and each is left in the JSON exactly as the runner
emitted it. **The JSON is preserved verbatim as the run produced it; the
corrections live in this authored layer, not by hand-editing the runtime
artifact.**

**5a. The budget ledger attributes every call to G3a.** It reads
`stage: "g3a"` and `reserved_by_stage: {"g3a": 55351929}` because the
budget's stage label was never advanced past G3a. The totals and cap
arithmetic are unaffected — the cap is charged against `reserved`, which
is stage-agnostic — but the per-stage breakdown is wrong.

**5b. The incident's `stage: "g2"` is the label of the pre-publish step,
not a claim that G2 is incomplete.** The stop happened *after* `run_g2`
returned and wrote the `g2_complete` checkpoint — at the pre-publication
reservation re-read, which the runner wraps in `_staged("g2", …)`, so a
failure there inherits the label `"g2"`. The same file's
`completed_gates.g2.status` is `"concordant"` over n = 1,072,696, and
`completed_gates_are_not_a_verdict` states the gates ran to their frozen
sample size. **A recovery audit must read `completed_gates`, not `stage`,
to decide what completed: G1 and G2 both did.** `stage` records only
where the run stopped.

**5c. The `g2_complete` checkpoint's `seed`/`rng_position` are the G1
stream, not G2's.** `Campaign.checkpoint` always serialises
`seed = o4b_g1_audit` (`40,000,401`) and `self.rng`, which is the G1
estimator's generator; `run_g2` draws from a *separate* local generator
seeded by `o4b_g2_leakage` (`40,000,411`) that is never stored on the
campaign and never serialised. So the RNG state in the `g2_complete`
record continues the **G1** stream. **A recovery audit must not read
`40,000,401` as the G2 seed:** G2's stream is `40,000,411`, and this
checkpoint does not carry its position. Both seeds are listed under the
incident's `seeds` block, so which is which is recoverable.

Defects 5a–5c are the runner's, distinct from §2's return-value bug but
of the same character — a value the runner writes without following its
meaning into the record. They are fixed with §2 in the recovery/fix
phase, each with a regression test.

## 6. Seeds, lineage, and what comes next

`40,000,401` / `40,000,411` move from FRESH to OBSERVED under their
functional names `o4b_g1_audit` / `o4b_g2_leakage`. They are spent: the
points were observed. They are **not** renamed `aborted`, because the
defect is in publication wiring and both gates completed with a status —
so, unlike O4, these streams may yet be the provenance of a verdict if
the recovery audit publishes one. `refs/o4b/reservation` is **retained,
not deleted**: deleting it would let the retired streams read as free
from a fresh clone. `refs/o4/reservation` (`c4da162`) and the unspent
`40,000,301` are unchanged.

Next, in order and each under its own approval:

1. this preservation record commits;
2. a recovery audit under a separate contract determines whether the
   preserved statistics can be published as an O4b result without any
   recomputation, and with what provenance;
3. the return-value defect (§2) and the runner-record provenance defects
   (§5a–5c) are fixed with regression tests that run the real `main` path
   end-to-end.

The retired streams may not be re-entered. Recovering the preserved
statistics is a reproduction of an observed run; a fresh O4b campaign
would need a new freeze with new scalars. This record makes no scientific
claim and grades neither S1 nor the oracle.
