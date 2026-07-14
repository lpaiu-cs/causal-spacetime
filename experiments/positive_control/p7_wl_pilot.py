"""P7 enhanced-sampling pilot: can Wang-Landau cross the action barrier?

This is a *feasibility probe*, not a preregistered experiment. It answers one
question and refuses to answer any other: at a given N, does a flat-in-S
walker complete round trips between the continuum and crystal ends of the
action window, and at what cost per round trip?

Why this has to run before N=600. The frozen N=600 grid showed
start-separated action gaps of 77--96 for beta >= 18; local Metropolis never
crosses them, so no beta_c is extractable. Wang-Landau is the proposed fix,
but a proposed fix that has never been shown to tunnel is just a longer way
to fail. The pilot sweeps N upward and reports round trips and cost, so the
N=600 production budget is set from measured tunneling times rather than
from hope.

Deliberately NOT done here:

- No beta_c, no G(beta), no phase claim. Round trips and cost only.
- No frozen constants. When the N=600 production protocol is frozen, it must
  be frozen in its own `docs/prereg/` file, with the bin count, window,
  flatness criterion, and seeds fixed before any production chain is run.
- No pooling across walkers that did not converge.

Usage:
    python p7_wl_pilot.py --n 60 --sweeps 400
    python p7_wl_pilot.py --n 150 --bins 60 --sweep-steps 20000
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

from causal_spacetime_lab.positive_control.multicanonical import (
    action_range,
    effective_sample_size,
    multicanonical_2d_order,
    reweight_to_beta,
    wang_landau_2d_order,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "positive_control"

# The P7 design axis holds eps*N fixed: the frozen N=600 grid runs eps=0.02, so
# eps*N = 12. A pilot at smaller N MUST rescale eps to keep this invariant --
# running N=60 at eps=0.02 gives eps*N = 1.2, a regime with no crystallization
# transition at all, so the walker has no barrier to cross and the run is
# vacuous. (This is not hypothetical: the first pilot run made exactly that
# mistake and reported a flat <S> across every beta.)
EPS_N = 12.0
REWEIGHT_BETAS = [0.0, 4.0, 8.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0, 28.0]

# The multicanonical run is the real acceptance test for ln_g, and it is a
# harsher one than Wang-Landau's own convergence flag. If ln_g were the true
# density of states, the production chain would random-walk flat across the
# whole action window. So: it must actually traverse it.
#
# These thresholds exist because a run passed a weaker check and was still
# worthless. An earlier pilot reported converged=True with a production chain
# that covered 4.7% of the window at 1.3% acceptance -- ln_f had reached its
# target while ln_g was still wrong. A spread-only guard waved it through.
MIN_MUCA_COVERAGE = 0.5
MIN_MUCA_ROUND_TRIPS = 2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=60, help="number of elements")
    parser.add_argument(
        "--eps",
        type=float,
        default=None,
        help="smearing; defaults to EPS_N/n so that eps*N is held fixed",
    )
    parser.add_argument(
        "--eps-n",
        type=float,
        default=EPS_N,
        help="the invariant eps*N of the P7 design axis (frozen grid uses 12)",
    )
    parser.add_argument("--bins", type=int, default=40)
    parser.add_argument("--sweep-steps", type=int, default=10_000)
    parser.add_argument("--max-sweeps", type=int, default=2000)
    parser.add_argument("--ln-f-final", type=float, default=1e-4)
    parser.add_argument("--production-steps", type=int, default=500_000)
    parser.add_argument("--seed", type=int, default=20260715)
    args = parser.parse_args()

    eps = args.eps if args.eps is not None else args.eps_n / args.n

    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    pi0 = rng.permutation(args.n)

    s_min, s_max = action_range(args.n, eps, seed=args.seed, probes=200)
    print(f"N={args.n} eps={eps:g} (eps*N={eps * args.n:g}) "
          f"window=[{s_min:.3f}, {s_max:.3f}] bins={args.bins}", flush=True)

    started = time.perf_counter()
    wl = wang_landau_2d_order(
        pi0=pi0,
        eps=eps,
        s_min=s_min,
        s_max=s_max,
        n_bins=args.bins,
        seed=args.seed,
        sweep_steps=args.sweep_steps,
        ln_f_final=args.ln_f_final,
        max_sweeps=args.max_sweeps,
    )
    wl_seconds = time.perf_counter() - started
    wl_moves = wl.sweeps * args.sweep_steps

    print(
        f"WL: converged={wl.converged} sweeps={wl.sweeps} "
        f"round_trips={wl.round_trips} acc={wl.acceptance:.3f} "
        f"ln_f={wl.final_ln_f:.2e} time={wl_seconds:.1f}s",
        flush=True,
    )
    # Edge pile-up check: action_range is a heuristic bracket, and a window that
    # clips reachable states rejects moves into them silently.
    if wl.visited[[0, -1]].any():
        print(
            "WARNING: the walker reached an edge bin. The action window may be "
            "truncating reachable states -- widen it before trusting ln_g.",
            flush=True,
        )

    # Refuse, loudly, rather than emit reweighted numbers from weights that were
    # still moving. An unconverged ln_g does not flatten the action landscape,
    # so the production chain never leaves the entropic bulk -- and then every
    # sample carries essentially the same S, every reweighting weight is equal,
    # and the Kish ESS comes out at its MAXIMUM. A failed run would print a full
    # beta table with "ok" beside every row. That trap is the reason for the
    # hard exit here.
    fatal = []
    if not wl.converged:
        fatal.append(
            f"Wang-Landau did not converge (ln_f={wl.final_ln_f:.2e}, target "
            f"{args.ln_f_final:.0e}): ln_g is not a density of states yet."
        )
    if wl.round_trips == 0:
        fatal.append(
            "Zero Wang-Landau round trips: the walker never crossed the action "
            "barrier, which is the one thing this sampler exists to do."
        )
    if fatal:
        for line in fatal:
            print(f"REFUSING TO REWEIGHT: {line}", flush=True)
        print(
            "\nNo beta table is produced. Raise --max-sweeps, widen the window, "
            "or reconsider the method at this N. Do not read the WL diagnostics "
            "above as a partial result.",
            flush=True,
        )
        raise SystemExit(1)

    production = multicanonical_2d_order(
        pi0=pi0,
        eps=eps,
        ln_g=wl.ln_g,
        bin_edges=wl.bin_edges,
        steps=args.production_steps,
        seed=args.seed + 1,
        sample_every=max(1, args.production_steps // 2000),
        burn_frac=0.2,
    )
    print(
        f"MUCA: samples={len(production.samples)} "
        f"round_trips={production.round_trips} acc={production.acceptance:.3f}",
        flush=True,
    )

    # Degeneracy guard. If the production chain barely moved in S, the beta
    # table below is meaningless no matter how healthy its ESS looks.
    sampled_actions = np.array([row["S"] for row in production.samples])
    action_spread = float(sampled_actions.max() - sampled_actions.min())
    window_span = s_max - s_min
    coverage = action_spread / window_span if window_span > 0 else 0.0
    print(
        f"MUCA action coverage: spread={action_spread:.3f} over a "
        f"{window_span:.3f} window ({coverage:.1%})",
        flush=True,
    )
    muca_fatal = []
    if coverage < MIN_MUCA_COVERAGE:
        muca_fatal.append(
            f"the production chain covered {coverage:.1%} of the action window "
            f"(need {MIN_MUCA_COVERAGE:.0%}). A correct ln_g makes this chain "
            f"flat in S; a chain that stays put means ln_g is still wrong, "
            f"whatever Wang-Landau's convergence flag says."
        )
    if production.round_trips < MIN_MUCA_ROUND_TRIPS:
        muca_fatal.append(
            f"the production chain completed {production.round_trips} round "
            f"trips (need {MIN_MUCA_ROUND_TRIPS}). Without traversals there is "
            f"no demonstrated overlap between the two basins."
        )
    if muca_fatal:
        for line in muca_fatal:
            print(f"REFUSING TO REWEIGHT: {line}", flush=True)
        print(
            "\nNo beta table is produced. The Wang-Landau stage needs a longer "
            "run, more bins, or a wider window.",
            flush=True,
        )
        raise SystemExit(1)

    rows = []
    for beta in REWEIGHT_BETAS:
        mean_s, _ = reweight_to_beta(production.samples, beta, "S")
        mean_n0, var_n0 = reweight_to_beta(production.samples, beta, "n0")
        mean_h, _ = reweight_to_beta(production.samples, beta, "height")
        ess = effective_sample_size(production.samples, beta)
        rows.append({
            "n_elements": args.n,
            "eps": eps,
            "beta": beta,
            "mean_S": mean_s,
            "mean_n0": mean_n0,
            "susceptibility_n0": var_n0 / args.n,
            "mean_height": mean_h,
            "reweight_ess": ess,
            # Low ESS means the flat-S run has no real support at this beta. The
            # row is reported, not silently dropped -- but it must not be read as
            # an estimate. Note this flag is only meaningful once the degeneracy
            # guard above has passed.
            "reweight_supported": bool(ess >= 20.0),
        })
        print(
            f"  beta={beta:5.1f} <S>={mean_s:+9.2f} <n0>={mean_n0:9.1f} "
            f"ESS={ess:8.1f} {'ok' if ess >= 20.0 else 'UNSUPPORTED'}",
            flush=True,
        )

    tag = f"n{args.n}_eps{eps:g}"
    with open(OUT / f"p7_wl_pilot_{tag}.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "n_elements": args.n,
        "eps": eps,
        "eps_n": eps * args.n,
        "window": [s_min, s_max],
        "bins": args.bins,
        "wl_converged": wl.converged,
        "wl_sweeps": wl.sweeps,
        "wl_moves": wl_moves,
        "wl_round_trips": wl.round_trips,
        "wl_acceptance": wl.acceptance,
        "wl_seconds": wl_seconds,
        "wl_moves_per_round_trip": (
            wl_moves / wl.round_trips if wl.round_trips else None
        ),
        "muca_round_trips": production.round_trips,
        "muca_acceptance": production.acceptance,
        "muca_samples": len(production.samples),
        "status": "feasibility probe; nothing frozen; no phase claim",
    }
    (OUT / f"p7_wl_pilot_{tag}_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"\nwrote {OUT / f'p7_wl_pilot_{tag}.csv'}", flush=True)


if __name__ == "__main__":
    main()
