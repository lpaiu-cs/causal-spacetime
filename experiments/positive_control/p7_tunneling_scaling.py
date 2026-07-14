"""P7: how does Wang-Landau tunneling time scale with N?

This measures the one number that decides whether the whole enhanced-sampling
approach reaches N=600, and it measures it before anything expensive is built
on top of it.

The N=60 pilot pinned the cost structure exactly. Wang-Landau's ln_f halves
only once the walker has traversed the action window, so after 6000 sweeps it
sat at ln_f = 4.88e-4 = 2^-11 with round_trips = 11: one halving per round
trip, precisely. Convergence to ln_f = 1e-5 therefore costs ~17 round trips,
and the only free parameter left is

    tau(N) = moves per round trip.

Everything follows from how tau grows:

- If tau ~ N^2, then N=600 costs roughly 100x the N=60 tunneling time, or ~7e9
  moves for a converged ln_g. With the Numba kernel (135 us/move) that is a few
  hundred core-hours -- heavy, but a window-split parallel run brings it to
  about a day. Worth building.
- If tau grows exponentially, N=600 is dead and no amount of kernel
  optimization saves it. The right response would be to abandon Wang-Landau
  for this ensemble, not to fund it harder.

The P4 transition is first-order-like. Taming the exponential barrier of a
first-order transition into a polynomial one is exactly what multicanonical
sampling is for -- but whether it succeeds on *this* ensemble is an empirical
question, and this script is the experiment that asks it.

eps*N is held at 12, the P7 design axis, so every N sits at the same point of
the phase diagram and the comparison is meaningful.

This is method development. Nothing is preregistered, and no beta_c, G(beta),
or phase claim is produced or implied.

Usage:
    python p7_tunneling_scaling.py --sizes 30 40 50 60
    python p7_tunneling_scaling.py --sizes 30 40 50 60 80 --target-round-trips 6
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
    wang_landau_2d_order,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "positive_control"

EPS_N = 12.0

# Convergence to ln_f = 1e-5 needs ~17 halvings, and the N=60 pilot showed the
# schedule spends exactly one round trip per halving.
HALVINGS_TO_CONVERGE = 17


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[30, 40, 50, 60])
    parser.add_argument("--eps-n", type=float, default=EPS_N)
    parser.add_argument(
        "--target-round-trips",
        type=int,
        default=6,
        help="stop each N once this many traversals are measured",
    )
    parser.add_argument(
        "--move-budget",
        type=int,
        default=60_000_000,
        help="per-N cap; an N that hits it reports a LOWER BOUND on tau",
    )
    parser.add_argument("--bins", type=int, default=60)
    parser.add_argument("--sweep-steps", type=int, default=8000)
    parser.add_argument("--seed", type=int, default=20260715)
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []

    for n in args.sizes:
        eps = args.eps_n / n
        rng = np.random.default_rng(args.seed + n)
        pi0 = rng.permutation(n)
        s_min, s_max = action_range(n, eps, seed=args.seed + n, probes=200)

        started = time.perf_counter()
        result = wang_landau_2d_order(
            pi0=pi0,
            eps=eps,
            s_min=s_min,
            s_max=s_max,
            n_bins=args.bins,
            seed=args.seed + n,
            sweep_steps=args.sweep_steps,
            max_sweeps=max(1, args.move_budget // args.sweep_steps),
            max_round_trips=args.target_round_trips,
            ln_f_final=1e-12,  # never converge; we are timing traversals only
        )
        elapsed = time.perf_counter() - started

        hit_budget = result.round_trips < args.target_round_trips
        if result.round_trips > 0:
            tau = result.moves / result.round_trips
        else:
            # No traversal at all inside the budget: tau is not measured, only
            # bounded below by the whole budget. Reporting moves/0 as a number
            # would invent a measurement that did not happen.
            tau = float("nan")

        rows.append({
            "n_elements": n,
            "eps": eps,
            "eps_n": eps * n,
            "window_low": s_min,
            "window_high": s_max,
            "moves": result.moves,
            "round_trips": result.round_trips,
            "tau_moves_per_round_trip": tau,
            # A run stopped by the move budget has not measured tau -- it has
            # only shown tau is at least this large. Mixing the two into one
            # column and fitting through them would manufacture an exponent.
            "tau_is_lower_bound": bool(hit_budget),
            "wl_acceptance": result.acceptance,
            "seconds": elapsed,
            "us_per_move": (
                1e6 * elapsed / result.moves if result.moves else float("nan")
            ),
            "projected_moves_to_converge": (
                tau * HALVINGS_TO_CONVERGE if np.isfinite(tau) else float("nan")
            ),
        })

        status = "LOWER BOUND (hit move budget)" if hit_budget else "measured"
        print(
            f"N={n:4d} eps={eps:.4f}  round_trips={result.round_trips:3d}  "
            f"moves={result.moves:,}  tau={tau:,.0f} moves/RT  [{status}]  "
            f"{elapsed:.0f}s",
            flush=True,
        )

    clean = [r for r in rows if not r["tau_is_lower_bound"] and r["round_trips"] > 0]
    exponent = None
    if len(clean) >= 2:
        # tau ~ N^alpha  =>  log tau = alpha log N + c
        log_n = np.log([r["n_elements"] for r in clean])
        log_tau = np.log([r["tau_moves_per_round_trip"] for r in clean])
        exponent, intercept = np.polyfit(log_n, log_tau, 1)
        exponent = float(exponent)

        print(f"\ntau ~ N^{exponent:.2f}  (fit over {len(clean)} measured sizes)")
        projected = float(
            np.exp(intercept) * 600.0**exponent * HALVINGS_TO_CONVERGE
        )
        print(
            f"extrapolated cost of a converged ln_g at N=600: "
            f"{projected:.3e} moves "
            f"({projected * 135e-6 / 3600:.0f} core-hours at the Numba kernel's "
            f"135 us/move)"
        )
        print(
            "\nThis extrapolation is a power-law fit over a short lever arm of "
            "small N. It is a planning number, not a result: if the true "
            "scaling is exponential, a power law will fit these few points and "
            "still be badly wrong at N=600."
        )
    else:
        print(
            "\nToo few measured sizes to fit an exponent. Raise --move-budget "
            "or lower --target-round-trips; do not fit through lower bounds."
        )

    with open(OUT / "p7_tunneling_scaling.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (OUT / "p7_tunneling_scaling_summary.json").write_text(
        json.dumps(
            {
                "eps_n": args.eps_n,
                "target_round_trips": args.target_round_trips,
                "move_budget": args.move_budget,
                "halvings_to_converge": HALVINGS_TO_CONVERGE,
                "fitted_exponent": exponent,
                "measured_sizes": [r["n_elements"] for r in clean],
                "rows": rows,
                "status": "method development; nothing frozen; no phase claim",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"\nwrote {OUT / 'p7_tunneling_scaling.csv'}", flush=True)


if __name__ == "__main__":
    main()
