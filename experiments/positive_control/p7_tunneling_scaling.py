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
import os

# Each (N, seed) run is an independent Wang-Landau walk, so the sweep is
# embarrassingly parallel across cores. But the inner loop calls small NumPy
# matmuls, and if their BLAS spawns threads while a dozen worker processes are
# already running, the machine oversubscribes and every worker slows down. Pin
# BLAS to one thread per worker; parallelism comes from processes, not threads.
# This must run before NumPy is imported (in every spawned child too, which
# re-executes this module), so it sits at the very top.
for _var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_var, "1")

import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

from causal_spacetime_lab.positive_control.multicanonical import (  # noqa: E402
    action_range,
    wang_landau_2d_order,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "positive_control"

EPS_N = 12.0

# Convergence to ln_f = 1e-5 needs ~17 halvings, and the N=60 pilot showed the
# schedule spends exactly one round trip per halving.
HALVINGS_TO_CONVERGE = 17


def _measure_one(task: dict) -> dict:
    """One independent tau measurement, sized to run in a worker process.

    Returns a plain dict (picklable) with the move count, round trips, and the
    tau estimate -- or a NaN tau and a lower-bound flag if the walker never
    completed the target number of traversals inside the move budget.
    """

    n = task["n"]
    eps = task["eps"]
    seed = task["seed"]

    rng = np.random.default_rng(seed)
    pi0 = rng.permutation(n)
    s_min, s_max = action_range(n, eps, seed=seed, probes=200)

    started = time.perf_counter()
    result = wang_landau_2d_order(
        pi0=pi0,
        eps=eps,
        s_min=s_min,
        s_max=s_max,
        n_bins=task["bins"],
        seed=seed,
        sweep_steps=task["sweep_steps"],
        max_sweeps=max(1, task["move_budget"] // task["sweep_steps"]),
        max_round_trips=task["target_round_trips"],
        ln_f_final=1e-12,  # never converge; we are timing traversals only
    )
    elapsed = time.perf_counter() - started

    hit_budget = result.round_trips < task["target_round_trips"]
    if result.round_trips > 0:
        tau = result.moves / result.round_trips
    else:
        # No traversal inside the budget: tau is not measured, only bounded
        # below. Reporting moves/0 as a number would invent a measurement.
        tau = float("nan")

    return {
        "n_elements": n,
        "repeat": task["repeat"],
        "seed": seed,
        "eps": eps,
        "eps_n": eps * n,
        "moves": result.moves,
        "round_trips": result.round_trips,
        "tau_moves_per_round_trip": tau,
        # A run stopped by the move budget has not measured tau -- it has only
        # shown tau is at least this large. Mixing the two into one column and
        # fitting through them would manufacture an exponent.
        "tau_is_lower_bound": bool(hit_budget),
        "wl_acceptance": result.acceptance,
        "seconds": elapsed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[30, 40, 50, 60])
    parser.add_argument("--eps-n", type=float, default=EPS_N)
    parser.add_argument(
        "--target-round-trips",
        type=int,
        default=25,
        help="traversals per run; round-trip times are heavy-tailed, so a "
        "handful of them does not measure tau",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="independent seeds per N, so tau gets an error bar rather than a "
        "single noisy point",
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
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 2) - 4),
        help="parallel worker processes; each (N, seed) run is independent",
    )
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)

    tasks = [
        {
            "n": n,
            "repeat": repeat,
            "seed": args.seed + 1000 * repeat + n,
            "eps": args.eps_n / n,
            "bins": args.bins,
            "sweep_steps": args.sweep_steps,
            "move_budget": args.move_budget,
            "target_round_trips": args.target_round_trips,
        }
        for n in args.sizes
        for repeat in range(args.repeats)
    ]

    print(
        f"dispatching {len(tasks)} runs over {args.workers} workers "
        f"(sizes {args.sizes}, {args.repeats} seeds each)",
        flush=True,
    )

    # Results arrive out of order -- the small-N runs finish in seconds while a
    # heavy-tailed large-N draw can take an hour -- so print each as it lands
    # and sort only at the end.
    rows = []
    wall_start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(_measure_one, task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            status = "LOWER BOUND" if row["tau_is_lower_bound"] else "measured"
            tau = row["tau_moves_per_round_trip"]
            done = len(rows)
            print(
                f"[{done:2d}/{len(tasks)}] N={row['n_elements']:4d} "
                f"rep={row['repeat']}  round_trips={row['round_trips']:3d}  "
                f"moves={row['moves']:,}  tau={tau:,.0f}  [{status}]  "
                f"({row['seconds']:.0f}s in-worker)",
                flush=True,
            )

    rows.sort(key=lambda r: (r["n_elements"], r["repeat"]))
    wall = time.perf_counter() - wall_start
    print(f"\nall runs done in {wall:.0f}s wall", flush=True)

    clean = [r for r in rows if not r["tau_is_lower_bound"] and r["round_trips"] > 0]

    # Per-N spread first. The previous run fit an exponent through single noisy
    # points and produced tau(N=50) > tau(N=60) -- a power law cannot be
    # non-monotonic, so the scatter was larger than the signal between adjacent
    # sizes and the fitted exponent was meaningless. Print the spread before any
    # exponent, so that failure mode is visible rather than averaged away.
    print("\nper-N tau (moves per round trip):", flush=True)
    by_size: dict[int, list[float]] = {}
    for row in clean:
        by_size.setdefault(row["n_elements"], []).append(
            row["tau_moves_per_round_trip"]
        )
    for n in sorted(by_size):
        taus = np.array(by_size[n])
        spread = taus.max() / taus.min() if taus.min() > 0 else float("inf")
        print(
            f"  N={n:4d}  n_runs={taus.size}  mean={taus.mean():,.0f}  "
            f"min={taus.min():,.0f}  max={taus.max():,.0f}  "
            f"max/min={spread:.1f}x",
            flush=True,
        )

    exponent = None
    interval = None
    fittable = {n: v for n, v in by_size.items() if len(v) >= 1}
    if len(fittable) >= 3:
        # tau ~ N^alpha  =>  log tau = alpha log N + c. Fit on every run, not on
        # per-N means, so the within-N scatter propagates into the uncertainty.
        log_n = np.log([r["n_elements"] for r in clean])
        log_tau = np.log([r["tau_moves_per_round_trip"] for r in clean])
        exponent, intercept = (float(v) for v in np.polyfit(log_n, log_tau, 1))

        # Bootstrap the exponent over runs. A tight point estimate from noisy
        # data is exactly what produced the bogus 2.95 before.
        boot_rng = np.random.default_rng(0)
        boots = []
        for _ in range(2000):
            pick = boot_rng.integers(0, log_n.size, log_n.size)
            if np.unique(log_n[pick]).size < 2:
                continue
            boots.append(np.polyfit(log_n[pick], log_tau[pick], 1)[0])
        if boots:
            interval = (
                float(np.percentile(boots, 5)),
                float(np.percentile(boots, 95)),
            )

        print(f"\ntau ~ N^{exponent:.2f}", flush=True)
        if interval:
            print(
                f"90% bootstrap interval on the exponent: "
                f"[{interval[0]:.2f}, {interval[1]:.2f}]",
                flush=True,
            )

        low, high = interval if interval else (exponent, exponent)
        for label, alpha in (("low", low), ("point", exponent), ("high", high)):
            projected = float(
                np.exp(intercept) * 600.0**alpha * HALVINGS_TO_CONVERGE
            )
            print(
                f"  N=600 converged ln_g, alpha={alpha:.2f} ({label}): "
                f"{projected:.2e} moves = "
                f"{projected * 135e-6 / 3600:,.0f} core-hours "
                f"(Numba kernel at 135 us/move)",
                flush=True,
            )

        print(
            "\nThese are planning numbers, not results. The lever arm is short "
            "(a factor of two in N) and the extrapolation to N=600 is a factor "
            "of ten beyond it. A power law fits a few small-N points perfectly "
            "well even when the true scaling is exponential -- in which case "
            "every cost above is an underestimate, and not by a little.",
            flush=True,
        )
    else:
        print(
            "\nToo few measured sizes to fit an exponent. Raise --move-budget "
            "or lower --target-round-trips; do not fit through lower bounds.",
            flush=True,
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
                "repeats": args.repeats,
                "move_budget": args.move_budget,
                "halvings_to_converge": HALVINGS_TO_CONVERGE,
                "fitted_exponent": exponent,
                "exponent_90pct_bootstrap": interval,
                "measured_sizes": sorted({r["n_elements"] for r in clean}),
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
