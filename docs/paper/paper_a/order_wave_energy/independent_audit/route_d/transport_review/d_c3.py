# ruff: noqa
# Referee's scratch script for ../transport_referee_report.md, archived verbatim
# apart from these three lines. Absolute paths point to the worktree it ran in.
"""C3 feasibility in 1+1 (FLAT2) at N = 32, 16 orders; exact formula check."""
import math, sys, pathlib, time
import numpy as np

AUD = pathlib.Path(r"E:/lab/causal-spacetime/.claude/worktrees/upbeat-newton-3e92fe/docs/paper/paper_a/order_wave_energy/independent_audit")
sys.path.insert(0, str(AUD / "route_d"))
sys.path.insert(0, str(AUD / "route_a"))
import flip_probe as fp  # noqa
import drift_probe as dp  # noqa

L = 1 / math.sqrt(2)
BINS_A = (0.0, 0.25, 0.5, 0.75, np.inf)
H3 = 1.4795  # 3+1 tip-to-tip time extent, for the scaled reading
BINS_B = tuple(b / H3 for b in BINS_A)


def uv(x):
    return (x[:, 0] + x[:, 1]) / math.sqrt(2), (x[:, 0] - x[:, 1]) / math.sqrt(2)


def exact_S(x0):
    n = len(x0)
    u, v = uv(x0)
    ru = np.argsort(np.argsort(u)); rv = np.argsort(np.argsort(v))
    pu = np.abs(ru[:, None] - ru[None]) / (n + 1)
    pv = np.abs(rv[:, None] - rv[None]) / (n + 1)
    return (pu * (1 - pu) + pv * (1 - pv)) / (n + 2)


def check_formula(rng, n=32, draws=4000):
    xa = fp.sample(fp.FLAT2, n, rng)[:n]
    u, v = uv(xa)
    ru = np.argsort(np.argsort(u)); rv = np.argsort(np.argsort(v))
    E = exact_S(xa)
    acc = np.zeros((n, n))
    for _ in range(draws):
        pos = []
        for _ in range(2):
            us = np.sort(rng.random(n) * L)[ru]
            vs = np.sort(rng.random(n) * L)[rv]
            pos.append(np.column_stack(((us + vs) / math.sqrt(2), (us - vs) / math.sqrt(2))))
        r = pos[0] - pos[1]
        acc += np.sum((r[:, None] - r[None]) ** 2, axis=2)
    acc /= draws
    up = np.triu(np.ones((n, n), bool), 1)
    ratio = acc[up].sum() / E[up].sum()
    far = up & (np.linalg.norm(xa[:, None] - xa[None], axis=2) >= 0.5)
    print(f"formula check N={n}: MC/exact all pairs {ratio:.4f}; far(d>=0.5) {acc[far].sum()/E[far].sum():.4f}; "
          f"max rel dev per pair {np.max(np.abs(acc[up]/E[up]-1)):.3f}")


def chain_states(xa, rng, times):
    n = len(xa) - 2
    x = xa.copy()
    idle = -np.ones(n, dtype=np.int64)
    tally = {"skipped": 0, "outside": 0}
    rel = fp.order_matrix(xa)
    states, now = {0: xa[:n].copy()}, 0
    import hit_probe_v3 as hv
    for target in sorted(set(times)):
        hv.run_block(x, fp.FLAT2, rng, target - now, idle, None, None, None, None, tally)
        now = target
        if not np.array_equal(fp.order_matrix(x), rel):
            raise RuntimeError("left the order cell")
        states[target] = x[:n].copy()
    return states, tally["skipped"]


def readout(states, lag, mask, E, x0, bins):
    xs = x0[mask]
    d = np.linalg.norm(xs[:, None] - xs[None], axis=2)
    n = len(xs)
    up = np.triu(np.ones((n, n), bool), 1)
    Em = E[np.ix_(mask, mask)]
    S = np.zeros((n, n)); I = np.zeros((n, n)); k = 0
    for s in (0, lag):
        if s + lag not in states:
            continue
        r = (states[s] - states[s + lag])[mask]
        S += np.sum((r[:, None] - r[None]) ** 2, axis=2)
        v = np.sum(r ** 2, axis=1); I += v[:, None] + v[None]
        k += 1
    S /= k; I /= k
    out = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        p = up & (d >= lo) & (d < hi)
        out.append((S[p].sum() / Em[p].sum() if p.any() else np.nan, S[p].mean() / I[p].mean() if p.any() else np.nan, p.sum()))
    return out


def main():
    rng = np.random.default_rng(4242)
    check_formula(rng)
    N, orders = 32, 16
    lags = [N * N // 4, N * N // 2, N * N, 2 * N * N, 4 * N * N]
    times = sorted({t for l in lags for t in (l, 2 * l)})
    res = {lag: [] for lag in lags}
    t0 = time.time()
    counts_pairs = []
    for o in range(orders):
        rr = np.random.default_rng(np.random.SeedSequence([99991, 2, N, o]))
        xa = fp.sample(fp.FLAT2, N, rr)
        rel = fp.order_matrix(xa)[:N, :N]
        interior = rel.any(axis=0) & rel.any(axis=1)
        E = exact_S(xa[:N])
        states, skipped = chain_states(xa, rr, times)
        allm = np.ones(N, bool)
        for lag in lags:
            res[lag].append((readout(states, lag, allm, E, xa[:N], BINS_A), readout(states, lag, interior, E, xa[:N], BINS_A),
                             readout(states, lag, allm, E, xa[:N], BINS_B)))
        counts_pairs.append((interior.sum(),))
    print(f"chains: {orders} orders at N={N}, {max(times)} sweeps, {time.time()-t0:.0f}s total")
    print("mean interior events", np.mean([c[0] for c in counts_pairs]))
    for lag in lags:
        for vi, label in ((0, "all, bins as 3+1"), (1, "interior, bins as 3+1"), (2, "all, bins/1.48")):
            rows = [r[vi] for r in res[lag]]
            Fn = np.array([r[0][0] for r in rows]); Ff = np.array([r[3][0] for r in rows])
            Xn = np.array([r[0][1] for r in rows]); Xf = np.array([r[3][1] for r in rows])
            Y = Xf / Xn
            npairs = np.array([[r[b][2] for b in range(4)] for r in rows]).mean(axis=0)
            def ci(a):
                a = a[np.isfinite(a)]
                return f"{a.mean():.3f} sd {a.std(ddof=1):.3f} (n={len(a)})"
            print(f"lag {lag/(N*N):.2f}N^2 [{label}] pairs/bin {np.round(npairs,1)} F_near {ci(Fn)} F_far {ci(Ff)} X_near {ci(Xn)} X_far {ci(Xf)} Y {ci(Y)}")


if __name__ == "__main__":
    main()
