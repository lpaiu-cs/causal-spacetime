# ruff: noqa
# Referee's scratch script for ../transport_referee_report.md, archived verbatim
# apart from these three lines. Absolute paths point to the worktree it ran in.
"""Mini C2 (chain-free): group element + Brownian field + noise, nonlinear fit
over interior events, against the linear-projection expectation."""
import math, sys, pathlib, time
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from b_chainfree import fp, dp, proj_complement, bin_stats, brownian_cov, EDGES, O  # noqa


def binned(r, xa):
    n = len(xa)
    d = np.linalg.norm(xa[:, None] - xa[None], axis=2)
    S = np.sum((r[:, None] - r[None]) ** 2, axis=2)
    v = np.sum(r ** 2, axis=1)
    I = v[:, None] + v[None]
    up = np.triu(np.ones((n, n), bool), 1)
    out = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        p = up & (d >= lo) & (d < hi)
        out.append((S[p].mean(), I[p].mean()))
    return out


def run(N, configs, draws, kappa, rng):
    rows = []
    for c in range(configs):
        xa = fp.sample(fp.BACKGROUND, N, rng)
        x = xa[:N]
        rel = fp.order_matrix(xa)[:N, :N]
        central = dp.central_mask(x)
        outer = ~central & rel.any(axis=0) & rel.any(axis=1)
        xi = x[central | outer]
        n = len(xi)
        cen_i = central[central | outer]
        M = proj_complement(xi)
        K4, k = brownian_cov(xi, kappa, N)
        L = np.linalg.cholesky(k + 1e-15 * np.eye(n))
        sigma = 10.0 / N / 2.0
        noise = 2 * sigma ** 2
        exp_tot = bin_stats(M @ (K4 + noise * np.eye(4 * n)) @ M.T, xi)
        exp_field = bin_stats(M @ K4 @ M.T, xi)
        for _ in range(draws):
            direction = rng.normal(size=7); direction /= np.linalg.norm(direction)
            scale = 0.1
            for _ in range(6):
                moved = dp.group_map(scale * direction, xi)
                rms = math.sqrt(np.mean(np.sum((moved - xi)[cen_i] ** 2, axis=1)))
                scale *= (3.2 / math.sqrt(N)) / rms
            truth = scale * direction
            gx = dp.group_map(truth, xi)
            tau = L @ rng.normal(size=(n, 4))
            na = sigma * rng.normal(size=xi.shape)
            nb = sigma * rng.normal(size=xi.shape)
            a = xi + na
            out = {}
            for label, b in (("with", gx + tau + nb), ("noise", gx + nb), ("field_only", gx + tau)):
                _, res, ok = dp.fit_group(a if label != "field_only" else xi, b)
                out[label] = binned(res, a if label != "field_only" else xi)
                out[label + "_ok"] = ok
            # also: field with no group element (projection at x_a is then exact in tangent)
            _, res0, _ = dp.fit_group(xi, xi + tau)
            out["field_nogroup"] = binned(res0, xi)
            rows.append(dict(exp_tot=exp_tot, exp_field=exp_field, **out))
    return rows


def summarise(N, kappa, rows):
    print(f"N={N} kappa={kappa} draws={len(rows)}")
    for bi, name in ((0, "near"), (3, "far")):
        ratio = np.array([r["with"][bi][0] / r["exp_tot"][bi][0] for r in rows])
        field_part = np.array([(r["with"][bi][0] - r["noise"][bi][0]) / r["exp_field"][bi][0] for r in rows])
        fo = np.array([r["field_only"][bi][0] / r["exp_field"][bi][0] for r in rows])
        f0 = np.array([r["field_nogroup"][bi][0] / r["exp_field"][bi][0] for r in rows])
        fieldshare = np.mean([r["exp_field"][bi][0] / r["exp_tot"][bi][0] for r in rows])
        Iratio = np.array([r["with"][bi][1] / r["exp_tot"][bi][1] for r in rows])
        print(f"  {name}: S rec/exp mean {ratio.mean():.4f} sd {ratio.std(ddof=1):.4f} | I rec/exp {Iratio.mean():.4f} sd {Iratio.std(ddof=1):.4f}"
              f" | field share of expected S {fieldshare:.3f}"
              f" | field part (with-noise)/exp mean {field_part.mean():.3f} sd {field_part.std(ddof=1):.3f}"
              f" | field only, group: {fo.mean():.3f} sd {fo.std(ddof=1):.3f}; no group: {f0.mean():.3f} sd {f0.std(ddof=1):.3f}")
    Y = np.array([(r["with"][3][0] / r["with"][3][1]) / (r["with"][0][0] / r["with"][0][1]) for r in rows])
    print(f"  Y per draw mean {Y.mean():.4f} sd {Y.std(ddof=1):.4f}; fits ok {all(r['with_ok'] for r in rows)}")


if __name__ == "__main__":
    rng = np.random.default_rng(777)
    t0 = time.time()
    for N, configs, draws in ((128, 12, 4), (256, 6, 4), (512, 3, 3)):
        for kappa in (1.0, 0.3):
            rows = run(N, configs, draws, kappa, rng)
            summarise(N, kappa, rows)
            print(f"  elapsed {time.time()-t0:.0f}s", flush=True)
