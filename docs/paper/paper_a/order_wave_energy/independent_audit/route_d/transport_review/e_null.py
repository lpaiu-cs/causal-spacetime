# ruff: noqa
# Referee's scratch script for ../transport_referee_report.md, archived verbatim
# apart from these three lines. Absolute paths point to the worktree it ran in.
"""Chain-free null models in 3+1: (1) cell-scale short-range correlation,
(2) N-independent position-dependent cell variance + composition drift."""
import math, sys, pathlib
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from b_chainfree import fp, dp, proj_complement, bin_stats, scale_of  # noqa


def fit_line(ns, ys):
    ns = np.asarray(ns, float); ys = np.asarray(ys, float)
    B = np.sum((ns - ns.mean()) * (ys - ys.mean())) / np.sum((ns - ns.mean()) ** 2)
    return ys.mean() - B * ns.mean(), B


def main():
    rng = np.random.default_rng(31337)
    res = {}
    for N, configs in ((128, 12), (256, 10), (512, 6)):
        acc = {}
        for c in range(configs):
            xa = fp.sample(fp.BACKGROUND, N, rng)
            x = xa[:N]
            rel = fp.order_matrix(xa)[:N, :N]
            central = dp.central_mask(x)
            outer = ~central & rel.any(axis=0) & rel.any(axis=1)
            inter = central | outer
            xi = x[inter]; n = len(xi)
            s = scale_of(xi)
            isout = outer[inter]
            M = proj_complement(xi)
            d = np.linalg.norm(xi[:, None] - xi[None], axis=2)
            noise = 2 * (10.0 / N / 2) ** 2 * np.eye(4 * n)
            for ell in (0.1, 0.2):
                k = (25.0 / N ** 2) * np.exp(-d ** 2 / (2 * ell ** 2))
                tot = M @ (np.kron(k, np.eye(4)) + noise) @ M.T
                b = bin_stats(tot, xi)
                X = [q[0] / q[1] for q in b]
                acc.setdefault(f"Y_cell_ell{ell}", []).append(X[3] / X[0])
                acc.setdefault(f"Xnear_cell_ell{ell}", []).append(X[0])
                acc.setdefault(f"Sfar_cell_ell{ell}", []).append(b[3][0] * N * N)
                acc.setdefault(f"Cnear_cell_ell{ell}", []).append((b[0][1] - b[0][0]) / 2 * N * N)
            for gamma in (0.5, 1.0, 1.5):
                var = (1.0 / (1 - np.minimum(s, 0.98))) ** gamma
                D = np.kron(np.diag(var), np.eye(4))
                b = bin_stats(M @ D @ M.T, xi)
                acc.setdefault(f"Sfar_g{gamma}", []).append(b[3][0])
                acc.setdefault(f"Ifar_g{gamma}", []).append(b[3][1])
                acc.setdefault(f"cen_g{gamma}", []).append(var[~isout].mean())
                acc.setdefault(f"out_g{gamma}", []).append(var[isout].mean())
            acc.setdefault("s_outer_mean", []).append(s[isout].mean())
            acc.setdefault("s_outer_gt95", []).append(np.mean(s[isout] > 0.95))
        res[N] = {k: np.mean(v) for k, v in acc.items()}
        print(f"N={N}: " + "; ".join(f"{k} {v:.4f}" for k, v in res[N].items()), flush=True)
    ns = sorted(res)
    for gamma in (0.5, 1.0, 1.5):
        # calibrate: central mean variance = 155 / N^2 (residual, both states) -> scale
        sf = [res[n][f"Sfar_g{gamma}"] / res[n][f"cen_g{gamma}"] * 155 for n in ns]
        A, B = fit_line(ns, sf)
        ratio = [res[n][f"out_g{gamma}"] / res[n][f"cen_g{gamma}"] for n in ns]
        print(f"gamma {gamma}: S_far N^2 (cell-scale only, calibrated) {np.round(sf,1)}; A {A:.1f} B {B:.3f}; "
              f"spurious kappa = sqrt(B/0.87) = {math.sqrt(max(B,0)/0.87):.3f}; outer/central variance {np.round(ratio,2)}")
    for ell in (0.1, 0.2):
        print(f"ell {ell}: Y {[round(res[n][f'Y_cell_ell{ell}'],3) for n in ns]}, X_near {[round(res[n][f'Xnear_cell_ell{ell}'],3) for n in ns]}, "
              f"S_far N^2 {[round(res[n][f'Sfar_cell_ell{ell}'],1) for n in ns]}, C_near N^2 {[round(res[n][f'Cnear_cell_ell{ell}'],1) for n in ns]}")


if __name__ == "__main__":
    main()
