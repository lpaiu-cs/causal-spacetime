# ruff: noqa
# Referee's scratch script for ../transport_referee_report.md, archived verbatim
# apart from these three lines. Absolute paths point to the worktree it ran in.
"""Causal character of interior pairs per bin (chain-free, 3+1)."""
import sys, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from b_chainfree import fp, dp, EDGES  # noqa

rng = np.random.default_rng(5)
for N, configs in ((128, 12), (512, 4)):
    tab = np.zeros((4, 3))
    for _ in range(configs):
        xa = fp.sample(fp.BACKGROUND, N, rng)
        x = xa[:N]
        rel = fp.order_matrix(xa)[:N, :N]
        central = dp.central_mask(x)
        outer = ~central & rel.any(axis=0) & rel.any(axis=1)
        xi = x[central | outer]; R = rel[np.ix_(central | outer, central | outer)]
        n = len(xi)
        dx = xi[None] - xi[:, None]
        d = np.linalg.norm(dx, axis=2)
        tau2 = dx[..., 0] ** 2 - np.sum(dx[..., 1:] ** 2, axis=2)
        nearnull = np.abs(tau2) < 0.1 * d ** 2  # |interval| small relative to Euclidean d^2
        related = R | R.T
        up = np.triu(np.ones((n, n), bool), 1)
        for b, (lo, hi) in enumerate(zip(EDGES[:-1], EDGES[1:])):
            p = up & (d >= lo) & (d < hi)
            tab[b] += [np.sum(p & related & ~nearnull), np.sum(p & nearnull), np.sum(p & ~related & ~nearnull)]
    tab /= tab.sum(axis=1, keepdims=True)
    print(f"N={N}: per bin fractions [timelike-not-near-null, near-null(|tau^2|<0.1 d^2), spacelike-not-near-null]")
    for b in range(4):
        print("  ", EDGES[b], np.round(tab[b], 3))
