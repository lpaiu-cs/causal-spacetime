# ruff: noqa
# Referee's scratch script for ../transport_referee_report.md, archived verbatim
# apart from these three lines. Absolute paths point to the worktree it ran in.
import sys, pathlib, math
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from b_chainfree import fp, dp, EDGES
rng = np.random.default_rng(8)
for N, configs in ((128, 12), (512, 4)):
    acc = [[] for _ in range(4)]; accd = [[] for _ in range(4)]
    for _ in range(configs):
        xa = fp.sample(fp.BACKGROUND, N, rng); x = xa[:N]
        rel = fp.order_matrix(xa)[:N, :N]
        c = dp.central_mask(x); o = ~c & rel.any(axis=0) & rel.any(axis=1)
        xi = x[c | o]; n = len(xi)
        m = (xi[:, 0] + xi[:, 1]) / math.sqrt(2)
        d = np.linalg.norm(xi[:, None] - xi[None], axis=2)
        dm = np.abs(m[:, None] - m[None])
        up = np.triu(np.ones((n, n), bool), 1)
        for b, (lo, hi) in enumerate(zip(EDGES[:-1], EDGES[1:])):
            p = up & (d >= lo) & (d < hi)
            acc[b].append(dm[p].mean()); accd[b].append(d[p].mean())
    print(N, "mean |dm| per bin", [round(np.mean(a), 3) for a in acc], "mean d per bin", [round(np.mean(a), 3) for a in accd])
