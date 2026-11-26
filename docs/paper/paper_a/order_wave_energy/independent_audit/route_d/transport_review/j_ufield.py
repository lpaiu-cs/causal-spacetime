# ruff: noqa
# Referee's scratch script for ../transport_referee_report.md, archived verbatim
# apart from these three lines. Absolute paths point to the worktree it ran in.
"""1-D count-like field along the wave coordinate: tau = W(m) e_m, E|W(m1)-W(m2)|^2 = 2|dm|/N (heuristic)."""
import sys, pathlib, math
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from b_chainfree import fp, dp, proj_complement, bin_stats, O
rng = np.random.default_rng(9)
e = np.array([1, 1, 0, 0]) / math.sqrt(2)
m0 = (O[0] + O[1]) / math.sqrt(2)
for N, configs in ((128, 12), (256, 8), (512, 4)):
    bc, sfar_pre, sfar_post, bc_iso = [], [], [], []
    for _ in range(configs):
        xa = fp.sample(fp.BACKGROUND, N, rng); x = xa[:N]
        rel = fp.order_matrix(xa)[:N, :N]
        c = dp.central_mask(x); o = ~c & rel.any(axis=0) & rel.any(axis=1)
        for mask, tag in ((c, "central"), (c | o, "interior")):
            xi = x[mask]; n = len(xi)
            m = (xi[:, 0] + xi[:, 1]) / math.sqrt(2) - m0
            # Brownian motion in m anchored at m0: cov = (k/2)(|m1|+|m2|-|m1-m2|) with variogram k|dm|, k = 2/N
            k = 2.0 / N
            km = (k / 2) * (np.abs(m)[:, None] + np.abs(m)[None] - np.abs(m[:, None] - m[None]))
            K = np.kron(km, np.outer(e, e))
            M = proj_complement(xi)
            post = M @ K @ M.T
            if tag == "central":
                C = post.reshape(n, 4, n, 4)
                wv = np.einsum("a,iaia,a->i", e, C, e)
                bc.append(wv.mean() * N * 4 * n / (4 * n - 7))
            else:
                b_post = bin_stats(post, xi, wave=True); b_pre = bin_stats(K, xi, wave=True)
                sfar_pre.append(b_pre[3][0] * N * N); sfar_post.append(b_post[3][0] * N * N)
    print(f"N={N}: central wave residual c_m^2 under heuristic 1-D u-field: {np.mean(bc):.4f} -> c_m = {math.sqrt(np.mean(bc)):.3f}; "
          f"interior far-bin S^m N^2 added: pre-fit {np.mean(sfar_pre):.1f}, post-fit {np.mean(sfar_post):.1f} (= {np.mean(sfar_post)/N:.3f} N)")
