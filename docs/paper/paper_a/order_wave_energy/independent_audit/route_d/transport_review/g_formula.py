# ruff: noqa
# Referee's scratch script for ../transport_referee_report.md, archived verbatim
# apart from these three lines. Absolute paths point to the worktree it ran in.
"""Direct-sampling check of E S_ij = [pu(1-pu)+pv(1-pv)]/(N+2) in FLAT2 (no chain)."""
import math, sys, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from d_c3 import fp, uv, exact_S, L  # noqa

rng = np.random.default_rng(11)
for n, draws in ((32, 200000), (64, 50000)):
    xa = fp.sample(fp.FLAT2, n, rng)[:n]
    u, v = uv(xa)
    ru = np.argsort(np.argsort(u)); rv = np.argsort(np.argsort(v))
    E = exact_S(xa)
    up = np.triu(np.ones((n, n), bool), 1)
    tot = np.zeros((n, n)); tot2 = 0.0
    batch = 2000
    sums = []
    for _ in range(draws // batch):
        us = np.sort(rng.random((2, batch, n)) * L, axis=2)[..., ru]
        vs = np.sort(rng.random((2, batch, n)) * L, axis=2)[..., rv]
        du = us[0] - us[1]; dv = vs[0] - vs[1]
        S = (du[:, :, None] - du[:, None, :]) ** 2 + (dv[:, :, None] - dv[:, None, :]) ** 2
        tot += S.sum(axis=0)
        sums.append(S[:, up].sum(axis=1) / E[up].sum())
    tot /= draws
    s = np.concatenate(sums)
    far = up & (np.linalg.norm(xa[:, None] - xa[None], axis=2) >= 0.5)
    print(f"N={n}, {draws} draws: sum S / sum E = {tot[up].sum()/E[up].sum():.4f} (SE {s.std()/math.sqrt(len(s)):.4f}); "
          f"far {tot[far].sum()/E[far].sum():.4f}; per-pair ratio range {np.min(tot[up]/E[up]):.3f}..{np.max(tot[up]/E[up]):.3f}; "
          f"all-event MSD/(1/(3(N+1))) {np.mean((us[0]-us[1])**2 + (vs[0]-vs[1])**2) * 3 * (n+1):.3f}")
