# ruff: noqa
# Referee's scratch script for ../transport_referee_report.md, archived verbatim
# apart from these three lines. Absolute paths point to the worktree it ran in.
import numpy as np
rng = np.random.default_rng(3)
for n in (32, 64):
    cnt = []
    for _ in range(4000):
        perm = rng.permutation(n)  # v-rank of the event at u-rank k
        # incomparable twins: adjacent in u-rank, v-ranks adjacent and decreasing
        cnt.append(np.sum(perm[1:] - perm[:-1] == -1))
    cnt = np.array(cnt)
    print(n, "mean incomparable twin pairs per order", cnt.mean().round(3), "orders with none", np.mean(cnt == 0).round(3))
