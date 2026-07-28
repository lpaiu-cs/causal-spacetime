"""Is the twin's nonzero contrast a real artifact or n=21 noise?
Same code path as the campaign (run_sample flat=True), high n,
design-check space 7600000+."""
import sys

import numpy as np

sys.path.insert(0, "experiments/positive_control")
import p13_tau_ell as P  # noqa: E402

NS = 150
ys = {}
for tau_c in P.RUNGS:
    vals = []
    for k in range(NS):
        r, ok = P.run_sample(tau_c, 7_600_000 + 100_000 * P.RUNGS.index(tau_c)
                             + 200 * k, flat=True)
        if ok:
            vals.append(r["y"])
    ys[tau_c] = np.array(vals)
    a = ys[tau_c]
    print(f"twin tau/ell={tau_c}: n {a.size} | mean y {a.mean():+.4f}"
          f" +/- {a.std(ddof=1)/np.sqrt(a.size):.4f} | s(y) {a.std(ddof=1):.4f}",
          flush=True)

b, t = ys[P.RUNGS[0]], ys[P.RUNGS[-1]]
d = t.mean() - b.mean()
se = np.sqrt(b.var(ddof=1)/b.size + t.var(ddof=1)/t.size)
print(f"\ntwin Delta (150/rung) = {d:+.4f} +/- {1.96*se:.4f}"
      f"   [campaign at n=21 read +0.0848, CI (0.032, 0.136)]")
# what does an n=21 subsample of these look like?
rng = np.random.default_rng(11)
sub = [float(rng.choice(t, 21, replace=False).mean()
             - rng.choice(b, 21, replace=False).mean()) for _ in range(4000)]
sub = np.array(sub)
print(f"n=21 sampling distribution of the twin contrast: "
      f"sd {sub.std():.4f}, P(|Delta| > 0.05) = {np.mean(np.abs(sub) > 0.05):.3f},"
      f" P(Delta > 0.0848) = {np.mean(sub > 0.0848):.3f}")
