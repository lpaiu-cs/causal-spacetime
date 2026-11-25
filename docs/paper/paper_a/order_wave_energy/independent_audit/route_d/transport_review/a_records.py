# ruff: noqa
# Referee's scratch script for ../transport_referee_report.md, archived verbatim
# apart from these three lines. Absolute paths point to the worktree it ran in.
"""Read existing Route D records (read-only) and recompute sizing inputs."""
import json, math, pathlib
import numpy as np
from scipy import stats

ROOT = pathlib.Path(r"E:/lab/causal-spacetime/.claude/worktrees/upbeat-newton-3e92fe/docs/paper/paper_a/order_wave_energy/independent_audit/route_d")

def load(p):
    return [json.loads(l) for l in open(p, encoding="utf-8")]

main = load(ROOT / "evidence/main/records.jsonl")
ext = load(ROOT / "evidence_extension/main/records.jsonl")
ctrl = load(ROOT / "evidence/control/records.jsonl")

# 1+1 control: per-order ratio at lag N^2 = 4096 (N=64)
vals = {}
secs = []
for r in ctrl:
    if not r["file"].startswith("chain"):
        continue
    secs.append(r["seconds"])
    for lag, pair in r["ratio_by_lag"].items():
        vals.setdefault(int(lag), []).append(np.mean(pair))
print("1+1 N=64 seconds per order", np.mean(secs), min(secs), max(secs))
for lag in sorted(vals):
    v = np.array(vals[lag])
    print(f"1+1 lag {lag/4096:.2f} N^2: mean {v.mean():.3f} sd {v.std(ddof=1):.3f} rel {v.std(ddof=1)/v.mean():.3f} n {len(v)}")

def per_order(records, lagfun):
    out = {}
    for r in records:
        if "meta" not in r:
            continue
        n = r["meta"]["n"]
        lag = lagfun(n)
        rows = [w for w in r["rows"] if w["lag"] == lag]
        if not rows:
            continue
        s1 = r["meta"]["s1_squared"]
        bins = np.array([[np.nan if b is None else b for b in w["pair_difference_bins"]] for w in rows])
        far = np.nanmean(bins[:, 2])
        near = np.nanmean(bins[:, 0])
        res = np.mean([w["central_residual"] for w in rows])
        resw = np.mean([w["central_residual_wave"] for w in rows])
        outer = np.mean([w["outer_residual"] for w in rows])
        out.setdefault(n, []).append(dict(far=far, near=near, s1=s1, res=res, resw=resw, outer=outer,
                                          nc=r["meta"]["central_events"], no=r["meta"]["outer_events"],
                                          sec=r["meta"]["seconds"]))
    return out

for label, lagf in (("N^2/2", lambda n: n * n // 2), ("N^2", lambda n: n * n), ("2N^2", lambda n: 2 * n * n)):
    d = per_order(main, lagf)
    d.update(per_order(ext, lagf) if label == "N^2/2" else {})
    print("lag", label)
    for n in sorted(d):
        rows = d[n]
        far = np.array([x["far"] for x in rows])
        s1 = np.array([x["s1"] for x in rows])
        print(f"  N={n} orders={len(rows)} far(>0.5) N^2: mean {np.mean(far*n*n):.1f} sd {np.std(far*n*n, ddof=1):.1f}; "
              f"far/s1^2 mean {np.mean(far/s1):.1f}; near N^2 {np.mean([x['near'] for x in rows])*n*n:.1f}; "
              f"res N^2 {np.mean([x['res'] for x in rows])*n*n:.1f} sd {np.std([x['res']*n*n for x in rows], ddof=1):.1f}; "
              f"resw N^2 {np.mean([x['resw'] for x in rows])*n*n:.1f}; outer N^2 {np.mean([x['outer'] for x in rows])*n*n:.1f} sd {np.std([x['outer']*n*n for x in rows], ddof=1):.1f}; "
              f"nc {np.mean([x['nc'] for x in rows]):.1f} no {np.mean([x['no'] for x in rows]):.1f} sec {np.mean([x['sec'] for x in rows]):.0f}")

# Projection of the draft: SE(B) for S_far N^2 = A + B N with 24,24,12 orders
def se_B(sizes, counts, sds):
    # OLS on per-order points with common variance estimated by pooled sd
    x = np.concatenate([[n] * c for n, c in zip(sizes, counts)]).astype(float)
    pooled = math.sqrt(sum((c - 1) * s * s for c, s in zip(counts, sds)) / (sum(counts) - len(sizes)))
    sxx = np.sum((x - x.mean()) ** 2)
    se_common = pooled / math.sqrt(sxx)
    # true SE of OLS slope under heteroscedastic variance
    w = (x - x.mean()) / sxx
    var = np.sum(w ** 2 * np.concatenate([[s * s] * c for s, c in zip(sds, counts)]))
    return se_common, math.sqrt(var), pooled

for sds in ((90, 90, 37), (150, 150, 37), (120, 120, 37)):
    print("sds", sds, "SE(B) common, hetero, pooled:", se_B((128, 256, 512), (24, 24, 12), sds))
