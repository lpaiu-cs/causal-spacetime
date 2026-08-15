"""O3' cost projection -- two independently measured curves, several
fit windows, one range. PROJECTIONS, NOT CERTIFICATIONS.

The O3' caps (900,000 calls / 24 h) are headroom above the range this
module computes, and a contract test recomputes the whole table, so
the freeze document's numbers cannot drift from the committed
artifacts they are fitted to.

Two data sources, fitted SEPARATELY (PI ruling):

  * the NEIGHBOR price ladder (`p14_oracle_price.json`): the four
    measured target crossings plus the final measured point of the
    same adaptive pass, on the (12, 18, 8.0)M neighbor anchors;
  * the executed O3 run's OWN refinement curve
    (`p14_o3_volume.json["curve"]`, 1,078 samples on the frozen
    (12, 18, 8.5)M anchors).

Model: calls ~ C0 * (1/ratio)^p, fitted by OLS on the log-log points
of each window. The exponent is NOT assumed; the v1 briefing's
p = 2 "cell ~ h^-2" scaling is retained only as the central
projection, bounded by the measured-window fits. The neighbor fits
are anchored analytically through the O3 executed point
(137,958 calls at ratio 0.009997) rather than through the neighbor
intercept, because the frozen diamond is a different (larger)
integrand; the O3-curve fits carry their own intercepts.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
_PRICE = _REPO / "docs" / "prereg" / "p14_oracle_price.json"
_O3 = _REPO / "docs" / "prereg" / "p14_o3_volume.json"

#: The executed O3 anchor point every neighbor fit is projected
#: through.
O3_CALLS, O3_RATIO, O3_WALL_H = 137_958, 0.009997, 3.72

#: O3-curve fit windows: keep samples with ratio <= bound (None =
#: the full curve).
O3_WINDOWS = (("full curve", None), ("ratio<=0.05", 0.05),
              ("ratio<=0.03", 0.03), ("ratio<=0.02", 0.02))

#: Neighbor-ladder fit windows: suffixes of the measured crossings.
LADDER_WINDOWS = (("all 5", 0), ("last 4", 1), ("last 3", 2),
                  ("last 2", 3))


def _fit_p(points: list[tuple[float, int]]) -> float:
    """OLS slope of ln(calls) on ln(1/ratio)."""

    xs = [math.log(1.0 / r) for r, _ in points]
    ys = [math.log(c) for _, c in points]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True))
    return sxy / sxx


def ladder_points() -> list[tuple[float, int]]:
    d = json.loads(_PRICE.read_text(encoding="utf-8"))
    pts = [(c["ratio"], c["calls"]) for c in d["measured_crossings"]]
    pts.append((d["final"]["ratio"], d["final"]["calls"]))
    return pts


def o3_curve_points() -> list[tuple[float, int]]:
    d = json.loads(_O3.read_text(encoding="utf-8"))
    return [(c["ratio"], c["calls"]) for c in d["curve"]]


def project(p: float, r_target: float,
            anchor: tuple[float, int] | None = None,
            intercept: float | None = None) -> float:
    """Projected calls at `r_target` under exponent `p`, either
    through an (ratio, calls) anchor or with a fitted intercept."""

    if anchor is not None:
        r0, c0 = anchor
        return c0 * (r0 / r_target) ** p
    return intercept * (1.0 / r_target) ** p


def table(r_target: float = 0.005) -> dict:
    """The whole projection table for one target, as data. The range
    is the min/max over every fitted window of both curves; the
    central projection is the v1 p = 2 scaling through the O3 point."""

    rows = []
    lad = ladder_points()
    for name, skip in LADDER_WINDOWS:
        p = _fit_p(lad[skip:])
        calls = project(p, r_target, anchor=(O3_RATIO, O3_CALLS))
        rows.append({"source": "neighbor ladder", "window": name,
                     "p": p, "calls": calls})
    curve = o3_curve_points()
    for name, bound in O3_WINDOWS:
        pts = curve if bound is None else [(r, c) for r, c in curve
                                           if r <= bound]
        xs = [math.log(1.0 / r) for r, _ in pts]
        ys = [math.log(c) for _, c in pts]
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        p = _fit_p(pts)
        c0 = math.exp(my - p * mx)
        calls = project(p, r_target, intercept=c0)
        rows.append({"source": "O3 own curve", "window": name,
                     "p": p, "calls": calls, "n_samples": n})
    central = project(2.0, r_target, anchor=(O3_RATIO, O3_CALLS))
    lo = min(r["calls"] for r in rows)
    hi = max(r["calls"] for r in rows)
    return {
        "r_target": r_target,
        "rows": rows,
        "central_projection_calls": central,
        "range_calls": [lo, hi],
        "wall_h_at_o3_rate": [lo * O3_WALL_H / O3_CALLS,
                              hi * O3_WALL_H / O3_CALLS],
        "grade": "projection, not certification",
    }


def main() -> None:
    t = table()
    print(f"O3' projection at r_target = {t['r_target']}:")
    for row in t["rows"]:
        extra = (f" (n={row['n_samples']})"
                 if "n_samples" in row else "")
        print(f"  {row['source']:16s} {row['window']:12s} "
              f"p={row['p']:.3f} -> ~{row['calls']:,.0f} calls{extra}")
    lo, hi = t["range_calls"]
    wl, wh = t["wall_h_at_o3_rate"]
    print(f"  central (p=2 through the O3 point): "
          f"~{t['central_projection_calls']:,.0f} calls")
    print(f"  RANGE: ~{lo:,.0f} - {hi:,.0f} calls "
          f"(~{wl:.1f} - {wh:.1f} h at the O3 call rate)")
    print(f"  frozen caps: 900,000 calls / 24 h "
          f"(+{(900_000 / hi - 1) * 100:.0f}% over the worst window); "
          f"{t['grade']}")


if __name__ == "__main__":
    main()
