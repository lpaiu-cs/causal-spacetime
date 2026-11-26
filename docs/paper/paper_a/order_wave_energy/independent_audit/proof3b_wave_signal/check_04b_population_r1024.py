# ruff: noqa: E501, UP031
"""Check 04b: one rung beyond the author's frozen population run (r = 1024).

The author's rows stop at r = 512 (0.53% from the limit).  This computes
the same reconstructed-coordinate coefficient at r = 1024 with the
referee code.  To stay near two minutes the probe-u rule is coarser than
in check 04 (56 panels of 8 Gauss nodes); the same rule scaled to r = 512
(28 panels) is first validated against the value check 04 obtains with
its finer rule, which is quoted as a constant.

Deterministic.  ASCII only.  Exit code 1 on any FAIL.  About 110 s.
"""

from __future__ import annotations

from popsignal import Probes, frozen_anchors, prediction
from wavegeom import Reporter

rep = Reporter()
Q, PHI = 1.0, 0.37
FINE_512 = complex(
    0.017393028, 0.008219749
)  # check 04 output, printed to 8 digits there

frozen = frozen_anchors(0.0)
coarse = Probes(Q, 512, PHI, u_panels=28, u_order=8)
out, _ = coarse.evaluate(frozen, [512.0])
v512 = out[0] * 512**2
print(f"r=512 with the coarse u rule: {v512.real:.9f}{v512.imag:+.9f}i")
rep.check(
    "coarse probe-u rule reproduces the fine-rule r=512 value",
    abs(v512 - FINE_512) / abs(FINE_512) < 2e-5,
    "rel diff %.2e" % (abs(v512 - FINE_512) / abs(FINE_512)),
)

big = Probes(Q, 1024, PHI, u_panels=56, u_order=8)
out, _ = big.evaluate(frozen, [1024.0])
v = out[0] * 1024**2
pred = prediction(big.wave, frozen)
err = abs(v - pred) / abs(pred)
print(
    f"r=1024: S/(lam^2 Ew) = {v.real:.9f}{v.imag:+.9f}i, limit {pred.real:.9f}{pred.imag:+.9f}i, rel err {err:.6f}"
)
rep.check(
    "r=1024 agrees with boxed limit (4) to 0.1%", err < 1e-3, f"rel err {err:.6f}"
)
rep.check(
    "error keeps falling beyond the author's last rung (0.00533 at r=512)",
    err < 0.00533 / 4,
)
rep.finish()
