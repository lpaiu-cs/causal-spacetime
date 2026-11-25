# The one-point strip bound, written out with an explicit constant

2026-11-26. **Proof note: the bound `TV(mu_lam, mu_0) <= C_1 lambda^2`
between the wave and reference one-point laws.**

- **Why this note exists.** [quantized_order_bound_r2.md](quantized_order_bound_r2.md)
  (argument 1, the lower end `4/3` of the bracket) imports this bound and
  does not prove it. That revision's fresh-context AI referee named it the
  step it trusts least
  ([independent_audit/argument1_revision/review.md](independent_audit/argument1_revision/review.md)).
  Every `alpha > 1` in argument 1 depends on it. The same bound is used in
  [coordinate_oracle_barrier.md](coordinate_oracle_barrier.md) Sections 1
  and 2 and in [full_order_resolution_probe.md](full_order_resolution_probe.md) (4).
- **What supported it before.** A six-line sketch in the audit referee's
  report ([proof1_lower_bound/report.md](independent_audit/proof1_lower_bound/report.md)
  Section 3.3), and two independent deterministic computations of
  `TV/lambda^2` at q = 0.25, 1 and 4, one phase, and `lambda` from 1/8 to
  1/128. The sketch was not re-checked by the second referee.
- **What this note does.** It proves the bound with an explicit constant
  `C_1(q)`, uniform in the phase and in `0 < lambda <= 1` (Proposition 1).
  The route differs from the earlier sketch: the transverse integral is
  done exactly in polar coordinates (Step 3), so no estimate of how far the
  section endpoints move is needed. Section 4 lists deterministic checks
  over a wider grid than before.
- **Authorship and review.** AI-written. Awaiting a fresh-context AI review.
- **What it is not.** It is not human-reviewed. It concerns the one-point
  laws only, not the order law. It uses the exact causal relation of the
  wave, which remains imported (Section 1).

## 1. Setting, and what is imported

Fix `0 < q < pi^2`, a phase `phi`, `0 < lambda <= 1` and the axial slack
`h > 0`; the model uses `h = 1`, and the result does not depend on `h`
(Remark 2). Put `s = sqrt(q)`, `a(u) = cos(s u)`, `f = a^-2`,
`f_max = sec^2(s/2)`, and, for `-1/2 < u < 1/2`,

\[
m_-(u)=\int_{-1/2}^{u}f,\qquad m_+(u)=\int_u^{1/2}f,\qquad M=m_-+m_+ .
\]

Write `t = 2 lambda s / sqrt(1 + 4 q lambda^2)`, so `0 < t < 1`. In the y
chart, `E_- = E_lambda(-1/2, u)` and `E_+ = E_lambda(u, 1/2)` are
symmetric, trace-free 2x2 matrices with eigenvalues `+/- e_-` and
`+/- e_+`, where

\[
e_\mp=t\,|G_\mp|,\qquad
G_-=e^{i\phi}\!\int_{-1/2}^{u}f(w)e^{iw/\lambda}dw,\quad
G_+=e^{i\phi}\!\int_u^{1/2}f(w)e^{iw/\lambda}dw .
\]

The tips are `p_- = (u, v, y) = (-1/2, 0, 0)` and `p_+ = (1/2, h, 0)`.
`D_lam` is the wave diamond `{X : p_- <=_lam X <=_lam p_+}` and `D_0` the
reference diamond, with the background relation. `mu_lam` and `mu_0` have
density proportional to `a(u)^2` in `du dv d^2y` on `D_lam` and `D_0`
(the physical factor `c` cancels on normalising). Write
`vol_a(S) = integral_S a^2 du dv d^2y`.

**Imported, not proved here.**
- The model's one-point laws as just stated
  ([full_order_resolution_probe.md](full_order_resolution_probe.md),
  [coordinate_oracle_barrier.md](coordinate_oracle_barrier.md) Section 1).
- The exact causal relation of the wave in the y chart
  ([quantized_order_bound_r2.md](quantized_order_bound_r2.md) Section 2;
  re-derived in the audit referee's report, Section 3.1): for
  `u_1 < u_2`, `X <=_lam Y` iff
  `2 Delta v >= Delta y^T (m(u_1,u_2) I + E_lambda(u_1,u_2))^-1 Delta y`.
  The background relation is the same with `E = 0`. No relation holds
  forward when `u_1 > u_2`.

**Used, and proved elsewhere in writing.** Lemma 2.1 of
[quantized_order_bound_r2.md](quantized_order_bound_r2.md): with
`C_E = 2 s (4 f_max - 2)`, for every segment, `e <= t m` and
`e <= t lambda (4 f_max - 2) <= C_E lambda^2`. In particular `e < m`.

## 2. Statement

Put `T_h = tan(s/2)` and

\[
I_c=\int_{-1/2}^{1/2}a^2=\tfrac12+\tfrac{\sin s}{2s},\qquad
I_s=\tfrac12-\tfrac{\sin s}{2s},
\]
\[
\widehat V=\int_{-1/2}^{1/2}a^2\,\frac{m_+m_-}{M}\,du=\frac{T_h^2I_c-I_s}{2sT_h},
\qquad
W=\int_{-1/2}^{1/2}a^2\,\frac{m_+^2+m_-^2}{M^2}\,du=\frac{I_s+T_h^2I_c}{2T_h^2}.
\]

**Proposition 1.** For every `0 < q < pi^2`, every phase, every `h > 0`
and every `0 < lambda <= 1`,

\[
\operatorname{TV}(\mu_\lambda,\mu_0)
\;\le\;\frac{1+t}{1-t}\,\frac{J_\lambda}{\widehat V}
\;\le\;C_1(q,t)\,\lambda^2
\;\le\;C_1(q)\,\lambda^2,
\tag{1}
\]

where

\[
J_\lambda=\int_{-1/2}^{1/2}a^2\,\frac{e_-m_+^2+e_+m_-^2}{M^2}\,du,\qquad
C_1(q,t)=\frac{1+t}{1-t}\,C_E\,\frac{W}{\widehat V},\qquad
C_1(q)=C_1(q,t_1),\quad t_1=\frac{2s}{\sqrt{1+4q}} .
\]

`C_1(q, t)` decreases to `C_1^0(q) = C_E W / Vhat` as `lambda -> 0`. So
`limsup TV / lambda^2 <= C_1^0(q)`.

| q | `f_max` | `C_E` | `W/Vhat` | `C_1^0` | `t_1` | `(1+t_1)/(1-t_1)` | `C_1(q)` |
|---|---|---|---|---|---|---|---|
| 0.25 | 1.065 | 2.261 | 3.820 | 8.636 | 0.7071 | 5.828 | 50.33 |
| 1 | 1.298 | 6.388 | 3.315 | 21.17 | 0.8944 | 17.94 | 379.9 |
| 4 | 3.426 | 46.81 | 1.754 | 82.09 | 0.9701 | 65.98 | 5417 |
| 9 | 199.9 | 4784 | 0.2147 | 1027 | 0.9864 | 146.0 | 1.50e5 |
| 9.8 | 3.248e4 | 8.134e5 | 0.01737 | 1.413e4 | 0.9875 | 158.8 | 2.24e6 |

The constants are uniform in the phase. They are not uniform as `q`
approaches `pi^2`, through `f_max`.

## 3. Proof

**Step 1 (2x2 algebra).** Let `E` be symmetric and trace-free with
eigenvalues `+/- e`, and `0 <= e < m`. Then `E^2 = e^2 I`, so
`(mI + E)^-1 = (mI - E)/(m^2 - e^2)`. For a unit vector `n`, put
`A = n^T (mI+E)^-1 n` and `B = 1/m`. Then

\[
A=\frac{m-n^TEn}{m^2-e^2}\in\Big[\frac1{m+e},\frac1{m-e}\Big],\qquad
A-B=\frac{e^2-m\,n^TEn}{m(m^2-e^2)},\qquad
|A-B|\le\frac{e(e+m)}{m(m^2-e^2)}=\frac{e}{m(m-e)} .
\tag{2}
\]

**Step 2 (sections).** Causal curves have `u` non-decreasing, so both
diamonds lie in `-1/2 <= u <= 1/2`. The slices `u = +/- 1/2` are null
sets. For `-1/2 < u < 1/2` and `X = (u, v, y)`, the imported relation gives

\[
X\in D_\lambda\iff
\tfrac12\,y^T\widetilde M_-^{-1}y\;\le\;v\;\le\;h-\tfrac12\,y^T\widetilde M_+^{-1}y,
\qquad \widetilde M_\mp=m_\mp I+E_\mp ,
\]

and the same with `E_mp = 0` for `D_0`. Write `y = r n(theta)`, with
`n = (cos theta, sin theta)` and `r >= 0`. Write `A_mp(theta)` for the
quantity `A` of Step 1 built from `(m_mp, E_mp, n)`, and `B_mp = 1/m_mp`.
Then the v-section of `D_lam` at `(u, y)` is
`[r^2 A_- / 2, h - r^2 A_+ / 2]`, and that of `D_0` is
`[r^2 B_- / 2, h - r^2 B_+ / 2]`. Their intersection is the section of
`D_lam cap D_0`:
`[r^2 max(A_-, B_-) / 2, h - r^2 max(A_+, B_+) / 2]`.

**Step 3 (polar reduction).** For `S > 0`,
`integral_0^inf (h - r^2 S / 2)^+ r dr = h^2 / (2 S)`: substitute
`x = r^2`. Put

\[
S_\lambda=A_++A_-,\qquad S_0=B_++B_-=\frac{M}{m_+m_-},\qquad
S_\cap=\max(A_+,B_+)+\max(A_-,B_-).
\]

All three are positive functions of `(u, theta)`. By Tonelli,

\[
\operatorname{vol}_a(D_\bullet)=\frac{h^2}{2}\int_{-1/2}^{1/2}a^2\!\int_0^{2\pi}\frac{d\theta}{S_\bullet}\,du,
\qquad \bullet\in\{\lambda,0,\cap\},
\]

with `D_cap = D_lam cap D_0`. Since `max(A, B) - A = (B - A)^+` and
`max(A, B) - B = (A - B)^+`,

\[
\Delta_\lambda:=\operatorname{vol}_a(D_\lambda\setminus D_0)
=\frac{h^2}{2}\iint a^2\,\frac{\sum_\mp(B_\mp-A_\mp)^+}{S_\lambda S_\cap},
\qquad
\Delta_0:=\operatorname{vol}_a(D_0\setminus D_\lambda)
=\frac{h^2}{2}\iint a^2\,\frac{\sum_\mp(A_\mp-B_\mp)^+}{S_0 S_\cap}.
\tag{3}
\]

Also `V_0 := vol_a(D_0) = pi h^2 Vhat`, because `1/S_0 = m_+ m_- / M`
does not depend on `theta`.

**Step 4 (TV in terms of the difference regions).** Write
`V_lam = vol_a(D_lam)`. The densities are `a^2 1_D / V`, so the overlap is
`integral min = vol_a(D_cap) / max(V_lam, V_0)`. Since
`max(V_lam, V_0) - vol_a(D_cap) = max(Delta_lam, Delta_0)`,

\[
\operatorname{TV}(\mu_\lambda,\mu_0)
=\frac{\max(\Delta_\lambda,\Delta_0)}{\max(V_\lambda,V_0)}
\le\frac{\max(\Delta_\lambda,\Delta_0)}{V_0}.
\tag{4}
\]

**Step 5 (pointwise bounds).** By (2) and `e_mp <= t m_mp`,

\[
(B_\mp-A_\mp)^+,\ (A_\mp-B_\mp)^+\ \le\ \frac{e_\mp}{m_\mp(m_\mp-e_\mp)}
\le\frac{1}{1-t}\,\frac{e_\mp}{m_\mp^2}.
\]

By (2) again, `A_mp >= 1/(m_mp + e_mp) >= B_mp / (1 + t)`, so
`S_lam >= S_0 / (1 + t)`. And `S_cap >= S_0`. Hence

\[
\frac{\sum_\mp(B_\mp-A_\mp)^+}{S_\lambda S_\cap}
\le\frac{1+t}{1-t}\,\frac{1}{S_0^2}\Big(\frac{e_-}{m_-^2}+\frac{e_+}{m_+^2}\Big)
=\frac{1+t}{1-t}\,\frac{e_-m_+^2+e_+m_-^2}{M^2},
\]

using `1/S_0 = m_+ m_- / M`. The same holds for the `Delta_0` integrand
without the factor `1 + t`. The right side does not depend on `theta`.

**Step 6 (integrate).** From (3) and Step 5,
`Delta_lam <= pi h^2 (1+t)/(1-t) J_lam` and
`Delta_0 <= pi h^2 J_lam / (1-t)`. With (4) and `V_0 = pi h^2 Vhat`, this is
the first inequality of (1). Lemma 2.1 gives `e_mp <= C_E lambda^2`, so
`J_lam <= C_E lambda^2 W`: the second inequality. `t` is increasing in
`lambda`, so `t <= t_1` for `lambda <= 1`: the third.

**Step 7 (closed forms).** `m_mp' = +/- a^-2` gives
`m_-(u) = (tan(s u) + T_h)/s`, `m_+(u) = (T_h - tan(s u))/s` and
`M = 2 T_h / s`. Then `a^2 m_+ m_- = (T_h^2 cos^2(s u) - sin^2(s u))/q` and
`a^2 (m_+^2 + m_-^2) = 2 (sin^2(s u) + T_h^2 cos^2(s u))/q`. Integrating
over `u` and dividing by `M` and `M^2` gives `Vhat` and `W`. []

**Remark 1 (uniformity).** The phase enters only through `G_mp`, and only
`|G_mp|` is used. The bound is therefore uniform in the phase. Near the
tips `m_mp -> 0`, and Step 5 has a factor `1/m_mp^2`. The factor
`1/S_0^2 = (m_+ m_- / M)^2` cancels it, so the integrand of `J_lam` is at
most `e_- + e_+ <= 2 C_E lambda^2` everywhere on the window.

**Remark 2 (`h`).** All of `Delta_lam`, `Delta_0`, `V_lam` and `V_0` scale as
`h^2`, so the TV does not depend on `h`.

**Remark 3 (what is not shown).** The note gives an upper bound only. That
`lambda^2` is the true order, rather than an overestimate, is suggested by
the numerics (Section 4), not proved. Argument 1 needs only the upper bound.

## 4. Deterministic checks

Script: [independent_audit/one_point_bound/checks.py](independent_audit/one_point_bound/checks.py).
No sampling. The u-integrals use composite 12-point Gauss panels of width
at most `lambda/6`, half the distance to the pole of `f`, and `1/32`. The
primitive `G` is integrated to every node. The `theta` integrals use 256
midpoints on a half period.

**Criteria, fixed before the full run.** One smoke case (q = 1,
lambda = 1/8, phase 0.3) was run first to debug the script.
- **K1:** at every node `(u, theta)` of every case:
  - `e <= t m`, and `e <= t lambda (4 f_max - 2)`;
  - `|A - B| <= e / (m (m - e))`;
  - `A >= 1/(m + e)`.

  All hold to relative slack `1e-9`.
- **K2:** the polar reduction (Step 3) against a Cartesian y-grid, at 5
  values of `u` in 2 cases, to relative `1e-3`.
- **K3:** the `theta` grid against the closed form
  `integral dtheta / S_lam = 2 pi / sqrt(P^2 - |Z|^2)`, for `V_lam - V_0`.
- **K4:** `V_0`, `W` and `Vhat` against their closed forms, to `1e-9`.
- **K5, the main check:** in every case,
  `TV <= B(lambda, phi) <= C_1(q,t) lambda^2 <= C_1(q) lambda^2`, and
  `Delta_lam` and `Delta_0` each lie below their bounds.
- **K6:** doubling the panels and the `theta` grid changes the TV by less
  than `1e-3` relative.
- **K7:** the previous referee's R8 values at phase 0.3 are reproduced to
  within 1%.

**Grid for K1, K3 and K5.**
- q in {0.25, 1, 4, 9, 9.8}. 9.8 is close to `pi^2`.
- 8 phases, `k pi / 4`.
- `lambda = 2^(-k/2)` for `k = 0..20`, that is from 1 down to 1/1024.

This makes 840 cases. The earlier computations covered three values of
q, one phase, and `lambda` from 1/8 to 1/128.

**Results.** To be filled in after the run.
