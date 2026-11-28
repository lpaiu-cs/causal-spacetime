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
- **What supported it before.**
  - A six-line sketch in the audit referee's report
    ([proof1_lower_bound/report.md](independent_audit/proof1_lower_bound/report.md)
    Section 3.3). The second referee did not re-check it.
  - Two computations of `TV/lambda^2` at q = 0.25, 1 and 4:
    - the audit referee's `check_c_one_point_law.py`: Monte Carlo in
      `(u, y)`, a random phase per case, `lambda` down to 1/128;
    - the second referee's R8: deterministic, phase 0.3, `lambda` down to
      1/64. It uses the same polar-reduction method as this note
      (corrected after review, finding F7).
- **What this note does.** It proves the bound with an explicit constant
  `C_1(q)`, uniform in the phase and in `0 < lambda <= 1` (Proposition 1).
  The route differs from the earlier sketch: the transverse integral is
  done exactly in polar coordinates (Step 3), so no estimate of how far the
  section endpoints move is needed. Section 4 lists deterministic checks
  over a wider grid than before. All of them pass.
- **Authorship and review.** AI-written. Reviewed by an AI referee in a
  fresh context:
  [independent_audit/one_point_bound/review.md](independent_audit/one_point_bound/review.md).
  - **Proof:** CORRECT WITH MINOR FIXES. Every step was re-derived, and
    no inequality or constant changed.
  - **Checks:** CORRECT WITH MINOR FIXES. The numbers were reproduced by
    an independent method. Several checks were weaker than this note
    claimed.

  Section 5 lists the corrections.
- **What it is not.** It is not human-reviewed. It concerns the one-point
  laws only, not the order law. It uses the exact causal relation of the
  wave. That relation was imported when this note was written, and has since
  been proved in [wave_causal_relation.md](wave_causal_relation.md).

## 1. Setting, and what is imported

Fix `0 < q < pi^2`, a phase `phi`, `0 < lambda <= 1` and the axial slack
`h > 0`; the model uses `h = 1`, and the result does not depend on `h`
(Remark 2). Put `s = sqrt(q)`, `a(u) = cos(s u)`, `f = a^-2`,
`f_max = sec^2(s/2)`, and, for `-1/2 < u < 1/2`, the quantities below.
*Notation warning:* here `s` is `sqrt(q)`. The probe, the oracle note and
both referee reports use `s = sinh(2 beta) = 2 lambda sqrt(q)` instead.
`t` means the same thing everywhere.

\[
m_-(u)=\int_{-1/2}^{u}f,\qquad m_+(u)=\int_u^{1/2}f,\qquad M=m_-+m_+ .
\]

Write `t = 2 lambda s / sqrt(1 + 4 q lambda^2)`, so `0 < t < 1`. In the y
chart, for `u_1 < u_2`,

\[
E_\lambda(u_1,u_2)=-t\int_{u_1}^{u_2}f(w)\,T(w/\lambda+\phi)\,dw,\qquad
T(\theta)=\begin{pmatrix}\cos\theta&\sin\theta\\ \sin\theta&-\cos\theta\end{pmatrix}
\]

([coordinate_oracle_barrier.md](coordinate_oracle_barrier.md) (1)).
`E_- = E_lambda(-1/2, u)` and `E_+ = E_lambda(u, 1/2)` are therefore
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

**Imported, not proved here.** (Status since review: the causal relation,
the matrix form of `E_lambda` and the volume density are now proved in
[wave_causal_relation.md](wave_causal_relation.md) (Theorem 2, Lemma 3).
That note was itself reviewed by a fresh-context AI referee. The model's
one-point laws remain a definition.)
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
- The matrix form of `E_lambda` displayed above, from the same sources.
  Proposition 1 uses only its eigenvalues `+/- t |G|`. Remark 4 also uses
  the matrix form.

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

**Step 4 (TV in terms of the difference regions).** The convention is
`TV(mu, nu) = sup_A |mu(A) - nu(A)| = 1 - integral min(p, p_0)`, with
values in `[0, 1]`, as in quantized_order_bound_r2.md. Write
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

**Remark 3 (what is not shown).** The note gives an upper bound only.
Argument 1 needs only the upper bound. That `lambda^2` is the true order,
rather than an overestimate, is not proved here. The numerics of Section 4
support it, and so does the referee's leading-order limit (review,
Section 4), which this note has not re-derived:

\[
\frac{\mathrm{TV}}{\lambda^2}\to L(q)
=\frac{2\sqrt q}{\pi\widehat V}\int a^2\,
\big\langle|f(u)-f_{\max}e^{i\chi}|\big\rangle_\chi\,
\frac{m_+^2+m_-^2}{M^2}\,du ,
\]

where `<.>_chi` averages over `chi`. Its values are 1.61982, 3.26922,
8.06357, 81.96 and 1124.5 at q = 0.25, 1, 4, 9 and 9.8. The first three
match Section 4's values at `lambda = 1/1024` to the printed digits.

**Remark 4 (the TV does not depend on the phase at all; added after the
run).** Let `R` be the rotation of the y-plane by `phi/2`. Since
`R T(theta) R^T = T(theta + phi)`, every `E_lambda(u_1, u_2)` at phase
`phi` equals `R E R^T` with `E` at phase 0. So the wave diamond at phase
`phi` is the image of the one at phase 0 under `(u, v, y) -> (u, v, R y)`.
`D_0` and the density `a^2` are invariant under that map, so the TV is the
same for every phase. Uniformity in the phase is therefore automatic.
Proposition 1 does not need this remark. It uses the matrix form of
`E_lambda` (Section 1), not only its eigenvalues.
- **The 8-phase grid of Section 4 does not test it** (corrected after
  review, finding F4). A phase `k pi / 4` rotates y by `k pi / 8`, which
  is a whole number of steps of the 256-point `theta` grid. So the
  discretised TV is phase-invariant by construction, and a wrong phase
  convention would pass as well.
- **The real test** is the referee's computation at the off-grid phase
  1.234, by a different method. It gave the phase-0 value to 7 digits.

**Remark 5 (a sharper constant, found in review, finding F12).** By (2),
`A - B` lies in `[-e/(m(m+e)), e/(m(m-e))]`. Step 5 bounds
`(B - A)^+` by the larger end, but the smaller one suffices:
`(B - A)^+ <= e/(m(m+e)) <= e/m^2`. This gives
`Delta_lam <= pi h^2 (1+t) J_lam`, while
`Delta_0 <= pi h^2 J_lam/(1-t)` is unchanged. Since `1 + t <= 1/(1-t)`,

\[
\mathrm{TV}(\mu_\lambda,\mu_0)\le\frac{C_E}{1-t}\,\frac{W}{\widehat V}\,\lambda^2 .
\]

That is `C_1(q, t)` divided by `1 + t`. The uniform constant
`C_E W / (Vhat (1 - t_1))` is 29.48, 200.5, 2749, 7.55e4 and 1.129e6 at
q = 0.25, 1, 4, 9 and 9.8. Proposition 1 keeps the reviewed constant, and
argument 1 needs only that some constant exists. The referee also notes
that the `theta` step loses a further factor of about `pi` at leading
order.

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

This makes 840 cases. Because of the phase invariance below, only 105
`(q, lambda)` pairs are distinct.

The earlier computations covered q = 0.25, 1 and 4. The audit referee's
was Monte Carlo with a random phase, down to `lambda = 1/128`. R8 was
deterministic at phase 0.3, down to 1/64.

**What each check actually tests** (added after review; the criteria above
are left as they were fixed before the run).
- **K5 cannot fail on its own** (F2). On the same nodes, Step 5 holds
  node by node once `e <= t m` (K1) and `|n^T E n| <= e` hold, and K5's
  chain then follows by summation. K5 therefore shows that Steps 4 to 6
  are implemented consistently. It is not evidence that the computed
  number is the TV. That evidence is K2, at 10 points, and the referee's
  independent computation (review, Section 4).
- **K1's last two parts cannot fail** (F3). `A - B` is evaluated from
  Step 1's own closed form, so they hold identically up to rounding. The
  inverse formula itself is tested by K2, and by the referee on 200,000
  random matrices.
- **K3 does not see the `max()` kinks** (F5). On each node the `S_cap`
  terms cancel in `Delta_lam - Delta_0`, so K3 tests only the smooth part
  of the `theta` integral. At the smallest `lambda`, its recorded
  `1.03e-6` is cancellation noise on the closed-form side.
- **The phase grid is no test** (F4; Remark 4).
- **K7 is not independent** (F7). R8 uses the same polar-reduction
  method, so K7 compares two implementations of one method.

**Results** (full run after commit `c4c3b78`; output:
[checks_output.txt](independent_audit/one_point_bound/checks_output.txt)).
**All 11 checks pass, and no case fails.**
- **K5.** The bound holds in all 840 cases. The ratio
  `TV / (C_1(q,t) lambda^2)` runs from `4.2e-7` (q = 9.8, `lambda = 1`)
  to 0.187 (q = 0.25, smallest `lambda`). So `C_1(q,t)` is **at least**
  5.3 times too large on this grid. As `lambda -> 0`, the overestimate
  tends to 5.3, 6.5, 10.2, 12.5 and 12.6 at q = 0.25, 1, 4, 9 and 9.8,
  using the referee's limits `L(q)`. The uniform `C_1(q)`, the constant
  quantized_order_bound_r2.md uses, overestimates at small `lambda` by
  about 31 times (q = 0.25) to about 2000 times (q = 9.8). (Corrected
  after review, finding F1: the first version said "at most 5 to 20
  times".)
- **K1.** `e / (t m)` and the ratio in (2) both reach 1.000000. These two
  inequalities are attained, on short segments and in the extreme
  direction. `e / (t lambda (4 f_max - 2))` reaches 0.942.
- **K2.** Agreement to `1.2e-6`.
- **K3.** Largest relative difference `1.03e-6`. It passes only through
  the absolute term `1e-12 V_0`, where `V_lam - V_0` is about `1e-10`
  of `V_0`. See "What each check actually tests".
- **K4.** `4.4e-16`.
- **K6.** Largest change `3.9e-6`, at the two values of `lambda` it
  tests, 1/8 and 1/256.
  - It does not cover q = 9.8 at `lambda = 1/64`. There the 256-point
    `theta` rule is off by `8.8e-6` relative, and doubling moves the
    value by `1.6e-5` (F6).
  - The converged value there is 120.0547; the output prints 120.06.
  - The bound is unaffected.
- **K7.** Agreement within `3.7e-4`.
- **The phase.** The 8 phases give the same TV in every `(q, lambda)`.
  This holds by construction of the grid (F4). The referee's off-grid
  phase 1.234 also gave the same value, to 7 digits.
- **The referee's independent computation** used:
  - a different `G` quadrature;
  - the adjugate inverse;
  - the ellipse formula for `V_lam`;
  - kink-aware `theta` quadrature;
  - `TV = (1/2) integral |p_lam - p_0|`.

  It reproduced the table below at 25 `(q, lambda)` points to within
  5-figure rounding. The exceptions are q = 9.8 at `lambda = 1/64` and
  1/256, off by about `1e-5`.
- **The referee's seeded random cases.** There were 24 of them, with q
  up to 9.8, `lambda` from `2^-8` to 1, and random phases. Proposition 1's
  chain held in all, with the largest ratio 0.116.

`TV / lambda^2`, the same for every phase (q = 0.25 at `lambda = 1/8`
corrected from 1.645 to 1.644 after review, F8):

| q | `lambda = 1` | 1/8 | 1/64 | 1/1024 | `C_1^0` | largest `TV/(C_1(q,t) lambda^2)` |
|---|---|---|---|---|---|---|
| 0.25 | 0.398 | 1.644 | 1.623 | 1.620 | 8.64 | 0.187 |
| 1 | 0.634 | 3.440 | 3.283 | 3.269 | 21.2 | 0.154 |
| 4 | 0.811 | 7.083 | 8.047 | 8.064 | 82.1 | 0.098 |
| 9 | 0.917 | 16.07 | 53.01 | 81.44 | 1027 | 0.078 |
| 9.8 | 0.947 | 22.54 | 120.1 | 765.7 | 14130 | 0.054 |

- **q = 0.25, 1 and 4.** `TV / lambda^2` settles to within `1e-4` by
  `lambda = 1/1024`.
- **q = 9 and 9.8.** It is still rising at `lambda = 1/1024`. Close to
  `pi^2` the distance from the window's ends to the pole of `f` is small
  (0.024 at q = 9, 0.0018 at q = 9.8). The `lambda^2` regime starts only
  once `lambda` is well below that distance. The bound holds throughout,
  and the proof does not need the regime to have started.
- **What this does not show.** These are checks at grid points. The
  statement for every `lambda` and every q rests on the proof of
  Section 3, not on the grid.

## 5. Review and the corrections made after it

[independent_audit/one_point_bound/review.md](independent_audit/one_point_bound/review.md).
An AI referee in a fresh context reviewed this note as committed at
`3fc1ffc`, together with the script and its output.
- It re-derived every step.
- It re-ran `checks.py`. The output matched the recorded one line for
  line, apart from the referee's own timing line.
- It computed the TV again by a different method.
- Its scripts and outputs are archived in
  [independent_audit/one_point_bound/referee/](independent_audit/one_point_bound/referee/).

**Verdicts.**
- **Proof: CORRECT WITH MINOR FIXES.** Steps 1 to 7 are correct. Proposition 1
  is the statement quantized_order_bound_r2.md needs: the same chart,
  tips, density and laws, uniform in phase and in `0 < lambda <= 1`.
  The referee found no mismatch.
- **Checks: CORRECT WITH MINOR FIXES.** The numbers are right, but
  several checks are weaker than this note first claimed.

**Corrected after the review.** No inequality, constant or statement of
Proposition 1 changed.

| finding | correction |
|---|---|
| F1 | Results, K5: the constant is *at least* 5.3 times too large, not "at most 5 to 20 times"; the overestimates of `C_1(q, t)` and `C_1(q)` are stated |
| F2 | K5 is stated to follow from K1 by summation; evidence that the number is the TV is K2 and the referee's computation |
| F3 | K1's last two parts are stated to hold by construction |
| F4 | Remark 4 and Results: the 8-phase grid is phase-invariant by construction; only 105 `(q, lambda)` pairs are distinct; the off-grid phase is the referee's |
| F5 | K3 is stated to be blind to the kinks, and its `1.03e-6` is cancellation noise |
| F6 | K6's coverage is stated; q = 9.8 at `lambda = 1/64` is 120.0547 |
| F7 | preamble and Section 4: the earlier support was Monte Carlo (audit referee) and R8 (same method, `lambda` down to 1/64); K7 is not independent |
| F8 | table: 1.645 becomes 1.644 |
| F9 | Section 1 states the matrix form of `E_lambda` and lists it as imported; Remark 4 says Proposition 1 does not need it |
| F10 | Section 1: a warning that `s` means `sqrt(q)` here and `sinh(2 beta)` elsewhere |
| F11 | Step 4 states the TV convention |
| F12 | Remark 5: the sharper constant `C_E W / (Vhat (1 - t))`; Proposition 1 keeps the reviewed one |

The referee's limit `L(q)` is added to Remark 3, as the referee's
derivation. This note has not re-derived it.

**The step the referee trusts least** is outside this note's own steps.
- **The step:** Step 2's reduction depends on the imported exact wave
  relation, in the form `mI + E` with `e = t |G|` and
  `t = tanh(2 beta)`. Every constant scales through `t` and `e <= t m`,
  and every check, the referee's included, assumes that form.
- **Its status:** re-derived once, in the audit referee's report
  Section 3.1, and not by this note's referee.
- **The note's own steps:** the referee found none it distrusts.

**Later the same day: that step is proved.** The exact wave relation in the
form `mI + E`, with `e = t |G|` and `t = tanh(2 beta)`, is now proved
in [wave_causal_relation.md](wave_causal_relation.md).
- **What it proves:** Theorem 2 gives the relation of any Rosen plane
  wave. Lemma 3 gives the circular family's form and the volume density.
  Lemma 4 gives the vacuum condition.
- **Its review:** a fresh-context AI referee graded it CORRECT WITH MINOR
  FIXES
  ([independent_audit/wave_relation/review.md](independent_audit/wave_relation/review.md)).
- **What remains here:** this note's imports are now only the model's
  definitions.

**Not checked by the referee:**
- the exact wave relation and the vacuum condition;
- Lemma 2.1 beyond its integration by parts;
- the rest of quantized_order_bound_r2.md;
- the earlier sketch;
- q between 9.8 and `pi^2`;
- `lambda > 1`;
- `kappa(q)`.