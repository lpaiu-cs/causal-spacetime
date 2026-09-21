# 진공 다이아몬드 부피 불일치: 직접 전개와 오차 추적

검토일: 2026-09-22. 기존 확인실험 이후의 결정론적 수학 검토다.
비교 대상은 [Wang, arXiv:1904.01034v2](https://arxiv.org/html/1904.01034v2)의
명시된 식 (74)–(79)이다. HTML과 arXiv PDF 12쪽에서 부호와 계수를 대조했다.
출판사 조판본이나 저자의 미공개 계산을 확인했다는 뜻은 아니다.

**결론은 두 개의 대수적 불일치다.** 식 (76)의 Lorentz 수축을 전개하면 식 (77)의
전기적 항 한 곳은 `+2`가 아니라 `−2`다. 또 인쇄된 식 (74),(77)을 그대로
적분해도 식 (79)의 자기적 항 계수는 나오지 않고 그 절반이 나온다.
두 차이는 기존의 평면파 부피 계수 차이 `11/18900`을 정확히 설명한다.
아래 차원 다항식은 명시된 중간식의 재적분 결과이며, 모든 진공 기하에서
생략항까지 독립적으로 재유도한 일반 부피 정리로 주장하지 않는다.

**1. 먼저 문헌식 없이 정확한 부피를 계산한다.**

기존 연구와 같은 metric, 끝점, 단위를 쓴다.

\[
 ds^2=2\,du\,dv+A(x^2-y^2)du^2+dx^2+dy^2,\qquad
 P=(0,\tfrac12,0,0),\quad Q=(1,-\tfrac12,0,0).
\]

끝점 사이 고유시간은 \(T=\sqrt2\), 반길이는 \(l=T/2=1/\sqrt2\)다.
이 좌표 단위에서 \(w=\sqrt A\)로 놓는다. \(0\le A<\pi^2\)인 정상 근방에서는
횡방향 경계값 문제가 첫 conjugate point 이전에 있다. 기존 실험 범위는
이보다 좁은 \(0\le A\le2\)다. 음의 A는 x와 y의 교환으로 환원된다.

u를 매개변수로 한 null 경로의 v 변화는 다음 작용으로 정해진다.

\[
 \Delta v=-\frac12\int
 \left[x'^2+y'^2+A(x^2-y^2)\right]du.
\]

Euler–Lagrange 방정식은 \(x''=Ax,\ y''=-Ay\)다. 한 끝점이 횡방향 원점이고
u 간격이 s이면, 경계값 해를 작용에 대입하여

\[
 S_s(x,y)=\frac w2\left[x^2\coth(ws)+y^2\cot(ws)\right]
\]

를 얻는다. 따라서 u 단면에서 허용되는 v의 폭은

\[
 \Delta v(u,x,y)=1-B_x(u)x^2-B_y(u)y^2,
\]
\[
 B_x=\frac w2[\coth(wu)+\coth(w(1-u))],\qquad
 B_y=\frac w2[\cot(wu)+\cot(w(1-u))].
\]

\(\sqrt{-g}=1\)이지만 적분 영역은 A에 의존한다. \(X=\sqrt{B_x}x\),
\(Y=\sqrt{B_y}y\)로 바꾸어 단위원판에서 적분하면

\[
 \int\!dx\,dy\,[1-B_xx^2-B_yy^2]_+
 =\frac{1}{\sqrt{B_xB_y}}\int_0^{2\pi}\!d\phi
   \int_0^1\rho(1-\rho^2)d\rho
 =\frac{\pi}{2\sqrt{B_xB_y}}.
\]

삼각함수 덧셈공식을 사용하면

\[
 V(A)=\frac\pi2\int_0^1S(u)du,
\quad
 S(u)=\frac2w\sqrt{
 \frac{\sin(wu)\sinh(wu)\sin(w(1-u))\sinh(w(1-u))}
 {\sin w\sinh w}}.
\]

여기서 직접 Taylor 전개한다.

\[
 \frac{\sin z\sinh z}{z^2}=1-\frac{z^4}{90}+O(z^8),
\]
\[
 S(u)=2u(1-u)\left[
 1+\frac{A^2}{180}\{1-u^4-(1-u)^4\}+O(A^4)\right].
\]

필요한 적분은 모두 다항식이다.

\[
 \int_0^1 2u(1-u)du=\frac13,\qquad
 \int_0^1 2u(1-u)\{1-u^4-(1-u)^4\}du=\frac5{21}.
\]
\[
 V(A)=\frac\pi2\left[\frac13+\frac{A^2}{756}+O(A^4)\right],
 \qquad
 \boxed{\frac{V(A)}{V(0)}=1+\frac{A^2}{252}+O(A^4)}.
\]

일반 반길이로 쓰면 \(V_0=2\pi l^4/3\),
\(V/V_0=1+A^2l^4/63+O(A^4l^8)\)다. 위의 \(l^4=1/4\)를 넣으면 같다.
이 계수는 곡선 맞춤이나 Monte Carlo 추정값이 아니다.

**2. 식 (76)에서 (77)로 넘어가는 부호를 성분별로 추적한다.**

시공간 차원을 d, 공간 차원을 \(m=d-1\)로 놓는다. 중간점의 직교 정상좌표에서
\(\eta=\operatorname{diag}(-1,1,\ldots,1)\), \(X^\mu=(t,rn^i)\),
\(n^in_i=1\)이다. [Wang의 식 (1)](https://arxiv.org/html/1904.01034v2#SII)
과 같은 정의를 쓴다.

\[
 E_{ij}=C_{0i0j},\quad H_{ijk}=C_{0ijk},\quad D_{ijkl}=C_{ijkl},
\qquad E^2=\sum E_{ij}^2,\ H^2=\sum H_{ijk}^2,\ D^2=\sum D_{ijkl}^2.
\]

식 (76)의 사차 수축을 \(K_{\alpha\beta}=C_{\alpha\mu\beta\nu}X^\mu X^\nu\),
\(Q=K_{\alpha\beta}K^{\alpha\beta}\)로 쓰면 \(\sqrt{-g}=1-Q/180+\cdots\)다.
Riemann 대칭성만 사용하여 성분을 전개한다.

\[
 K_{00}=r^2E_{ij}n^in^j,
\]
\[
 K_{0i}=-trE_{ij}n^j+r^2H_{jik}n^jn^k,
\]
\[
 K_{ij}=t^2E_{ij}-tr(H_{ijk}+H_{jik})n^k
             +r^2D_{ikjl}n^kn^l.
\]

시간 지수를 한 번 올릴 때 음수가 붙으므로

\[
 \boxed{Q=K_{00}^2-2\sum_iK_{0i}^2+\sum_{ij}K_{ij}^2}.
\]

그러므로 \(-2K_{0i}^2\)가 주는 \(t^2r^2\) 항은 반드시

\[
 -2t^2r^2 E_{ik}E^k{}_j n^in^j
\]

다. 식 (77)에는 이 위치가 \(+2E_{ik}E^k{}_j\)로 인쇄되어 있다.
각도에 대해 홀수인 항은 원문에서도 생략되므로, 정확한 비교는

\[
 \boxed{Q_{77}-Q_{76,\mathrm{even}}
        =4t^2r^2(E^2)_{ij}n^in^j}
\]

다. 여기서 \((E^2)_{ij}=E_{ik}E^k{}_j\)는 행렬이고, 지수 없는 \(E^2\)는 그 trace다.
\(C\mapsto-C\)라는 Riemann 부호 규약 변경은 양변의 곡률 제곱을 바꾸지 못한다.

**3. 각도 평균에서는 어느 항이 상쇄되는지 확인한다.**

\(\langle\cdot\rangle\)는 단위 \(S^{m-1}\)의 정규화된 평균이다.

\[
 \langle n^in^j\rangle=\frac{\delta^{ij}}m,
\quad
 \langle n^in^jn^kn^l\rangle
 =\frac{\delta^{ij}\delta^{kl}+\delta^{ik}\delta^{jl}
                  +\delta^{il}\delta^{jk}}{m(m+2)}.
\]

\(e(n)=E_{ij}n^in^j\), \(e_2(n)=(E^2)_{ij}n^in^j\),
\(h_a(n)=H_{iaj}n^in^j\), \(F_{ij}(n)=D_{ikjl}n^kn^l\)로 둔다.
trace-free 조건과 Bianchi 항등식으로

\[
 \langle e^2\rangle=\frac{2E^2}{m(m+2)},\quad
 \langle e_2\rangle=\frac{E^2}m,\quad
 \langle h_ah_a\rangle=\frac{3H^2}{2m(m+2)},
\]
\[
 \langle F_{ij}F_{ij}\rangle=\frac{E^2+\frac32D^2}{m(m+2)},\quad
 E_{ij}\langle F_{ij}\rangle=\frac{E^2}m,
\]
\[
 \left\langle[(H_{ijk}+H_{jik})n^k]^2\right\rangle
 =\frac{3H^2}m.
\]

특히 마지막 식의 \(3H^2\)는 두 제곱의 \(2H^2\)와 교차항의 \(H^2\)를 더한 값이다.
\(H_{ijk}+H_{jki}+H_{kij}=0\)에 \(H_{ijk}\)를 곱해 합하면
\(\sum H_{ijk}H_{jki}=-H^2/2\), 따라서
\(\sum H_{ijk}H_{jik}=H^2/2\)다. 사차 D 수축의 세 pairing은 각각
\(E^2,D^2,D^2/2\)를 준다. 위 H 사차 평균은 원문의 식 (11d)에 적힌
\(3H^2/2\)와도 일치한다.

이제 \(Q\)의 \(t^2r^2\) 전기적 항은 \(+2E^2/m-2E^2/m=0\)으로 상쇄되어

\[
 \boxed{\langle Q\rangle
 =E^2t^4+\frac{3H^2}{m}t^2r^2
 +\frac{3E^2-3H^2+\frac32D^2}{m(m+2)}r^4}.
\]

식 (77)을 그대로 쓰면 여기에 \(4E^2t^2r^2/m\)가 추가된다.
따라서 부호 문제는 최종 부피의 전기적 계수까지 영향을 준다.

**4. 경계 이동과 부피 요소를 각각 적분한다.**

이 단계는 [원문의 식 (74),(75)](https://arxiv.org/html/1904.01034v2#SVI)가
채택한 적분용 경계식과 위·아래 cone의 환원을 입력으로 사용한다.
\(r_0=l-t\), \(r^+=r_0(1+b_1+b_2)+\cdots\)로 쓰면 필요한 항은

\[
 b_1=\frac{l^2e}{6},\qquad
 b_2=\frac{l^4e^2}{24}-\frac{l^3t e_2}{45}
            +\frac{l^2r_0^2}{90}(e_2+h_ah_a-e^2).
\]

여기서는 \(b_1\) 자체가 곡률 일차, \(b_2\)가 이차다.
반지름 적분을 먼저 하면 \((r^+)^m/m\)이므로 경계의 이차 기여는

\[
 \delta V_b=2\Omega_{d-2}\int_0^l r_0^m
 \left[\langle b_2\rangle+\frac{m-1}{2}\langle b_1^2\rangle\right]dt.
\]

\(\langle b_1\rangle=0\)이지만 \(\langle b_1^2\rangle\)는 남는다.
이를 포함하면 대괄호는

\[
 \frac{l^4E^2}{36m}-\frac{l^3tE^2}{45m}
 +\frac{l^2r_0^2(mE^2+\frac32H^2)}{90m(m+2)}.
\]

부피 요소의 곡률 이차 기여는 이미 이차이므로, 그 항의 경계에는 \(r_0\)를 쓰면 된다.

\[
 \delta V_g=-\frac{\Omega_{d-2}}{90}
    \int_0^l dt\int_0^{l-t}dr\,r^{m-1}\langle Q\rangle.
\]

시간 적분에 쓰는 공식은

\[
 \int_0^l t^a(l-t)^b dt
 =l^{a+b+1}\frac{\Gamma(a+1)\Gamma(b+1)}{\Gamma(a+b+2)}.
\]

\(\delta V_b=2\Omega_{d-2}l^{m+5}\mathcal B\),
\(\delta V_g=-\Omega_{d-2}l^{m+5}\mathcal M/90\)로 정의하면

\[
 \mathcal B=
 \frac{E^2}{36m(m+1)}-\frac{E^2}{45m(m+1)(m+2)}
 +\frac{mE^2+\frac32H^2}{90m(m+2)(m+3)},
\]
\[
 \begin{aligned}
 \mathcal M={}&\frac{24E^2}{m(m+1)(m+2)(m+3)(m+4)(m+5)}\\
 &+\frac{6H^2}{m(m+2)(m+3)(m+4)(m+5)}\\
 &+\frac{3E^2-3H^2+\frac32D^2}{m(m+2)(m+4)(m+5)}.
 \end{aligned}
\]

특히 자기적 항만 남기면 계산이 짧아진다.

\[
 \frac{\delta V_{b,H}}{\Omega_{d-2}l^{m+5}}
 =\frac{H^2}{30m(m+2)(m+3)},
\]
\[
 \frac{\delta V_{g,H}}{\Omega_{d-2}l^{m+5}}
 =\frac{(m+1)H^2}{30m(m+2)(m+3)(m+4)(m+5)}.
\]

이를 더하고 아래 공통분모로 통분하면 \(3(d+2)(d+6)H^2\)가 된다.
이 과정은 전기적 부호를 고치든 그대로 두든 같다.

**5. 각 단계의 분자를 비교한다.**

\[
 \mathcal D_d=90(d-1)(d+1)(d+2)(d+3)(d+4),\qquad
 \delta V=\Omega_{d-2}l^{d+4}\frac{P}{\mathcal D_d}.
\]

\(2\mathcal B-\mathcal M/90\)을 정리한 식 (74),(76)의 재적분 분자는

\[
 \boxed{P_{74,76}=(7d^3+58d^2+146d+116)E^2
 +3(d+2)(d+6)H^2-\tfrac32(d+2)D^2}.
\]

식 (77)의 추가 \(4E^2t^2r^2/m\)를 적분하면

\[
 \Delta\mathcal M_{77}
 =\frac{8E^2}{m(m+2)(m+3)(m+4)(m+5)},
 \qquad P_{74,77}=P_{74,76}-8E^2.
\]

따라서 비교는 다음과 같다. [원문 식 (79)](https://arxiv.org/html/1904.01034v2#SVI)
의 세 계수 중 전기적 계수는 인쇄된 식 (77)을 적분한 값과 맞고,
자기적 계수는 맞지 않는다.

| 경로 | \(E^2\) 계수 | \(H^2\) 계수 | \(D^2\) 계수 |
|---|---|---|---|
| 식 (74),(76) 재적분 | \(7d^3+58d^2+146d+116\) | \(3(d+2)(d+6)\) | \(-3(d+2)/2\) |
| 인쇄된 식 (74),(77) 재적분 | \(7d^3+58d^2+146d+108\) | \(3(d+2)(d+6)\) | \(-3(d+2)/2\) |
| 인쇄된 최종 식 (79) | \(7d^3+58d^2+146d+108\) | \(6(d+2)(d+6)\) | \(-3(d+2)/2\) |

따라서 두 불일치는 서로 별개다.

\[
 P_{79}-P_{74,76}=-8E^2+3(d+2)(d+6)H^2.
\]

두 지수 자기적 텐서 B를 쓰는 4차원 관례에서는 \(H^2=2B^2\)다.
그러나 이 비교는 처음부터 끝까지 원문의 세 지수 H를 쓰고, 원문의 식 (11d)도
같은 H 평균을 준다. H와 B 사이의 표기 변경으로 위의 단계별 불일치를 없앨 수 없다.
저자의 내부 계산에서 왜 이 계수가 발생했는지는 이 비교만으로 단정하지 않는다.

**6. 평면파에서 기존 차이가 정확히 복원된다.**

첫 metric에서 직접 계산하면 \(R_{uxux}=-A\), \(R_{uyuy}=A\),
\(R_{\mu\nu}=0\), \(\nabla R_{\alpha\beta\gamma\delta}=0\)다.
축의 직교 시간벡터 \(e_0=(\partial_u-\partial_v)/\sqrt2\)를 택하면

\[
 E=\operatorname{diag}(-A/2,A/2,0),\qquad
 E^2=A^2/2,\quad H^2=A^2,\quad D^2=2A^2.
\]

\(d=4\)에서 \(\mathcal D_4=453600\), \(V_0=\Omega_2l^4/6\),
\(l^4=1/4\)이므로 정규화는

\[
 c=\frac{V/V_0-1}{A^2}\bigg|_{A\to0}
 =\frac32\frac{P/A^2}{453600}.
\]

| 경로 | \(P/A^2\) | 상대 부피의 \(A^2\) 계수 c |
|---|---:|---:|
| 정확 단면 적분, 독립 정상좌표 계산, 식 (74),(76) | 1200 | \(1/252=300/75600\) |
| 식 (74),(77)를 그대로 적분 | 1196 | \(299/75600\) |
| 식 (79)에 직접 대입 | 1376 | \(43/9450=344/75600\) |

\[
 \boxed{\frac{43}{9450}-\frac1{252}
 =\underbrace{-\frac1{75600}}_{\text{전기적 부호}}
 +\underbrace{\frac1{1680}}_{\text{자기적 계수 두 배}}
 =\frac{11}{18900}}.
\]

따라서 기존 평면파의 \(1/252\)를 유지한다. 단면 적분과 정상좌표 부피 계산이
같은 답을 주므로 좌표나 전체·반 고유시간의 혼동으로 설명되지 않는다.
이 family는 \(\nabla R=0\)이므로 누락된 곡률 미분항을 추가해서 차이를
해소할 수도 없다. 이 판단은 시공간 4-부피에 대한 것이며 최대 공간 초곡면의
3-부피와 비교한 것이 아니다.

**검산과 기록의 범위.**

[check_diamond_geometry.py](check_diamond_geometry.py)의 실행 가능한 검산을 확장했다.
기존 metric·정상좌표 계산에 더하여, 전기적 강도 e와 자기적 강도 b를 독립적으로
변화시키는 4차원 대수적 Weyl 텐서도 사용했다. 이는 실제 새 진공해의 구성은 아니다.
이 텐서는 Bianchi 항등식과 Ricci trace가 0인 것을 확인한 뒤 수축했다.

\[
 E=e\operatorname{diag}(1,2,-3),\qquad
 B=b\begin{pmatrix}0&1&0\\1&0&2\\0&2&0\end{pmatrix},
\quad (E^2,H^2,D^2)=(14e^2,20b^2,56e^2).
\]

인쇄된 식 (77)의 모든 짝수 항을 별도로 입력하여, 식 (76)과의 차이가
\(4e^2r^2t^2(n_x^2+4n_y^2+9n_z^2)\)뿐임을 기호식으로 확인했다.
일반 m의 beta 적분으로 위 분자 다항식과 두 계수 차이도 검사한다.
SymPy 1.14.0 / Python 3.11.9에서 모든 assertion과 Ruff 검사를 통과했다.
출력과 script의 LF-normalized SHA256은 [geometry_audit.json](geometry_audit.json)에 있다.

```powershell
python docs/paper/paper_a/window_response/check_diamond_geometry.py
```

계산 경로는 여러 개지만 작성·검산 주체는 같은 AI다. 독립 연구자의 검토나
정정문 수용이 완료된 결과로 표현하지 않는다. 현재 원고 PDF와 기존 9월 22일
bundle은 이 상세 추적 전의 버전이며, 이 문서와 검산 출력이 후속 검토 기록이다.
기존 동결 예측, 확인실험 자료와 정식 판정은 바꾸지 않았다.

후속 [일반화 유도](vacuum_generalization_derivation.md)는 이 문서가 입력으로
삼은 경계를 Jacobi 방정식과 지오데식 작용에서 직접 다시 구하고, 곡률 미분항과
두 null 경계의 접합부를 검토했다. 별도의 Kasner·비상수 평면파 계산 및
우주상수 확장도 그 기록에 있다. 여기의 조건부 중간식 추적과 후속 일반 유도를
구분해 읽어야 한다.
