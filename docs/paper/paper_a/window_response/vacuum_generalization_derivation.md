# 진공 다이아몬드 부피의 일반화: 국소 점근식과 독립 검산

2026-09-22. 상태: 일반 차원의 국소 전개를 직접 유도하고 아래 기호식·특수해
검산을 완료했다. 외부 동료검토 전의 연구 기록이다. 표본 확인실험이나 새
사전등록이 아니며 기존 동결 결과를 변경하지 않는다.

**확보한 범위는 매끄러운 진공 시공간의 충분히 작은 다이아몬드다.**
앞선 평면파의 `1/252`가 모든 기하의 동일한 수치가 되는 것은 아니다.
그 계수를 포함하는 아래의 Weyl 불변량 공식이 일반화된다. 큰 다이아몬드,
특이점, 전역 인과구조까지 한 유한 차수 식으로 기술한다는 주장은 하지 않는다.

**1. 명제와 정규화.**

차원 \(d\ge4\)인 매끄러운 시간지향 Lorentz 시공간에서 O의 열린 근방이
\(R_{ab}=0\)이라고 하자. \(g(U,U)=-1\)인 미래지향 벡터 U를 O에 고정하고
\(p_\ell=\exp_O(-\ell U)\), \(q_\ell=\exp_O(\ell U)\)로 놓는다.
다이아몬드 \(I^+(p_\ell)\cap I^-(q_\ell)\)가 인과적으로 볼록한 정상 근방 안에
들어가고 빛의 경계에 conjugate point가 없을 만큼 \(\ell>0\)가 작다고 가정한다.
\(2\ell\)은 끝점 사이의 고유시간이며 근사의 극한에서 O와 U를 고정한다.

O의 직교틀에서 \(e_0=U\), \(\eta=\mathrm{diag}(-1,1,\ldots,1)\)로 두고

\[
 E_{ij}=C_{0i0j},\qquad H_{ijk}=C_{0ijk},\qquad D_{ijkl}=C_{ijkl},
\]

각 제곱은 공간 지수를 모두 수축한 양의 제곱합으로 정의한다.
\(V_0=2\Omega_{d-2}\ell^d/[d(d-1)]\)일 때 직접 유도한 결과는

\[
 \boxed{\frac{V(\ell)}{V_0(\ell)}
 =1+\mathcal C_d\ell^4+O(\ell^5)},
\]
\[
 \boxed{\mathcal C_d=
 \frac{d\{(7d^2+44d+58)E^2+3(d+6)H^2-\frac32D^2\}}
 {180(d+1)(d+3)(d+4)}}. \tag{G1}
\]

이는 앞선 중간식 재적분 다항식에서 공통인자 \(d+2\)를 약분한 식이다.
이번에는 경계식을 가정하는 대신 계량의 Jacobi 전개에서 다시 유도했다.
나머지는 고정된 기하·O·U에서의 국소 점근항이며, 모든 진공해에 공통인
유한 \(\ell\)의 오차 상한이나 허용 크기를 주지 않는다.

4차원에서는 \(D^2=4E^2\), \(H^2=2B^2\)이므로

\[
 \boxed{\frac{V(\ell)}{2\pi\ell^4/3}
 =1+\frac{17E^2+3B^2}{630}\ell^4+O(\ell^5)}. \tag{G2}
\]

여기서 B는 두 지수 자기적 Weyl 텐서다. \(T=2\ell\) 표기로 바꾸면
\(\ell^4/630=T^4/10080\)이다. 전자기장을 뜻하는 E와 B가 아니다.

**2. Jacobi 방정식에서 계량과 부피 요소를 얻는다.**

원점 O의 정상좌표 벡터를 Y라고 하자. 아래 \(K_j(Y)\)는
\((\nabla_Y^jR)_{aYbY}\)의 한 지수를 올린 행렬이며 차수는 \(j+2\)다.
방사 지오데식 \(s\mapsto\exp_O(sY)\)를 따라 평행이동한 틀에서
Jacobi 행렬은 \(J(0)=0,\ J'(0)=I\)를 만족하고

\[
 J''+(K_0+sK_1+\tfrac12s^2K_2+\cdots)J=0.
\]

계수별로 두 번 적분하면

\[
 J(s)=sI-\frac{s^3}{6}K_0-\frac{s^4}{12}K_1
 +s^5\left(\frac{K_0^2}{120}-\frac{K_2}{40}\right)+O(|Y|^5).
\]

\(g(Y)=J(1)^T\eta J(1)\)이므로, 내린 지수로

\[
 g_{ab}(Y)=\eta_{ab}-\frac13(K_0)_{ab}-\frac16(K_1)_{ab}
 -\frac1{20}(K_2)_{ab}+\frac2{45}(K_0)_{ac}(K_0)^c{}_b+O(|Y|^5). \tag{G3}
\]

이 식은 수정된 [Brewin arXiv:0903.2087v3, 식 (11.1)](https://arxiv.org/html/0903.2087v3)
과도 일치한다. 해당 논문은 구버전의 고차 계산 오류를 명시적으로 정정했다.
이 기록은 구버전의 arc-length 식을 그대로 전제하지 않는다.

진공에서는 \(\mathrm{tr}K_0=\mathrm{tr}K_1=\mathrm{tr}K_2=0\)이다.
\(\log\det J\)의 곡률 제곱 계수는 \(1/120-1/72=-1/180\)이므로

\[
 \sqrt{-g(Y)}=1-\frac1{180}\mathrm{tr}(K_0^2)+O(|Y|^5). \tag{G4}
\]

여기서 미분항이 없어지는 이유는 \(R_{ab}\)가 **열린 근방 전체에서** 0이기 때문이다.
한 점에서만 Ricci가 0이라는 조건으로는 이 명제가 성립하지 않는다.

**3. 지오데식 작용으로 null 경계를 직접 유도한다.**

위쪽 끝점 \(q=\ell U\), 일반 점 \(X=tU+rn\), \(\Delta=X-q\)를 잡는다.
\(n\cdot n=1\)이며 시간 성분은 0이다. 매개변수 구간이 [0,1]인 지오데식의
에너지 작용 \(\mathscr S=\int g_{ab}\dot Y^a\dot Y^b ds\)는 signed squared
distance, 즉 Synge 함수의 \(2\sigma\)다. null 경계는 \(\mathscr S=0\)이다.

평탄 직선은 \(Y_0(s)=q+s\Delta\)다. 공변 벡터

\[
 P_a=R_{a q\Delta q},\qquad W_a=R_{a\Delta\Delta q}
\]

를 정의하면 곡률 일차의 경로 변위는
\(Z^a(s)=s(s-1)W^a/3\)이다. 이는 양 끝에서 0이며 선형화한 지오데식
방정식 \(Z''{}^a=2W^a/3\)을 만족한다.
직선 위 계량 사차항과 경로의 정상성에 따른 보정을 합하면

\[
 \mathscr S_{R^2}
 =\frac2{45}\int_0^1(P+sW)^2ds-\int_0^1 Z'^2ds
 =\boxed{\frac{2P^2+2P\cdot W-W^2}{45}}. \tag{G5}
\]

두 번째 적분은 \(W^2/27\)이다. Lorentz 부호를 포함한 수축이며
양의 Euclidean norm으로 바꿀 수 없다.

\(e=E_{ij}n^in^j\), \(e_2=(E^2)_{ij}n^in^j\),
\(h_a=H_{iaj}n^in^j\), \(e\!h=(E_{ij}n^j)h_i\)로 쓰면

\[
 \mathscr S_0=-(\ell-t)^2+r^2,\qquad
 \mathscr S_R=-\frac{\ell^2r^2}{3}e,
\]
\[
 \mathscr S_{R^2}=\frac{\ell^2r^2}{45}
 \left[(4\ell t-t^2-\ell^2)e_2+2r(t-2\ell)e\!h
                         +r^2(e^2-h_ah_a)\right]. \tag{G6}
\]

미분항은 (G3)의 나머지 항을 직선 위에서 적분하면 된다. 아래에서
\(e_{;a}=n^in^j\nabla_a C_{0i0j}\)이며, 미분은 완전한 Weyl 텐서에 작용한다.
별도의 가속 관측자장이나 회전하는 틀의 성분을 미분하는 것이 아니다.
\(\nabla_n=n^k\nabla_k\)로 두면

\[
 \mathscr S_{\nabla R}=-\frac{\ell^2r^2}{12}
                  [ (\ell+t)\nabla_0+r\nabla_n]e,
\]
\[
 \mathscr S_{\nabla^2R}=-\frac{\ell^2r^2}{60}
 \left[(\ell^2+\ell t+t^2)\nabla_0^2
 +r(\ell+2t)\nabla_{(0}\nabla_{n)}+r^2\nabla_n^2\right]e. \tag{G7}
\]

곡률 미분에 의한 추가 경로 변위는 작용의 6차항에 기여하지 않는다.
그 최초 교차항은 7차 이상이다. 괄호 속 미분은 순서를 보존하며
임의로 교환하지 않는다.

이제 \(r_0=\ell-t\), \(r^+=r_0(1+b_2+b_3+b_4)+O(\ell^6)\)로 쓴다.
이번 문서의 첨자 2,3,4는 **\(\ell\)의 상대 차수**이며, 앞선 문서의 곡률
일차·이차 표기와 다르다. null 방정식의 계수를 순서대로 풀면

\[
 b_2=\frac{\ell^2e}{6},\qquad
 b_3=\frac{\ell^2}{24}[(\ell+t)\nabla_0+r_0\nabla_n]e,
\]
\[
 \begin{aligned}
 b_{4,R^2}={}&\frac{\ell^4e^2}{24}
 +\frac{\ell^2}{90}\left[(t^2-4\ell t+\ell^2)e_2
       -2r_0(t-2\ell)e\!h+r_0^2(h_ah_a-e^2)\right],\\
 b_{4,\nabla^2R}={}&\frac{\ell^2}{120}
 \left[(\ell^2+\ell t+t^2)\nabla_0^2
 +r_0(\ell+2t)\nabla_{(0}\nabla_{n)}+r_0^2\nabla_n^2\right]e.
 \end{aligned} \tag{G8}
\]

이 계산은 물리적 경계의 구면대칭을 가정하지 않는다.
각 방향 n마다 달라지는 경계를 유지하며, 다음 단계에서 필요한 평균을 취한다.

**4. 생략했던 미분항과 두 cone의 접합부를 닫는다.**

공간 차원 \(m=d-1\)에서 구면 평균은

\[
 \langle n^in^j\rangle=\delta^{ij}/m,\qquad
 \langle n^in^jn^kn^l\rangle
 =\frac{\delta^{ij}\delta^{kl}+\delta^{ik}\delta^{jl}+\delta^{il}\delta^{jk}}
 {m(m+2)}.
\]

\(\langle b_2\rangle=0\)은 E의 trace가 0이기 때문이다.
\(b_3\)의 시간 미분항도 같은 trace의 미분이며 공간 미분항은 n에 대해 홀수다.
\(b_{4,\nabla^2R}\)의 시간·시간 항은 trace, 시간·공간 항은 홀수 평균으로
없어진다. 남는 공간·공간 항의 세 수축은

\[
 \delta^{ij}\delta^{kl}\nabla_k\nabla_l C_{0i0j},\quad
 \sum_{ij}\nabla_i\nabla_j C_{0i0j},\quad
 \sum_{ij}\nabla_j\nabla_i C_{0i0j}.
\]

첫 항은 Ricci trace의 이차 미분이라 0이다. 뒤 두 항은 각각
안쪽 미분이 Weyl 텐서의 발산이다. 진공의 contracted Bianchi 항등식
\(\nabla^aC_{abcd}=0\)으로 각각 0이다. 공간합을 시공간합으로 보충할 때
추가되는 시간 항은 Weyl의 반대칭인 한 쌍에 U가 두 번 들어가므로 0이다.
따라서 **미분을 교환해 곡률 제곱 항을 만들어낼 필요 없이**

\[
 \boxed{\langle b_3\rangle=\langle b_{4,\nabla^2R}\rangle=0}. \tag{G9}
\]

이 논증은 \(\nabla C=0\) 또는 \(\nabla\nabla C=0\)을 가정하지 않는다.

위·아래 cone이 만나는 시각을 \(\tau(n)\)라 하자. 일차 곡률 경계는 두 cone에서
같으므로 \(\tau=O(\ell^4)\)이다. 실제 접합면 대신 t=0에서 적분을 나누면
두 적분의 차이는 접합부에서 서로 상쇄된다. 남는 오차는

\[
 O(\ell^{d-2}\tau^2)=O(\ell^{d+6}),
\]

이므로 관심 있는 \(\ell^{d+4}\) 계수에 영향을 주지 않는다.
접합 시각의 평균만 보고 선형 항을 버리는 것보다 강한, 두 영역의 차이에 대한 추정이다.

**5. 남는 항을 적분한다.**

\(Q=R_{aXbX}R^{a}{}_{X}{}^{b}{}_{X}\)에 대해

\[
 \langle Q\rangle=E^2t^4+\frac{3H^2}{m}t^2r^2
 +\frac{3E^2-3H^2+\frac32D^2}{m(m+2)}r^4.
\]

경계에서 \(e\!h\)는 홀수이고,
\(\langle h_ah_a\rangle=3H^2/[2m(m+2)]\)다. 반지름 적분의 비선형 항도
포함하면

\[
 \left\langle b_{4,R^2}+\frac{m-1}{2}b_2^2\right\rangle
 =\frac{\ell^4E^2}{36m}-\frac{\ell^3tE^2}{45m}
 +\frac{\ell^2(\ell-t)^2(mE^2+\frac32H^2)}{90m(m+2)}.
\]

따라서

\[
 \begin{aligned}
 \delta V={}&2\Omega_{d-2}\int_0^\ell(\ell-t)^m
       \left\langle b_{4,R^2}+\frac{m-1}{2}b_2^2\right\rangle dt\\
 &-\frac{\Omega_{d-2}}{90}\int_0^\ell dt\int_0^{\ell-t}
                          r^{m-1}\langle Q\rangle dr.
 \end{aligned} \tag{G10}
\]

\(\int_0^\ell t^a(\ell-t)^b dt
=\ell^{a+b+1}\Gamma(a+1)\Gamma(b+1)/\Gamma(a+b+2)\)로 적분하고
\(V_0\)로 나누면 (G1)이 나온다. 미분항과 경계 변형까지 포함한 유도이며
특정 Petrov 유형, 정적 기하, 평면파, 물리적 구면대칭을 가정하지 않았다.

**6. 원래 Kasner 좌표에서 별도로 계산한다.**

위 증명과 다른 계산 경로로 실제 진공해

\[
 ds^2=-d\tau^2+\sum_{i=1}^3\tau^{2p_i}dx_i^2,
 \qquad\sum_i p_i=\sum_i p_i^2=1
\]

를 사용했다. 끝점은 \((\tau,\mathbf x)=(1\pm\ell,0)\)이다.
곡률이 시간에 따라 변하므로 \(\nabla C=0\) 특수 경우의 검산이 아니다.
\(\sqrt{-g}=\tau\)이고 x=0은 고유시간으로 매개화한 지오데식이다.

\(\tau=1\)에서 출발한 미래 null 면을 \(\tau=F(\mathbf x)\)라 하면

\[
 \sum_i F^{-2p_i}(\partial_iF)^2=1.
\]

\(r^2=\sum x_i^2\), \(z_i=x_i^2/r^2\), \(z_3=1-z_1-z_2\)로 두고
\(F=1+\sum_{k=1}^5r^kf_k(z_1,z_2)+O(r^6)\), \(f_1=1\)을 대입한다.
\(r^kf\)의 i방향 미분은

\[
 \partial_i(r^kf)=r^{k-1}n_i
 \left[kf+2\partial_{z_i}f-2\sum_jz_j\partial_{z_j}f\right].
\]

위 식에서 f는 \(z_1,z_2\)의 함수로 연장하고 \(\partial_{z_3}f=0\)으로 둔다.
eikonal 식의 \(r^{k-1}\) 항에 새로 들어오는 \(f_k\)의 계수는 \(2k\)다.
따라서 앞서 정해진 잔차를 \(-2k\)로 나누어 모든 \(f_2,\ldots,f_5\)를 얻는다.
부피의 예상 계수를 이 재귀에 입력하지 않는다.

시작 시각이 \(\tau_e\)인 면은 Kasner의 자기유사성으로
\(T(\mathbf x)=\tau_e F(\tau_e^{p_i-1}x_i)\)다.
과거 면에는 r의 부호를 반대로 적용한다. \(\mathbf x=\ell\rho n\)일 때
아래쪽 도달시각을 \(T_-=1+\sum_{j=1}^5\ell^ja_j(\rho,n)\)라 쓰면
위쪽 면은 \(T_+(\ell)=T_-(-\ell)\)이며 \(a_1=\rho-1\)이다.

시간 적분을 먼저 수행하면

\[
 \int_{T_-}^{T_+}\tau d\tau
 =\tfrac12(T_+^2-T_-^2)
 =\ell H_0+\ell^3H_2+\ell^5H_4+O(\ell^7),
\]
\[
 H_0=2(1-\rho),\quad H_2=-2(a_3+a_1a_2),\quad
 H_4=-2(a_5+a_2a_3+a_1a_4).
\]

교차 반지름도 이동한다. \(\rho_e=1+\ell^2b+O(\ell^4)\),
\(b=-a_3(1,n)\)이므로 상대 부피의 계수는

\[
 c_2=6\left\langle\int_0^1\rho^2H_2d\rho\right\rangle,\qquad
 c_4=6\left\langle\int_0^1\rho^2H_4d\rho+b^2\right\rangle. \tag{G11}
\]

끝의 \(b^2\)는 경계 이동이 만드는 항으로 생략하지 않았다. 이 방법은
정상좌표 곡률 전개나 Wang의 경계식을 사용하지 않는다.

| Kasner 지수 | 직접 계산 \(c_2\) | 직접 계산 \(c_4\) | (G2)의 예측 |
|---|---:|---:|---:|
| \((1,0,0)\), 평탄 Milne 대조군 | 0 | 0 | 0 |
| \((-1/3,2/3,2/3)\) | 0 | \(68/8505\) | \(68/8505\) |
| \((-2/7,3/7,6/7)\) | 0 | \(68/12005\) | \(68/12005\) |

모두 유리수 기호식으로 일치한다. Kasner의 \(B=0\),
\(E_{ii}=p_i(1-p_i)\)로부터 두 곡선 경우의 \(E^2\)는 각각
\(8/27\), \(72/343\)이다. 각 주축의 null 도달시각
\([1+(1-p_i)r]^{1/(1-p_i)}\) 또는 \(p_i=1\)일 때 \(e^r\)와도 비교했다.
Milne 대조군은 전체 각도의 정확한 평탄 null 면과 추가로 일치한다.

**7. 비상수 평면파와 고차원의 독립 성분.**

비상수 profile \(A(u)=a_0+a_1u+a_2u^2\)를 가진
\(ds^2=2du\,dv+A(u)(x^2-y^2)du^2+dx^2+dy^2\)도 정확한 진공해다.
끝점을 \((u,v)=(\pm L,\mp L)\)에 놓으면 \(\ell=\sqrt2L\)이다.
횡방향 Jacobi 전달함수는 \(D_\pm''=\pm A(u)D_\pm\)이고
한 끝점에서 \(D=0,D'=1\)을 만족한다.

\(u=Lx\)로 축척한 전달함수를 Volterra 적분으로 전개하면
\(D_+D_-\)의 \(L^2,L^3\) 항은 0이며 \(L^4\) 항은
\(-a_0^2(x-y)^6/90\)이다. \(a_1\)과 \(a_2\)에 선형인 항이 소거된다.
세 전달함수의 Wronskian 관계로 diamond 단면을 적분하여

\[
 \frac V{V_0}=1+\frac{a_0^2}{63}\ell^4+O(\ell^5)
\]

를 얻었다. profile의 일차·이차 변화가 있어도 (G2)와 일치한다.
\(a_0=0\)인 경우 이 선도 계수는 0이며, 이것만으로 유한 영역 전체가
평탄하다고 결론 내릴 수 없다.

4·5·6차원에서는 전기적 성분, 자기적 성분, 혼합 성분 및 가능한 차원에서
순수 공간 Weyl 성분을 독립적으로 넣은 **대수적 곡률 텐서** 11개를 검산했다.
각 텐서의 Bianchi 항등식과 Ricci trace를 먼저 검사하고, (G5)에서 경계를
계산한 뒤 각도·반지름·시간을 적분했다. 모두 (G1)과 일치했다.
이 대수적 텐서를 새로운 시공간 해 11개로 세지 않는다.

고차원의 D 항에는 별도의 해석도 있다. \(ds^2=-dt^2+h\)이고 h가
Ricci-flat Riemannian metric이면 E=H=0이며 단면은 h의 지오데식 공이다.
그 공의 Jacobi 부피 요소를 적분하면

\[
 \frac{V_{\rm ball}(r)}{\omega_mr^m}
 =1-\frac{D^2r^4}{120(m+2)(m+4)}+O(r^5).
\]

이를 \(r=\ell-|t|\)로 한 번 더 적분한 결과는
\(-dD^2/[120(d+1)(d+3)(d+4)]\)로, (G1)의 D 계수와 같다.
따라서 4차원에서의 선도 부피 증가를 모든 차원으로 확대하면 안 된다.

**8. 우주상수를 포함하는 진공의 확장.**

\(R_{ab}=\lambda g_{ab}\)에서 \(\lambda\)가 상수이고
\(\kappa=\lambda/(d-1)=2\Lambda/[(d-1)(d-2)]\)라고 하자.
곡률은 \(R_{abcd}=C_{abcd}+\kappa(g_{ac}g_{bd}-g_{ad}g_{bc})\)로 분해된다.
Weyl의 발산은 여전히 0이고
Ricci의 미분은 0이므로 앞의 미분항 소거가 유지된다.
\(\kappa C\)에 선형인 평균 scalar는 Weyl의 trace로 환원되어 0이다.
그러므로 Weyl 제곱 항은 (G1)의 \(\mathcal C_d\)이며, 나머지는 상수곡률에서
정할 수 있다.

상수곡률 static patch의 null 반지름은
\(r(s)=\tanh(\sqrt\kappa s)/\sqrt\kappa\)이고, \(\kappa<0\)에서는
해석적 연속으로 tan 식이 된다. 정확 부피는

\[
 V_\kappa=\frac{2\Omega_{d-2}}{d-1}\int_0^\ell r(s)^{d-1}ds.
\]

직접 Taylor 적분한 결과, Einstein 진공의 국소 전개는

\[
 \boxed{\frac V{V_0}=1-\frac{d(d-1)}{3(d+2)}\kappa\ell^2
 +\left[\frac{d(d-1)(5d+2)}{90(d+4)}\kappa^2+\mathcal C_d\right]\ell^4
 +O(\ell^5)}. \tag{G12}
\]

4차원에서는

\[
 \frac V{V_0}=1-\frac{2\Lambda}{9}\ell^2
 +\left[\frac{11\Lambda^2}{270}+\frac{17E^2+3B^2}{630}\right]\ell^4
 +O(\ell^5).
\]

여기서 E,H,D는 **Weyl** 성분이다. 전체 Riemann 성분을 대신 넣으면
우주상수 항을 중복 계산한다. 알려진 4차원 de Sitter 정확식
([Gibbons–Solodukhin, 식 (21)](https://arxiv.org/html/hep-th/0703098))과도 맞는다.
4·5차원에서 상수 \(\kappa\)와 혼합 Weyl 텐서를 동시에 넣고 Ricci 제곱 및
부피 요소의 경계 이동 항까지 직접 계산한 두 추가 검산도 통과했다.

**9. 무엇이 끝났고 무엇이 남았는가.**

| 질문 | 이 기록의 결론 |
|---|---|
| \(\nabla C=0\)인 평면파에서만 성립하는가? | 아니다. 미분항을 명시적으로 유지한 일반 유도와 비상수 곡률의 독립 계산이 있다. |
| 임의의 매끄러운 Ricci-flat 진공 근방에 적용되는가? | 명제 1의 국소성과 작은 크기 조건에서 (G1)이 유도된다. |
| 우주상수가 있는 진공도 포함하는가? | Einstein 조건에서 (G12)로 확장된다. |
| 유한 크기와 오차율까지 인증했는가? | 아니다. 특정 기하의 나머지 상한과 허용 크기는 별도 문제다. |
| Schwarzschild·Kerr를 각각 직접 적분했는가? | 아직 아니다. 매끄러운 진공 근방은 일반 유도의 적용 대상이지만, 별도 exact-coordinate 검산 대상은 위 Kasner와 평면파다. |
| 임의의 전역 인과구조·caustic·특이점을 포함하는가? | 아니다. 정상 근방의 국소 정리다. |
| 원문 수정식의 최초성과 외부 판정이 확정됐는가? | 아니다. 문헌 검색은 제한적이며 독립 연구자의 검토는 아직 없다. |

남은 과학적 위험은 일반 유도 자체의 외부 검토, 구체적 유한 크기에서의 나머지
상한, 그리고 문헌상 동일 계수·정정의 존재 여부다. 이 기록은 예제 일치만으로
보편성을 주장하지 않는다. 보편성의 근거는 (G3)–(G10)의 차원 독립 유도이며,
실행 가능한 검산은 그 유도를 반박하거나 전사 오류를 찾는 별도의 수단이다.

[check_vacuum_generalization.py](check_vacuum_generalization.py)는 원래 metric의
Kasner eikonal, 변화하는 평면파의 전달함수, 대수적 곡률과 정확 상수곡률식을
유리수·기호식으로 계산한다. 코드에서 요구하는 계수와 직접 계산을 비교하며
새 표본이나 과거 확인 자료를 읽지 않는다. 기록은
[vacuum_generalization_audit.json](vacuum_generalization_audit.json)에 있고,
소스의 LF-normalized SHA256, Python/SymPy 버전과 실제 실행 시간이 포함된다.

```powershell
python docs/paper/paper_a/window_response/check_vacuum_generalization.py
```

SymPy 1.14.0 / Python 3.11.9에서 검산했다. 모든 assertion과 Ruff 검사가
통과했다. 작성·검산 주체는 같은 AI이며 독립 좌표 계산과 독립 인간 검토를
구분한다. 현재 원고 PDF와 기존 Git bundle에는 이 후속 유도가 아직 반영되지 않았다.

직접 비교한 문헌은 Wang arXiv:1904.01034v2의 식 (74)–(79),
Brewin arXiv:0903.2087v3의 (11.1),(11.3),(11.21),
Gibbons–Solodukhin hep-th/0703098의 de Sitter 정확식이다.
Wang 논문의 최종 계수 불일치의 상세 추적은
[앞선 기록](volume_discrepancy_derivation.md)에 보존했다.

후속 [전역 검토](global_diamond_exact_cases.md)는 큰 다이아몬드, 위상, 초점 및
특이점 접근을 다룬다. 국소 곡률 jet가 모두 같아도 전역 부피가 달라지는
평탄한 반례와, 실제 매끄러운 Kasner 다이아몬드의 유리수 하한 반례를 얻었다.
따라서 이 문서의 국소 점근식을 임의 크기의 정확식으로 확대할 수 없다.
후속 기록에는 조건부 정확 적분식과 부피 발산 경계도 있다.
