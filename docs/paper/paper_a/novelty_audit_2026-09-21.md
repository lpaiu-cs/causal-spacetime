# Paper A 보완연구 G0 — 차별성 감사

검토일: 2026-09-21. 저장소 기준: `b27a8dfd0b3f37fa9b01ffaa850c071b0c1c029c`.
대상: [보완연구 계획](supplementary_research_plan_2026-09-21.md)의 평면파 관계 비율 연구.
종류: 문헌·기존 결과 대조와 적용 조건 검토. 새 캠페인 결과나 사전등록이 아니다.

## 1. 판정: 축소 진행 — 차별성 확정은 보류

**G0 문헌 조사는 완료했지만, 본 캠페인으로 넘어갈 만큼 새 기여가 확정되지는 않았다.**
G1 중 응답식 유도를 먼저 검토할 근거는 남는다. 3,840-poset 파일럿과
19,200-poset 확인 캠페인은 아직 시작하지 않는다.

넓은 기여 후보들은 이미 선행연구 또는 이 저장소가 다뤘다. 남겨둘 질문은 다음이다.

> 4차원 상수 진공 평면파에서, 명시한 관측 영역에 뿌린 한 인과집합의
> 관계 비율은 조석 세기에 어떻게 반응하는가? 영역별 첫 비영 응답 계수와
> 공유 사건 공분산을 계산하여, 새 세기·N 조건의 판별 성능을 예측할 수 있는가?

이는 **새로운 시공간 이론, 최초 Weyl 검출, 일반적인 곡률 복원**의 주장이 아니다.
알려진 인과구조를 대상으로 관측 영역이 검출 성능을 어떻게 바꾸는지 설명하는
제한된 방법 연구다. 같은 단어 조합을 찾지 못했다는 것만으로 최초성을 주장하지 않는다.

| 기여 후보 | 이번 판정 | 근거 |
|---|---|---|
| 인과순서에서 곡률 정보를 얻는다 | 새 기여에서 제외 | chain 곡률 보정, 그래프 관측량, Weyl 관련 선행연구 |
| 유한 표본·잡음·대칭성이 측정을 제한한다 | 새 기여에서 제외 | Roy 등, Kambor–Nomaan X, Ghosh 등 |
| 관계쌍은 독립 표본이 아니다 | 표준 방법으로 사용 | U-statistic 이론, Roy 등의 중첩 chain 계산, 기존 P14 |
| 대칭 영역에서는 선형 응답이 상쇄될 수 있다 | 이것만으로는 부족 | P14 설계와 P3-E의 기존 설명·탐색 |
| 더 많은 N·진폭에서 높은 AUC를 보인다 | 이것만으로는 부족 | 기존 양성 결과의 범위 확대이며 원인 설명을 대신하지 못함 |
| 정확한 모형에서 영역별 응답 계수·오차를 계산하고 미관측 성능을 예측한다 | **조건부로 유지할 후보** | 아래 비교 범위에서 동일한 계산·검증을 찾지 못함; 아직 우리가 계산한 결과도 아님 |

판정의 확실성은 비대칭이다. **위의 넓은 최초성 주장을 버려야 한다는 근거는 강하다.**
좁은 질문의 비중복과 학술적 중요성은 응답식이 실제로 얻어지고 어떤 현상을
설명하는지에 달려 있다. 이것은 편집자의 구체적 심사평을 복원한 결론이 아니다.

## 2. 가장 먼저 겹치는 선행연구는 우리 저장소다

외부 문헌보다 먼저 [P14 설계](../../prereg/p14_weyl_curvature.md)와
[P3-E 결과](../../prereg/p14_probe_p3e_results.md)를 대조했다.

- P14 §§4.4–4.5는 trace-free 횡방향 응답의 상쇄와 부피의 이차 응답을 이미 논의한다.
  §5.2는 공유 점들을 가진 통계량을 U-statistic으로 다룬다.
- P3-E에는 고정된 square slab, 평균 N=300의 진폭 탐색이 있다. 순관계 변화의
  작은-A 유효 기울기는 1.977, bootstrap 95% CI [1.827, 2.173]이다.
  이는 새 계수 유도나 독립적인 점근 법칙 검증으로 다시 이름 붙일 결과가 아니다.
- P3-E의 운영점 간 부호 변화에서는 영역과 세기 등이 함께 바뀐다.
  같은 영역·같은 N에서 세기만 바꾼 영점 통과로 해석하지 않는다.
- 기존 Class C의 diamond 관측량은 좌표에서 구한 eligibility label도 사용한다
  (P14 §4.6.2). 단일 팔이라는 사실과 순서 정보만 받는다는 주장은 다르다.

따라서 이번 연구는 **기존 그림을 더 촘촘히 그리는 것**으로 끝나면 중단한다.
P14 POSITIVE, S4 CONFIRMED, S5 DETECTED 및 S6의 동결 판정은 바꾸지 않는다.
그 결과가 성립한다는 사실과 새 논문이 필요한 기여를 갖췄다는 판단은 별개다.

## 3. 주장 단위의 문헌 비교표

`f_N`은 N개 원소 중 관계가 있는 무순서 쌍의 비율이다.
우리 후보는 주어진 모형·관측자·영역에 조건부인 검출 문제다. 보편적인
curvature invariant 추정 문제와 같은 목표라고 놓고 성능을 비교하지 않는다.

| 1차 출처·확인 위치 | 실제 입력·관측량·영역 | 이미 답한 부분 | 우리 후보와 남는 차이 |
|---|---|---|---|
| [Bombelli, 2000](https://arxiv.org/html/gr-qc/0002053), §§II–VI | 주어진 Lorentzian 영역에서 고정 n의 poset 분포 `P_n(C\|G)` | 유한 표본으로 기하를 통계적으로 구별하는 틀 | 새 통계적 기하 틀을 제안할 수는 없음. 특정 진공 평면파의 영역별 2점 확률 응답 계산은 별도 문제 |
| [Roy–Sinha–Surya, 2013](https://arxiv.org/html/1212.0631), 곡률 전개 및 §6, 식 (62)–(66) | 작은 Alexandrov interval의 chain abundance; R, R00의 선도 보정 | chain으로 곡률 추정; 겹치는 chain의 공분산과 밀도에 따른 요동 차수 | Ricci 보정을 반복하는 연구가 아님. Ricci-flat 모형의 `f_N` 계수와 지정 영역의 유한-N 성능은 별도 계산 |
| [Kambor–Nomaan X, 2020](https://arxiv.org/html/2007.03835), §§3–4 | Minkowski·de Sitter·FLRW의 chain 통계와 비다양체 대조군 | 곡률·시간·차원 추정의 수치적 안정성과 오차 | 단순한 관계 비율만으로 일반 기하를 식별할 수 있다는 주장 불가. higher-chain 대조군은 기존 방법으로 취급 |
| [Wang, 2019](https://arxiv.org/html/1904.01034), §§III–V | 세 종류의 작은 causal diamond의 연속체 부피·면적 | 진공 곡률의 고차 기하 응답과 영역 정의의 중요성 | 부피나 면적의 계수를 관계 비율 계수에 대입할 수 없음. 두 점의 인과관계 적분이 필요 |
| [Alfyorov–Shnyukov, Paper 7](https://inspirehep.net/files/32388edfda552e0b4cfb6697bda66d7b), §§2.3–2.4, 5.6, 11.7 | 같은 좌표 점에서 만든 flat/curved 두 Hasse diagram의 변화와 위치별 strata; electric Weyl 제곱에 대한 CJ | 평면파 Weyl 응답과 정량적 연결식을 직접 제안 | 본문 정의의 CJ는 한 poset만 받는 `f_N`과 입력이 다름. broad Weyl 최초성은 충돌하지만 같은 추정량의 중복은 아님 |
| [Ghosh v4, 2026-08-14](https://zenodo.org/records/21931074), §§2–3, 6–6.4, 8 | 국소 count 평균의 비용 가정; chain 계열의 bulk 응답 채널 | 유한 표본 비용·조건부 채널 제약을 이미 연구 | 우리 공유-pair 통계에 가정을 확인해야 함. 측도 변화와 순서 변화의 정보 경로를 구분해야 하며, §4와 부록의 적용 조건 문제가 남음 |
| [Eichhorn–Mack–Le–Wagner, 2026 v1](https://arxiv.org/html/2605.27514v1), §§II.1.1, III, V | Hasse link degree·graph spectrum·interval abundance 등; 아홉 구성 계열 | 저비용 그래프 관측량으로 여러 기하/비다양체 계열 구별 | 연구한 manifoldlike 계열은 conformally flat로 Weyl=0. 진공 평면파 응답 계수는 이 논문의 검증 대상 밖 |
| [Ahmed–Bombelli, APS 2026](https://meetings-archive.aps.org/smt/2026/apr-r93/1/), 발표 초록 | plane-wave causal set에서 geodesic deviation을 이용한 Weyl 존재 탐색 | 같은 물리적 주제의 선행 발표가 존재 | 전문을 확보하지 못했으므로 입력·정량 성능의 비중복은 미확정. “최초 검출” 근거로는 사용할 수 없음 |
| [Hoeffding, 1948](https://doi.org/10.1214/aoms/1177730196), U-statistic 이론 | iid 점들로 만드는 대칭 kernel 평균 | 중첩 index를 고려한 분산 분석의 표준 틀 | 새로울 수 있는 것은 우리 kernel의 계수 계산과 적용 결과 |

보조 확인: [Braun, 2025](https://arxiv.org/abs/2507.01907)의 재구성 정리를
유한-N 스칼라 하나의 추정 보장으로 읽지 않는다.
[Barton–Borza–Röhrig, 2026](https://arxiv.org/abs/2606.04910)는 초록 범위에서
optimal transport를 통한 timelike Ricci 관측량 연구로 분류했다. 이 논문 전체의
Weyl 고차항을 감사한 것은 아니다. 이 두 자료는 직접 중복 판정의 근거로 삼지 않았다.

## 4. 겹친다는 판단과 겹치지 않는다는 판단 모두 좁혀야 한다

### Alfyorov–Shnyukov: 직접 비교하되 정보 입력을 맞춘다

확인한 PDF pp.6–7의 CJ는 같은 좌표 점에서 만든 두 order의 Hasse 경로 점수 차이,
proper time·radial position·Hasse depth의 strata를 사용한다. 우리의 주 판별기는
검사 대상 한 poset의 관계와 N을 받도록 계획됐다. 두 방법의 주장을 같은
정보 조건으로 놓을 수 없다. 해당 PDF의 연결식에는 미증명 가정이 명시돼 있고,
§3.1의 술어 설명에는 섭동 근사가 들어간다. 수치 baseline으로 옮기려면 먼저
버전과 술어 규약을 일치시켜야 한다. [확인한 원문](https://inspirehep.net/files/32388edfda552e0b4cfb6697bda66d7b)

그렇다고 그 연구를 무관하다고 분류할 수도 없다. 곡률 응답 계수와 유한 표본을
연결한다는 목표에서 가장 가까운 비교 대상이다. 우리의 입력 요구가 적다는 차이는
유용할 수 있으나, 목표도 검출로 더 좁다. 이것만으로 더 나은 추정기라고 하지 않는다.

### Ghosh: 조건부 비용과 전체 order의 불가능성을 구분한다

v4 §§2–3의 sufficient cost에는 H1–H3가 있다. 특히 H2의 비양의 inter-site
공분산은 우리 관계쌍에 자동 적용되지 않는다. §6.3의 정확한 bulk kernel 결과도
명시한 tensor 구조를 전제로 하며, §6.4는 ordering/light-cone 경계의 완전한
곡률 전개를 마치지 않는다. 이를 모든 Weyl 검출의 no-go로 읽지 않는다.
§6의 position-space 정보 상한을 우리 모형으로 옮기는 데에는 아래 부록의
매개변수 의존 mapping 문제가 있다. [v4 원문·기록](https://zenodo.org/records/21931074)

우리 통계에서 공유 사건 공분산은 독립인 X2, X3를 X1에 조건화하면
`Cov(h12,h13)=Var(E[h12|X1]) >= 0`이다. 이것은 서로 다른 관측량을 다루는
H2의 반박이 아니라, 같은 비용식을 그대로 쓸 수 없다는 이유다. 또한 한 매개변수
평면파의 flat/curved 판별은 일반 곡률 불변량 여러 개를 동시에 복원하는 문제보다 좁다.

### 대칭성과 유한-N 법칙: 도구와 발견을 분리한다

`Var(f_N) = 2p(1-p)/[N(N-1)] + 4(N-2)zeta/[N(N-1)]`은 표준 항등식이다.
Roy 등의 §6에도 점을 공유하는 chain의 요동 차수가 이미 나온다.
`zeta>0`일 때 표준편차가 N의 제곱근에 반비례한다는 것, 따라서 신호가
`c lambda^k`이면 검출 경계가 조건부로 `N^(-1/(2k))` 꼴이라는 것을
새 물리 법칙으로 제시하지 않는다. [Roy–Sinha–Surya 원문](https://arxiv.org/html/1212.0631)

미분 가능성, c의 비영 여부, zeta의 퇴화 여부, 유한 표본 꼬리가 실제로 이 모형에서
어떤지는 계산할 가치가 있다. 대칭성만으로 c의 값이나 근사 유효 범위가 나오지는 않는다.

## 5. 차별성 후보를 실제 결과로 만들 최소 요건

G1의 첫 작업은 새 runner 제작이 아니라 아래 수학적 산출물이다.

1. **명시적인 응답 계수.** R+/R0/R-에 대해 signed A와 축 교환의 관계를 쓰고,
   `p_R(A)-p_R(0)`의 첫 비영 항을 계산한다. 부호·크기·적용 범위를 제시한다.
   짝함수라는 관찰이나 P3-E 기울기 재인용만으로 완료하지 않는다.
2. **계수가 생기는 원인.** 인과 경계 이동, 표본 측도, 관측 영역 이동을 구분한다.
   고정 Brinkmann box는 측도가 같지만 order가 바뀐다. D는 영역·정규화까지 바뀐다.
   이것이 관측 영역의 선택 때문에 생기는 신호인지 물리적으로 해석해야 한다.
3. **잡음 계수와 검증 가능한 예측.** p와 zeta를 2점·3점 적분으로 계산할 수 있는지
   확인하고 수치 오차를 분리한다. 이미 관측한 AUC를 fit한 재설명을 넘어서,
   아직 보지 않은 세기·N에 대한 분포나 검출률을 예측할 수 있어야 한다.
4. **효용 또는 한계의 명확한 추가.** 관측 영역을 달리했을 때 필요한 N이 왜 변하는지,
   혹은 특정 영역에서 선도 민감도가 왜 사라지는지 설명해야 한다. 정확한 술어를
   사용한다거나 사전등록했다는 절차만으로 중요성을 대신하지 않는다.

현실적인 첫 범위는 고정된 box 세 개다. D는 경계 이동을 포함하는 별도 검산이며,
그 유도가 예산에 맞지 않으면 D까지 해결했다는 넓은 결론을 포기한다.
정확한 해석식이 불가능해도 오차가 관리된 수치 적분은 후보가 될 수 있으나,
그 경우 결과의 이름도 수치 방법/예측 검증으로 좁힌다.

“한 poset만 사용한다”는 문장에도 조건을 붙인다. 통계량 f 자체는 순서와 N으로
계산하지만, 관측자·box·diamond tips를 정해 표본을 모으는 절차는 연속체 규약이다.
그 규약을 알고 보정한 검출기이며, 임의의 인과집합에서 관측 영역까지 순서만으로
찾아내는 알고리즘을 제안하는 것은 아니다. D를 큰 poset에서 추출할 때 추가 labels를
사용한다면 입력 표에 포함한다.

| G1 첫 유도의 결과 | 후속 결정 |
|---|---|
| 짝함수·표준 분산식·기존 기울기 재현만 얻음 | 차별성 부족. 이 연구의 본 캠페인 중단 |
| 기존 논문에 같은 영역·관측량·응답 계수가 있음 | 우선권 주장 삭제. 실질적인 추가 질문이 없으면 중단 |
| 새 계수 또는 오차가 관리된 응답식과 구체적인 미관측 예측을 얻음 | G0를 다시 평가한 뒤 G1 검산·G2 가격 산정 후보로 진행 |
| 계수가 예상과 달리 0이거나 비해석적 항/경계 문제가 있음 | 예상 지수를 고집하지 않고 원인을 검토. 이것만으로 출판 가능한 음성 결과라 선언하지 않음 |
| box에는 성립하나 D에는 성립하지 않음 | 관측 프로토콜에 제한된 결론. 부적합한 D 결과를 제외하고 보편성 주장하지 않음 |

첫 이론 검토의 작업 한도는 **2–3일의 집중 작업을 제안**한다. 비용 약속이나
자동 실행 승인이 아니다. 종료 시 계수식·오차·남은 장애를 남기며, 대규모 표본 생성으로
유도의 빈자리를 메우지 않는다. 정확도 주장은 기존 95% 양측 CI 규칙을 유지하고,
확인 캠페인의 표본 수는 기존 90% joint power 기준으로 다시 산정한다.

## 6. 검색 범위와 검토 한계

2026-09-21에 arXiv, INSPIRE의 저자 원문, Zenodo v4 원문, APS 발표 기록을 대조했다.
주요 검색 묶음은 `causal set / plane wave / ordering fraction`,
`curvature / U-statistic / boundary`, `Weyl / pair / probability`,
최근 저자명·정확한 제목이었다. 알려진 비교 논문의 관련 본문으로 역추적했다.
비관련 검색 결과와 일반 설명 사이트는 판정 근거에서 제외했다.

이는 범위를 정한 1차 출처 감사이며, 전체 문헌의 비존재 증명이나 모든 논문의
완전한 검증은 아니다. Ahmed–Bombelli의 전문은 이 검색에서 확보하지 못했다.
Hoeffding은 표준 이론의 서지 대조이며 이번에 전체 원문을 재검증한 것은 아니다.
각 preprint의 동료심사 통과 여부도 확인했다고 하지 않는다.

두 근접 preprint의 정확한 확인본은 다음과 같다. PDF 페이지는 뷰어의 1-based 순서다.

| 확인본 | 버전·분량 | SHA-256 |
|---|---|---|
| [Alfyorov–Shnyukov PDF](https://inspirehep.net/files/32388edfda552e0b4cfb6697bda66d7b) | 표지 March 2026, 37쪽; 표지 DOI 10.5281/zenodo.19354508 | `d571888a4c082768d14a979f4be8cd4e57be6b41fc91984b4767616eb2069d09` |
| [Ghosh v4 PDF](https://zenodo.org/records/21931074/files/Sample-Complexity-of-Causet-Curvature_V4.pdf?download=1) | v4 기록 2026-08-14, 31쪽 | `e9cb5d56d56429d05313df4234222c2cb1f18f6446b1be43dd67472f98a51266` |

가져온 시각은 manifest에 UTC로 기록했다. 첫 PDF를 기존 priority 문서의 다른 Zenodo
record와 동일한 버전이라고 가정하지 않았다. 로컬 확인본·텍스트·중요 페이지 렌더·manifest는
`outputs/literature_audit_2026-09-21/`에 있으며 Git 추적 대상은 아니다.
후속 논문 재현 묶음에는 사용한 버전과 적법한 보존 방법을 별도로 고정해야 한다.

## 부록 A. 정보 상한의 적용 조건을 확인하는 두 점 예

이 부록은 우리의 수학적 적용 조건 점검이다. Ghosh 전체 연구를 반박하거나,
새로운 정보이론 정리를 주장하는 결과가 아니다. §6의 position-to-order
data-processing 논리를 이 평면파족에 그대로 적용할 수 있는지만 확인한다.

통상적인 data processing은 매개변수 A에 의존하지 않는 공통 channel K를 전제로 한다.
좌표 표본 X가 주어졌을 때 인과관계를 정하려면 계량도 필요하므로 여기서는
`C_A = F_A(X)`다. 모든 A에서 `X ~ P`가 같아도 F_A가 달라지면 C_A의 분포는
달라질 수 있다. 입력을 `(X,A)`로 확장하면 A 자체가 들어가므로, 그 정보량을
좌표 표본의 정보량만으로 대신할 수도 없다.

구체적으로 `g_A = 2du dv + A(x²-y²)du² + dx²+dy²`에서 A=0,1을 비교한다.
고정된 convex box를

```text
R = [0,1]_u × [0,1]_v × [-0.1,0.1]_x × [0.9,1.1]_y
```

로 잡는다. 두 계량 모두 `sqrt(-det g_A)=1`이므로 이 R에서 균일하게 뽑은 점의
좌표 분포는 같다. R 전체에서 `x²-y² <= -0.8`이다. 따라서 flat future-causal
두 점을 잇는 직선은 R 안에 머물고, 그 접벡터의 g1 norm은 g0 norm보다 작거나 같다.
flat의 인과쌍은 curved에서도 유지된다.

반면 다음 내부의 두 점을 보자. 좌표 순서는 `(u,v,x,y)`다.

```text
p = (0.25, 0.25, 0, 1)
q = (0.75, 0.30, 0, 1)
g0(q-p,q-p) = +0.05
g1(q-p,q-p) = -0.20   (이 직선에서는 x,y가 일정)
```

flat에서는 spacelike지만 curved에서는 future-timelike 직선이 존재한다.
엄격한 부등식은 주변의 열린 점쌍 집합에서도 유지된다. 따라서 균일한 독립 두 점이
만드는 **unlabelled 두 원소 poset**의 chain 확률은 A=1에서 엄격히 더 크다.
좌표 표본의 법칙은 동일하지만 order 법칙은 다르다. 그러므로 raw coordinates만을
입력으로 하는 매개변수 독립 channel을 두 모형에 공통으로 둘 수 없다.
유효한 data-processing 부등식의 위반이 아니라, 그 전제의 누락을 보여주는 예다.

이 예는 scalar curvature invariant 복원의 표본복잡도 지수를 결정하지 않는다.
알려진 parameter-independent causal map 안에서 측도만 바뀌는 문제의 하한도
부정하지 않는다. 또한 R은 외부에서 정의한 관측 영역이므로 관측 영역을 모르는
검출기로의 일반화는 없다.

기존 정확 술어와 독립적인 직선 norm 계산으로 한 쌍을 대조했다.
flat margin은 약 -0.05, curved margin은 약 +0.205341921221036이며,
기존 술어의 오차 구간은 둘 다 0을 제외한다. 실행 기록은
`outputs/literature_audit_2026-09-21/order_channel_witness.json`이다.
아래 검산은 저장소 루트에서 재실행할 수 있다. 통계 표본 생성은 없다.

```python
import math
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path("experiments/positive_control").resolve()))
from p14_plane_wave import Slab, PlaneWaveGeometry, causal_relation, transverse_cost

p = np.array([0.25, 0.25, 0.0, 1.0])
q = np.array([0.75, 0.30, 0.0, 1.0])
# 이 enclosing slab는 술어의 conjugate-point guard용이며 위 R의 대체가 아니다.
slab = Slab(1.0, 2.0, 2.0, 4.0)
assert causal_relation(PlaneWaveGeometry(slab, 0.0), p, q).related is False
assert causal_relation(PlaneWaveGeometry(slab, 1.0), p, q).related is True
du, dv = q[0] - p[0], q[1] - p[1]
assert 2 * du * dv > 0 > 2 * du * dv - du**2
assert 0.1**2 - 0.9**2 < 0
assert math.isclose(transverse_cost(du, 0, 0, 1, 1, 1), -math.tan(du/2), abs_tol=1e-14)
```
