# P14 O4 — G3 재설계 (G3a / G3b) 설계 문서

상태: **설계 확정, 크기 미동결.** 이 문서는 G3 를 무엇으로 대체할지와
그 결정을 위해 census 가 무엇을 세어야 하는지를 고정한다. 표본 수·CP
상계·적격률 규칙은 census 이후 사전등록 재개방에서 동결한다.

**결과 없음.** census 실행은 별도 승인 사안이다.

선행 기록: `p14_o4_incident.md` (abort), `p14_o4_replay_diagnostic.md`
(replay 진단 프로토콜).

---

## 1. 무엇이 고장났는가 — 하나의 결함군

동결 G3 는 **적응하지 않는 여유(non-adaptive margin)** 를 세 군데에서
썼다. 세 번째는 재설계를 준비하다 발견했다.

1. probe 를 경계에서 고정 `10⁻⁶` 만큼 떼어 놓았는데, 그 거리가 국소
   `err` 에 적응하지 않았다. `tol = 10⁻⁸` 대비 두 자릿수 여유를 둔 것이
   **`err` 대비 여유를 전혀 보장하지 않는다** — `err` 는 인증된 격납이
   아니라 Gauss–Legendre 정지 휴리스틱이기 때문이다.
2. 삼치 술어를 boolean 으로만 읽어 **정당한 `None`** 을 실패로 취급했다.
3. window 를 `T(θ)` 로 만들고 판정은 `t_min(dpsi)` 로 이루어지는데,
   `dpsi = acos(cos θ)` 의 각도 복원 오차가 **ulp 단위로 유계가 아니다**
   (§4.1). 몇 ulp 짜리 여유로는 덮이지 않는다.

재설계는 셋을 각각 닫는다. 1·3 은 G3b 에서, 2 는 G3a 에서.

---

## 2. G3a — wrapper 의미론 계약 (선행 전제조건)

### 2.1 지위와 실행 시점

G3a 는 **별도의 과학 stage 가 아니라 같은 stage 의 선행 전제조건**이다.
같은 `flight_time` 을 공유하는 wrapper 의미론·좌표 변환·삼치 분기의
계약 검사이지 독립적인 solver 검증이 아니다.

실행 순서를 이렇게 고정한다:

1. **exact-freeze checkout 에서 G3a 를 먼저 수행한다.**
2. 실패하면 **fresh 시드를 건드리기 전에** stage 를 blocked/invalid 로
   종료한다.
3. 통과한 경우에만 G1 / G2 를 실행한다.
4. G1 의 수락 점에서 G3b 적격 cluster 를 구성해 검사한다.

이것이 O4 abort 의 직접적 교훈이다. 동결 러너는 wrapper 계약을 12 시간
치 표본을 소비한 **뒤에** 처음 만졌다. G3a 를 preflight 로 올리면 같은
결함이 시드를 쓰기 전에 잡힌다.

### 2.2 검사 항목

각 검사점에서 `(t_min, err)` 를 **술어가 만들 `dpsi` 로** 먼저 구하고,
그로부터 `dt` 를 지어 술어에 넣는다. 답이 기하 직관이 아니라 **구성상**
정해진다.

| 행 | 구성 | 예상 | 근거 |
|---|---|---|---|
| A | `dt = t_min + err + η` | `True` | `\|dt−t_min\| > err`, `dt > t_min` |
| B | `dt = t_min − err − η`, **`dt ≥ 0` 고정 사례에서만** | `False` | `\|dt−t_min\| > err`, `dt < t_min` |
| C | `dt = t_min` | **`None`** | `\|dt−t_min\| = 0 ≤ err` — 정당한 undecided |
| D | `dt < 0` | `False` | 술어의 단락. 솔버를 호출하지 않는다 |

**행 C 가 통과라는 것이 재설계의 핵심이다.** 동결 G3 는 어떤 `None`
이든 abort 시켰다. 삼치는 설계이므로 계약은 "`None` 이 없다"가 아니라
"`None` 이 나와야 할 때만 나온다"여야 한다.

**행 B 와 행 D 는 겹치면 안 된다.** B 는 `dt ≥ 0` 이 보장되는 고정
사례에서만 실행한다 — `t_min − err − η < 0` 인 점에서 B 를 돌리면 술어가
단락으로 `False` 를 반환해 우연히 통과하고, 정작 검사하려던 **하단
error-band 바깥 판정 경로**는 한 번도 실행되지 않는다. 그러면 B 는 D 의
중복일 뿐이다. `t_min − err − η < 0` 인 점은 B 에서 **construction
-unavailable** 로 분류하고 그 개수를 보고한다.

### 2.3 η 의 적용

행 A·B 의 `dt` 도 §4 의 **실현 여유 검사**를 통과해야 한다. 즉 구성 후
술어와 동일한 방식으로 재계산해 `|dt − t_min| − err ≥ η` 를 확인한다.

---

## 3. G3b — 결정 가능성 계약

### 3.1 좌표 — 전부 술어 자신의 것

각 후보 공간점 `(r, θ)` 에서:

```
dpsi   = acos(clamp(sin p_θ sin q_θ cos(p_φ−q_φ) + cos p_θ cos q_θ))
t₁, e₁ = flight_time(R_IN,  r,     dpsi)
t₂, e₂ = flight_time(r,     R_OUT, dpsi)
lo = t₁                hi = Δt − t₂
L  = hi − lo
```

`t₁·t₂·e₁·e₂·L` 은 **모두 술어 좌표 `dpsi` 에서 재계산한 값**이며,
`_ell` 이 `θ` 로 계산한 `T₁·T₂·L_S1` 과 구별된다. 후자는 표본을
고르는 데만 쓰이고, 계약 대수에는 들어가지 않는다. 이유는 §4.1.

### 3.2 적격 규칙

```
W_robust = L − e₁ − e₂ − 2η
적격  ⟺  W_robust > 0        (strict)
```

### 3.3 probe 세 개 — 판정 가능함이 증명된다

`η > 0`, `W_robust > 0` 아래에서 (실현 여유 검사는 §4):

**(i) determined-inside** — robust window 의 중점
`t_x = ½((lo + e₁ + η) + (hi − e₂ − η))`

- `p→x`: `dt − t₁ = t_x − t₁ ≥ e₁ + η > e₁` → 판정, `True`
- `x→q`: `(Δt − t_x) − t₂ ≥ e₂ + η > e₂` → 판정, `True`
- ⟹ `a and b = True`

두 leg 의 여유의 합은 `L − e₁ − e₂ = W_robust + 2η` 이고, 중점에서 각
leg 가 그 절반 `≥ η` 를 갖는다.

**(ii) determined-outside-above** — `t_x = hi + e₂ + η`

- `x→q`: `dt = Δt − t_x = t₂ − e₂ − η`, `\|dt − t₂\| = e₂ + η > e₂`
  → 판정, `dt < t₂` 이므로 `False`
- `p→x`: `\|dt − t₁\| = L + e₂ + η`. 적격에서 `L > e₁ + e₂ + 2η` 이므로
  `L + e₂ + η > e₁ + 2e₂ + 3η > e₁` → 판정, `True`
- ⟹ `a and b = False`

**(iii) determined-outside-below** — `t_x = lo − e₁ − η`

- `p→x`: `\|dt − t₁\| = e₁ + η > e₁` → 판정, `dt < t₁` 이므로 `False`
- `x→q`: `dt = Δt − t_x = Δt − lo + e₁ + η`, `\|dt − t₂\| = L + e₁ + η`
  → (ii) 와 같은 논거로 `> e₂` → 판정, `True`
- ⟹ `a and b = False`

동결 G3 는 (ii) 만 있었고 그나마 offset 이 `e₂` 에 적응하지 않았다.
**(iii) 은 아래쪽 경계를 처음으로 검사한다** — 상단과 다른 leg 를 통해
다른 시간 산술 경로를 지나므로, census 를 보고 포함 여부를 정할 이유가
없다. 지금 포함한다.

### 3.4 (iii) 에서 `t_x < 0` 인 경우 — 세 결과를 분리한다

`t_x = lo − e₁ − η < 0` 이면 술어는 `dt < 0` 단락으로 `False` 를 반환해
**우연히 원하는 답**을 낸다. 그것은 하단 경계 solver 검사가 아니다.
따라서 다음을 각각 다른 결과로 기록한다:

| 결과 | 조건 | 뜻 |
|---|---|---|
| `lower-boundary-solver-probe` | `t_x ≥ 0`, 실현 여유 통과 | (iii) 이 의도한 검사 |
| `negative-dt-short-circuit` | `t_x < 0` | 단락 경로. G3a 행 D 가 이미 검사한다 |
| `construction-unavailable` | 실현 여유 검사 실패 | §4.3 |

`negative-dt-short-circuit` 은 **G3b 의 통과 증거로 세지 않는다.**

### 3.5 무적격의 처분

`W_robust ≤ 0` 은 mismatch 도 실패도 아닌 **availability outcome** 이다.

적격 점만으로 rate 를 만들면 분포가 `G1 measure | L_S1>0 ∧ W_robust>0`
로 바뀐다 — 자기선택된 부분집합이다. 따라서 G3b 의 rate 는 그 조건부
법칙 위의 진술로 명시하고, **availability rate 를 함께 게시한다.**

**census 에는 gate 를 붙이지 않는다** (임계 근거가 아직 없다). 다만 새
O4 freeze 전에 다음 중 **하나를 반드시 사전 고정**한다:

- 고정 G1 표본에서의 **최소 적격률**, 또는
- 고정된 적격 cluster 수를 얻을 때까지 스캔하는 **규칙**(스캔 상한
  포함), 또는
- 적격 부족 시 `INVALID` / `INCONCLUSIVE` 로 종료하는 **규칙**.

**이전 stress set 의 availability 는 설계 입력일 뿐이며, 새 캠페인의
통과 증거로 재사용하지 않는다.**

---

## 4. η — 실현된 최소 판정 여유

### 4.1 반드시 없애야 할 항 — 각도 복원 오차

술어는 각도를 `dpsi = acos(cos θ)` 로 되찾는다. 이 복원의 오차는 **ulp
단위로 유계가 아니고 `1/θ` 로 커진다**(실측, 예측 상계
`½·ulp(1)/sin θ` 와 일치):

| θ | `\|acos(cos θ) − θ\|` | θ 의 ulp 배수 |
|---|---|---|
| 10⁻⁵ | 4.14e-13 | 2.4×10⁸ |
| 10⁻⁴ | 2.62e-13 | 1.9×10⁷ |
| 10⁻³ | 7.83e-15 | 3.6×10⁴ |
| 10⁻² | 1.44e-15 | 8.3×10² |
| 0.0626 | 4.16e-17 | 3 |

샘플러는 `cos θ` 를 균등하게 뽑으므로 작은 θ 가 드물지만 배제되지
않는다. `T(θ)` 로 window 를 만들고 `t_min(dpsi)` 로 판정하면 그 불일치가
최대 ~10⁻¹² 급이며, 이는 §1 의 결함 1 과 **같은 계열**이다.

**설계 결정: §3.1 대로 window·probe·적격을 전부 술어 좌표에서 만든다.**
그러면 이 항이 항등적으로 0 이 된다. 그 크기는 census 가 진단으로
측정해 기록에 남긴다(§5.7) — 설계를 census 로 정당화하는 것이 아니라,
설계가 없앤 항의 규모를 남기는 것이다.

### 4.2 η 의 정의 — addend 가 아니라 검사되는 하한

**η = 10⁻¹² (frozen).**

η 는 단순한 addend 가 아니라 **실현된 최소 판정 여유**로 정의한다.
후보 `t_x` 를 만든 뒤, **wrapper 와 동일한 방식으로** `dt`·`dpsi`·
`t_min`·`err` 를 재계산하여

```
realized_margin = |dt − t_min| − err  ≥  η
```

인지 **실제로 검사한다.** 이렇게 하면 "연산 3~4회의 반올림이 약
5×10⁻¹⁵ 이므로 η = 10⁻¹² 면 충분하다"는 **휴리스틱에 인증을 맡기지
않아도 된다** — 여유는 가정되지 않고 확인된다.

η = 10⁻¹² 를 고른 근거는 두 가지뿐이며, 둘 다 census 와 무관하다:

- `ulp(Δt) = 1.78×10⁻¹⁵` 규모의 반올림보다 충분히 크다(약 560배).
- 관측된 `err` 규모(`e₂ = 2.19×10⁻⁶`)보다 **6자릿수 아래**라 적격을
  좌우하지 않는다. η 는 부동소수점 분리라는 제 일만 하고 과학에
  개입하지 않는다.

**η 는 `err₂ ≥ 10⁻⁶` 의 빈도로 정하지 않는다.** census 는 η 별
적격률·비용을 보여주는 **탐색 입력**이지 η 의 정당성 근거가 아니다.

### 4.3 부족할 때 — nextafter 조정, 아니면 construction-unavailable

`realized_margin < η` 이면:

- **(ii)·(iii)** 단측 probe: `t_x` 를 바깥 방향으로 `math.nextafter` 로
  한 ulp 씩 밀며 재검사한다. 바깥으로 미는 것은 해당 leg 의 여유를
  단조 증가시킨다.
- **(i)** 양측 probe: 두 leg 의 여유가 상충하므로, 부족한 쪽으로 밀며
  **두 leg 를 모두** 재검사한다. 여유의 합이 `W_robust + 2η` 로 고정돼
  있으므로 `W_robust` 가 반올림보다 클 때만 성공한다.

**`MAX_NUDGES = 64` (frozen).** 이 값도 census 이전에 정해야 한다 —
cap 이 같은 cluster 를 유효 probe 로 만들지 `construction-unavailable`
로 만들지 가르므로, 나중에 정하면 **census 가 어떤 cluster 를 셀지
고르게 된다.** η 와 같은 이유로 여기서 동결한다.

근거는 반올림뿐이다. 후보는 명목상 여유가 정확히 η 인 자리에 놓이고,
그것을 만드는 산술(`t_x` 2~3연산, `dt` 1연산, 차 1연산)이 각각 최대
`½·ulp(8.5) = 8.9×10⁻¹⁶` 만큼 반올림하므로 실현 여유는 최대 약
`3.6×10⁻¹⁵` 부족할 수 있다. 한 걸음은 `dt` 를 probe 시각의 1 ulp 만큼
움직이며 최악이 `ulp(1.0) = 2.2×10⁻¹⁶` 이므로 약 17 걸음이면 충분하고,
64 는 약 3.7 배의 여유다.

**세는 규칙.** 걸음 수는 `(cluster, probe)` 단위로 센다. 64 걸음 안에
검사를 통과하지 못하면 그 점의 그 probe 를
**`construction-unavailable`** 로 분류한다. 이는 mismatch 가 아니며
**mismatch 와 분리해서** 집계한다. 실제로 쓴 걸음 수도 함께 기록한다 —
routinely 여러 걸음이 필요한 실행은 위 반올림 논거가 틀렸다는 신호이며,
그것이 숨겨져서는 안 된다.

### 4.4 전제 — `flight_time` 의 결정론

같은 인자에 대해 비트 동일한 결과를 반환해야 한다. 계약 테스트로 못
박는다.

### 4.5 `err` 를 써도 되는가

G3b 는 `err` 로 probe 를 **배치**한다. 이는 정당하다 — G3a/G3b 는
instrumentation 계약이지 추정량이 아니다. 금지된 것은 `err` 로 `V_S1`
의 구간을 만드는 것이며(O4 freeze review P1-5), 그 금지는 그대로
유지된다. 두 용법의 구분은 여기서 명시적이다.

---

## 5. Census 충분통계 스키마 (전체 실행 전 동결)

census 는 **replay** 다. 동결 stress 집합을 그대로 재생하며, 새 probe 를
실행하지 않는다. 아래 양들은 이미 재유도하는 `flight_time` 호출에서
파생되므로 **추가 solver 비용이 없다.**

모든 비율의 분모는 `clusters_probed` 이며 함께 게시한다.

### 5.1 `(probe, leg, outcome)` 완전 분리

`probe ∈ {midpoint, outside}` × `leg ∈ {p→x, x→q}` ×
`outcome ∈ {true, false, undecided}` 의 12 칸을 각각 센다.
**어떤 축으로도 합산하지 않는다.** 이것이 현행 스키마의 결함
(`outside-in-error-band` 가 두 leg 를 합산)을 직접 고친다.

### 5.2 `outside x→q` 의 `err₂ ≥ 10⁻⁶` 개수와 교차검증

- `nominal`: `e₂ ≥ 10⁻⁶` (동결 offset 과의 비교, **`≥`**)
- `realized`: `outside · x→q` 의 undecided 수 (= `e₂ ≥ |dt − t_min|`)
- `disagreements`: 둘이 갈리는 cluster 수와 **첫 좌표**

두 조건은 `|dt − t_min|` 이 `10⁻⁶` 와 정확히 같지 않기 때문에 원리적으로
갈릴 수 있다(반올림 규모). 이 차이 자체가 §1 결함군의 또 다른 표식이므로
**같다고 단언하지 않고 세어서 보고한다.** 계약 테스트는 세 필드가 모두
존재하고 `disagreements` 가 좌표와 함께 기록됨을 요구한다.

### 5.3 사전 동결 요약 — `e₁`, `e₂`, `L`, `L − e₁ − e₂`

모두 **술어 좌표 `dpsi` 에서 재계산한 값**이다.

- **histogram.** bin edge 를 실행 전에 이 문서와 코드에 동시에 못
  박는다. `e₁`·`e₂` 는 log 격자(십진 decade 를 3 등분, `10⁻¹⁶`–`10⁻²`),
  `L`·`L − e₁ − e₂` 는 선형 격자(`0`–`L_max_ub`, 32 등분). 각각
  underflow / overflow bin 을 둔다.
  **bin 경계는 `[lo, hi)` 로 왼쪽 닫힘·오른쪽 열림**이며, 마지막 bin 만
  `[lo, hi]`.
- **min / max.**
- **분위수.** 확률 격자 `{0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999}`,
  **interpolation 은 `numpy.quantile(..., method="linear")`** 로
  고정한다.

평균·분산은 내지 않는다 — 설계에 필요한 것은 꼬리이고, 아티팩트의
result-shaped 키 금지와도 충돌한다.

### 5.4 η 후보 격자별 `W_robust > 0` 개수

격자를 **정확히** 이 목록으로 고정한다:

```
η ∈ {0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6}
```

각 η 에 대해 적격 개수와 적격률. **비교는 strict:
`W_robust > 0`**, 즉 `L − e₁ − e₂ − 2η > 0`.

η = 10⁻¹² 가 η = 0 과 사실상 같은 적격률을 준다는 것이 §4.2 둘째 근거의
실증이 된다(정당성의 근거가 아니라 실증이다).

### 5.5 boolean mismatch — 분모는 fully-decided probe 만

`bool(a and b) ≠ want` 를 **두 leg 모두 판정된 probe** 에 한해 세고, 그
분모(`fully_decided_probes`)를 probe 종류별로 함께 낸다.

**`construction-unavailable` 은 mismatch 와 분리해서 집계한다**(§4.3).

다만 census 에는 `construction-unavailable` 의 직접 대응물이 없다 —
census 는 replay 이므로 **자기 probe 를 구성하지 않는다.** 가장 가까운
census 양은 §5.8 의 하단 probe 도달 실패 개수이며, 아티팩트가 이 사실을
문장으로 싣는다. 실현 여유 검사의 실패는 G3b 가 실제로 probe 를 구성할
때 발생하며, 그때 mismatch 와 분리해 집계한다.

### 5.6 범주별 첫 발생 좌표

각 범주의 첫 좌표를 전 좌표(§`p14_o4_replay_diagnostic.md` §4)와 함께
기록한다. 분모도 함께.

### 5.7 (진단 전용) `T(θ)` 대 `t_min(dpsi)` 의 leg 별 차이

§4.1 이 제거한 항의 규모. histogram + max + 분위수. **gate 가 아니다.**

### 5.8 (진단 전용) 하단 probe 도달 실패 — 적격과 **결합해서** 센다

probe (iii) 의 `lo − e₁ − η ≥ 0` 위반 빈도. probe coverage 정보이며
**gate 가 아니다.** §3.4 의 `negative-dt-short-circuit` 이 얼마나 흔할지를
미리 알려준다.

**주변합만으로는 답이 나오지 않는다.** probe (iii) 은 `W_robust > 0` 인
**적격 cluster 에서만** 실행되므로, coverage 질문은 "적격 집합 안에서
얼마나 도달하지 못하는가" 이다. `eligible_clusters` 와 전체
`lower_probe_t_x_negative` 두 주변합만 남기면, 도달 실패 cluster 대부분이
애초에 무적격인 경우에 coverage 를 크게 **과장**하고, 총수와 함께 두어도
결합분포가 복원되지 않는다.

따라서 η 별로 교집합
**`eligible_and_lower_probe_t_x_negative`** 를 함께 센다. 이 한 칸에
두 주변합과 분모를 더하면 2×2 가 결정된다.

---

## 6. census 이후에야 동결할 수 있는 것

- G3b 의 cluster 수와 CP 상계 (적격률이 유효 표본을 바꾼다)
- §3.5 의 세 규칙 중 어느 것을 쓸지와 그 임계값
- G3a 의 검사점 수와 그 표본 추출 방식
- 총 비용과 상한

**η 와 `MAX_NUDGES` 는 여기 없다.** 둘 다 census 없이 §4 에서
반올림 논거만으로 정해졌고, census 는 그 값을 바꾸지 않는다. 두 값이
사후 동결 대상이 아님을 명시해 두는 이유는, 그렇지 않으면 census 를 본
뒤에 고를 여지가 남기 때문이다.

---

## 7. 순서

1. ~~G3a/G3b 설계 메모와 census 스키마~~ — 이 문서
2. **진단 스키마 보강 PR (결과 없음) 및 리뷰** — 현재 단계
3. clean checkout 에서 100,000-cluster census — **별도 실행 승인**
4. census 를 탐색 입력으로 §6 항목 동결 → 사전등록 재개방
5. 새 freeze·새 G1/G2 시드로 O4 전면 재실행 — **실행 승인 또 별도**

O4 는 여전히 **무판정(verdict = null)** 이며, 이 문서의 어떤 내용도 그
상태를 바꾸지 않는다.
