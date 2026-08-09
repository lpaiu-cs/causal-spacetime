# P14 사전등록 — v0.4 (준비 커밋; 최종 freeze 아님)

Status: **STAGE COMPLETE — POSITIVE.** 설계 v0.4 승인(in-session
4라운드 + R-5) → 준비 커밋(PR #49) → ledger 단일화 P′=`51875a2` →
clean checkout P′에서 preflight 인증 → 최종 freeze F=`b858b08`
(PR #51, S1 가격 첨부) → 실행 승인 후 **정확한 F checkout**에서
캠페인 단일 실행 → 이 결과 커밋 R. 판정: **C1 confirmed ∧ C2
confirmed** (아래 결과 절).

Runner: `experiments/positive_control/p14_prereg.py`. 이 문서의
상수는 러너의 동결 상수와 테스트로 동기화된다.

## §0. 입력 계보

p14_weyl_curvature.md v0.8의 상태줄("NO STOCHASTIC PROBE ... HAS
RUN")을 이어 쓰지 않는다. 시작점:

| 입력 | 출처 | 값 |
|---|---|---|
| P14 체인 완결 | **PR #38–#48**, 최종 merge **`1025c50`** | foundation(#38)·P0·P1·P2·P3-E·P3-C·P4 |
| 운영점 선정 | P3-E (`p14_probe_p3e_results.json`) | aniso-a1.0 primary |
| 분리 확인 (unpaired) | P3-C confirmed 블록 (시드 20260851/52) | s=11.3489, AUC 하한 0.999232, BA=1.0 |
| C1 효과 (paired) | **P3-E 원시 paired 표본** (aniso-a1.0, 250쌍) | paired pilot 평균 0.05018773 |
| Paired variance | P4 + 직접 pair-bootstrap | var(Δ3a)=2.40086e-5, 95% CI [2.01316e-5, 2.79232e-5], 상한 SD 0.00528424 |
| 검정 기계 | P3-C 동결분 | 3분기, DeLong+disjoint-pair 경계, 마진 3종, preflight 프로토콜 |

House-rule 선행 출처(체인 외부): 구간 규칙(§6, AGENTS.md)과 P12의
seed-window helper — PR #37은 P12 Stage B이며 P14 계보가 아니다.

## §1. 동결 주장과 운영점

**운영점:** aniso-a1.0 — w=1.0, slab (Δu,Δv,Δx,Δy)=(1.0,1.0,2.0,6.0),
E[N]=300, a=1.0<π, profile `A(u)(x²−y²)` 상수-A slab. 러너는 P3-C의
`POINT`·`E_N`을 identity로 계승한다.

**C1 (paired ensemble 평균):** "aniso-a1.0에서 paired ensemble의 평균
global relation-fraction 차이 θ_Δ = E[f_A − f_0]의 95% CI가 +ε_Δ 위에
있다." 검정량은 관계 비율 **순변화**의 앙상블 평균이다 — 개별
점집합의 이동량도, gained+lost 총이동도 아니다.

**C2 (classifier 앙상블 판별, 독립 preregistered replication):** "단일
poset의 global relation fraction만을 입력받는 동결 classifier가
curved와 flat ensemble을 동결된 s/AUC/BA 규칙으로 분리한다." P3-C
확인의 독립 replication이며, 실패해도 P3-C 기록을 소급 취소하지 않는다.

**비주장(계승):** Petrov type N 한정, Weyl 정량 복원 아님, null은
밀도·계기 스코프 진술.

## §2. 통계량과 paired/unpaired 역할

| | C1 | C2 |
|---|---|---|
| 설계 | paired — 같은 sprinkling을 A=1.0/A=0으로 두 번 읽음 (§4.1) | unpaired — 독립 스트림 두 팔 |
| 통계량 | Δ̄ (Δ3a_i = f_A,i − f_0,i 평균), Student-t 95% CI | s, AUC, BA (P3-C 동결 함수 재사용) |
| 복제 단위 | sprinkling | sprinkling |

## §3. 판정 문법

**C1 마진(동결):** ε_Δ = 0.08061025 × 0.0044395871 = **0.0003578762**
— P3-C confirmed 블록에 anchored된 **operational margin**(물리적 마진
아님). Stage SD 재표준화 대안 기각(R-1).

**Stage 판정:**
- **POSITIVE = C1 confirmed ∧ C2 confirmed** — conjunction이므로
  power도 joint event로 인증(§4).
- **Joint equivalence 판정은 정의하지 않는다** — 현재 n에서 C1
  equivalence 0.961 × C2 equivalence 0.91045 ≈ 0.875 < 0.90.
- 비양성 결과는 주장별 3분기 판정을 **각자의 동결 문장**으로 병기,
  어떤 조합에서도 합성하지 않는다:
  - C1 confirmed: "paired ensemble 평균 이동이 ε_Δ를 넘는다" /
    equivalent: "이동은 ε_Δ 안이다" / inconclusive: 그대로 보고.
  - C2 confirmed: "P3-C 분리를 독립적으로 재현했다." /
    inconclusive: "이번 preregistered replication은 재확인하지
    못했다; 기존 P3-C 확인 기록과 병기한다." / equivalent: "이번
    replication은 equivalence를 지지해 기존 P3-C와 충돌한다."
    (conflicting replication — 소급 취소 없음)

## §4. Power — 90%, exact CP 하한 인증

보수적 입력: sd_Δ는 paired variance **자체의** pair-bootstrap 95%
상한 **0.00528424** (gain-하한 나눗셈은 분자도 추정량이라 폐기);
C1 effect bootstrap 출처는 **P3-E paired 표본으로 고정**.

**표본수 동결:** n_C1 = **3000** (equivalence branch 구속: 필요량
≈2834), n_C2 = **4800/arm** (P3-C 인증 크기 identity). 상향 조항
없음 — 인증 실패 시 stage **blocked/review-reopened**.

**반복수 동결:**

| branch | B | root entropy | 구조 |
|---|---|---|---|
| joint effect | 4,000 | [781, 60] | replicate당 spawn(3) → (c1, c2_curved, c2_flat); C1 pair-bootstrap n=3000 + C2 독립 재표집 n=4800/arm, 두 파이프라인 모두 실행, `C1 confirmed AND C2 confirmed` 빈도 |
| C1 centered-null | 20,000 | [781, 61] | replicate당 spawn(1); 중심화 Δ 표본 bootstrap, equivalence 빈도 |
| C2 null | 20,000 | [781, 62] | replicate당 spawn(2); flat 표본 양팔 독립 재표집, equivalence 빈도 |

각 branch의 pass 빈도의 exact CP 95% 하한 ≥ 0.90. 어느 인증이
실패해도 **B나 n을 자동으로 올리지 않는다** — 설계 리뷰 재개만이
상수를 움직일 수 있고, 움직이면 새 동결 커밋이다.

독립 audit 참고치(n=3000, 고정 절대마진): effect 20000/20000,
centered-equivalence 19220/20000 = 0.961. 본 인증은 위 자체 고정
시드·B로 재실행된다.

## §5. 시드 — 정확 고정, 이중 검사

**Bootstrap 스트림 의미(동결, v0.4):** `[781, k]`는 **entropy
벡터**다(`np.random.default_rng([781, k])` 관례 계승) — spawn key가
아니며, 두 구성은 다른 스트림이다. 구조:

```python
root = np.random.SeedSequence([781, 60])   # entropy vector
rep_nodes = root.spawn(B)                  # replicate당 1 노드
c1_ss, c2_curved_ss, c2_flat_ss = rep_node.spawn(3)   # joint
```

충돌 없음의 증명은 전체 구조적 key 유일성: `(root_entropy,
spawn_key)` 전수 열거로 joint 12,000개 / C1 null 20,000개 / C2 null
40,000개 전부 유일 + 신규 root [781,60/61/62]가 기존
[781,41]/[781,50]/[781,51]과 상이함을 테스트가 단언. 초기 state
pin은 재현성 테스트로 병행(충돌 증명의 대체 아님).

**실행 시드(정확히 3개, 예비 없음):** c1_paired **40000061**,
c2_curved **40000071**, c2_flat **40000072**.

> **R-5 (closed).** v0.3이 승인한 20260861/71/72는 기계적 검사에서
> **P13 campaign v3 ledger range [15000000..21503999] 내부**로
> 판정되어 사용 불가. 같은 이유로 과거 P14 날짜형 시드(2026xxxx)는
> P13 v3의 장부상 예약·소비 envelope 안에 있어 **house freshness를
> 인증할 수 없었다** — 다만 P13 v3의 가능한 실제 RNG 시드 58,880개
> 전수 재구성과 날짜형 시드의 교집합은 0이므로, 정확히 같은 정수
> RNG 스트림이 재사용됐다는 증거는 없으며 **기존 probe 결과를 소급
> 강등하지 않는다** (거부 테스트가 증명하는 것은 ledger overlap이지
> 스트림 재사용이 아님을 테스트 자체에 명시). Stage 시드는 P12의
> **할당 decade [30000000..39999999] 전체**를 예약 범위로 취급한
> 뒤(사용 상한 33263999가 아니라 문서화된 할당 경계; P12는 이후
> 공간을 40000000+로 선언) 그 위의 40000061/71/72로 확정(끝자리
> 유지). 회귀 테스트: 20260861 거부(P13 v3), 34000061 거부(P12
> 할당 decade), 40000061/71/72 수용.

검사는 두 종류 모두 기계적으로: ① scalar 실행 시드 — P14 소각
집합(∪ 20260851/52) 및 house 정수 ledger(P11–P13 spent ranges +
P12 할당 decade)와 서로소(`seed_windows.assert_point_seeds_fresh`;
P12의 helper는 공용 순수 함수로 추출되었고 P12 wrapper 유지);
② SeedSequence 경로 — 위 구조적 유일성 열거.

**중단 복구:** 같은 시드로 처음부터 재실행. Ambiguity 단언 실패 시
시드 교체 금지 — stage 중단, 조사, 정정 기록 후에만 재개(정정이
시드를 소각하면 새 동결 커밋이 새 시드를 핀).

## §6. 이탈 정책

마진·방향·통계량·운영점·n_C1·n_C2·B 3종·시드 전부 동결 후 불변.
인증 실패·ambiguity 위반·기타 이탈 필요 시 stage는
blocked/review-reopened — 자동 조정 없음. 완전분리 AUC는 동결
disjoint-pair 경계 규칙; 사후 규칙 변경은 블록 강등 사유.

## §7. S1 첨부 (최종 freeze 선행조건)

"S1은 해당 solver·patch·N에서 인과 술어 구성요소의 비용을 기록한다.
감당 가능하더라도 Schwarzschild 부피·캠페인 경로는 별도
미해결이며(설계 문서 §8.1: 다이아몬드 부피 oracle 부재), 감당
불가능하면 해당 solver·도메인·예산 경로만 보류한다." — S1은 일반
Weyl 경로를 열지도 닫지도 못하며, scope 문장은 비용 기록 이상을
싣지 않는다. S1 문장은 freeze manifest에 기록된다.

## §8. 실행 순서와 code_version

1. **준비 커밋 `P`** — 러너·테스트·helper 추출·이 문서. 아티팩트
   없음(아티팩트 기대 테스트는 명시적 skip). **최종 freeze 아님.**
2. **S1 병행 완료** (§7 문장의 값 확보).
3. clean checkout `P`에서 **preflight 실행** → preflight 아티팩트에
   **`code_version = P`** 기록.
4. **최종 freeze 커밋 `F`** — 인증 아티팩트 + freeze manifest
   (preflight digest, `preflight_code_version = P`, S1 문장, 동결
   n·시드). "freeze"라는 단어는 여기서만 쓴다.
5. clean checkout `F`에서 **캠페인 실행** (별도 실행 승인 후) →
   결과 아티팩트에 **`code_version = F`** 기록. 캠페인 모드는 3중
   기계 게이트를 통과해야만 실행된다: manifest 존재, manifest의
   preflight digest 일치, 그리고 **preflight 아티팩트에 기록된
   실행 관련 소스 digest들(`_FROZEN_SOURCES`: 러너·P3-C/P3-E/P2/P1
   파이프라인·geometry·seed helper)이 현재 파일과 일치**. 마지막
   게이트가 P→F 사이 코드 변경을 차단한다 — F는 아티팩트 추가만
   가능하며, 코드가 바뀌면 preflight부터 다시다.
6. **후속 결과 커밋 `R`.**

**Ancestry 계약:** 테스트가 `P ≺ F ≺ R`을 git으로 단언
(`git merge-base --is-ancestor`). 상수 동일성 테스트는 값의 일치만
보이고 실행 커밋을 증명하지 못하므로 병행이지 대체가 아니다. 모든
확률적 모드는 dirty worktree에서 실행을 거부한다.

## 결과 (커밋 R; `code_version = F = b858b08`)

실행 전 점검: HEAD == F 정확 일치, `git status --porcelain` 공백,
결과 아티팩트 부재 — 셋 다 확인 후 동결 `campaign` 모드 1회 실행
(2026-08-09). 아티팩트: `p14_prereg_results.json`.

| 주장 | 판정 | 수치 (95% CI) |
|---|---|---|
| C1 (paired, n=3000, 시드 40000061) | **confirmed** | Δ̄ = 0.05029290 [0.05010457, 0.05048123]; 하한이 +ε_Δ = 0.0003579의 140배 |
| C2 (replication, n=4800/arm, 시드 40000071/72) | **confirmed** | s = 11.1985 [11.0351, 11.3618]; AUC = 1.0 [0.999232, 1.0]; BA = 1.0 [0.9859, 1.0141] |

**Stage: POSITIVE** (C1 confirmed ∧ C2 confirmed — joint 사전
인증분). 동결 문장 그대로: C1 "paired ensemble 평균 이동이 ε_Δ를
넘는다"; C2 **"P3-C 분리를 독립적으로 재현했다."** 네 relation
census(C1 양읽기 + C2 양팔) 모두 ambiguous = escalated = 0,
아티팩트에 기록·단언. C2는 완전분리를 재현했고 AUC 구간은 동결
disjoint-pair 경계 규칙이다. 검증: 지표·판정은 저장 원시 표본에서
재계산되고(테스트), 시드 스트림은 prefix로 재현되며(slow), git
ancestry `P′ ≺ F ≺ R`가 테스트로 단언된다.

라이선스 범위(§1 그대로): aniso-a1.0의 고정 slab·밀도·profile —
paired 앙상블 평균 이동의 존재(C1)와 단일-poset 3a classifier의
앙상블 판별(C2). 더 넓은 주장 없음; 일반화 경로의 가격은 S1
문장(freeze manifest)이 운반한다.

## 리뷰 이력

- v0.1 미승인: claim estimand 초과, joint power 부재, S1 범위 과대,
  분산 상한 유도 오류, 계보 오류.
- v0.2 부분: B 미동결, C2 부분문장 소급 약화, freeze 순서 충돌.
- v0.3 승인 + 2건: bootstrap 시드 의미 미동결, preflight
  code_version 부재.
- v0.4: 위 2건 델타 반영. R-5 발견(승인 시드의 P13 v3 ledger
  overlap) → 조건부 승인으로 40000061/71/72 확정, P12 할당 decade
  전체를 예약 취급, "스트림 재사용" 표현은 "ledger freshness 인증
  불가"로 정정(소급 강등 없음).
