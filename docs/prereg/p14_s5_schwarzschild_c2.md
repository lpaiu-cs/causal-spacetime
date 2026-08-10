# S5 사전등록 — Schwarzschild C2-unpaired 판별

상태: **freeze 대상** (인세션 설계 리뷰 3라운드로 확정; 결과는 이후
커밋에만 등장한다). 선행: S3 탐색·S4 C1 확증, C2 feasibility 감사
(`p14_c2_feasibility_audit.json`, 보수 S3 앵커).

## 1. 목적·허용 문장·주장 클래스

S4의 type-D 결과를 C1 counterfactual에서 **단일-poset 판별(C2-급)**로
승격한다. 허용 문장(양성 시)은 §8의 DETECTED 동결 문장이며, 프로그램
내부 진술로서 동결 좌표·도메인에 결부된다.

## 2. 마진 선언 (선언 순서가 규율이다)

**탐색 결과를 설계 입력으로 사용하되, AUC 0.60(ε_AUC = 0.10)을
예산·선택된 n·향후 확인 데이터와 독립적으로 동결한 "최소한 실질적으로
의미 있는 single-poset 판별력"으로 정의한다** — anchor-independent가
아니라 **confirmation-data-independent**다. 이 아래의 판별력은 §1의
허용 문장을 유용한 의미로 지지하지 못한다. 표본수는 이 선언 이후
power가 결정했다(§6). ε_BA = 0.10 (BA 문턱 0.60 = 0.5 + ε_BA).

## 3. 동결 프로토콜

- 두 **독립** 팔: 곡면팔(신선 Schwarzschild 스프링클, S1 술어 tol
  1e-8→1e-10, undecided-never-silent-False), 평탄팔(신선 평탄 스프링클,
  정확 chord 술어). 공통 측도 샘플러(ρ·r²sinθ), reading당
  N ~ Poisson(300), S1 동결 도메인(M=1, r∈[10,20], 캡 1.0, T=40).
- 분류 점수 = per-reading global relation fraction f. **방향 동결**:
  곡면이 평탄 아래(S3·S4 실측 부호), AUC = P(f_curved < f_flat) + ½동률.
- 곡면팔 identified 처리: undecided가 만든 (lower, upper) 두 계열
  **모두에서 같은 분기가 나올 때만** 그 분기를 채택하고, 불일치는
  INCONCLUSIVE(**identified-agreement 규칙**). ambiguity 0이면 두 계열이
  일치해 통상 규칙으로 환원된다.

## 4. Primary gate — AUC, DeLong, 4분기 (엄격 부등호)

| 조건 | 판정 |
| --- | --- |
| AUC CI95 하한 > 0.60 | **DETECTED** |
| CI95 전체 ⊂ (0.40, 0.60) | **EQUIVALENT-AT-MARGIN** (판별력 부재의 단언이 아님 — §8 문장에 명시) |
| CI95 상단 < 0.40 | **DIRECTION-REVERSED** (역방향 강결과의 INCONCLUSIVE 은닉 금지) |
| 그 외 | **INCONCLUSIVE** |

각 분기의 동결 문장은 §8에 있으며, 러너 `SENTENCES`가 정본이다.

**DeLong 동결 명세**: 동률은 midrank(ψ = 1(c<f) + ½·1(c=f)), 분산 =
S10/m + S01/n (placement 성분의 표본분산, ddof=1), CI = Â ±
1.959964·SE, **clipping 없음**(계산값 그대로 엄격 비교). brute-force
전수 쌍 AUC와의 일치를 계약 테스트가 검증한다.

## 5. Secondary gate — out-of-sample BA (결정론적 CI, stage 판정 불참)

- 분할 동결: reading 인덱스 앞 절반 = TRAIN, 뒤 절반 = TEST.
- 문턱 = TRAIN 두 팔 평균의 midpoint(TRAIN에서만); 방향 동결(곡면 <
  문턱); **문턱 동률은 평탄측 배정**; 곡면 TEST는 불리한 upper 계열로
  평가.
- CI: 팔별 TEST 정확도에 **exact Clopper–Pearson 구간(측면당 α=0.0125,
  양측 97.5% 커버리지; Bonferroni 결합 95%)**, BA 구간 = 두 경계의 평균.
  내부 bootstrap 없음 — 따라서 power 인증이 단층으로 성립한다.
- Gate: BA 결합 하한 > 0.60. **자체 effect-power 90% 인증 의무**를
  지며(§6), stage 판정에는 불참(별도 문장만). 표준화 d는 기술량, gate
  없음.

## 6. Power 인증과 표본수 선택 (동결)

- **선택 규칙(사전 선언)**: "{300, 600, 1000} 중 필수 3행(AUC effect,
  AUC equivalence, BA effect)의 CP95 하한이 모두 ≥ 0.90인 **최소 n**을
  동결한다. 세 후보가 모두 실패하면 새 n을 추가하지 않고
  **blocked/review-reopened**."
- **재표집 소스(동결)**: null = **S3 curved marginal 경험분포에서 양팔
  독립 재표집, 무스케일**(AUC는 rank 통계량이라 양팔 공통 양-스케일이
  결과를 바꾸지 않는다 — 스케일에 보수성 주장 없음). effect = S3
  flat/curved marginal을 팔별 독립 재표집하되, **각 팔의 잔차(자기 평균
  중심)를 그 팔의 χ² one-sided 95% SD-상한 인수로 스케일 후 원래 팔
  평균에 재중심** — 산포는 넓어지고(검출에 보수적) 평균 간격은 정확히
  보존된다.
- 스트림: [781, 4860, n, j] (후보 n의 행 j; design lineage), B =
  20,000/행, 판정은 replicate마다 실제 gate 함수로 직접 평가.
- **결과(테스트가 재계산 pin)**: n=300에서 AUC effect **20000/20000**
  (CP95 하한 0.999816), AUC equivalence **19564/20000 = 97.8%** (CP95
  하한 0.976080), BA effect **20000/20000** (0.999816) — **선택 n/arm =
  300** (최소-n 규칙; 600/1000 행은 규칙상 불필요).
- **부정 대조(동결, 테스트 실행)**: 무효과 쌍에서 DETECTED 불발화,
  실효과 쌍에서 EQUIVALENT-AT-MARGIN 불발화, 역방향 강효과에서
  **DIRECTION-REVERSED 도달**(INCONCLUSIVE 아님), 무효과에서 BA gate
  불통과, 경계 strict 규칙.

## 7. 시드·실행 규율

- 신선 스칼라 2개: **곡면팔 40_000_251, 평탄팔 40_000_261**
  (`probe_seed_ledger` fresh 등록, entry에서 단언; smoke는 관측된
  40_000_221 스트림만).
- Freeze 커밋(이 문서 + 시드 + 러너 + 인증 + 테스트, 결과 없음) →
  정확 freeze 체크아웃 단일 실행(entry/exit git-state + 콘텐츠-주소
  manifest 양단 검증, dirty 거부) → 결과 커밋: 아티팩트 + 두 시드
  FRESH→OBSERVED + 러너 replay 전환 + **executed-freeze manifest
  스냅샷 즉시 보존**(S4/PR #59 교훈의 선반영). replay 출력은 전용
  경로만 소유하며 fresh 아티팩트를 대체할 수 없다.
- 중단 시 동일-시드 재시작. 예상 비용: 곡면팔 ≈ 2.9 h, 평탄팔 < 1분.
- freeze 후 마진·게이트·n 불변; 실행 실패·설계 모순 발견 시
  blocked/review-reopened + 새 freeze 커밋.

## 8. 동결 문장

정본은 러너 `SENTENCES`이며, 아래 인용은 그것과 **바이트 동일**하다
(줄바꿈 없이 한 줄; 계약 테스트가 exact containment로 pin). DETECTED
문장이 §1의 허용 문장이다.

- DETECTED: "동결된 Schwarzschild 도메인·밀도에서, 단일 causal set의 global relation fraction은 flat/Schwarzschild 앙상블을 우연 수준보다 판별하는 정보를 운반한다 (AUC CI95 하한 > 0.60, 프로그램 내부 진술)."
- EQUIVALENT-AT-MARGIN: "AUC가 chance의 동결 ±0.10 무시가능성 밴드 안으로 해상됐다 — 판별력 부재의 단언이 아니다."
- DIRECTION-REVERSED: "판별이 동결 방향과 반대로 해상됐다 — 동결 방향 규칙 위반, stage 실패로 기록하고 원인을 재검토한다."
- INCONCLUSIVE: "어느 분기 조건도 충족되지 않아 판정을 유보한다."
- BA 통과: "out-of-sample balanced accuracy의 결합 95% 하한이 0.60을 넘어, 학습-외 판별이 확인됐다 (secondary)."
- BA 불통과: "out-of-sample balanced accuracy는 0.60 문턱을 넘지 못했다 (secondary; stage 판정에 불참)."

## 9. 범위 제외

M-일반성(단일 M=1), 부피 오라클(예측-정박 없음; [TO CERTIFY] 이후
아크), W2/W3, 완전분리 주장(AUC 상한 게이트 없음 — "불완전 분리"는
한계 문장).
