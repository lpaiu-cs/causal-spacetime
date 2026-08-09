# S4 사전등록 — Schwarzschild C1-paired 확증

상태: **freeze 대상 초안** (인세션 설계 리뷰 완료 후 freeze 커밋으로 동결;
결과는 이후 커밋에만 등장한다). 선행: S2 설계 메모(vault), S3 탐색
(`p14_s3_probe_results.json`, PR #56 merge `cb520b0`).

## 1. 목적과 주장 클래스

S3 탐색이 실측한 paired 이동을 사전등록 게이트로 확증한다. 주장 클래스는
**C1 counterfactual sensitivity**(같은 점집합의 이중 census)로 한정한다 —
C2 단일-poset 판별이 아니며, 모든 문장은 프로그램 내부 진술이고 아래 동결
좌표·도메인에 결부된다.

## 2. 동결 프로토콜

S3 탐색과 동일(변경 없음):

- 도메인: Schwarzschild 좌표, M=1, r ∈ [10, 20] 외부 셸, 극관각 캡 반각
  1.0(pairwise Δψ ≤ 2.0), t ∈ [0, 40]
- 측도: 공통 강도 ρ·r²sinθ (질량 무관 — 측도 항등성), reading당
  N ~ Poisson(300)
- 술어: 곡면 = S1 solver tol 1e-8, undecided 시 1e-10 에스컬레이션;
  평탄 = 정확 chord. undecided-is-never-a-silent-False — 잔여 undecided는
  identified 구간이 흡수한다(별도 정지 규칙 없음)
- 추정량: θ_Δ = E[f_M − f_0], reading 단위 paired
- 동결 방향: 음(−) — Shapiro 지연 기제, S3 실측 부호

## 3. 블록

- **S4 블록**: 신규 n = 300 readings, 신선 시드 **40_000_241**
  (`probe_seed_ledger`에 fresh 등록; smoke는 관측된 40_000_221 스트림만)
- **비교 블록(동결 데이터)**: S3 공식 아티팩트
  `docs/prereg/p14_s3_probe_results.json` — n₃ = 300,
  Δ̄₃ = −0.0360870056, s₃ = 0.0012108375 (per-reading 원자료 보존)

## 4. Gate 구조

모든 게이트는 **identified 양**으로 평가하며(§6의 정확 공식), 모든 비교는
**엄격 부등호**다 — 정확히 경계에 닿은 값은 통과하지 않는다.

- **Gate A (primary — C1 검출)**: S4의 identified CI95 상단 < −ε_det
- **Gate B (secondary — 정량 재현)**: 독립 두-표본 차이 θ_S4 − θ_S3의
  identified Welch CI95로 삼분 판정:
  - **REPLICATED**: CI 전체가 (−ε_rep, +ε_rep) 안
  - **DISCORDANT**: CI 전체가 밴드 밖 (하단 > +ε_rep 또는 상단 < −ε_rep)
  - **B-INCONCLUSIVE**: 그 외 (경계에 걸침)

**stage 종합 판정(완전 분할)**:

| 조건 | 판정 |
| --- | --- |
| A ∧ B=REPLICATED | **CONFIRMED** |
| A ∧ B≠REPLICATED | **DETECTED-NOT-REPLICATED** |
| S4 identified CI95 ⊂ (−ε_det, +ε_det) | **NO-DETECTION** |
| 나머지 | **INCONCLUSIVE** |

(A 통과 시 equivalence 조건과는 상호배타적이므로 위 분할은 전 경우를
정확히 한 번씩 덮는다.)

**포함관계(동결 명시)**: Gate B=REPLICATED는 Gate A 통과를 강제한다 —
Welch 반폭 h_W는 S4 단독 반폭 h₄보다 크므로 B=REPLICATED 시
Δ̄₄ + h₄ < Δ̄₃ + ε_rep = −0.034887 < −ε_det = −0.0036. 따라서
**CONFIRMED = A∧B이지만 현재 마진에서는 B=REPLICATED와 동치**다. 두
판정은 별도 문장으로 보고하며, effect-side joint power는 Gate B power와
동일하므로 별도의 joint-power 주장은 두지 않는다.

## 5. 마진 (동결 후 불변)

freeze 전에는 n을 포함해 전 항목이 재심의 가능하다. **freeze 후에는
마진·게이트·n 모두 불변**이며, 실행 실패나 설계 모순이 발견되면 결과
해석이 아니라 **blocked/review-reopened 처리 후 새 freeze 커밋**이
필요하다.

- **ε_det = 0.0036** — 앵커 |Δ̄₃| = 0.0360870의 **약 9.98%(≈약 10%)**.
  운영 문장: "탐색 효과의 최소 약 10%를 유지하는 이동만 검출로 인정한다."
- **ε_rep = 0.0012** — S3 single-reading SD 0.0012108375에서 유도한
  **반올림 운영 마진**. 프로토콜·밀도(E[N]=300, 이 도메인) 의존 기준이며
  불변량이 아니다.

## 6. identified 구간의 정확 공식

ambiguity가 0이면 아래의 lower/upper 계열이 일치하고 통상 CI로 환원된다.
각 끝점은 97.5% one-sided coverage이므로 바깥 끝점 결합은 union bound로
보수적 95% 결합구간이다.

- S4 자체 (Gate A·NO-DETECTION):
  [ Δ̄₄,lower − t_{n₄−1} · s₄,lower/√n₄ ,
    Δ̄₄,upper + t_{n₄−1} · s₄,upper/√n₄ ],
  여기서 lower/upper는 undecided 쌍을 각각 무관계/관계로 계상한 per-reading
  계열, t는 양측 95%(one-sided 97.5%) Student-t 임계값.
- Welch 차이 (Gate B):
  [ (Δ̄₄,lower − Δ̄₃) − t_{ν_lo} · SE_lo ,
    (Δ̄₄,upper − Δ̄₃) + t_{ν_hi} · SE_hi ],
  SE_b = √( s₄,b²/n₄ + s₃²/n₃ ),
  ν_b = ( s₄,b²/n₄ + s₃²/n₃ )² / ( (s₄,b²/n₄)²/(n₄−1) + (s₃²/n₃)²/(n₃−1) )
  (Welch–Satterthwaite; b ∈ {lower, upper}). 비정수 ν는 **floor(ν)**로
  내림해 임계값을 취한다(작은 df → 큰 t → 보수적).

## 7. Power 인증 (동결 n = 300)

- 목표: **각 branch 90% 이상**
- 보수 분산 입력: **SD_cons = 0.001298615** — S3 공식 SD의 χ² one-sided
  95% 상한(df = 299), 양 블록 모두 적용
- 사전 고정 alternative: Gate A·null은 각각 θ = Δ̄₃, θ = 0; Gate B는 주
  alternative θ_S4 = θ_S3와 **assurance 영역 |θ_S4 − θ_S3| ≤ κ = 3×10⁻⁴**
  (ε_rep의 1/4)에서의 최악 power
- **재표집 규칙(동결)**: 소스 분포 = S3 per-reading 원자료를 중심화 후
  SD_cons로 스케일한 경험분포. **Gate B는 S4 블록(n=300)과 S3 블록
  (n=300)을 모두 독립 재표집**한다(관측 S3를 고정하는 조건부 power가
  아니라 replication power). RNG: [781, 4840]은 probe 체인 관례의
  **entropy vector**(`np.random.default_rng([781, 4840, j])` 의미론,
  spawn key 아님)이고, 표의 행 j마다 child stream [781, 4840, j] 하나를
  쓰며, replicate 안에서 **S4 블록을 먼저, S3 블록을 다음에** 뽑는다.
  design-lineage 부계열이며 캠페인 스트림이 아니다.
- **재표집 인증** (B = 20,000/행):

  | 행 j | Branch | alternative | power | CP95 하한 |
  | --- | --- | --- | --- | --- |
  | 0 | Gate A | θ = Δ̄₃ | 20000/20000 = 1.00000 | 0.9998156 |
  | 1 | Gate B | θ_S4 = θ_S3 | 20000/20000 = 1.00000 | 0.9998156 |
  | 2 | Gate B | θ_S4 = θ_S3 − κ | 20000/20000 = 1.00000 | 0.9998156 |
  | 3 | Gate B | θ_S4 = θ_S3 + κ | 20000/20000 = 1.00000 | 0.9998156 |
  | 4 | NO-DETECTION | θ = 0 | 20000/20000 = 1.00000 | 0.9998156 |

  Clopper–Pearson exact 95% 하한 0.9998156 ≥ 0.90 — 전 분기 목표 인증.
- 해석 교차검증(표준화 여유, SD_cons 기준): Gate B 최악 z ≈ 6.52,
  Gate A z ≈ 431.3, null z ≈ 46.0.
- **부정 대조(동결, 계약 테스트가 실행)** — 게이트가 항상 통과하는
  결함을 배제하는 반증가능성 검사. child stream [781, 4840, 99]의 구성
  블록으로: (i) θ = −ε_det(경계)에서 Gate A 불통과, (ii) θ = −ε_det/2
  (문턱 안쪽)에서 Gate A 불통과, (iii) 차이 = +3ε_rep에서 Gate B
  ≠REPLICATED, (iv) θ = +2ε_det에서 NO-DETECTION 불통과, (v) A-only 구성
  (θ_S4 = Δ̄₃ + 3ε_rep): Gate A 통과 ∧ Gate B ≠REPLICATED — 종합 판정
  DETECTED-NOT-REPLICATED 경로 도달 확인. 경계 규칙: 모든 비교는 엄격
  부등호이며 정확 경계값은 불통과(§4).
- 위 수치·대조 전부를 계약 테스트가 저장 원자료에서 재계산해 pin한다.

## 8. 실행 규율

- **Freeze-commit 분리**: 이 문서 + 시드 등록 + 러너 + 테스트가 freeze
  커밋(결과 없음) → 정확히 그 freeze 내용의 클린 체크아웃에서 단일 실행
  → 결과는 이후 커밋
- **콘텐츠-주소 freeze 검증**: 프로토콜 표면 8개 파일 — 이 문서, 러너
  `s4_schwarzschild_c1.py`, `s3_schwarzschild_probe.py`,
  `s1_schwarzschild_cost.py`, `p14_probe_p2.py`, `seed_windows.py`,
  `probe_seed_ledger.py`, S3 비교 아티팩트
  `p14_s3_probe_results.json` — 의 raw SHA-256을
  `p14_s4_freeze_manifest.json`에 동결하고, 러너가 **entry에서 검증,
  exit에서 재검증**한다. clean이지만 freeze 이후의 커밋(드리프트된
  프로토콜·비교 데이터)에서 일회성 신선 시드를 소모하는 경로를
  차단하며, 커밋 SHA가 아닌 내용 주소라 merge 커밋에도 생존한다.
  manifest 자체의 무결성은 freeze 커밋 리뷰가 정박한다. (해시 대상
  경로는 전부 `.gitattributes`로 LF 핀 — 워킹트리 bytes == blob bytes.)
- 실행 계보: entry/exit git-state 봉인(cwd=repo, check=True),
  `run_kind: fresh_observation`; 결과 커밋에서 40_000_241을 OBSERVED로
  이동하고 러너를 replay 경로로 전환
- 중단 시 동일-시드 재시작, 시드 교체 불가
- 예상 비용: ≈ 2.8 h (S3 실측 기준)

## 9. 동결 문장

- **Gate A (검출)**: "동결된 Schwarzschild 좌표·도메인(M=1, r∈[10,20],
  극관각 캡 1.0, T=40)의 공통 측도 위에서, paired 앙상블 평균 이동의
  identified CI95가 동결 방향으로 검출문턱 ε_det = 0.0036을 초과했다 —
  유한밀도 인과 census가 type D 진공 곡률의 빛원뿔 변형을 C1-급으로
  검출했다(프로그램 내부 진술)."
- **Gate B REPLICATED**: "S4 블록과 S3 탐색 블록의 독립 두-표본 차이의
  identified Welch CI95가 ±ε_rep = ±0.0012 안에 들어, 탐색 효과가
  정량적으로 재현됐다."
- **Gate B DISCORDANT**: "차이의 identified Welch CI95가 ±ε_rep 밖에
  전부 놓여, S4 효과 크기가 탐색 효과와 불일치한다."
- **Gate B INCONCLUSIVE**: "차이의 identified Welch CI95가 ±ε_rep
  경계에 걸쳐 재현 판정을 유보한다."
- **DETECTED-NOT-REPLICATED (종합)**: "검출문턱은 초과했으나 탐색 효과의
  정량 재현에는 이르지 못했다 — 검출 판정은 유효하되 효과 크기는 재심의
  대상이다."
- **NO-DETECTION**: "동결 조건과 검출문턱 ε_det에서 평균 이동이 0과
  동등한 것으로 판정됐다."
- **INCONCLUSIVE**: "어느 분기 조건도 충족되지 않아 판정을 유보한다."

## 10. 범위 제외

M-사다리(M ≠ 1), C2-unpaired 판별, 부피 오라클 예측-정박은 이 사전등록의
범위 밖이다(각각 이후 탐색 / [TO DESIGN] / [TO DERIVE]).
