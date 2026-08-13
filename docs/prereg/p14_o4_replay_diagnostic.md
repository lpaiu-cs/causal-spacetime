# P14 O4 — G3 abort의 replay 전용 진단 (프로토콜)

동결 `1eb9461` 에서 실행된 O4 캠페인은 G3 에서 fail-closed abort 했고
과학적 verdict 를 발행하지 않았다 (`p14_o4_incident.json`,
`p14_o4_incident.md`). 그 러너는 **실패 좌표를 기록하지 않았다** — v1 의
관측성 결함이다. 이 문서는 그 좌표를 사후에 읽어내기 위한 진단의
프로토콜이며, 진단 자체는 `experiments/oracle/o4_replay_diagnostic.py`
에 있다.

이 문서에는 **결과가 없다.** 진단 실행은 별도 승인 사안이고, 산출물은
승인 후 별도 커밋에서 `p14_o4_replay_diagnostic.json` 으로 게시된다.

## 1. 이것이 아닌 것

**replay 는 재현이지 관측이 아니다.** 여기서 나오는 어떤 값도 과학적
결과로 승격될 수 없다. 구체적으로:

- G1 통계량을 계산하지 않는다. G2 는 아예 재실행하지 않는다.
- 어떤 게이트도 status 를 갖지 않고, verdict 는 존재하지 않는다.
- 아티팩트는 `run_kind: "replay"`, `replay_of: "PR #70 incident"` 로
  표시되며 결과·게이트 스키마와 **키가 겹치지 않는다.** 계약 테스트가
  `mean`·`var`·`half_width`·`v_s1`·`identified`·`leak`·`cp_upper`·
  `concordant`·`discordant`·`power`·`band`·`verdict`·`gate`·`status`·
  `estimate` 를 포함하는 키가 하나라도 있으면 실패시킨다.
- **이 진단의 결과로 G3 규칙을 바꾸지 않는다.** 재설계(G3a/G3b)는
  사전등록을 다시 여는 별도 절차다. 진단은 재설계의 입력이지 근거가
  아니다.

## 2. 시드와 예약 경계

| 대상 | 처분 |
|---|---|
| `o4_aborted_g1` (40,000,281) | `replay_scalar` 로만 재생. 은퇴 상태 유지 |
| `o4_aborted_g2` (40,000,291) | 건드리지 않음 — G2 는 재실행 대상이 아니다 |
| `40,000,301` | 미할당 상태 그대로 보존 |
| `refs/o4/reservation` | 읽지도 쓰지도 않음. 영구 보존 |
| `FRESH_PROBE_SCALARS` | 비어 있어야 하며, 진단은 여기서 뽑지 않는다 |

진단은 `assert_fresh_scalar` 를 호출하지 않는다. 예약 권위에 닿는 함수
(`reserve_remote` · `remote_reservation` · `reservation_authority` ·
`probe_reservation_namespace`)와 `subprocess` 자체를 **AST 수준에서**
사용하지 않음이 테스트로 강제된다. 산문에서 예약을 언급하는 것과
실제로 ref 를 건드리는 것을 구분하기 위해 문자열 검색이 아니라 구문
트리를 본다.

## 3. 무엇을 어디까지 재현하는가

동결 시점의 코드로 재생해야 재현이다. 진단은 **executed snapshot**
(`p14_o4_executed_freeze_manifest.json`)을 기준으로 자기 자신을
검증한다 — 현재 manifest 는 이미 움직였기 때문이다(abort 기록이 캠페인
스칼라를 은퇴시키면서 원장 다이제스트가 바뀌었다).

- snapshot 의 모든 파일이 바이트 단위로 일치해야 한다.
- 예외는 `probe_seed_ledger.py` **하나**뿐이며, 그 근거는 "이 계산에
  기여하는 것이 시드 하나뿐"이라는 것이다. 그래서 그 시드 값을 진단
  모듈에 못 박고 원장과 대조한다 — 면제하되 검증한다.
- 환경 잠금(python·gmpy2·MPFR·GMP·numpy)도 동일하게 대조한다.

재현 범위는 **원래 G3 stress point 생성에 필요한 G1 prefix 까지**다.
G1 의 추첨 루프(같은 generator, 같은 chunk 경계, 같은 수락 조건
`L_S1 > 0`)를 그대로 밟되, 누적 통계는 하나도 만들지 않고, 요청한
클러스터 수를 채우면 즉시 멈춘다.

**interleaving 은 근사가 아니다.** 동결 러너는 G1 안에서 stress 목록을
쌓고 G3 에서 따로 훑지만, 한 클러스터의 4회 probe 호출은 그 클러스터
자신의 `(r, θ, T₁, T₂)` 만의 결정론적 함수다. 따라서 수락 즉시 probe
해도 `causal_relation` 에 들어가는 입력은 동결 순서와 바이트 단위로
같다. 이 성질 덕분에 `--clusters` 를 작게 잡아 앞쪽 실패만 싸게 볼 수
있다.

## 4. 무엇을 기록하는가

동결 G3 는 첫 undecided 에서 죽는다. **진단은 죽지 않고 계속 간다** —
동결 러너가 원리적으로 볼 수 없었던 것이 바로 그 뒤이기 때문이다.

첫 undecided 에 대해:

- cluster index, 좌표 `(r, θ)`, `T₁`, `T₂`, window `[lo, hi]` 와 `L`
- probe 종류(`midpoint` / `outside`)와 leg(`p→x` / `x→q`)
- `dt`, `t_min`, `err`, `|dt − t_min|`, 판정 여유 `|dt − t_min| − err`
- `dpsi` — 술어가 `acos(cos θ)` 로 되찾는 각도. `_ell` 에 넘어간 `θ`
  와 마지막 ulp 에서 다를 수 있으므로 **술어가 실제로 쓴 값**을 쓴다.

원인별로도 첫 발생 지점을 같은 형식으로 남기고, 고정 stress 집합 전체를
훑은 경우 원인별 개수를 보고한다:

| 원인 | 조건 |
|---|---|
| `midpoint-in-error-band` | midpoint probe 의 어느 leg 가 undecided |
| `outside-in-error-band` | outside probe 의 어느 leg 가 undecided |
| `boolean-mismatch` | 두 leg 모두 판정됐으나 `bool(a and b) ≠ want` |

**이 개수는 진단 빈도이지 추정량이 아니다.** 은퇴한 스트림에서 재현한
하나의 고정 stress 집합 위의 도수이며, 어떤 구간도 붙지 않는다.
아티팩트가 이 문장을 직접 싣고 계약 테스트가 그것을 강제한다.

읽기가 술어를 오독하지 않았음도 함께 검증한다. 각 leg 에서 술어가 쓰는
것과 같은 `flight_time` 호출로 `t_min`·`err` 를 재유도하고, 그로부터
예측되는 삼치값을 `causal_relation` 의 실제 반환값과 대조해 불일치
개수를 센다. 이 값이 0 이 아니면 **술어가 아니라 진단이 틀린 것**이다.

## 5. 진단이 겨냥하는 항등식

사건 기록이 실행 없이 대수적으로 닫은 표다. `lo = T₁`,
`hi = Δt − T₂`, `L = hi − lo` 로 두면:

| probe | leg | `\|dt − T_min\|` | undecided 조건 |
|---|---|---|---|
| midpoint | `p→x` | `L/2` | `L/2 ≤ err₁` |
| midpoint | `x→q` | `L/2` | `L/2 ≤ err₂` |
| outside `hi+1e-6` | `p→x` | `L + 1e-6` | `L + 1e-6 ≤ err₁` |
| **outside `hi+1e-6`** | **`x→q`** | **`1e-6` 정확히** | **`err₂ ≥ 1e-6`** |

마지막 행에 기하가 들어가지 않는다는 것이 요점이다. 진단은 이 항등식을
계약 테스트로 못 박는다(실제 클러스터에서
`|dt − t_min| == 1e-6 ± 1e-12`, midpoint 는 `L/2`).

여기서 미리 말해 둘 것이 하나 있다. offset `1e-6` 은 `tol = 1e-8` 보다
두 자릿수 위로 잡힌 값이지만, 솔버가 보고하는 `err` 는 **`tol` 로
상계되지 않는다** — `err` 는 인증된 격납이 아니라 Gauss–Legendre 정지
휴리스틱이기 때문이며, 이는 O4 freeze 가 `err` 를 추론에서 배제한
바로 그 이유다. 따라서 마지막 행의 조건이 실제로 발동하는지는 기하가
아니라 솔버의 오차 보고가 결정한다. 진단은 그 사실을 측정할 뿐, 이
문서에서 결론내지 않는다.

## 6. 실행

```
python experiments/oracle/o4_replay_diagnostic.py --clusters 2000
python experiments/oracle/o4_replay_diagnostic.py \
    --out docs/prereg/p14_o4_replay_diagnostic.json
```

`--out` 없이는 아무것도 쓰지 않는다. `--out` 은 no-clobber 이며 첫
게시가 확정이다. `--clusters` 는 동결 `g3_clusters = 100,000` 을 넘을
수 없다 — stress 집합은 동결돼 있고 진단이 그 밖으로 나갈 수 없다.

## 7. 이 진단이 닫지 않는 것

- G3 재설계. `[T₁+err₁+η, Δt−T₂−err₂−η]` 적격 규칙과 G3a/G3b 분리는
  사전등록 재개방에서 다룬다.
- O4 재실행. 새 freeze 와 새 G1/G2 스칼라가 필요하고, 실행 승인은
  언제나 동결 승인과 별도다.
- 러너의 관측성 결함. 향후 러너는 모든 abort 경로에서 write-once
  incident artifact 를 남겨야 하고, G3 는 실패 좌표를 스스로 기록해야
  한다. 이 진단은 그 부재를 사후에 메우는 임시 수단이지 그 요건을
  대체하지 않는다.
