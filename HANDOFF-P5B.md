# HANDOFF: P5-B 확정 실행 마무리 (2026-07-10 저녁 작성)

## 상황 요약

작업 위치: **worktree** `/Users/lpaiu/vs/lab/causal-spacetime/.claude/worktrees/serene-bouman-86979f`,
브랜치 `restart/m17-baseline` (최신 커밋 `859ad98`). 오늘 P4 확정 완료(PASS,
`c3363c3`), P5 동결 완료(`546bc63`). 전체 맥락은 auto-memory
(`causal-spacetime-review-findings.md`)에 있음.

**P5의 질문**: 작용 가중 2D-orders 앙상블의 연속체 상(β<β_c) 평형 샘플이
동결된 order-intrinsic 기하 판별기(P3 프로토콜, N=600 재교정 10/10)를
통과하는가? 결정 대조(β=32)는 차단되는가?

## 지금 돌아가는 것 (세션 독립, nohup detached)

PID 13221–13227, 7개 프로세스 (시작 ~19:00, **5–7시간 소요 예상**):

| 체인 | 내용 | 로그 |
|---|---|---|
| β=2 × 시드 100/101/102 | 3M스텝 → 평형 구성 3개씩 판별기 판정 | `outputs/positive_control/p5b_logs/b2_s{100,101,102}.log` |
| β=8 × 시드 100/101/102 | 〃 | `p5b_logs/b8_s{100,101,102}.log` |
| β=32 (대조, 이분 시작) 시드 100–101 | 1M스텝, 구조차단 기대 | `p5b_logs/b32_s100-101.log` |

진행 확인: `pgrep -fl "stage b --betas"` (7개), 로그 tail.
완료 판정: `outputs/positive_control/`에 CSV 7개 생성 —
`p5_stage_b_2_s100-100.csv` … `p5_stage_b_8_s102-102.csv`, `p5_stage_b_32_s100-101.csv`.

## 완료 후 할 일 (순서대로)

worktree 루트에서:

```bash
# 1. verdict 집계 (동결 게이트·기대와 대조)
cd experiments/positive_control
PYTHONPATH=../../src /Users/lpaiu/vs/lab/causal-spacetime/.venv/bin/python \
  p5_two_orders_emergence.py --stage verdict
cd ../..

# 2. 산출물 동결 디렉토리로 복사
cp outputs/positive_control/p5_stage_b_decision_registry.json docs/prereg/frozen/
head -1 outputs/positive_control/p5_stage_b_2_s100-100.csv > docs/prereg/frozen/p5_stage_b_all.csv
for f in outputs/positive_control/p5_stage_b_*_s*.csv; do tail -n +2 "$f" >> docs/prereg/frozen/p5_stage_b_all.csv; done
```

3. `docs/prereg/p5_two_orders_emergence.md` §7(Confirmatory outcome)에 결과
   기록 — **게이트·규칙 변경 금지**, 사실만. 기대 충족 여부:
   - β=2, β=8 각각: 9구성 중 게이트(heldout≤0.10 ∧ null_gap≥0.10 ∧
     truth≤0.40) 통과 ≥75%(≥7개)면 "pass" 기대 충족.
   - β=32: 통과 0개면 "block" 기대 충족(구조차단도 비통과로 집계).
   - `all_expectations_met` true → **P5 PASS** = "검증된 계기로 판정한
     기하 창발" 대표 결과. 아깝게 미달이면(예: β=8 null_gap 근소 미달)
     그대로 정량 경계로 보고 — 사후 게이트 이동 없음.
4. 커밋(관행: `P5: confirmatory outcome — ...`, Co-Authored-By 포함),
   auto-memory의 P5 항목을 결과로 갱신.
5. (선택) Paper B에 P3/P4/P5 창발 절 초안 — memory의 "전체 재시작 프로그램
   요약" 참조.

## 사고 대응

- 체인이 죽어 있으면(CSV 없음 + 프로세스 없음): 해당 샤드만 재실행 —
  위 표의 β/시드로 `--stage b --betas <β> --seeds <s>-<s>`. 결정론적이라
  안전하게 재현됨.
- β=32 대조가 h≈25 부근 층상 상태로 나와 6체인이 뽑히는 경우: 판별기
  게이트가 판정하므로 그대로 기록(구조차단이 아니어도 게이트 차단이면
  기대 충족).

## 파일 지도

- 동결 상수: `docs/prereg/frozen/p5_test_constants.json` (게이트 0.10/0.10/0.40,
  기대, 시드, 스텝)
- 사전등록: `docs/prereg/p5_two_orders_emergence.md`
- 실험 코드: `experiments/positive_control/p5_two_orders_emergence.py`
  (P3의 `analyze_order` 재사용)
- 모듈: `src/causal_spacetime_lab/positive_control/{action,two_orders}.py`
  (테스트 `tests/test_action.py`, `tests/test_two_orders.py`)
