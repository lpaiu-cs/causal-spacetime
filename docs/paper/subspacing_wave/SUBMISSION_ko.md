# Sub-discreteness gravitational-wave detector 논문 투고 준비 기록 (초안)

## 패키지 생성 방법

```powershell
python docs/paper/subspacing_wave/make_figures.py
python docs/paper/subspacing_wave/build_submission.py outputs/submissions/subspacing-wave-draft
```

## Step 4 완료 후 최종 점검 체크리스트

1. `note.tex` 내 `[PENDING ...]` 마커 전부 실제 검증 JSON 수치로 교체 (Step 3 유한 $N$ 크로스오버 + Step 4 사전등록 확증 수치·95% 양측 CI·교차검증 감사 결과).
2. 초록 단어 수 `150-250` 확인 (`build_submission.py` 자동 검증), 순수 ASCII 유지 (`pdfLaTeX` 서버 호환).
3. `Fig1.pdf`, `Fig2.pdf`, `Fig3.pdf` 생성 및 시각적 검수 (8 pt 산세리프, 폭 119 mm, 흑백 인쇄 대비 모양/채움 구분).
4. 독립 검증기(`check_subspacing_confirmatory.py`) 실행 결과·커밋 SHA·Clean worktree 기록을 Data and Code Availability에 명시.
5. `cover_letter.md`의 `[PENDING ...]` 요약 갱신 후 `cover_letter.pdf`·`manuscript.pdf`·`manuscript_source.zip`·`abstract.txt` 패키지 재빌드.
