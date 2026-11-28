# Sub-discreteness gravitational-wave detector 논문 투고 준비 기록

## 패키지 생성 방법

```powershell
python docs/paper/subspacing_wave/make_figures.py
python docs/paper/subspacing_wave/build_submission.py outputs/submissions/subspacing-wave-final
```

## 최종 점검 완료 내역 (Step 1--4)

1. `note.tex` 전 절(`Abstract`, `Sect. 4.2`, `Sect. 5.1`, `Sect. 5.2`, `Sect. 6`, `Statements and Declarations`) 완성 및 `S1`(`36e551a`, 16/18 통과 및 $\omega \le 32$ 저주파 잔여 꼬리 진단) -> `S2`(`4192369 -> 022166c`, 20/20 게이트 전수 CONFIRMED, $N \le 2^{23}$) 검증 JSON 수치 반영.
2. 초록 단어 수 `150-250` 및 순수 ASCII 계약 (`build_submission.py`) 자동 검증 완료.
3. `Fig1.pdf`, `Fig2.pdf`, `Fig3.pdf` vector PDF 생성 및 시각적 검수 완료 (임시 `Fig2_continuum_preview.pdf` 제거).
4. 독립 clean-room 검증기(`check_subspacing_confirmatory.py`, `check_subspacing_s2_confirmatory.py`) 통과 결과(`audit_receipt.json`)와 `P1->F1->R1` / `P2->F2->R2` 커밋 계보를 Data and Code Availability에 명시.
5. `cover_letter.md`, `cover_letter.pdf`, `manuscript.pdf`, `manuscript_source.zip`, `abstract.txt` 패키지 빌드 완료.
