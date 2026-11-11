# GRG 재투고 준비 기록 (2026-11-12)

대상: *General Relativity and Gravitation*, Research Article, 구독 경로.
원고는 2026-09-16에 거절된 Paper A의 수정본이 아니라 새 원고다.

## 패키지 만들기

    python docs/paper/weyl_note/build_submission.py outputs/submissions/grg-2026-11-12

| 파일 | 용도 |
|---|---|
| `manuscript.pdf` | 원고 10쪽. 패키지의 `source/`에서 직접 컴파일한 것 |
| `manuscript_source.zip` | `main.tex`와 `Fig1.pdf`. 편집 가능한 원본(필수 제출물) |
| `cover_letter.pdf`, `cover_letter.md` | 커버레터 1쪽 |
| `abstract.txt` | 온라인 양식용 초록, 245단어(허용 150–250) |

## 투고 전에 본인이 해야 하는 일

1. **분기 공개.** 원고의 데이터 가용성 문장과 커버레터는 코드와 기록이
   GitHub에 공개돼 있다고 적는다. 커밋 `0ff9802`, `c14d4a2`, `a54ac97`이
   공개 저장소에서 열리는지 확인한다. 열리지 않으면 그 문장은 거짓이다.
2. **원고 통독.** 4.3절의 "The author has examined the results and the text
   and takes full responsibility for them"은 본인이 실제로 읽은 뒤에만 참이다.
3. **각주 1.** Wang의 식 (79)를 재현하지 못했다는 문장이다. 직접 확인하지
   않을 것이면 그 문장만 지운다(계수 자체와 수치 검증 문장은 남긴다).
4. **두 논문의 인용 문장.** Madsen 정리 4.9의 오차율과 결론부, Ibrahim–Bombelli
   6절의 식 (19)–(20), 상자 부피와 두께 법칙을 원문과 한 번 대조한다. 이
   문장들은 arXiv HTML을 자동 요약 도구로 읽어 확인한 것이다.
5. **저자 정보.** 이름·소속·이메일 표기, ORCID가 있으면 포털에 입력한다.
6. **비용.** 오픈 액세스 등 유료 선택 항목은 고르지 않는다. Research Square
   프리프린트 공개 여부는 본인이 정한다.
7. **서버 PDF.** 포털이 `manuscript_source.zip`에서 만든 PDF를 열어 표·그림·
   참고문헌이 온전한지 확인한다. 서버는 pdfLaTeX를 쓰고, 로컬 빌드는
   Tectonic(XeTeX 계열)이어서 서버 변환은 여기서 시험하지 못했다.
8. **추천 심사자(선택).** 포털에서 물으면 본인이 정한다.

## 포털 입력란에 붙여 넣을 문구

- **Keywords:** Causal sets; Weyl tensor; Curvature estimation; Plane gravitational waves; Poisson sprinkling
- **Funding:** No funding was received for conducting this study.
- **Competing interests:** The author has no competing interests to declare that are relevant to the content of this article.
- **Author contribution:** J.K. is the sole author: J.K. conceived and directed the study, examined the results and wrote the manuscript, with the assistance of the AI tools described in Sect. 4.3.
- **Data availability:** The source code, the pre-specified protocols, the summaries of all runs and the records of the checks are available at https://github.com/lpaiu-cs/causal-spacetime (directories experiments/positive_control and docs/paper/paper_a/intrinsic_weyl). The order archives can be regenerated from the recorded seeds and are available from the author on reasonable request.

## 이번 검토에서 고친 것

**인용**

- 제목 오기 2건: Brightwell–Gregory와 Benincasa–Dowker의 제목에 없는 관사 "The"를 삭제.
- 설명 부족 1건: Alfyorov–Shnyukov가 pp-wave 표본에서 수치 실험을 했다는 사실을 추가.
- 빠진 인용 7건 추가: Zeeman 1964, Szekeres 1965, Kambor–Nomaan 2021,
  Barton–Borza–Röhrig 2026, Ghosh 2026, NumPy, SciPy.
- DOI가 있는 23건은 Crossref 기록과, arXiv·Zenodo·CERN·MIT 문헌 10건은 각
  기록과 저자·제목·권호·쪽을 대조했다.
- 참고문헌을 학술지 서식(저자: 제목. 학술지 권, 쪽 (연도). DOI 링크)과 첫 인용
  순서로 정리했다.

**표현**

- 내부 용어 삭제: frozen, gate, verdict, "committed before data" → pre-specified
  criterion, outcome 등.
- 제목의 "an executed", "Draft for author review", "that verdict stands",
  "not a contest", "desktop computer" 수정.
- 초록에서 인용 번호를 없애고(규정) 245단어로 줄임.

**삭제한 내용**

- 포함률 실험이 한 번 중단됐다는 실행 기록 문장(저장소 기록에는 남아 있음).
- 신뢰공이 5억–6억 사건에서 참값을 놓친다는 외삽.
- 사후에 추가한 둘째 통계량(평탄 적합 대비 카이제곱 감소).
- 우리 쪽의 이전 외삽이 틀렸다는 내부 경과.

**규정 대응**

- 키워드 5개, "Statements and Declarations" 구획(자금, 이해관계, 자료·코드, 저자 기여).
- AI 사용 내용을 방법 절(4.3)에 기술.
- 소제목을 번호가 있는 3단계 이내 구조로 변경.
- 그림: 그림 안 제목 삭제, 패널 문자 a/b, 폭 119 mm, 글자 8 pt, 색 대신 채움·모양으로 구분.
- 원본은 ASCII만 쓰고 pdfLaTeX 분기를 넣음.
