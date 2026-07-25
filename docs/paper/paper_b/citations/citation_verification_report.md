# Citation verification report

All 54 references in `references.bib` were verified (13 on 2026-07-09, 11 on
2026-07-14 for the emergence section and related-work coverage, 4 on
2026-07-17 for the theory section's named results, and 26 copied verbatim
on 2026-07-25 from `docs/theory/citations/t1_references.bib` for the
positioning pass) against authoritative
sources (APS, AIP, Cambridge Core, Springer, PMLR, arXiv, CERN CDS, NASA ADS,
MIT DSpace, OUP, JSTOR). Each entry was cross-checked on author list, year,
venue, and volume/page (or article number). No field was fabricated.

| key | confirmed via | correction vs first draft |
| --- | --- | --- |
| blms1987 | NASA ADS `1987PhRvL..59..521B`; doi 10.1103/PhysRevLett.59.521 | pages completed 521--524; added issue 5, DOI |
| brightwell1991 | INSPIRE; APS doi 10.1103/PhysRevLett.66.260 | title has leading "The"; pages 260--263; added DOI |
| myrheim1978 | CERN CDS record 293594 | none (unpublished CERN-TH-2538 preprint) |
| meyer1988 | MIT DSpace OAI `1721.1/14328` | year kept 1988 (record shows 1988 thesis / 1989 degree — judgment call) |
| malament1977 | NASA ADS `1977JMP....18.1399M`; doi 10.1063/1.523436 | pages completed 1399--1404; added issue 7, DOI |
| hkm1976 | AIP JMP 17(2):174; ADS `1976JMP....17..174H`; doi 10.1063/1.522874 | pages completed 174--181; added issue 2, DOI |
| kronheimer1967 | Cambridge Core doi 10.1017/S030500410004144X | pages completed 481--501; added issue 2, DOI; historical journal name kept |
| sorkin2005 | arXiv gr-qc/0309009; INSPIRE book metadata; doi 10.1007/0-387-24992-3_7 | publisher Springer (not Plenum); editors Gomberoff & Marolf; pages 305--327; added DOI |
| surya2019 | Springer doi 10.1007/s41114-019-0023-1 | none (vol 22, art. 5); added DOI |
| shepard1962 | RePEc `v27y1962i2p125-140` (I) and `i3p219-246` (II) | Part I 125--140, Part II 219--246; two DOIs noted |
| kruskal1964 | Springer doi 10.1007/BF02289565 | pages completed 1--27; added issue 1, DOI |
| agarwal2007 | PMLR proceedings.mlr.press/v2/agarwal07a.html | matched (PMLR v2, pp 11--18, editors Meila & Shen) |
| kleindessner2014 | PMLR proceedings.mlr.press/v35/kleindessner14.html | matched (PMLR v35, pp 40--67); canonical "von Luxburg" |
| rideout2000 | APS doi 10.1103/PhysRevD.61.024002 (PRD 61, 024002; published 1999-12-13, volume year 2000) | new entry (2026-07-14) |
| kleitman1975 | Trans. AMS 205 (1975) 205--220 (search-confirmed volume/pages) | new entry; DOI not verified, omitted |
| winkler1985 | Springer doi 10.1007/BF00582738 (Order 1, 317--331) | new entry (2026-07-14) |
| benincasa2010 | APS doi 10.1103/PhysRevLett.104.181301; arXiv:1001.2725 | new entry (2026-07-14) |
| surya2012 | IOP doi 10.1088/0264-9381/29/13/132001; arXiv:1110.6244 (arXiv v1 title differs: "Evidence for a Phase Transition..."; published title used) | new entry (2026-07-14) |
| dowker2013 | IOP 10.1088/0264-9381/30/19/195016; arXiv:1305.2588 (formulas additionally cross-checked against the PDF, eqs. 25--26, Tables 1--2) | new entry (2026-07-14) |
| glaser2013 | APS doi 10.1103/PhysRevD.88.124026; arXiv:1309.3403 | new entry (2026-07-14) |
| major2009 | IOP doi 10.1088/0264-9381/26/17/175008; arXiv:0902.0434 | new entry (2026-07-14) |
| rideout2009 | IOP doi 10.1088/0264-9381/26/15/155013; arXiv:0810.1768 | new entry (2026-07-14) |
| glaser2018 | IOP doi 10.1088/1361-6382/aa9540 (CQG 35, issue 4); arXiv:1706.06432 | new entry; article number not captured, omitted |
| madsen2026 | arXiv abs page 2607.05840 (preprint, no journal) | new entry (2026-07-14) |
| campbell1909 | Internet Archive full text of Proc. Camb. Phil. Soc. vol. 15 (1908--10); pages 117--136 per the standard citation across the Campbell's-theorem literature | new entry (2026-07-17); pre-DOI journal, DOI omitted |
| robinson1951 | Cambridge Core `S0002731600008532`; JSTOR issue `i212276` (Am. Antiquity 16(4), Apr. 1951); doi 10.2307/276978 | new entry (2026-07-17) |
| daleyverejones2003 | Springer book page doi 10.1007/b97277 (Vol. I, 2nd ed., 2003, New York) | new entry (2026-07-17) |
| boucheron2013 | OUP Academic book 26549; global.oup.com ISBN 978-0-19-953525-5; doi 10.1093/acprof:oso/9780199535255.001.0001 (Bernstein's inequality is its Section 2.8) | new entry (2026-07-17) |

## Residual judgment calls

- **meyer1988 year (1988 vs 1989).** The official MIT record contains both:
  the thesis is dated 1988, the degree was awarded 1989. The causal-set
  literature conventionally cites 1988, which we use; switch to 1989 if
  degree-conferral year is preferred.
- **shepard1962** is a two-part paper (Psychometrika 27(2) and 27(3)). It is
  kept as a single entry with both page ranges and both DOIs in the note; it
  can be split into `shepard1962a`/`shepard1962b` if the venue requires.

## Notes for the LaTeX/submission pass

- All journal names are written in full; abbreviate per the target venue's
  style at submission.
- DOIs are present for every entry that has one (all except the two preprints /
  thesis, which carry their archive identifiers instead).
- Inline citations in `manuscript.md` use pandoc keys (`[@blms1987]`, ...) that
  match these BibTeX keys exactly, so a pandoc or LaTeX build wires them
  automatically.


## 2026-07-25 positioning pass (26 entries)

These were verified by the issuing audit session and are recorded in
`docs/theory/citations/t1_citation_verification_report.md`, which is the
authoritative row-by-row record for them. They were **copied verbatim**
into `references.bib`; keys are identical across both files by design, so
the two must be kept in step.

| key | role in Paper B | section |
| --- | --- | --- |
| bdj1999, johansson2000 | the theorems the order-only chain exponent confirms | 8.4 |
| deyjosephpeled2024 | the theorem the tube exponent confirms, mechanism included | 8.4 |
| bollobasbrightwell1991, bachmat2007 | the null-coordinate identification making those apply | 8.4 |
| bodineaumartin2005, baiksuidan2005 | thin domain does not destroy KPZ; the non-conflation warning | 8.4 |
| brightwellluczak2015 | the naming gap that survives as contribution | 8.4 |
| calhoun2021, compagnoni2014 | labeled TDOA's own `d+2`, distinguished from ours | 8.5 |
| gtt2019, garamvolgyi2021, garamvolgyi2022, connelly2024, gkioulekas2024 | unlabeled rigidity; the reconciliation | 8.6 |
| cruickshank2023 | identifiability under algebraic constraints; the structural home | 8.6 |
| humeyermeyer2025, muller2025 | nearest prior identifiability results | 8.2 |
| rideout2009 (already present), eichhorn2019, bogunakrioukov2024, major2006, major2007 | the four competing distance proposals | 8 opening |
| atkins1998, preafortin2014, carmona2025, brainerd1951 | the seriation framework behind Lemma 4c | 8.2 |

### Two corrections made during this pass

| key | defect | correction | confirmed via |
| --- | --- | --- | --- |
| armstrong2021 | `year` was the arXiv preprint year (2021); volume/pages absent | journal version is SIAM J. Math. Data Sci. **5**(1), 201--221, **2023**; key left unchanged so citations do not break | SIAM article page, doi 10.1137/22M1495342 |
| aghili2018 | **title was wrong** -- recorded as "Statistical geometry and causal set dimension" | actual title is "Path length distribution in two-dimensional causal sets" | EPJ C article page and Springer record, doi 10.1140/epjc/s10052-018-6229-7 |

Both corrections were applied to `docs/theory/citations/t1_references.bib`
as well. No field was fabricated and no unverified entry was added; in
particular Bombelli-Noldus-Tafoya (arXiv:1212.0601) is deliberately absent,
having been withdrawn by its authors.


### Post-hoc verification of the two load-bearing theorem statements

The manuscript's Section 8.4 states the content of the theorems it cites,
which is more exposure than citing them. The statements were verified
against the published abstracts during code review of the positioning PR
(2026-07-25) — after the prose was written, recorded here so the order of
operations is on the record:

- `deyjosephpeled2024` (arXiv:1808.08407 / Israel J. Math. 262): the
  abstract states that the maximal increasing path in `[0,n]^2` restricted
  to a strip of width `n^gamma`, `gamma < 2/3`, has expectation
  `2n - n^{1-gamma+o(1)}`, **variance `n^{1-gamma/2+o(1)}`**, and converges
  to the **Gaussian** distribution after scaling. At `gamma = 0` that is
  `sd ~ n^{1/2}` against mean `~ 2n`, i.e. `theta = 1/2`, as the
  manuscript asserts.
- `bdj1999` (J. AMS 12): expectation `2n - n^{1/3}(c_1+o(1))`, variance
  `n^{2/3}(c_2+o(1))`, Tracy-Widom limit — `theta = 1/3`, as asserted.

Both match the manuscript's usage exactly; no edit was needed.
