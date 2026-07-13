# Citation verification report

All 24 references in `references.bib` were verified (13 on 2026-07-09, 11 on
2026-07-14 for the emergence section and related-work coverage) against
authoritative sources (APS, AIP, Cambridge Core, Springer, PMLR, arXiv, CERN
CDS, NASA ADS, MIT DSpace). Each entry was cross-checked on author list, year,
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
