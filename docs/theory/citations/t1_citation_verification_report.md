# T1 citation verification report

Companion to `t1_references.bib`. Verified 2026-07-25 by targeted literature
search against arXiv, publisher pages, DOI records, JSTOR, Project Euclid,
INSPIRE-HEP and NASA ADS. Same rule as the Paper B report: **no field is
fabricated**; a field that could not be confirmed is omitted from the BibTeX
entry and recorded below.

Scope note: this bibliography was assembled to answer a specific question —
*what is already known that bears on the T1 claims* — not to be exhaustive.
Sections 6-8 (rigidity, localization, KPZ) are new territory for this
repository and are the ones with residual risk.

## 1. Confidence tiers

| tier | meaning | count |
| --- | --- | --- |
| **A** | author list, title, venue, volume/number, pages and DOI (or arXiv ID) all confirmed | 34 |
| **B** | identity confirmed; one bibliographic field omitted as unconfirmed | 14 |
| **C** | preprint with no journal version — arXiv ID is the identifier | 7 |

No entry is below tier C. Entries carried over from
`docs/paper/paper_b/citations/references.bib` retain that file's verification
(2026-07-09 / 07-14 / 07-17) and were not re-checked here; they are marked
`[also in paper_b]` in the `.bib`.

## 2. Per-entry residuals (tier B and C)

| key | what is missing | why it was left out |
| --- | --- | --- |
| `bollobasbrightwell1991` | DOI | JSTOR/AMS returned 403. Volume 324(1), pp. 59--72 confirmed via the JSTOR stable record `2001495`. The AMS DOI pattern is predictable but was not confirmed, so it is omitted rather than constructed. |
| `atkins1998` | DOI | SIAM publisher page returned 403. Volume 28(1), pp. 297--310 confirmed from several independent secondary listings. |
| `asimowroth1978` | DOI | Trans. AMS 245 (1978) 279--289 confirmed; DOI not retrieved. A Part II in *J. Math. Anal. Appl.* 68 (1979) is commonly cited but was **not** independently confirmed and is therefore absent. |
| `schaurobinson1987` | DOI | IEEE record not reached; volume/number/pages confirmed from consistent secondary citation. |
| `connelly2024` | volume, pages | Online-first in *Combinatorica* at time of check; DOI confirmed. |
| `huangdokmanic2021` | pages | IEEE TSP vol. 69 and DOI confirmed; page range not retrieved. |
| `baiksuidan2005` | volume, pages | IMRN 2005; arXiv ID confirmed. |
| `carmona2024siam`, `carmona2025` | volume, pages | Both confirmed at DOI / ScienceDirect record level; issue assignment not retrieved. |
| `armstrong2021` | volume, pages, year of record | SIAM J. Math. Data Sci.; DOI and arXiv ID confirmed. The BibTeX `year` is the arXiv year, not the journal year — **correct this before submission.** |
| `aghili2018` | title wording | Confirmed at DOI-record level (Eur. Phys. J. C 78, 744, 2018) with the author list; the exact published title was not read off the article page. **Re-check the title string before citing.** |
| `eren2004`, `thrun2005`, `pollefeysnister2008`, `kuang2013` | DOI | Conference proceedings; identity and pages confirmed where shown, DOIs not retrieved. |
| `ferranti2021` | — | Tier A (DOI + arXiv both confirmed); listed here only because the specific minimal-configuration tables inside it were **not** extracted (PDF text extraction failed). Do not cite specific `(m, n)` pairs from it without re-reading the published version. |
| `myrheim1978` | — | Unpublished CERN preprint CERN-TH-2538; INSPIRE record confirms title/author/number/year. |
| `meyer1988` | — | MIT DSpace handle `1721.1/14328` confirmed at search level only. The Paper B report's 1988-vs-1989 judgment call carries over unchanged. |
| `madsen2026`, `muller2025`, `sivaramakrishnan2025`, `cruickshank2023`, `brightwellluczak2015`, `calhoun2021` | journal reference | Preprints with no published version found. |

## 3. Entries deliberately excluded

- **Bombelli, Noldus & Tafoya, "Lorentzian manifolds and causal sets as partially ordered measure spaces" (arXiv:1212.0601).** Withdrawn by its authors. Do not cite.
- **The uDGP survey pair** (Duxbury, Granlund, Gujarathi, Juhás & Billinge, *Discrete Applied Mathematics*; Billinge, Duxbury, Gonçalves, Lavor & Mucherino, *4OR*, DOI `10.1007/s10288-016-0314-2`). Both are real and relevant as background for the "unassigned distance geometry problem" name, but volume/pages were not confirmed for either, and neither is load-bearing for any T1 claim. Add them only if the paper needs the uDGP-survey framing, and verify first.
- **Nixon, Owen & Power, "Rigidity of frameworks supported on surfaces"** and the two DCG follow-ups on generic global rigidity on surfaces. The *program* is cited in the positioning document as the nearest structural analogue, but the individual author lists on the two DCG papers could not be confirmed, so no BibTeX entry was created. Verify before citing any of them individually.
- **Compagnoni & Notari, "TDOA-based localization in two dimensions: the bifurcation curve"** (*Fundamenta Informaticae*, DOI `10.3233/FI-2014-1118`, arXiv:1402.1530). Real, and directly on the `R = d+1` two-solution locus, but the arXiv version carries template-placeholder volume/pages. Verify against the published version if the bifurcation curve is cited.

## 4. Known verification gaps — read before claiming priority

Two areas were searched but **not exhausted**, and both bear directly on
novelty:

1. **The TOA/TDOA sensor self-calibration corpus** (Lund group: Kuang,
   Åström, Oskarsson, Burgess; plus Ferranti et al., El Badawy et al., Cao
   et al.). This literature solves exactly the "unknown receivers, unknown
   sources" problem, but publishes *minimal solvers for small fixed
   `(d, #receivers, #senders)`* rather than general-`d` identifiability
   theorems. It is large, partly conference-only, patent-adjacent and poorly
   indexed. Estimated residual risk that a general-`d` count matching
   `(R = d+2, n = d^2+3d+1)` already exists there: **~40%.**
   *Recommended action: hand-check this corpus before any priority claim.*

2. **The causal-set literature on chain-length fluctuations.**
   `brightwellluczak2015` was full-text searched and contains zero
   occurrences of "Tracy-Widom", "KPZ", "Kardar", "last-passage" or
   "Hammersley", which is good evidence for the naming gap — but it is a
   2015 survey. A 2015-2026 sweep of causal-set chain-statistics papers was
   started, not completed.

## 5. Notes for a submission pass

- Journal names are written in full; abbreviate per the target venue.
- `armstrong2021`'s `year` field is the preprint year and must be corrected
  to the journal year at submission.
- `aghili2018`'s title string must be re-read off the article page.
- Overlap with `docs/paper/paper_b/citations/references.bib` is intentional
  and keys are kept identical, so the two bibliographies can be merged
  without collision if the theory work is folded back into Paper B.
