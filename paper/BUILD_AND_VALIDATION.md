# Build and validation record

**Date:** 2026-08-13

## 1. Scope and evidence boundary

This record covers the final corrected bounded manuscript build. The first builder's `BLOCKED` outcome was caused only by a split-root dispatch contract; it was not a scientific or toolchain failure.

The manuscript is an honest computational-exploration report. The source package was not independently rerun. Static source-to-claim correspondence, matching manifests and hashes, a clean manuscript build, and visual review do **not** prove or disprove four-dimensional Kakeya, establish a new dimension bound, establish novelty or current-best status, or constitute an independent reproduction.

During preparation, source assets and the target repository remained read-only. No package was installed, no research computation was rerun, and no journal or arXiv submission occurred.

## 2. Prerequisite evidence gate

Before prose, the source-first review read the two frozen matrices, all five literature-packet files, the scientific assessment, all 39 source-package files, and then the manuscript package. The frozen matrix hashes were:

```text
92959c933457c120cd9b5bd25bf41cdae4ff84737b3efa1927cdfc44ec05c683  CLAIM_EVIDENCE_MATRIX.md
35dd8c394a1677954d517e57507cb9a8046f62ef1f2b70aa4db6cd06e39e6fe6  INTERNAL_EVIDENCE_MANIFEST.md
```

The retained literature matrix supplied every external citation key. The internal matrix supplied each quantitative claim and its exact artifact. The independent source-first manuscript review concluded `PASS_WITH_MINOR_NOTES`: no material overclaim, unsupported number, citation-class error, privacy leak, broken intended-layout link, or visual blocker was found.

The review's sole correction was applied across the TeX, Chinese summary, HTML, and validator: the deterministic open-path interface screen sampled, per reversal-inequivalent backbone, 24/160 pair interfaces, 40/640 triple interfaces, and 16/64 structured-full interfaces; reversal-symmetric backbones used smaller antiparallel samples of 16/160, 24/640, and 8/64. These are bounded sampled screens, not exhaustive no-go theorems.

## 3. Toolchain and clean build

Installed tools used without installation:

- Latexmk 4.85;
- pdfTeX 3.141592653-2.6-1.40.26, TeX Live 2024;
- BibTeX 0.99d;
- Poppler `pdfinfo`, `pdffonts`, `pdftotext`, and `pdftoppm` 26.04.0;
- Python 3 and its standard HTML parser;
- local image-review utilities for raster inspection.

Final clean rebuild command:

```sh
latexmk -gg -pdf -interaction=nonstopmode \
  -halt-on-error -file-line-error main.tex > build-clean.log 2>&1
```

Result: exit 0. The final PDF is 10 pages and 305,493 bytes. `main.log` contained:

```text
overfull=0
underfull=0
hyperref_warnings=0
unresolved citations/references or rerun markers=0
```

`pdfinfo` reported a US Letter PDF 1.5 with `Suspects: no`, `JavaScript: no`, and `Encrypted: no`. `pdffonts` listed 17 font rows; every font was embedded and Unicode-mapped. Layout-preserving extracted text contained 458 lines, 3,372 words, and 28,520 bytes, with no private path, worker ID, `??`, or `undefined` placeholder.

## 4. Deterministic editorial and correspondence suite

Command shape (the source argument points to the read-only unpacked public source package):

```sh
python3 validation_checks.py --project <private-build-project-root> \
  --source <read-only-unpacked-source-root>
```

Result: exit 0 and `VALIDATION_CHECKS_PASS`.

```text
PASS package-presence: manuscript, bibliography, Chinese summary, HTML, PDF, and final log
PASS frozen-matrix-hashes: 2/2 match supplied SHA-256 values
PASS decisive-source-hashes: 6/6 internal-manifest rows match
PASS source-payload-hashes: 38/38 manifest entries match
PASS citations: 6 cited keys = 6 bibliography keys; every key occurs in a retained verified row
PASS final-reference-log: no unresolved citation/reference or rerun marker
PASS public-privacy-and-omission: no private path/worker ID; Boltzmann--Grad layer absent
PASS claim lint: no affirmative proof/disproof/new-bound/novelty/current-best/exhaustiveness/reproduction claim; mandatory boundaries present
PASS screen-regression-correspondence: named source certificates match all reported C4/N counts, widths, ranks, versions, and open-path sampling fractions
PASS exact-correspondence: C1/C2/C3 values anchored in named frozen artifacts; all public asset links allowlisted
PASS cross-format claim parity: key exact counts/margins/rank and PSITIP/GLOP limits occur in TeX, Chinese summary, and HTML
PASS HTML-static: 16 IDs, 36 links, 14 mapped local assets, zero external resources/broken fragments
PASS draft-text-diff-check: no trailing whitespace, conflict marker, or missing final newline
VALIDATION_CHECKS_PASS
```

The script performs editorial and static correspondence checks only. It does **not** execute the research searches or certificate checkers. Its principal assertions are:

- both frozen matrix hashes match;
- all six decisive-source hashes and all 38 payload hashes match;
- all six citations are present in the retained verified-literature matrix;
- the final LaTeX log has no unresolved references;
- the C1/C2/C3 exact values and all reported C4/N counts, beam widths, rank ratios, deficits, solver versions, and sampled-interface fractions agree with named frozen source artifacts;
- every public asset link maps to one of 14 allowlisted repository assets;
- the public text contains no local absolute path, private framework path, or worker ID;
- mandatory non-claims are present and unsupported proof/disproof/new-bound/novelty/current-best/exhaustiveness/reproduction language is absent;
- the unsupported Boltzmann--Grad layer is absent;
- key exact claims and limitations agree across TeX, Chinese summary, and HTML;
- public draft text has no trailing whitespace, merge marker, or missing final newline.

## 5. HTML and PDF visual review

`paper.html` is a self-contained HTML5 document. All CSS is inline; it imports no external script, stylesheet, font, image, iframe, or other rendering resource. Static parsing found 16 IDs, 36 links, 14 mapped local assets, zero external resources, duplicate IDs, or broken fragments.

A real browser render was inspected after the sampling correction. Every section, reference, footer, and the expanded disclosure row was visible. The only HTML note was blank tail below the document in the tall capture viewport; this is a screenshot-framing issue, not hidden content.

All 10 PDF pages were rasterized and independently inspected after the correction. There was no clipping, overlap, broken mathematics, missing glyph, blank page, or unreadable content. Minor cosmetic notes only:

- page 6: dense justified table text;
- page 9: stretched spacing around a long URL;
- page 10: journal-title hyphenation.

These notes do not change the scientific or evidence boundary and are not publication blockers.

## 6. Public artifact hashes

```text
af85daca490e96a1217d3bf17c6c6e9f743aafa4168bebd50a1ecceaee489b7e  main.tex
c06419021f9361eafb1735fd5d54ca2851b02dafb424dadb15d9d42edcffb1f5  references.bib
433929652fa5f488bf748992dabdd7dd915a4556bf57a4ea2c977f8625ee6794  SUMMARY.zh-CN.md
848a72706a52613c4a41c43031f96197415a7c81b4555695a8a5ef0d8b6c9eff  paper.html
473a714e71b3e01f5162c34365f5ff1399f40e7b53d9ef33fb5f892585165a10  main.pdf
cb449e4245b09a263e7b8f9ebce4fbc964918a702ba7eb994c21e054290a268b  validation_checks.py
```

## 7. Remaining limits

- The source package has not been independently rerun.
- PSITIP/GLOP remains floating-point output without an exact rational dual certificate.
- Beam-null and sampled-interface results remain bounded screens, not no-go theorems.
- The worldline soft-mode-to-grain and collision-compilation/global-tail interfaces remain unresolved.
- Literature novelty and current-best status were not established.
- The 2026 v1 preprint and any problem-status sentence require a fresh literature-status check before any future formal submission.
- This repository publication is a transparent research record, not journal/arXiv submission or peer-reviewed acceptance.

## 8. 2026-08-13 LingTai tool-use disclosure follow-up

Runyuan explicitly required the paper to state that LingTai was used to write it. This follow-up adds one matched disclosure to `main.tex`, `paper.html`, and `SUMMARY.zh-CN.md`, plus an exact cross-format validator gate. The wording states that LingTai AI assisted evidence organization, source-to-claim checking, drafting, consistency validation, and document build, while Runyuan Wang retains sole authorship, scientific judgment, and final responsibility.

The scientific source package and both frozen editorial matrices were not changed or rerun. The validator-only private literature view had been removed during post-publication cleanup; it was faithfully reconstructed from the already accepted six-row public bibliography and the exact manuscript claim contexts, adding no new citation or scientific claim.

Clean rebuild and validation results:

- `latexmk -gg -pdf -interaction=nonstopmode -halt-on-error -file-line-error main.tex`: exit 0;
- full `validation_checks.py --project ... --source ...`: `VALIDATION_CHECKS_PASS`;
- frozen matrices 2/2, decisive source hashes 6/6, source payload hashes 38/38, citations 6/6: PASS;
- LingTai disclosure parity across TeX, HTML, and Chinese summary: PASS;
- HTML static checks: 17 IDs, 38 links, 14 mapped local assets, zero external resources or broken fragments;
- final PDF: 10 pages, 306,172 bytes, no suspects, JavaScript, encryption, or form;
- final LaTeX log: zero overfull/underfull boxes, hyperref warnings, unresolved citations/references, or rerun markers;
- PDF fonts: 17/17 embedded, subset, and Unicode;
- extracted PDF text: 616 lines, 3,378 words, 23,308 bytes;
- visual review: PDF page 9 PASS with the disclosure integrated cleanly between Conclusion and References; full 1440×9000 HTML capture PASS with the disclosure visible, readable, and professionally placed before References. Page 10 retains the prior cosmetic reference spacing/hyphenation and blank-tail notes; the disclosure correctly appears on page 9.

Final follow-up hashes before commit:

```text
48af59572b4a11cc4660778f6969006e289154f620703d8438229f9c65c51864  main.pdf
a39908384b5d47723c62f1cc7a17c20a737f309687f911cc880c836b296e388a  main.tex
0362f0ff13639906a8d9fb4181942c904b0a06248c9f060dd06a45f4542f1798  paper.html
cce0e8f6e4de77f24a5f011367e43c54827d1f26e47356b76510bbfcf4539522  SUMMARY.zh-CN.md
73028b8c9d946d5124af3892f39c2ffdb7e82a990cb3f60a9de30cb21d8110cd  validation_checks.py
```

No journal/arXiv submission, release, Pages, issue, pull request, new experiment, literature-status update, or scientific-claim expansion is part of this follow-up.
