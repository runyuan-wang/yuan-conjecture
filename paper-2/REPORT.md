# Second-paper draft report

## Status

A complete local second-paper draft package has been prepared for parent review. The manuscript is self-contained at the mathematical level: it defines the elastic rectangle operator and canonical normalization, proves the explicit bounded-span theorem, separates exact finite certificates from floating-point gap evidence, formulates the remaining nonlocal problem as open, and states the Kakeya/non-novelty boundaries repeatedly and consistently.

No source package, repository, public record, runtime configuration, credential, or user data was modified. Nothing was committed, pushed, published, submitted, released, or communicated externally.

## Required artifacts complete

- `PAPER_PLAN.md`: thesis, evidence taxonomy, section logic, contribution boundary, forbidden inferences, citation strategy, and acceptance gates.
- `SOURCE_MATRIX.md`: strong claim → exact local source or already verified public citation, with source tier and disallowed inference.
- `main.tex`: complete scholarly manuscript, full proof, finite evidence table, validation description, limitations, open problem, and tool-assistance disclosure.
- `SUMMARY.zh-CN.md`: plain-Chinese theorem explanation, finite-evidence split, open boundary, and next step.
- `REPORT.md`: this report.

Additional local artifacts include `references.bib`, `validation_checks.py`, `main.pdf`, LaTeX/PDF/lint logs, extracted PDF text, and a fresh `theorem_validation_rerun.json`.

## Source review completed

The drafting pass read the substantive first-paper/publication manuscript, claim/evidence and internal-evidence matrices, verified-literature matrix, build/validation record, publication status/receipt, and Chinese summary under `../yuan-conjecture-paper-20260813/`. It then read the current repository manuscript/bibliography and complete update report, validation, builder brief, patch evidence, receipt and final validation under `../yuan-conjecture-latest-experiment-update-20260814/`.

For the accepted sparse elastic rectangle package, the pass read the parent final review and final progress, report, validation, mathematical theorem paper, README, literature boundary, input receipt, failure log, mechanical preflight, complete core/generator/search/exact-validator/theorem-validator scripts, independent raw-gap script/result, bounded-span rerun, attack accounting, certificate comparison, all four complete fresh exact certificates (including full minor index payloads), and the raw candidate/result fields needed for every manuscript claim. `validation_checks.py` independently traverses every row of all four raw candidate arrays and rechecks their hashes, geometry, five exact mode residuals, degrees, and listed odd minors.

Only citations already verified in the first-paper package are used: Katz–Tao 1999 and Katz–Tao 2002, solely for broad Kakeya-project context. No search result or snippet was used as scholarly evidence, and no additional citation was introduced.

## Scientific boundary preserved

The only infinite-family rigorous conclusion is:

`gamma_m <= (L^2/2) sqrt(D(2m+1)^3/G_m) = L^2 sqrt(45D/32) m^-2(1+o(1))`

under positive degrees `1<=d_v<=D` and uniformly bounded Euclidean side lengths `|p|,|q|<=L`. This excludes fixed-size local/fixed-block-glued uniform expanders under the canonical normalization.

The draft does **not** infer anything from that theorem about nonlocal families whose spans grow with `m`. The unrestricted nonlocal proposition is explicitly open. The four rational-kernel statements and degree/count data are exact but finite; the four normalized gaps are floating-point and finite. Exact kernel dimension is explicitly not treated as a uniform gap. There is no Kakeya extraction, theorem, proof/disproof, dimension bound, global novelty, current-best, priority, literature-exhaustiveness, search-exhaustiveness, or attack-exhaustiveness claim.

## Validation performed

All commands were run from this draft directory with existing local tools only; no dependency was installed.

### 1. Fresh bounded-span validation

Command:

```text
PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="../yuan-conjecture-collision-rectangle-expander-20260814/.deps:../yuan-conjecture-collision-rectangle-expander-20260814/collision_rectangle_expander" python3 ../yuan-conjecture-collision-rectangle-expander-20260814/collision_rectangle_expander/validate_bounded_span_theorem.py --results ../yuan-conjecture-collision-rectangle-expander-20260814/collision_rectangle_expander/results --output theorem_validation_rerun.json
```

Result: PASS. Exact moments and five orthogonality identities passed for `m=1,…,10`; the theorem inequality chain passed on four block candidates. Standard error was empty. The fresh JSON has `all_checks_passed=true`.

### 2. Deterministic source and claim validator

Command:

```text
PYTHONDONTWRITEBYTECODE=1 python3 validation_checks.py --source ../yuan-conjecture-collision-rectangle-expander-20260814 --root .
```

Checks include:

- every one of the 3,794 designated candidate rows;
- canonical row geometry, orthogonal sides, five integer mode residuals, uniqueness, and domain-separated SHA-256;
- exact degree ranges and zero-isolated-vertex condition;
- all four fresh certificates and independent reconstruction of each listed odd minor over `F_2`;
- rational rank/kernel fields and two large-prime screens;
- independent gap rerun values, warnings and residual thresholds;
- exact moment/orthogonality arithmetic for `m=1,…,10` and the parent theorem rerun;
- artifact presence, synchronization of finite values, citation allowlist/bibliography equality, brace/document structure, private-path leakage, claim-boundary language, conflict markers, and trailing whitespace.

Result: PASS: `26,803` deterministic assertions and four candidate notes; standard error is empty. The exact record is `validation_checks.stdout.txt`.

### 3. LaTeX build

Command:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error main.tex
```

Result: PASS. BibTeX completed, citations resolved, and `main.pdf` was produced. `latexmk.stderr.txt` retains the expected first-pass notice that `main.bbl` did not yet exist; LaTeXmk then generated it and completed normally. A scan of the final `main.log` found no LaTeX/package warnings, undefined references, overfull boxes, or underfull boxes. The PDF has 8 letter-sized pages, is unencrypted, has no JavaScript, and `pdfinfo` reports no suspects.

### 4. PDF readback and fonts

Commands:

```text
pdftotext main.pdf main.pdftotext.txt
pdfinfo main.pdf
pdffonts main.pdf
```

Result: PASS. Extracted text is nonempty and contains the title, unrestricted-nonlocal boundary, and bibliography text. Every listed font is embedded, subsetted and Unicode mapped.

### 5. ChkTeX

The installed `chktex` executable ran but exited `2` because its bundled regular-expression configuration is incompatible with the host regex engine (`repetition-operator operand invalid`). It nevertheless emitted only style suggestions (mainly math grouping and nonbreaking-space advice), not a compile or semantic error. This tool-specific failure is preserved in `chktex.stderr.txt` and `chktex.stdout.txt`; it was not concealed or treated as a pass. The authoritative LaTeX build and log scan passed.

## Unresolved gaps

### Citation and novelty

- No systematic literature/priority review specific to sparse elastic rectangle frames or associated Poincaré inequalities has been completed.
- The draft intentionally makes no novelty or priority claim.
- Before submission, an updated official-source citation audit and a systematic subject-specific literature review are required.

### Technical

- The unrestricted nonlocal bounded-degree proposition remains open.
- No recursive nonlocal family with a uniform lower gap is constructed.
- No obstruction is proved when side lengths grow with `m`.
- The four spectral gaps are not interval-certified or exact algebraic values.
- No Kakeya extraction, regularization, bounded-congestion lemma, or geometric transfer is supplied.
- The fresh finite data end at `m=5`; increasing scale alone would not resolve the structural question.

### Review

- This is an AI-assisted local draft and has not received independent expert mathematical or editorial review.
- The theorem proof should be checked line by line by the parent and then by an independent specialist before any external use.

## Exact recommended next step

Do **not** prioritize another scale-only `m` run. The next mathematical task should target the unrestricted nonlocal gap directly: either (a) give an explicit recursive bounded-degree rectangle family with side lengths allowed to grow and prove a scale-free Poincaré lower bound after quotienting the five collision modes, or (b) construct a new test function/decomposition that extends the upper-bound obstruction to a clearly stated nonlocal subclass. In parallel, before submission, conduct a systematic official-source literature/priority review and obtain independent mathematical review of Theorem 3.1 and its normalization.

## Draft disposition

All required manuscript artifacts are present, the accepted scientific boundary is preserved, the deterministic source checks and local LaTeX/PDF checks pass, and the unresolved research/novelty gaps are disclosed rather than papered over. The package is ready for parent review as a draft, not for publication or submission.

SECOND_PAPER_DRAFT_READY_FOR_PARENT_REVIEW
