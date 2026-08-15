# Build and validation

## Frozen artifacts

- `main.pdf`: 8 letter-sized pages, 259,336 bytes, SHA-256 `9eda214d86520015c6017ce0bebccd41de38c687aa377f6dfbe2259600185eb9`
- `main.tex`: SHA-256 `911a21ff8138e336a0355c7e57cf70613fb9c197b5c3147f09e2ba56ed207403`
- `validation_checks.py`: SHA-256 `6ce9983736d499231cbd92c77534ed11c4d1e71e32779c08db4fd713ba47faed`
- `assets/collision_rectangle_expander_assets.zip`: 2,001,758 bytes, 200 regular files, SHA-256 `cd4b859a8dfb51ac07427748a9ba852e5b90c2b166e75c9c04999600d829bf2a`

The experiment ZIP is preserved byte-for-byte. The six JSON files under `parent-reruns/` are the exact certificate, gap and bounded-span artifacts consumed by the manuscript validator.

## Deterministic manuscript/source validator

From `paper-2/`:

```bash
rm -rf .repro
mkdir -p .repro/PARENT_RERUNS/certificates
unzip -q assets/collision_rectangle_expander_assets.zip -d .repro
cp parent-reruns/independent_gap_recompute.json .repro/PARENT_RERUNS/
cp parent-reruns/bounded_span_validation.json .repro/PARENT_RERUNS/
cp parent-reruns/certificates/*.json .repro/PARENT_RERUNS/certificates/
PYTHONDONTWRITEBYTECODE=1 python3 validation_checks.py --source .repro --root .
```

Accepted parent result: **PASS — 26,818 deterministic assertions and four candidate notes; empty standard error**.

The validator checks every designated candidate row, canonical rectangle geometry, five exact collision-mode residuals, domain-separated hashes, degree ranges, four exact rational-kernel certificates and listed odd minors over `F_2`, two large-prime screens, finite gap values/residuals, bounded-span arithmetic, citations, synchronized manuscript values, claim-boundary language, structure and path/privacy hygiene.

## Independent bounded-span arithmetic rerun

This optional rerun requires Python with NumPy available:

```bash
PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH=".repro/collision_rectangle_expander" \
python3 .repro/collision_rectangle_expander/validate_bounded_span_theorem.py \
  --results .repro/collision_rectangle_expander/results \
  --output bounded_span_validation.local.json
```

Accepted result: `all_checks_passed=true`; exact moment/orthogonality checks pass for `m=1,…,10`, and the theorem inequality chain passes on four block candidates.

## Paper build

With a standard LaTeX distribution:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error main.tex
```

Accepted build checks:

- BibTeX and all citations resolved;
- no unresolved references, overfull boxes or underfull boxes in the final log;
- PDF is unencrypted, has no JavaScript or reported suspects;
- 20 fonts are embedded, subsetted and Unicode mapped;
- all eight pages passed visual review.

## Independent review and parent acceptance

One corrected source-first independent Sol review read the proof, source scripts, raw evidence and certificate package before reading the parent review. It reported `INDEPENDENT_SOL_PASS` with no MATERIAL or MINOR manuscript findings on the original manuscript. The later disclosure-only amendment changed no proof, evidence or scientific claim; parent reran the 26,818-assertion validator, exact theorem arithmetic, PDF metadata/log/font gates and all-eight-page visual review on the amended build.

## Limits

This package does not prove the unrestricted nonlocal proposition, interval-certify the four numerical gaps, establish a Kakeya result, or claim novelty/priority/current-best/exhaustiveness. A systematic subject-specific literature/priority review and independent human mathematical review are still required before formal submission.
