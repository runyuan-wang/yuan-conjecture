# Build and validation

## Frozen artifacts

- `main.pdf`: 9 letter-sized pages, 268,723 bytes, SHA-256 `9de879a158a7a99bc9df7bfbb1a8541e559390c3738087647d13c2f1903626c7`
- `main.tex`: SHA-256 `04a541e03a64407bbe574c540643752be0606a373fdec408f31d8e08b42c731d`
- `validation_checks.py`: SHA-256 `b28b00c214425c730250b09813981d651d479f92b316b598c97bcd89253d28f0`
- `assets/collision_rectangle_expander_assets.zip`: 2,073,624 bytes, 200 regular files, SHA-256 `e0b045bb43f5c37ea4a52f27644dcb8b964db25b4782ced47cc4cba5a5ced4a9`

The experiment ZIP was rebuilt after repairing the `basic_statistics.min_incidence_degree` generator bug in `rectangle_core.py`; row arrays and candidate SHA-256 identifiers are unchanged. The six JSON files under `parent-reruns/` are the exact certificate, gap and bounded-span artifacts consumed by the manuscript validator.

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

Current revision result: **PASS** when run against the rebuilt `.repro` source package; the four repaired stored minima are `4,3,3,3` and all zero-degree counts are `0`.

The validator checks every designated candidate row, canonical rectangle geometry, five exact collision-mode residuals, domain-separated hashes, degree ranges, stored degree metadata against recomputed incidence degrees, four exact rational-kernel certificates and listed odd minors over `F_2`, two large-prime screens, finite gap values/residuals, weighted-span arithmetic, citations, synchronized manuscript values, claim-boundary language, structure and path/privacy hygiene.

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
- no unresolved references in the final log; underfull box warnings remain in the compact literature-boundary table and are cosmetic;
- PDF is unencrypted, has no JavaScript or reported suspects;
- 20 fonts are embedded, subsetted and Unicode mapped;
- `main.pdf` was rebuilt locally; no visual review beyond successful compilation was performed in this revision pass.

## Independent review and parent acceptance

Earlier parent/independent review records applied to the pre-revision draft. This revision supersedes their validation counts with the current fresh-unpack validator result: `PASS 26823 deterministic assertions; 4 candidate notes`. No new independent visual review was performed in this narrow documentation-consistency repair.

## Limits

This package does not prove the unrestricted nonlocal proposition, interval-certify the four numerical gaps, or claim novelty/priority/current-best/exhaustiveness. A systematic subject-specific literature/priority review and independent human mathematical review are still required before formal submission.
