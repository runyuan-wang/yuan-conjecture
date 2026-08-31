# Build and validation

## Frozen author-reviewed artifacts

- `main.pdf`: 10 pages, 336,655 bytes, SHA-256 `7c2b63a6241acc3efb93781ebf522afa493bb68ba3953f51f55b4a73ee0cf081`
- `main.tex`: 30,971 bytes, SHA-256 `dc805326178f57bda5a37f436d73cf37ab2e06d8e2409c1bcd5b1dbac0f0bb91`
- `assets/legal_collision_movie_2026-08-29/`: the six original companion files plus their original SHA-256 manifest, preserved byte-for-byte

The manuscript PDF and LaTeX source are byte-identical to the parent-accepted author-review artifacts. The companion manifest checks the compilation audit, compiler source, certificate JSON, independent verifier source/result and master-progress record. The two preserved historical Markdown records intentionally retain their original double-space Markdown hard breaks; the publication diff check treats those byte-locked lines as a documented exception while requiring every newly authored file to pass whitespace checks.

## Static publication validator

From `paper-3/`:

```bash
python3 validate_publication.py
```

This checks:

- frozen PDF/LaTeX hashes;
- every file in the original companion SHA-256 manifest;
- the exact seven-family inventory;
- key ranks, left-null dimensions, feasibility classifications and exhaustive finite optima;
- required scientific-boundary language in the LaTeX manuscript;
- absence of local paths, agent run IDs, mailbox identifiers and internal review markers from the public tree.

The static validator uses only the Python standard library.

## Executable numerical rerun

The two programs require Python with NumPy. The accepted rerun used Python 3.13.14 and NumPy 2.5.2. Floating-point values may vary slightly across BLAS/LAPACK and NumPy builds, so the acceptance gate is categorical agreement on status, ranks, feasibility, finite subset optima and scaling conclusions—not byte equality of regenerated JSON.

Run:

```bash
python3 validate_publication.py --run
```

Equivalent direct commands are:

```bash
python3 assets/legal_collision_movie_2026-08-29/legal_collision_movie_compiler.py \
  --output compiler.local.json
python3 assets/legal_collision_movie_2026-08-29/verify_legal_collision_movie.py \
  --output verifier.local.json
```

The independent verifier does not import the compiler and uses a different variable ordering. Its frozen result reports `status: PASS`.

## Paper build

With a standard LaTeX distribution:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error main.tex
```

The accepted PDF is unencrypted and passed deterministic content checks plus visual review of representative pages 1, 6 and 10.

## Interpretation limits

A successful validator run reproduces the finite model calculations and publication integrity only. It does not prove a full positive-radius hard-sphere trajectory, a Kakeya theorem, an entropy exponent, a universal sparse obstruction, a Paper-1/Paper-2 bridge, novelty or priority.
