# Global First-Order Compatibility Obstructions for Finite-Radius Contact Wiring on Rational Worldline Networks

**Author:** Runyuan Wang（王润圆，昆明医科大学硕士）

**Status:** author-reviewed computational/mathematical research record; not a formal submission

- [Paper (PDF)](main.pdf)
- [LaTeX source](main.tex)
- [中文摘要与边界](SUMMARY.zh-CN.md)
- [Claim-to-evidence matrix](SOURCE_MATRIX.md)
- [Build and validation](BUILD_AND_VALIDATION.md)
- [SHA-256 manifest](SHA256SUMS.txt)
- [Legal-collision-movie companion package](assets/legal_collision_movie_2026-08-29/)

## Main result

This paper studies a finite first-order compatibility problem for exact rational worldlines. Pairwise concurrence events are interpreted locally as equal-mass central collisions, and finite-radius contact corrections are compiled into a global linear matching system.

The audited planar `K6` weave is wireable despite a 20-dimensional left nullspace. In contrast, the rational ruled `K4,4` network and a determinant-one shear of it have stable nonzero compatibility residuals. Exhaustive finite subset search retains at most `12/16` contacts without a degree cap or with cap 4, and `8/16` with cap 2. The independent verifier reproduces the ranks, feasibility classifications and finite subset optima with a different variable ordering.

The `12/16` result is a **finite first-order contact-wiring skeleton**. It is not a complete positive-radius hard-sphere trajectory.

## Reproducible package

The companion directory preserves the supplied compilation audit, compiler, certificate JSON, independent verifier, independent result JSON, project progress record and original SHA-256 manifest. The compiler uses exact rational incidence geometry plus float64 SVD/least squares for the contact-wiring system. The independent verifier does not import the compiler and uses a different variable ordering.

Run `python3 validate_publication.py` for static integrity and claim checks. Add `--run` to execute both numerical programs when NumPy is available. See [BUILD_AND_VALIDATION.md](BUILD_AND_VALIDATION.md).

## Scientific boundaries

This repository does **not** claim:

- a theorem about arbitrary hard-sphere dynamics;
- a proof or disproof of the four-dimensional Kakeya conjecture;
- a Kakeya dimension bound or entropy exponent;
- a universal sparse obstruction;
- a proof bridge from Paper 1 or Paper 2;
- novelty, priority, current-best status or an exhaustive search.

The passage from a finite first-order wiring obstruction to a genuine positive-radius collision movie—and from there to any broader geometric statement—remains open.

## AI-production disclosure

The manuscript prose and structure were generated with **LingTai AI**. Runyuan Wang supplied the conjectural agenda and experimental direction, reviewed the final research record, and retains final scientific judgment and responsibility. LingTai AI assisted with source organization, mathematical/computational checking, code validation and document build. Claims must stand on the displayed derivations and reproducible evidence rather than on AI authorship. Independent human mathematical and literature/priority review remains necessary before formal submission.
