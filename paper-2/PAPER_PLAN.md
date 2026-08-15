# Paper plan — bounded-span elastic rectangle frames

## Working identity

- **Working title:** *Bounded-Span Obstructions for Sparse Elastic Rectangle Frames on Velocity Grids*
- **Author:** Runyuan Wang
- **Document status:** local draft for parent review; not submitted or published
- **Central thesis:** On the cubic velocity grids `V_m={-m,…,m}^3`, the canonical degree-normalized elastic-rectangle operator cannot have a scale-independent quotient gap when both incidence degree and Euclidean side length are uniformly bounded. The proof gives an explicit upper bound of order `m^-2`.

## Evidence classes used throughout

| Label | Meaning | Permitted conclusion |
|---|---|---|
| `[R]` | Rigorous, self-contained derivation | May be stated as a theorem within its explicit hypotheses |
| `[E]` | Exact finite certificate over integers/rationals | May establish exact facts only for the listed finite matrices |
| `[N]` | Finite floating-point computation | May report reproducible values and residuals, but no asymptotic inference |
| `[O]` | Open proposition | Must be stated as unresolved |

The manuscript will not let `[E]` or `[N]` imply an infinite-family theorem.

## Precise theorem boundary `[R]`

Let `V_m={-m,…,m}^3`, and let `R_m` be elastic rectangles whose four vertices lie in `V_m`. Let `C_m` be the signed rectangle matrix, `d_v` the incidence degree, `D_m=diag(d_v)`, and

`A_m=(1/2) C_m D_m^{-1/2}`.

Assume every vertex has positive degree, `1 <= d_v <= D`, and every rectangle can be written with orthogonal side vectors `p,q` satisfying `|p|,|q| <= L`. If `K_m=span{1,x,y,z,|v|^2}` and

`gamma_m=min{||A_m x||_2 : ||x||_2=1, x ⟂ D_m^{1/2}K_m}`,

then

`gamma_m <= (L^2/2) sqrt(D(2m+1)^3/G_m)`,

where

`G_m=sum_{v∈V_m}(x^2-y^2)^2=2S_0(S_0S_4-S_2^2)`,

`S_j=sum_{t=-m}^m t^j`. Consequently,

`gamma_m <= L^2 sqrt(45D/32) m^-2(1+o(1))`.

**Permitted corollary:** any family with uniform `D` and `L` has gap tending to zero; this excludes fixed-size local/block-glued constructions satisfying these hypotheses from being uniform expanders under this normalization.

**Not permitted:** any conclusion about bounded-degree constructions whose rectangle spans grow with `m`.

## Proof architecture

1. Parameterize an elastic rectangle as `a, a+p, a+q, a+p+q` with `p·q=0` and define its second-difference row.
2. Verify that each row annihilates the five collision modes.
3. Choose `g(x,y,z)=x^2-y^2`; by cubic symmetry, `g` is unweighted-orthogonal to all five modes.
4. Project `g` onto `K_m` in the degree-weighted norm and use `f=g-h`. Then `D_m^{1/2}f` lies in the quotient space and every rectangle has `Δf=Δg`.
5. Compute `Δg=2(p_xq_x-p_yq_y)` and bound `|Δg|<=2L^2`.
6. Use `4|R_m|=sum_v d_v<=D|V_m|` for the numerator.
7. Use `d_v>=1` and the unweighted orthogonality of `g` to obtain `||f||_{d}^2>=G_m` for the denominator.
8. Evaluate `G_m` exactly and take its asymptotics.

## Section logic

1. **Introduction and scope.** Explain the project motivation and state the evidence taxonomy before any result. Cite only already verified Kakeya context; explicitly deny a Kakeya transfer.
2. **Elastic rectangle operators.** Give the rectangle parameterization, signed matrix, degree normalization, invariant space and quotient gap.
3. **Bounded-span obstruction.** State and prove the explicit theorem and local-family corollary.
4. **Canonical finite model and certificates.** Describe deterministic pool deduplication and the exact rational-kernel certificate logic, with no claim that exact kernel controls the gap.
5. **Finite full-pool search evidence.** Report the four residual-greedy matrices `m=2,3,4,5`. Separate exact integer/rational columns from floating gaps.
6. **Validation and reproducibility.** Identify source scripts, fresh parent reruns, independent raw-JSON recomputation, theorem arithmetic checks, and manuscript-local checks.
7. **Open problem and limitations.** State the unrestricted nonlocal proposition and the missing asymptotic lower bound/proof. Record that no systematic novelty review or exhaustiveness claim is present.
8. **Conclusion.** Restate only the bounded theorem and the remaining open problem.

## Contribution boundary

### Included

- A self-contained proof of the bounded-degree, bounded-side-length normalized-gap upper bound `[R]`.
- The exact implication for uniformly local or fixed-block-glued families `[R]`.
- Four exact finite rational-kernel statements and exact degree/count data `[E]`.
- Four finite normalized gap estimates, with independent raw-JSON recomputation `[N]`.
- A precise formulation of the unrestricted nonlocal proposition as open `[O]`.

### Explicitly excluded

- A Kakeya theorem, dimension estimate, counterexample extraction, or reduction.
- A proof or disproof of the unrestricted bounded-degree rectangle-expander proposition.
- An asymptotic lower bound inferred from `m=2,…,5`.
- Any implication from `ker_Q(C_m)=K_m` to a uniform spectral gap.
- Any global novelty, priority, current-best, literature-exhaustiveness, or attack-exhaustiveness claim.
- Any use of search snippets, unverified future publications, or the 2026 preprint as settled evidence.

## Citation strategy

Only citations already verified in the first-paper package are eligible. The manuscript uses Katz–Tao (1999, 2002) solely for broad Kakeya-project context and does not use them to support the rectangle theorem. No additional public citation is introduced. The mathematical result itself is proved in full. Before any submission, a separate systematic literature and priority review remains mandatory.

## Draft acceptance gates

- All five required artifacts exist and agree on the theorem/open boundary.
- Every strong claim maps to an exact local source or a verified bibliographic entry in `SOURCE_MATRIX.md`.
- Citations in `main.tex` match `references.bib` and the verified allowlist.
- A standard-library source/claim checker validates all four candidates, certificates, parent gap rerun and theorem rerun.
- LaTeX compiles locally if a TeX engine is available; logs and PDF checks are recorded.
- English and Chinese disclosures state exactly that all manuscript prose and structure were generated by LingTai AI under Runyuan Wang's direction; LingTai also performed source, mathematical/computational and code-validation assistance; the original conjectural/experimental direction is Runyuan's; evidence rather than AI status supports claims; independent human review is mandatory; and Runyuan retains final scientific judgment and responsibility.
- The final report contains exactly one authorized verdict.
