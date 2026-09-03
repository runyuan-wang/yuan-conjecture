# Paper plan — weighted-span elastic rectangle frames

## Working identity

- **Working title:** *Weighted-Span Obstructions for Sparse Elastic Rectangle Frames on Velocity Grids*
- **Author:** Runyuan Wang
- **Document status:** local draft for parent review; not submitted or published
- **Central thesis:** On the cubic velocity grids `V_m={-m,…,m}^3`, the canonical degree-normalized elastic-rectangle operator obeys a weighted-span test. Uniformly bounded degree plus sublinear side lengths gives `gamma_m=O_D((L_m/m)^2)`, so a scale-independent gap requires quantitatively macroscopic rectangles on a positive fraction of rows.

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

Assume every vertex has positive degree. If `K_m=span{1,x,y,z,|v|^2}` and

`gamma_m=min{||A_m x||_2 : ||x||_2=1, x ⟂ D_m^{1/2}K_m}`,

then

`gamma_m^2 <= G_m^-1 sum_R (p_x q_x-p_y q_y)^2 <= G_m^-1 sum_R |p_R|^2 |q_R|^2`,

where

`G_m=sum_{v∈V_m}(x^2-y^2)^2=2S_0(S_0S_4-S_2^2)`,

`S_j=sum_{t=-m}^m t^j`. Consequently, if `1<=d_v<=D` and `|p_R|,|q_R|<=L_m`,

`gamma_m <= O_D((L_m/m)^2)`.

**Permitted corollary:** any family with uniform `D` and `L_m=o(m)` has gap tending to zero; this excludes fixed-size local, fixed-block-glued and all sublinear-span constructions satisfying these hypotheses from being uniform expanders under this normalization.

**Macroscopic necessary condition:** if `gamma_m>=c>0` with constants `c,D`, then at least `eta=min(1,c^2/(540D))` of all rectangles have both side lengths at least `rho*m`, where `rho=min(1,c/sqrt(135D))`. This is necessary, not sufficient.

## Proof architecture

1. Parameterize an elastic rectangle as `a, a+p, a+q, a+p+q` with `p·q=0` and define its second-difference row.
2. Verify that each row annihilates the five collision modes.
3. Choose `g(x,y,z)=x^2-y^2`; by cubic symmetry, `g` is unweighted-orthogonal to all five modes.
4. Project `g` onto `K_m` in the degree-weighted norm and use `f=g-h`. Then `D_m^{1/2}f` lies in the quotient space and every rectangle has `Δf=Δg`.
5. Compute `Δg=2(p_xq_x-p_yq_y)`, so the `1/2` normalization cancels the factor `2`.
6. Use `d_v>=1` and the unweighted orthogonality of `g` to obtain `||f||_{d}^2>=G_m` for the denominator.
7. Apply Cauchy-Schwarz for the weighted product bound and then count rows using `4|R_m|=sum_v d_v<=D|V_m|`.
8. Derive the sublinear-span and macroscopic-span corollaries.

## Section logic

1. **Introduction and scope.** Explain the independent rectangle-frame problem and state the evidence taxonomy before any result.
2. **Elastic rectangle operators.** Give the rectangle parameterization, signed matrix, degree normalization, invariant space and quotient gap.
3. **Weighted-span obstruction.** State and prove the proposition, sublinear-span corollary, and macroscopic necessary condition.
4. **Canonical finite model and certificates.** Describe deterministic pool deduplication and the exact rational-kernel certificate logic, with no claim that exact kernel controls the gap.
5. **Finite full-pool search evidence.** Report the four residual-greedy matrices `m=2,3,4,5`. Separate exact integer/rational columns from floating gaps.
6. **Validation and reproducibility.** Identify source scripts, fresh parent reruns, independent raw-JSON recomputation, theorem arithmetic checks, and manuscript-local checks.
7. **Open problem and limitations.** State the unrestricted nonlocal proposition and the missing asymptotic lower bound/proof. Record that no systematic novelty review or exhaustiveness claim is present.
8. **Conclusion.** Restate only the weighted-span theorem and the remaining macroscopic open problem.

## Contribution boundary

### Included

- A self-contained proof of the weighted-span normalized-gap upper bound `[R]`.
- The exact implication for uniformly local, fixed-block-glued, and sublinear-span families `[R]`.
- A quantitative macroscopic-span necessary condition `[R]`.
- Four exact finite rational-kernel statements and exact degree/count data `[E]`.
- Four finite normalized gap estimates, with independent raw-JSON recomputation `[N]`.
- A precise formulation of the unrestricted nonlocal proposition as open `[O]`.

### Explicitly excluded

- A claim about a relation to any other research route.
- A proof or disproof of the remaining macroscopic bounded-degree rectangle-expander proposition.
- An asymptotic lower bound inferred from `m=2,…,5`.
- Any implication from `ker_Q(C_m)=K_m` to a uniform spectral gap.
- Any global novelty, priority, current-best, literature-exhaustiveness, or attack-exhaustiveness claim.
- Any use of search snippets, unverified future publications, or the 2026 preprint as settled evidence.

## Citation strategy

The manuscript uses targeted literature-boundary checks for DVM/normal invariants and locality vs spectral gap. No novelty or priority claim is made. Before any submission, a separate systematic literature and priority review remains mandatory.

## Draft acceptance gates

- All five required artifacts exist and agree on the theorem/open boundary.
- Every strong claim maps to an exact local source or a verified bibliographic entry in `SOURCE_MATRIX.md`.
- Citations in `main.tex` match `references.bib` and the verified allowlist.
- A standard-library source/claim checker validates all four candidates, certificates, parent gap rerun and theorem rerun.
- LaTeX compiles locally if a TeX engine is available; logs and PDF checks are recorded.
- English and Chinese disclosures state exactly that all manuscript prose and structure were generated by LingTai AI under Runyuan Wang's direction; LingTai also performed source, mathematical/computational and code-validation assistance; the original conjectural/experimental direction is Runyuan's; evidence rather than AI status supports claims; independent human review is mandatory; and Runyuan retains final scientific judgment and responsibility.
- The final report contains exactly one authorized verdict.
