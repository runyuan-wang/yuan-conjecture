#!/usr/bin/env python3
"""Search a probability-legal heterogeneous tertiary-fold family.

This is the next experiment after ``protein_tertiary_fold_search.py``.  The
earlier screen glued two copies of the *same* five-atom colored path.  Here the
first local motif is the strongest path near miss from that screen,

    X --X-- X --X-- X --Y-- X,

while the second motif ranges over every rooted, four-coloured five-vertex
star or fork, modulo rooted coloured-graph isomorphism.  Thus the two local
motifs have genuinely different topology, analogous only at a structural
level to different secondary motifs meeting at a common tertiary core.

Probability legality
--------------------
Each coloured tree is generated from a base atom U=(X,Y) by four successive
conditional copies over one of X, Y, X+Y, X+2Y.  Every vertex therefore still
has marginal law U.  Two differently distributed tree molecules can then be
coupled conditionally independently over either

* one allowed projection of a selected core atom (cost at most M), or
* the full core atom U, represented by the independent projections X and Y
  (cost at most H(U) <= H(X)+H(Y) <= 2M).

Consequently the collision cost is 9 or 10.  For the meaningful target
17/10=1.7, the encoding budget is respectively 8 or 7.  These constructions
do not require the full local motifs to have the same joint distribution.

Before search, a dimension gate gives a cheap necessary condition.  If an
encoding h has rank b and determines the target X-Y on k anchor atoms, adding
the full anchor atoms can increase rank by at most k further dimensions.
Since the quotient dimensions are 11 and 10, every candidate needs at least
three anchors.  This does not rule the family out, so the screen is warranted.

The first pass uses a narrow exact-rank beam over the entire finite family.
Only its strongest folds are rerun with a wider beam.  Any reported positive
candidate is checked directly by exact modular row-space arithmetic.  A zero
result is a deterministic screened-family negative, not a no-go theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations, product
from time import perf_counter
from typing import Iterable, Sequence

from eight_atom_cycle_search import Edge, Observation, P, SLOPES, quotient_model
from protein_tertiary_fold_search import (
    BeamState,
    constraint_signature,
    diverse_beam_search,
    exact_subspace_search,
    five_path_edges,
    format_observations,
)


N = 10
TARGET_TOTAL_COST = 17
COARSE_BEAM_WIDTH = 2
REFINE_BEAM_WIDTH = 32
REFINE_COUNT = 64
EXACT_AUDIT_COUNT = 3
SEED_BACKBONE = ("X", "X", "X", "Y")
LABELS = tuple(SLOPES)

TOPOLOGIES: dict[str, tuple[tuple[int, int], ...]] = {
    "star": ((0, 1), (0, 2), (0, 3), (0, 4)),
    "fork": ((0, 1), (1, 2), (1, 3), (3, 4)),
}


@dataclass(frozen=True)
class RootedMotif:
    topology: str
    root: int
    edges: tuple[Edge, ...]


@dataclass(frozen=True)
class HeterogeneousFold:
    seed_root: int
    motif: RootedMotif
    glue: tuple[str, ...]


@dataclass(frozen=True)
class FoldResult:
    fold: HeterogeneousFold
    collision_cost: int
    encoding_budget: int
    quotient_dimension: int
    encoding: tuple[Observation, ...]
    anchors: tuple[int, ...]
    reconstruction_rank: int
    feasible: bool
    beam_states: int


def canonical_rooted_coloured_graph(
    edges: Sequence[tuple[int, int]],
    colours: Sequence[str],
    root: int,
) -> tuple[int, tuple[Edge, ...]]:
    """Canonicalize a rooted coloured five-vertex tree over all relabellings."""
    best: tuple[int, tuple[Edge, ...]] | None = None
    for relabel in permutations(range(5)):
        mapped: list[Edge] = []
        for (left, right), colour in zip(edges, colours):
            a, b = sorted((relabel[left], relabel[right]))
            mapped.append((a, b, colour))
        candidate = (relabel[root], tuple(sorted(mapped)))
        if best is None or candidate < best:
            best = candidate
    assert best is not None
    return best


def rooted_motifs() -> tuple[RootedMotif, ...]:
    motifs: list[RootedMotif] = []
    for topology, edges in TOPOLOGIES.items():
        canonical = {
            canonical_rooted_coloured_graph(edges, colours, root)
            for colours in product(LABELS, repeat=4)
            for root in range(5)
        }
        for root, coloured_edges in sorted(canonical):
            motifs.append(RootedMotif(topology, root, coloured_edges))
    # These counts are useful guards against accidental search-space changes.
    counts = {
        topology: sum(motif.topology == topology for motif in motifs)
        for topology in TOPOLOGIES
    }
    assert counts == {"star": 115, "fork": 736}
    return tuple(motifs)


def fold_edges(fold: HeterogeneousFold) -> list[Edge]:
    edges = five_path_edges(SEED_BACKBONE, 0)
    edges.extend(
        (left + 5, right + 5, colour)
        for left, right, colour in fold.motif.edges
    )
    edges.extend(
        (fold.seed_root, fold.motif.root + 5, colour)
        for colour in fold.glue
    )
    return edges


def generate_folds() -> Iterable[HeterogeneousFold]:
    # Any two independent allowed projections identify the whole core atom.
    # X,Y is the canonical representative of that rank-two rowspace.
    glues = tuple((label,) for label in LABELS) + (("X", "Y"),)
    for motif in rooted_motifs():
        for seed_root in range(5):
            for glue in glues:
                yield HeterogeneousFold(seed_root, motif, glue)


def evaluate_fold(fold: HeterogeneousFold, width: int) -> FoldResult:
    glue_rank = len(fold.glue)
    collision_cost = 8 + glue_rank
    encoding_budget = TARGET_TOTAL_COST - collision_cost
    model = quotient_model(N, fold_edges(fold))
    assert model.quotient_dimension == 12 - glue_rank
    best, visited = diverse_beam_search(model, encoding_budget, width=width)
    encoding = tuple(model.observation_labels[index] for index in best.indices)
    return FoldResult(
        fold=fold,
        collision_cost=collision_cost,
        encoding_budget=encoding_budget,
        quotient_dimension=model.quotient_dimension,
        encoding=encoding,
        anchors=best.check.anchors,
        reconstruction_rank=best.check.reconstruction_rank,
        feasible=best.check.feasible,
        beam_states=visited,
    )


def score(result: FoldResult) -> tuple[int, int, int, int]:
    return (
        int(result.feasible),
        result.reconstruction_rank - result.quotient_dimension,
        len(result.anchors),
        -result.collision_cost,
    )


def print_result(prefix: str, result: FoldResult) -> None:
    print(prefix)
    print(
        f"  motif/root pair = {result.fold.motif.topology}/"
        f"{result.fold.seed_root}->{result.fold.motif.root}"
    )
    print(f"  motif edges = {result.fold.motif.edges}")
    print(f"  core glue = {result.fold.glue}")
    print(
        f"  collision cost / encoding budget = "
        f"{result.collision_cost}/{result.encoding_budget}"
    )
    print(f"  encoding = {format_observations(result.encoding)}")
    print(f"  anchors = {result.anchors}")
    print(
        f"  reconstruction rank = {result.reconstruction_rank}/"
        f"{result.quotient_dimension}; feasible={result.feasible}"
    )


def coarse_screen() -> tuple[list[FoldResult], int]:
    results: list[FoldResult] = []
    seen: set[tuple[int, tuple]] = set()
    beam_states = 0
    started = perf_counter()
    for fold in generate_folds():
        glue_rank = len(fold.glue)
        signature = (glue_rank, constraint_signature(fold_edges(fold)))
        if signature in seen:
            continue
        seen.add(signature)
        result = evaluate_fold(fold, COARSE_BEAM_WIDTH)
        results.append(result)
        beam_states += result.beam_states
        if result.feasible:
            break
        if len(results) % 1000 == 0:
            best = max(results, key=score)
            print(
                f"coarse unique folds={len(results)}; "
                f"best={best.reconstruction_rank}/{best.quotient_dimension}; "
                f"anchors={len(best.anchors)}",
                flush=True,
            )
    print(f"coarse distinct folds = {len(results)}")
    print(f"coarse beam rowspaces = {beam_states}")
    print(f"coarse candidates = {sum(result.feasible for result in results)}")
    print(f"coarse elapsed seconds = {perf_counter() - started:.3f}")
    return results, beam_states


def refine(results: Sequence[FoldResult]) -> list[FoldResult]:
    selected = sorted(results, key=score, reverse=True)[:REFINE_COUNT]
    refined: list[FoldResult] = []
    started = perf_counter()
    for result in selected:
        rerun = evaluate_fold(result.fold, REFINE_BEAM_WIDTH)
        refined.append(rerun)
        if rerun.feasible:
            break
    print(f"refined folds = {len(refined)}")
    print(f"refined candidates = {sum(result.feasible for result in refined)}")
    print(
        f"refined beam rowspaces = "
        f"{sum(result.beam_states for result in refined)}"
    )
    print(f"refine elapsed seconds = {perf_counter() - started:.3f}")
    return sorted(refined, key=score, reverse=True)


def exact_audit(results: Sequence[FoldResult]) -> None:
    print(f"\n[exact rowspace audit of strongest {EXACT_AUDIT_COUNT} folds]")
    for index, result in enumerate(results[:EXACT_AUDIT_COUNT], start=1):
        model = quotient_model(N, fold_edges(result.fold))
        exact = exact_subspace_search(model, result.encoding_budget)
        print(
            f"audit {index}: rowspaces={exact.rowspaces_tested}; "
            f"feasible={exact.feasible}; best="
            f"{exact.best_rank}/{exact.quotient_dimension}; "
            f"anchors={exact.best_anchors}",
            flush=True,
        )


def main() -> None:
    print("Yuanjiang heterogeneous-tertiary-fold search")
    print(f"rank prime = {P}")
    print("target gate: 1.674733895 < 17/10 = 1.7 < 7/4")
    print("legality gate: differently distributed trees glue only over one atom")
    print("dimension gate: every candidate requires at least three anchors")
    motifs = rooted_motifs()
    print(
        f"rooted motifs = {len(motifs)} "
        f"(star=115, fork=736); raw folds={len(motifs) * 5 * 5}"
    )
    coarse, _ = coarse_screen()
    if not coarse:
        print("no folds generated")
        return
    if any(result.feasible for result in coarse):
        refined = sorted(coarse, key=score, reverse=True)
    else:
        refined = refine(coarse)
    print("\n[strongest heterogeneous folds]")
    for index, result in enumerate(refined[:5], start=1):
        print_result(f"rank {index}:", result)
    exact_audit(refined)
    if refined[0].feasible:
        print("RESULT: a 17/10 linear certificate candidate was found")
    else:
        print(
            "RESULT: no 17/10 candidate in this deterministic heterogeneous "
            "star/fork screen; this is not a global no-go theorem"
        )


if __name__ == "__main__":
    main()
