#!/usr/bin/env python3
"""Search Möbius-style twists beyond the aligned eight-atom copy search.

The script reuses the exact rank engine in ``eight_atom_cycle_search.py`` and
keeps two logically different searches separate.

PROVED-SAFE SEAMS
-----------------
Take two copies of any cost-7 four-atom seed.  Glue an arbitrary atom i in the
first block to an arbitrary atom j in the second block over either:

* one input projection (cost at most M), or
* the full atom, represented by two independent projections (cost at most 2M).

This coupling is always legitimate: every atom has the same U=(X,Y) marginal,
so the two seam variables have identical distributions and the standard gluing
/ copy lemma applies.  The previous search only allowed i=j.

FORMAL FULL-RUNG DIAGNOSTIC
---------------------------
Also connect all four atoms through a permutation pi, one colored projection
per rung.  This is the graph-theoretic analogue of twisting one boundary before
closing a four-rung ladder.  It is a candidate generator, not automatically an
entropy proof: the two four-coordinate boundary tuples need not have identical
joint distributions.  Any candidate from this mode must pass a separate
copy/gluing validity audit before it can imply an inequality.

For either mode, a strict improvement on 7/4 requires accounted total cost at
most 13 on eight atoms.  Rank computations are exact over Q by the same
Hadamard argument documented in ``eight_atom_cycle_search.py``.
"""

from __future__ import annotations

from contextlib import redirect_stdout
from dataclasses import dataclass
from io import StringIO
from itertools import permutations, product
from time import perf_counter
from typing import Iterable, Sequence

from eight_atom_cycle_search import (
    Edge,
    FourAtomSchema,
    GlobalNearMiss,
    Observation,
    RREF,
    SLOPES,
    constraint_rows,
    format_observations,
    four_atom_edges,
    general_four_atom_search,
    quotient_model,
    search_encodings,
)


@dataclass(frozen=True)
class Seed:
    schema: FourAtomSchema
    edges: tuple[Edge, ...]


@dataclass(frozen=True)
class TwistCandidate:
    seed_index: int
    mode: str
    seam: str
    encoding: tuple[Observation, ...]
    anchors: tuple[int, ...]
    collision_cost: int
    total_cost: int
    permutation: tuple[int, ...] | None = None
    rung_labels: tuple[str, ...] | None = None


def load_seeds() -> tuple[Seed, ...]:
    # Suppress the already-recorded four-atom report while reusing its exact
    # enumeration rather than reimplementing that wheel.
    with redirect_stdout(StringIO()):
        feasible = general_four_atom_search()
    unique: dict[tuple[Edge, ...], FourAtomSchema] = {}
    for schema in feasible:
        edges = tuple(sorted(four_atom_edges(schema.first_copy, schema.second_copy)))
        unique.setdefault(edges, schema)
    assert len(unique) == 16
    return tuple(Seed(schema=schema, edges=edges) for edges, schema in unique.items())


def doubled_seed(seed: Sequence[Edge]) -> list[Edge]:
    edges = list(seed)
    edges.extend((left + 4, right + 4, label) for left, right, label in seed)
    return edges


def constraint_signature(edges: Sequence[Edge]) -> tuple[tuple[int, tuple[int, ...]], ...]:
    basis = RREF(16, constraint_rows(8, edges))
    return tuple((pivot, tuple(row)) for pivot, row in sorted(basis.rows.items()))


def update_near_miss(
    current: GlobalNearMiss | None,
    candidate: GlobalNearMiss,
) -> GlobalNearMiss:
    if current is None or (
        candidate.reconstruction_rank - candidate.quotient_dimension,
        len(candidate.anchors),
        -candidate.total_cost,
    ) > (
        current.reconstruction_rank - current.quotient_dimension,
        len(current.anchors),
        -current.total_cost,
    ):
        return candidate
    return current


def safe_seam_search(seeds: Sequence[Seed]) -> tuple[list[TwistCandidate], int, int]:
    """Exhaust every universally valid single-atom twisted seam."""
    candidates: list[TwistCandidate] = []
    graph_count = 0
    encoding_count = 0
    seen: dict[tuple[int, tuple[tuple[int, tuple[int, ...]], ...]], None] = {}
    near_miss: GlobalNearMiss | None = None
    started = perf_counter()

    print("[proved-safe arbitrary-atom seams]")
    for seed_index, seed in enumerate(seeds, start=1):
        base = doubled_seed(seed.edges)
        variants: list[tuple[str, int, list[Edge]]] = [("no seam", 0, base)]
        for left in range(4):
            for right in range(4):
                for label in SLOPES:
                    variants.append(
                        (
                            f"{label}[{left}]={label}[{right + 4}]",
                            1,
                            base + [(left, right + 4, label)],
                        )
                    )
                # Any two distinct input projections determine U.  X and Y
                # are a canonical representation of the full-atom seam.
                variants.append(
                    (
                        f"U[{left}]=U[{right + 4}]",
                        2,
                        base
                        + [
                            (left, right + 4, "X"),
                            (left, right + 4, "Y"),
                        ],
                    )
                )

        seed_graphs = 0
        seed_encodings = 0
        for seam, seam_cost, edges in variants:
            signature = constraint_signature(edges)
            key = (seam_cost, signature)
            if key in seen:
                continue
            seen[key] = None
            graph_count += 1
            seed_graphs += 1
            model = quotient_model(8, edges)
            result = search_encodings(model, 5 - seam_cost)
            encoding_count += result.tested
            seed_encodings += result.tested

            near_miss = update_near_miss(
                near_miss,
                GlobalNearMiss(
                    reconstruction_rank=result.best_rank,
                    quotient_dimension=result.quotient_dimension,
                    cross_components=(),
                    encoding=result.best_encoding,
                    anchors=result.best_anchors,
                    total_cost=8 + seam_cost + len(result.best_encoding),
                ),
            )
            if result.feasible:
                candidates.append(
                    TwistCandidate(
                        seed_index=seed_index,
                        mode="proved-safe seam",
                        seam=seam,
                        encoding=result.encoding,
                        anchors=result.anchors,
                        collision_cost=8 + seam_cost,
                        total_cost=8 + seam_cost + len(result.encoding),
                    )
                )
        print(
            f"seed {seed_index:02d}/{len(seeds)}: "
            f"new graphs/encodings={seed_graphs}/{seed_encodings}; "
            f"candidates so far={len(candidates)}",
            flush=True,
        )

    elapsed = perf_counter() - started
    print(f"distinct cost-aware graphs = {graph_count}")
    print(f"encodings tested = {encoding_count}")
    print(f"strict candidates = {len(candidates)}")
    print(f"elapsed seconds = {elapsed:.3f}")
    if near_miss is not None:
        print(
            f"best rejected reconstruction rank = "
            f"{near_miss.reconstruction_rank}/{near_miss.quotient_dimension}; "
            f"anchors={near_miss.anchors}; accounted cost={near_miss.total_cost}"
        )
    return candidates, graph_count, encoding_count


def full_rung_edges(
    seed: Sequence[Edge],
    permutation: Sequence[int],
    labels: Sequence[str],
) -> list[Edge]:
    edges = doubled_seed(seed)
    edges.extend(
        (source, 4 + permutation[source], labels[source]) for source in range(4)
    )
    return edges


def formal_full_rung_search(
    seeds: Sequence[Seed],
) -> tuple[list[TwistCandidate], int, int]:
    """Search all four-rung permutation twists; candidates need a gluing audit."""
    candidates: list[TwistCandidate] = []
    raw_graphs = 0
    encoding_count = 0
    unique: dict[
        tuple[tuple[int, tuple[int, ...]], ...],
        tuple[int, tuple[int, ...], tuple[str, ...], list[Edge]],
    ] = {}
    started = perf_counter()
    labels_universe = tuple(SLOPES)

    print("\n[formal four-rung Möbius/permutation diagnostic]")
    for seed_index, seed in enumerate(seeds, start=1):
        for permutation in permutations(range(4)):
            for labels in product(labels_universe, repeat=4):
                raw_graphs += 1
                edges = full_rung_edges(seed.edges, permutation, labels)
                signature = constraint_signature(edges)
                unique.setdefault(
                    signature,
                    (seed_index, tuple(permutation), tuple(labels), edges),
                )

    print(f"raw labeled graphs = {raw_graphs}")
    print(f"distinct fixed-label constraint spaces = {len(unique)}")
    near_miss: GlobalNearMiss | None = None
    for index, (_, payload) in enumerate(unique.items(), start=1):
        seed_index, permutation, labels, edges = payload
        model = quotient_model(8, edges)
        result = search_encodings(model, 1)
        encoding_count += result.tested
        near_miss = update_near_miss(
            near_miss,
            GlobalNearMiss(
                reconstruction_rank=result.best_rank,
                quotient_dimension=result.quotient_dimension,
                cross_components=(),
                encoding=result.best_encoding,
                anchors=result.best_anchors,
                total_cost=12 + len(result.best_encoding),
            ),
        )
        if result.feasible:
            seam = ", ".join(
                f"{labels[source]}[{source}]="
                f"{labels[source]}[{4 + permutation[source]}]"
                for source in range(4)
            )
            candidates.append(
                TwistCandidate(
                    seed_index=seed_index,
                    mode="formal four-rung diagnostic",
                    seam=seam,
                    encoding=result.encoding,
                    anchors=result.anchors,
                    collision_cost=12,
                    total_cost=12 + len(result.encoding),
                    permutation=permutation,
                    rung_labels=labels,
                )
            )
        if index % 5000 == 0:
            print(
                f"tested {index}/{len(unique)} distinct spaces; "
                f"candidates={len(candidates)}",
                flush=True,
            )

    elapsed = perf_counter() - started
    print(f"encodings tested = {encoding_count}")
    print(f"formal rank candidates = {len(candidates)}")
    print(f"elapsed seconds = {elapsed:.3f}")
    if near_miss is not None:
        print(
            f"best rejected reconstruction rank = "
            f"{near_miss.reconstruction_rank}/{near_miss.quotient_dimension}; "
            f"anchors={near_miss.anchors}; accounted cost={near_miss.total_cost}"
        )
    return candidates, len(unique), encoding_count


def print_candidates(title: str, candidates: Iterable[TwistCandidate]) -> None:
    items = list(candidates)
    print(f"\n[{title}: {len(items)}]")
    for index, candidate in enumerate(items[:20], start=1):
        print(f"candidate {index}:")
        print(f"  seed = {candidate.seed_index}; mode = {candidate.mode}")
        print(f"  seam = {candidate.seam}")
        if candidate.permutation is not None:
            print(f"  permutation = {candidate.permutation}")
        print(f"  encoding = {format_observations(candidate.encoding)}")
        print(f"  anchors = {candidate.anchors}")
        print(
            f"  total cost = {candidate.total_cost}; "
            f"ratio = {candidate.total_cost}/8 = {candidate.total_cost/8:.6f}"
        )
    if len(items) > 20:
        print(f"  ... {len(items) - 20} additional candidates omitted")


def main() -> None:
    print("Yuanjiang Möbius-twist collision search")
    seeds = load_seeds()
    print(f"cost-7 four-atom seeds loaded = {len(seeds)}")
    safe, _, _ = safe_seam_search(seeds)
    print_candidates("proved-safe strict candidates", safe)
    formal, _, _ = formal_full_rung_search(seeds)
    print_candidates("formal candidates requiring gluing validity", formal)


if __name__ == "__main__":
    main()
