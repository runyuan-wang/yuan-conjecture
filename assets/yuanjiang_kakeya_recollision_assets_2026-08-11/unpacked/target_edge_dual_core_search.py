#!/usr/bin/env python3
"""Go/no-go screen for a target-conditioned dual-core molecule.

Literature boundary
-------------------
Conditionally independent copying relative to a difference such as X-Y is a
standard additive-combinatorics device (for example Tao, 2010).  This script
does not claim that ingredient as new.  It tests one project-specific use:

* start with the recorded five-atom path and five-atom star;
* attach one new leaf to each by a Z=X-Y conditional-copy edge;
* match the two resulting Z-coloured edge laws over two input projections.

There are 12 atoms, 8 local input-coloured edges, 2 Z-coloured edges and a
rank-r input interface.  If an input-projection encoding has m components and
the usual anchor/reconstruction certificate holds, entropy bookkeeping gives

    (12-k)H(U) + (k-2)H(Z) <= (8+r+m)M.

For k>=2, H(U)>=H(Z), r=2 and m=7 this implies 10H(Z)<=17M.  Thus two target
edges add two atoms without changing the effective 17/10 gate.

The family is intentionally tiny: all 5x5 attachment positions, both endpoint
orientations and all 4x4 ordered input-projection interfaces (800 raw and 800
distinct rowspaces).  A deterministic width-16 beam is only a mechanism gate.
If its best deficit is not smaller than the preceding 2-dimensional deficit,
the script stops without expensive exact enumeration and makes no family-wide
no-go claim.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Sequence

from eight_atom_cycle_search import (
    Edge,
    Observation,
    QuotientModel,
    RREF,
    SLOPES,
    TARGET,
    projection_row,
)
from protein_tertiary_fold_search import (
    basis_signature,
    diverse_beam_search,
    format_observations,
)


N = 12
C = "C=X+Y"
D = "D=X+2Y"
Z_LABEL = "Z=X-Y"
INPUT_LABELS = tuple(SLOPES)
EXTENDED_SLOPES = {**SLOPES, Z_LABEL: TARGET}
BEAM_WIDTH = 16
PREVIOUS_BEST_DEFICIT = 2

SEED_EDGES: tuple[Edge, ...] = (
    (0, 1, "X"),
    (1, 2, "X"),
    (2, 3, "X"),
    (3, 4, "Y"),
)
STAR_EDGES: tuple[Edge, ...] = (
    (6, 7, C),
    (6, 8, C),
    (6, 9, D),
    (6, 10, "Y"),
)


def extended_collision_row(n: int, edge: Edge) -> list[int]:
    left, right, label = edge
    coefficients = EXTENDED_SLOPES[label]
    row = projection_row(n, left, coefficients)
    other = projection_row(n, right, coefficients)
    return [a - b for a, b in zip(row, other)]


def constraint_rows(edges: Sequence[Edge]) -> list[list[int]]:
    return [extended_collision_row(N, edge) for edge in edges]


def input_observation_model(edges: Sequence[Edge]) -> QuotientModel:
    constraints = RREF(2 * N, constraint_rows(edges))
    free_columns = [
        column for column in range(2 * N) if column not in constraints.rows
    ]

    def quotient(row: Sequence[int]) -> tuple[int, ...]:
        remainder = constraints.reduce(row)
        return tuple(remainder[column] for column in free_columns)

    # Z is allowed as a conditional-copy edge but not as an M-priced encoding
    # coordinate.  Only the four input projections enter this observation set.
    representatives: dict[tuple[int, ...], Observation] = {}
    for atom in range(N):
        for label, coefficients in SLOPES.items():
            vector = quotient(projection_row(N, atom, coefficients))
            if any(vector):
                representatives.setdefault(vector, (atom, label))

    observation_rows = tuple(representatives)
    return QuotientModel(
        n=N,
        constraint_rank=constraints.rank,
        quotient_dimension=len(free_columns),
        observation_labels=tuple(
            representatives[row] for row in observation_rows
        ),
        observation_rows=observation_rows,
        target_rows=tuple(
            quotient(projection_row(N, atom, TARGET)) for atom in range(N)
        ),
        x_rows=tuple(
            quotient(projection_row(N, atom, (1, 0))) for atom in range(N)
        ),
        y_rows=tuple(
            quotient(projection_row(N, atom, (0, 1))) for atom in range(N)
        ),
    )


def constraint_signature(edges: Sequence[Edge]) -> tuple:
    return basis_signature(RREF(2 * N, constraint_rows(edges)))


@dataclass(frozen=True)
class ScreenResult:
    seed_attachment: int
    star_attachment: int
    reverse: bool
    first_projection: str
    second_projection: str
    interface_rank: int
    encoding_budget: int
    reconstruction_rank: int
    quotient_dimension: int
    anchors: tuple[int, ...]
    encoding: tuple[Observation, ...]
    beam_states: int
    feasible: bool

    @property
    def deficit(self) -> int:
        return self.quotient_dimension - self.reconstruction_rank


def score(result: ScreenResult) -> tuple[int, int, int]:
    return (int(result.feasible), -result.deficit, len(result.anchors))


def run_screen() -> tuple[list[ScreenResult], int, int]:
    raw = 0
    seen: set[tuple] = set()
    results: list[ScreenResult] = []

    for seed_attachment in range(5):
        for star_local_attachment in range(5):
            star_attachment = 6 + star_local_attachment
            local_edges = (
                *SEED_EDGES,
                (seed_attachment, 5, Z_LABEL),
                *STAR_EDGES,
                (star_attachment, 11, Z_LABEL),
            )
            local_rank = RREF(2 * N, constraint_rows(local_edges)).rank
            assert local_rank == 10

            for reverse in (False, True):
                star_first, star_second = (
                    (star_attachment, 11)
                    if not reverse
                    else (11, star_attachment)
                )
                for first_projection in INPUT_LABELS:
                    for second_projection in INPUT_LABELS:
                        raw += 1
                        edges = (
                            *local_edges,
                            (
                                seed_attachment,
                                star_first,
                                first_projection,
                            ),
                            (5, star_second, second_projection),
                        )
                        signature = constraint_signature(edges)
                        if signature in seen:
                            continue
                        seen.add(signature)

                        model = input_observation_model(edges)
                        interface_rank = model.constraint_rank - local_rank
                        if interface_rank not in (1, 2):
                            continue
                        # 8 local input edges + r interface edges + m encoding
                        # components must total at most 17.
                        encoding_budget = 9 - interface_rank
                        best, visited = diverse_beam_search(
                            model, encoding_budget, width=BEAM_WIDTH
                        )
                        encoding = tuple(
                            model.observation_labels[index]
                            for index in best.indices
                        )
                        results.append(
                            ScreenResult(
                                seed_attachment=seed_attachment,
                                star_attachment=star_local_attachment,
                                reverse=reverse,
                                first_projection=first_projection,
                                second_projection=second_projection,
                                interface_rank=interface_rank,
                                encoding_budget=encoding_budget,
                                reconstruction_rank=best.check.reconstruction_rank,
                                quotient_dimension=model.quotient_dimension,
                                anchors=best.check.anchors,
                                encoding=encoding,
                                beam_states=visited,
                                feasible=best.check.feasible,
                            )
                        )
    return results, raw, len(seen)


def main() -> None:
    results, raw, distinct = run_screen()
    assert raw == 800
    assert distinct == 800
    assert len(results) == 800
    results.sort(key=score, reverse=True)

    candidates = sum(result.feasible for result in results)
    total_states = sum(result.beam_states for result in results)
    deficit_histogram = Counter(result.deficit for result in results)
    interface_histogram = Counter(result.interface_rank for result in results)
    best = results[0]

    assert candidates == 0
    assert best.reconstruction_rank == 8
    assert best.quotient_dimension == 12
    assert best.deficit == 4
    assert best.deficit >= PREVIOUS_BEST_DEFICIT

    print("Yuanjiang target-edge dual-core gate")
    print(f"raw/distinct/evaluated = {raw}/{distinct}/{len(results)}")
    print(f"beam width / total states = {BEAM_WIDTH}/{total_states}")
    print(f"interface-rank histogram = {dict(sorted(interface_histogram.items()))}")
    print(f"deficit histogram = {dict(sorted(deficit_histogram.items()))}")
    print(f"candidates = {candidates}")
    print(
        f"best rank/quotient/deficit = {best.reconstruction_rank}/"
        f"{best.quotient_dimension}/{best.deficit}"
    )
    print(f"best anchors = {best.anchors}")
    print(
        "best attachment/orientation/interface = "
        f"{best.seed_attachment}/{best.star_attachment}/"
        f"{'reverse' if best.reverse else 'parallel'}/"
        f"({best.first_projection},{best.second_projection})"
    )
    print(f"best encoding = {format_observations(best.encoding)}")
    print(
        "RESULT: mechanism gate failed; target edges increased the best "
        "deficit from 2 to 4, so exact enumeration was not authorized"
    )


if __name__ == "__main__":
    main()
