#!/usr/bin/env python3
"""Independent rational-arithmetic audit for the collision-cycle search.

This file deliberately does not import ``eight_atom_cycle_search``.  It uses
``fractions.Fraction`` rather than finite-field arithmetic and independently
rebuilds the two bounded searches that matter most for the certificate:

* every four-atom, two-layer schema of accounted cost at most 7; and
* every extra-copy tuple W for the original Katz--Tao seed at total cost at
  most 13.

It also checks the reported global near-miss from the broader 16-seed search.
The broader search itself is exact over Q by the Hadamard bound documented in
``eight_atom_cycle_search.py``; this slower script is a separate implementation
intended to catch coding or criterion mistakes.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from time import perf_counter
from typing import Iterable, Sequence


SLOPES: dict[str, tuple[int, int]] = {
    "X": (1, 0),
    "Y": (0, 1),
    "C=X+Y": (1, 1),
    "D=X+2Y": (1, 2),
}
TARGET = (1, -1)

Edge = tuple[int, int, str]
Observation = tuple[int, str]
Vector = tuple[Fraction, ...]


class RationalBasis:
    """Small reduced row-space basis over Q."""

    def __init__(self, width: int, rows: Iterable[Sequence[int | Fraction]] = ()):
        self.width = width
        self.rows: dict[int, list[Fraction]] = {}
        for row in rows:
            self.add(row)

    def copy(self) -> "RationalBasis":
        result = RationalBasis(self.width)
        result.rows = {pivot: row.copy() for pivot, row in self.rows.items()}
        return result

    @property
    def rank(self) -> int:
        return len(self.rows)

    def reduce(self, vector: Sequence[int | Fraction]) -> list[Fraction]:
        work = [Fraction(value) for value in vector]
        for pivot in sorted(self.rows):
            coefficient = work[pivot]
            if coefficient:
                row = self.rows[pivot]
                work = [a - coefficient * b for a, b in zip(work, row)]
        return work

    def contains(self, vector: Sequence[int | Fraction]) -> bool:
        return not any(self.reduce(vector))

    def add(self, vector: Sequence[int | Fraction]) -> bool:
        work = self.reduce(vector)
        pivot = next((index for index, value in enumerate(work) if value), None)
        if pivot is None:
            return False
        scale = work[pivot]
        work = [value / scale for value in work]
        for old_pivot, old_row in list(self.rows.items()):
            coefficient = old_row[pivot]
            if coefficient:
                self.rows[old_pivot] = [
                    a - coefficient * b for a, b in zip(old_row, work)
                ]
        self.rows[pivot] = work
        return True


def projection_row(
    n: int, atom: int, coefficients: tuple[int, int]
) -> tuple[int, ...]:
    row = [0] * (2 * n)
    row[2 * atom], row[2 * atom + 1] = coefficients
    return tuple(row)


def collision_row(n: int, left: int, right: int, label: str) -> tuple[int, ...]:
    lhs = projection_row(n, left, SLOPES[label])
    rhs = projection_row(n, right, SLOPES[label])
    return tuple(a - b for a, b in zip(lhs, rhs))


@dataclass(frozen=True)
class Model:
    n: int
    constraint_rank: int
    quotient_dimension: int
    observations: tuple[Observation, ...]
    observation_rows: tuple[Vector, ...]
    targets: tuple[Vector, ...]
    xs: tuple[Vector, ...]
    ys: tuple[Vector, ...]


def make_model(n: int, edges: Sequence[Edge]) -> Model:
    constraints = RationalBasis(
        2 * n,
        (collision_row(n, left, right, label) for left, right, label in edges),
    )

    def quotient(vector: Sequence[int]) -> Vector:
        return tuple(constraints.reduce(vector))

    representatives: dict[Vector, Observation] = {}
    for atom in range(n):
        for label, coefficients in SLOPES.items():
            vector = quotient(projection_row(n, atom, coefficients))
            if any(vector):
                representatives.setdefault(vector, (atom, label))

    return Model(
        n=n,
        constraint_rank=constraints.rank,
        quotient_dimension=2 * n - constraints.rank,
        observations=tuple(representatives.values()),
        observation_rows=tuple(representatives),
        targets=tuple(
            quotient(projection_row(n, atom, TARGET)) for atom in range(n)
        ),
        xs=tuple(
            quotient(projection_row(n, atom, (1, 0))) for atom in range(n)
        ),
        ys=tuple(
            quotient(projection_row(n, atom, (0, 1))) for atom in range(n)
        ),
    )


@dataclass(frozen=True)
class SearchResult:
    feasible: bool
    tested: int
    encoding: tuple[Observation, ...]
    anchors: tuple[int, ...]
    best_rank: int


def check_encoding(model: Model, indices: Sequence[int]) -> tuple[bool, int, tuple[int, ...]]:
    basis = RationalBasis(
        2 * model.n, (model.observation_rows[index] for index in indices)
    )
    anchors = tuple(
        atom for atom, target in enumerate(model.targets) if basis.contains(target)
    )
    reconstruction = basis.copy()
    for atom in anchors:
        reconstruction.add(model.xs[atom])
        reconstruction.add(model.ys[atom])
    return reconstruction.rank == model.quotient_dimension, reconstruction.rank, anchors


def search_encodings(model: Model, max_size: int) -> SearchResult:
    tested = 0
    best_rank = -1
    count = len(model.observation_rows)
    for size in range(max_size + 1):
        for indices in combinations(range(count), size):
            tested += 1
            feasible, rank, anchors = check_encoding(model, indices)
            best_rank = max(best_rank, rank)
            if feasible:
                return SearchResult(
                    True,
                    tested,
                    tuple(model.observations[index] for index in indices),
                    anchors,
                    best_rank,
                )
    return SearchResult(False, tested, (), (), best_rank)


def four_atom_edges(
    first_copy: Sequence[str], second_copy: Sequence[Observation]
) -> tuple[Edge, ...]:
    edges: list[Edge] = []
    for label in first_copy:
        edges.extend(((0, 1, label), (2, 3, label)))
    edges.extend((atom, atom + 2, label) for atom, label in second_copy)
    return tuple(edges)


def audit_all_four_atom_schemas() -> tuple[tuple[Edge, ...], ...]:
    first_universe = tuple(SLOPES)
    second_universe = tuple(
        (atom, label) for atom in range(2) for label in SLOPES
    )
    schema_count = 0
    encoding_count = 0
    feasible_edges: dict[tuple[Edge, ...], None] = {}
    strict_count = 0
    started = perf_counter()

    for first_size in range(4):
        for first in combinations(first_universe, first_size):
            max_second_size = 7 - 2 * first_size
            for second_size in range(min(8, max_second_size) + 1):
                collision_cost = 2 * first_size + second_size
                for second in combinations(second_universe, second_size):
                    schema_count += 1
                    edges = four_atom_edges(first, second)
                    result = search_encodings(
                        make_model(4, edges), 7 - collision_cost
                    )
                    encoding_count += result.tested
                    if result.feasible:
                        feasible_edges[tuple(sorted(edges))] = None
                        if collision_cost + len(result.encoding) <= 6:
                            strict_count += 1

    elapsed = perf_counter() - started
    assert schema_count == 1725
    assert encoding_count == 313762
    assert len(feasible_edges) == 16
    assert strict_count == 0
    print("[independent Q audit: all four-atom schemas]")
    print(f"schemas/encodings = {schema_count}/{encoding_count}")
    print(f"cost-7 seeds = {len(feasible_edges)}; cost<=6 seeds = {strict_count}")
    print(f"elapsed seconds = {elapsed:.3f}")
    return tuple(feasible_edges)


KT4_EDGES: tuple[Edge, ...] = (
    (0, 1, "X"),
    (2, 3, "X"),
    (0, 2, "D=X+2Y"),
    (1, 3, "Y"),
)


def doubled_edges(seed: Sequence[Edge], cross: Sequence[Observation]) -> tuple[Edge, ...]:
    return tuple(seed) + tuple(
        (left + 4, right + 4, label) for left, right, label in seed
    ) + tuple((atom, atom + 4, label) for atom, label in cross)


def audit_all_kt_extra_copies() -> None:
    universe = tuple((atom, label) for atom in range(4) for label in SLOPES)
    graphs = 0
    encodings = 0
    started = perf_counter()
    for size in range(6):
        for cross in combinations(universe, size):
            graphs += 1
            result = search_encodings(
                make_model(8, doubled_edges(KT4_EDGES, cross)), 5 - size
            )
            encodings += result.tested
            assert not result.feasible

    elapsed = perf_counter() - started
    assert graphs == 6885
    assert encodings == 558271
    print("\n[independent Q audit: Katz--Tao seed plus every extra-copy W]")
    print(f"graphs/encodings = {graphs}/{encodings}")
    print("strict cost<=13 candidates = 0")
    print(f"elapsed seconds = {elapsed:.3f}")


def exact_global_near_miss_check() -> None:
    seed_one = four_atom_edges(
        ("X",), ((0, "Y"), (1, "C=X+Y"))
    )
    model = make_model(8, doubled_edges(seed_one, ()))
    selected = {
        observation: index for index, observation in enumerate(model.observations)
    }
    encoding = (
        (0, "X"),
        (0, "Y"),
        (4, "C=X+Y"),
        (6, "D=X+2Y"),
        (7, "Y"),
    )
    feasible, rank, anchors = check_encoding(
        model, tuple(selected[observation] for observation in encoding)
    )
    assert not feasible
    assert model.quotient_dimension == 8
    assert rank == 6
    assert anchors == (0, 5)
    print("\n[independent Q audit: reported global near-miss]")
    print("quotient reconstruction rank = 6/8; anchors = (0, 5)")
    print("two independent degrees of freedom remain; candidate rejected")


def main() -> None:
    print("Yuanjiang collision-cycle certificate verifier (exact Q arithmetic)")
    audit_all_four_atom_schemas()
    audit_all_kt_extra_copies()
    exact_global_near_miss_check()
    print("\nALL INDEPENDENT RATIONAL CHECKS PASSED")


if __name__ == "__main__":
    main()
