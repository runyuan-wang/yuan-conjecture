#!/usr/bin/env python3
"""Exact anatomy of the two missing dimensions in the best dual-core fold.

The strongest certificate in ``dual_core_tertiary_certificate_2026-08-12.json``
has quotient dimension 10 and reconstruction rank 8.  This script does not
rerun the 15,936-fold search.  It takes that one recorded near miss, computes
the nullspace over Q, and identifies the two invisible motions explicitly.

It then performs only the pre-authorized cheap test suggested by that anatomy:
rewire one C-coloured root--leaf edge as a leaf--leaf edge in each of the four
input colours.  A width-256 deterministic beam is used only as a go/no-go
gate.  Since none improves the 8/10 rank, no expensive exact enumeration is
run or claimed.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
from math import gcd, lcm
from typing import Sequence

from eight_atom_cycle_search import (
    SLOPES,
    TARGET,
    collision_row,
    projection_row,
    quotient_model,
)
from protein_tertiary_fold_search import diverse_beam_search, format_observations


N = 10
C = "C=X+Y"
D = "D=X+2Y"

STRONGEST_EDGES = (
    (0, 1, "X"),
    (1, 2, "X"),
    (2, 3, "X"),
    (3, 4, "Y"),
    (5, 6, C),
    (5, 7, C),
    (5, 8, D),
    (5, 9, "Y"),
    (3, 5, C),
    (4, 9, "X"),
)

STRONGEST_ENCODING = (
    (0, "X"),
    (4, C),
    (5, D),
    (0, "Y"),
    (1, "Y"),
    (2, "Y"),
    (8, "X"),
)

ANCHORS = (0, 1, 2, 8, 9)


def rational_rref_nullspace(
    rows: Sequence[Sequence[int]], ncols: int
) -> tuple[tuple[int, ...], list[list[Fraction]]]:
    """Return pivot columns and a basis of the exact Q-nullspace."""
    matrix = [[Fraction(value) for value in row] for row in rows]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(ncols):
        found = next(
            (
                row
                for row in range(pivot_row, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if found is None:
            continue
        matrix[pivot_row], matrix[found] = matrix[found], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                left - scale * right
                for left, right in zip(matrix[row], matrix[pivot_row])
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break

    free = [column for column in range(ncols) if column not in pivots]
    nullspace: list[list[Fraction]] = []
    for free_column in free:
        vector = [Fraction(0) for _ in range(ncols)]
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -matrix[row][free_column]
        nullspace.append(vector)
    return tuple(pivots), nullspace


def primitive_integer_vector(vector: Sequence[Fraction]) -> tuple[int, ...]:
    denominator = reduce(lcm, (value.denominator for value in vector), 1)
    integers = [int(value * denominator) for value in vector]
    divisor = reduce(gcd, (abs(value) for value in integers if value), 0) or 1
    integers = [value // divisor for value in integers]
    first = next((value for value in integers if value), 1)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def dot(left: Sequence[int], right: Sequence[int]) -> int:
    return sum(a * b for a, b in zip(left, right))


def exact_residual_modes() -> tuple[tuple[int, ...], tuple[int, ...]]:
    rows = [collision_row(N, *edge) for edge in STRONGEST_EDGES]
    rows.extend(
        projection_row(N, atom, SLOPES[label])
        for atom, label in STRONGEST_ENCODING
    )
    for atom in ANCHORS:
        rows.append(projection_row(N, atom, (1, 0)))
        rows.append(projection_row(N, atom, (0, 1)))

    pivots, rational_basis = rational_rref_nullspace(rows, 2 * N)
    basis = tuple(primitive_integer_vector(vector) for vector in rational_basis)
    assert len(rows) == 27
    assert len(pivots) == 18
    assert len(basis) == 2

    expected_first = tuple(
        coordinate
        for atom in range(N)
        for coordinate in ((1, -1) if atom == 6 else (0, 0))
    )
    expected_second = tuple(
        coordinate
        for atom in range(N)
        for coordinate in ((1, -1) if atom == 7 else (0, 0))
    )
    assert basis == (expected_first, expected_second)
    return basis  # type: ignore[return-value]


def print_residual_signatures(basis: Sequence[Sequence[int]]) -> None:
    for atom in range(N):
        entries: list[str] = []
        for label, coefficients in (*SLOPES.items(), ("Z=X-Y", TARGET)):
            row = projection_row(N, atom, coefficients)
            signature = tuple(dot(row, mode) for mode in basis)
            if any(signature):
                entries.append(f"{label}:{signature}")
        if entries:
            print(f"atom {atom} residual signatures = {', '.join(entries)}")


def rewiring_gate() -> None:
    print("\n[leaf-to-leaf rewiring gate]")
    for bridge in SLOPES:
        edges = (
            (0, 1, "X"),
            (1, 2, "X"),
            (2, 3, "X"),
            (3, 4, "Y"),
            (5, 6, C),
            (6, 7, bridge),
            (5, 8, D),
            (5, 9, "Y"),
            (3, 5, C),
            (4, 9, "X"),
        )
        model = quotient_model(N, edges)
        best, visited = diverse_beam_search(model, max_size=7, width=256)
        encoding = tuple(
            model.observation_labels[index] for index in best.indices
        )
        assert model.quotient_dimension == 10
        assert best.check.reconstruction_rank == 8
        assert not best.check.feasible
        print(
            f"bridge={bridge}; rank={best.check.reconstruction_rank}/10; "
            f"anchors={best.check.anchors}; states={visited}; "
            f"encoding={format_observations(encoding)}"
        )


def main() -> None:
    basis = exact_residual_modes()
    print("strongest fold exact matrix rows/rank/nullity = 27/18/2")
    for index, mode in enumerate(basis, start=1):
        support = {
            atom: (mode[2 * atom], mode[2 * atom + 1])
            for atom in range(N)
            if mode[2 * atom] or mode[2 * atom + 1]
        }
        print(f"invisible mode {index} = {support}")
    print_residual_signatures(basis)
    rewiring_gate()
    print(
        "RESULT: the two exact obstructions are independent X-Y fibre slides "
        "at sibling leaves 6 and 7; simple rewiring did not reduce the deficit"
    )


if __name__ == "__main__":
    main()
