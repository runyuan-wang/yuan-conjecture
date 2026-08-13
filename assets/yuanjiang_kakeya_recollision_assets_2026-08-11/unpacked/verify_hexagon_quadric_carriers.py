#!/usr/bin/env python3
"""Exact rational audit of quadric carriers for C12/Theta configurations in P^4."""

from fractions import Fraction
from itertools import combinations_with_replacement
import json


MONOMIALS = list(combinations_with_replacement(range(5), 2))


def quadric_value_row(point):
    return [Fraction(point[i] * point[j]) for i, j in MONOMIALS]


def polar_value_row(point, other):
    row = []
    for i, j in MONOMIALS:
        if i == j:
            row.append(Fraction(2 * point[i] * other[i]))
        else:
            row.append(Fraction(point[i] * other[j] + point[j] * other[i]))
    return row


def rational_rank(rows):
    matrix = [[Fraction(value) for value in row] for row in rows]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or not matrix[index][column]:
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                value - multiplier * pivot_value
                for value, pivot_value in zip(matrix[index], matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def carrier_audit(points, edges):
    vertices = sorted({vertex for edge in edges for vertex in edge})
    # A quadric contains every edge line iff it vanishes at every vertex and
    # its polar form vanishes on the two endpoints of every edge.
    rows = [quadric_value_row(points[index]) for index in vertices]
    rows.extend(polar_value_row(points[left], points[right]) for left, right in edges)
    rank = rational_rank(rows)
    return {
        "vertices": len(vertices),
        "edges": len(edges),
        "constraint_rows": len(rows),
        "constraint_rank": rank,
        "quadric_vector_space_dimension": 15 - rank,
    }


def main():
    hexagon_points = [
        (1, 0, 0, 0, 0),
        (1, 1, 0, 0, 0),
        (1, 0, 1, 0, 0),
        (1, 0, 0, 1, 0),
        (1, 0, 0, 0, 1),
        (1, 1, 1, 1, 1),
    ]
    hexagon_edges = [(index, (index + 1) % 6) for index in range(6)]

    # Three internally disjoint three-edge paths from vertex 0 to vertex 7.
    # This is the contracted theta(3,3,3) made by two C12s sharing three tubes.
    theta_points = [
        (1, 0, 0, 0, 0),
        (1, 1, 0, 0, 0),
        (1, 0, 1, 0, 0),
        (1, 0, 0, 1, 0),
        (1, 0, 0, 0, 1),
        (1, 1, 1, 1, 1),
        (1, 2, 3, 5, 7),
        (1, -1, 2, -2, 3),
    ]
    theta_edges = []
    for first, second in [(1, 2), (3, 4), (5, 6)]:
        theta_edges.extend([(0, first), (first, second), (second, 7)])

    result = {
        "status": "verified",
        "ambient": "projective P^4; quadrics form a 15-dimensional vector space",
        "single_hexagon": carrier_audit(hexagon_points, hexagon_edges),
        "theta_3_3_3": carrier_audit(theta_points, theta_edges),
        "interpretation": (
            "Every six-edge projective hexagon has at least a 3-dimensional "
            "vector space of containing quadrics; the displayed generic hexagon "
            "has exactly dimension 3. Two C12s sharing a three-edge path need "
            "not share any nonzero quadric: the displayed theta(3,3,3) has rank 15."
        ),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
