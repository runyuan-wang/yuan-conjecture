#!/usr/bin/env python3
"""Enumerate small R^4 parallel-rigidity gadgets.

The exact certificate uses integer moment-curve coordinates and rational
Gaussian elimination.  A normalized floating-point matrix is also used to
record a conditioning diagnostic.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import numpy as np


def points(n: int) -> np.ndarray:
    # Deterministic integer moment-curve points.  The first coordinate differs
    # on every edge, which lets us use three exact wedge equations per edge.
    t = np.arange(1, n + 1, dtype=float)
    return np.column_stack((t, t**2, t**3, t**4))


def rational_rank(rows: list[list[int]]) -> int:
    """Exact row rank over Q using fraction Gaussian elimination."""
    a = [[Fraction(x) for x in row] for row in rows]
    m = len(a)
    n = len(a[0]) if m else 0
    rank = 0
    for col in range(n):
        pivot = next((r for r in range(rank, m) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        pv = a[rank][col]
        a[rank] = [x / pv for x in a[rank]]
        for r in range(m):
            if r != rank and a[r][col]:
                q = a[r][col]
                a[r] = [x - q * y for x, y in zip(a[r], a[rank])]
        rank += 1
        if rank == m:
            break
    return rank


def perp_rows(direction: np.ndarray) -> np.ndarray:
    """Return a 3 x 4 orthonormal row basis of direction^perp."""
    u = direction / np.linalg.norm(direction)
    # SVD of a 1x4 matrix: the last three right singular vectors span u^perp.
    _, _, vh = np.linalg.svd(u.reshape(1, 4), full_matrices=True)
    return vh[1:, :]


def rigidity_rank(n: int, edges: list[tuple[int, int]]) -> dict:
    p = points(n)
    mat = np.zeros((3 * len(edges), 4 * n))
    exact_rows: list[list[int]] = []
    for k, (i, j) in enumerate(edges):
        direction = p[j] - p[i]
        q = perp_rows(direction)
        mat[3 * k : 3 * k + 3, 4 * i : 4 * i + 4] = q
        mat[3 * k : 3 * k + 3, 4 * j : 4 * j + 4] = -q
        d = [int(x) for x in direction]
        for coordinate in range(1, 4):
            row = [0] * (4 * n)
            # d0*(h_j,k-h_i,k)-dk*(h_j,0-h_i,0)=0.
            row[4 * i] = d[coordinate]
            row[4 * i + coordinate] = -d[0]
            row[4 * j] = -d[coordinate]
            row[4 * j + coordinate] = d[0]
            exact_rows.append(row)
    singular = np.linalg.svd(mat, compute_uv=False)
    tol = 1e-9 * singular[0] if singular.size else 1e-9
    numeric_rank = int(np.sum(singular > tol))
    rank = rational_rank(exact_rows)
    if rank != numeric_rank:
        raise RuntimeError(f"exact rank {rank} != numeric rank {numeric_rank}")
    return {
        "vertices": n,
        "edges": len(edges),
        "circuit_rank": len(edges) - n + 1,
        "rank": rank,
        "rank_is_exact_over_Q": True,
        "nullity": 4 * n - rank,
        "extra_flex_beyond_translation_scale": 4 * n - rank - 5,
        "parallel_rigid": rank == 4 * n - 5,
        "smallest_nonzero_singular": float(singular[rank - 1]) if rank else 0.0,
    }


def cycle(k: int) -> tuple[int, list[tuple[int, int]]]:
    return k, [(i, (i + 1) % k) for i in range(k)]


def theta(a: int, b: int, c: int) -> tuple[int, list[tuple[int, int]]]:
    """Three internally disjoint A--B paths of lengths a,b,c."""
    edges: list[tuple[int, int]] = []
    nxt = 2
    for length in (a, b, c):
        path = [0]
        path.extend(range(nxt, nxt + length - 1))
        nxt += length - 1
        path.append(1)
        edges.extend(zip(path, path[1:]))
    return nxt, edges


def main() -> None:
    out: dict[str, object] = {
        "ambient_dimension": 4,
        "coordinate_model": "integer moment curve (t,t^2,t^3,t^4)",
        "cycles": {},
        "theta": {},
        "chorded_C6": {},
    }
    for k in range(3, 11):
        n, e = cycle(k)
        out["cycles"][f"C{k}"] = rigidity_rank(n, e)
    for r in range(1, 6):
        # Union of two C6s sharing a consecutive r-edge path.
        n, e = theta(r, 6 - r, 6 - r)
        out["theta"][f"shared_{r}_edges"] = rigidity_rank(n, e)
    n, base = cycle(6)
    for chord in ((0, 2), (0, 3), (0, 4)):
        out["chorded_C6"][f"chord_{chord[0]}_{chord[1]}"] = rigidity_rank(
            n, base + [chord]
        )

    target = Path("parallel_rigidity_gadget_certificate.json")
    target.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
