#!/usr/bin/env python3
"""Two-slice compression and C4 statistics for spacetime line families.

At times s<u, each line gives an edge between its spatial cell at s and its
spatial cell at u.  Four-cycle counts measure repeated 2-by-2 grids of
trajectories.  In exact geometry every such C4 obeys the direction relation

    v11 - v12 - v21 + v22 = 0.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
import math
import numpy as np

from cross_time_recollision import (
    Family,
    make_book,
    make_random,
    make_regulus,
    make_smooth_fold,
    make_star,
    make_two_bush_grid,
    make_worst_affine,
)


def cells_at(family: Family, t: float, cell: float) -> tuple[np.ndarray, int]:
    x = family.a + t * family.v
    keys = np.floor(x / cell + 1e-10).astype(np.int64)
    _, inv = np.unique(keys, axis=0, return_inverse=True)
    return inv, int(inv.max() + 1)


def bipartite_stats(family: Family, s: float, u: float, cell: float) -> dict[str, int | float]:
    left, na = cells_at(family, s, cell)
    right, nb = cells_at(family, u, cell)
    edge_mult: dict[tuple[int, int], int] = defaultdict(int)
    for a, b in zip(left.tolist(), right.tolist()):
        edge_mult[(a, b)] += 1

    adj: dict[int, set[int]] = defaultdict(set)
    for a, b in edge_mult:
        adj[a].add(b)

    common_pair_count: dict[tuple[int, int], int] = defaultdict(int)
    for neighbours in adj.values():
        for b1, b2 in combinations(sorted(neighbours), 2):
            common_pair_count[(b1, b2)] += 1
    c4 = sum(k * (k - 1) // 2 for k in common_pair_count.values())

    unique_edges = len(edge_mult)
    parallel_edges = family.n - unique_edges
    return {
        "s": s,
        "u": u,
        "A": na,
        "B": nb,
        "lines": family.n,
        "unique_edges": unique_edges,
        "edge_excess": parallel_edges,
        "c4": int(c4),
    }


def best_pair(family: Family, cell: float, extra_times: list[float] | None = None) -> dict[str, int | float]:
    times = list(np.linspace(0.0, 1.0, 101))
    if extra_times:
        times.extend(extra_times)
    times = sorted(set(round(float(t), 12) for t in times))
    cache = {t: cells_at(family, t, cell)[1] for t in times}
    candidates = []
    for i, s in enumerate(times):
        for u in times[i + 1 :]:
            if u - s >= 0.2:
                candidates.append((cache[s] * cache[u], cache[s] + cache[u], s, u))
    _, _, s, u = min(candidates)
    return bipartite_stats(family, s, u, cell)


def main() -> None:
    roots = [
        0.5 - 0.5 / math.sqrt(2.0),
        0.5,
        0.5 + 0.5 / math.sqrt(2.0),
    ]
    families = [
        (make_random(7), []),
        (make_star(7), [0.5]),
        (make_worst_affine(7), roots),
        (make_smooth_fold(7), []),
        (make_two_bush_grid(16), [0.27, 0.73]),
        (make_book(10, 18), []),
        (make_regulus(55), []),
    ]
    for cell in (0.04, 0.02, 0.01):
        print(f"cell={cell:g}")
        print(
            "family                              s      u      A      B   "
            "edges excess             C4"
        )
        for family, extra in families:
            z = best_pair(family, cell, extra)
            print(
                f"{family.name:34s} {z['s']:5.3f} {z['u']:5.3f} "
                f"{z['A']:6d} {z['B']:6d} {z['unique_edges']:7d} "
                f"{z['edge_excess']:6d} {z['c4']:14d}"
            )
        print()


if __name__ == "__main__":
    main()
